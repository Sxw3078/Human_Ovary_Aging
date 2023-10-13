
library(optparse)
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
	    make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-b", "--RMbE"), help="remove batch effect by function IntegrateData , default %default",default=TRUE),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
        make_option(c("-r", "--mt_regress"), help="regree the mt gene in the SCTransform, default %default",default=FALSE),
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),

        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),

        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default=FALSE)

)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

################################


################################

if(! is.null(opt$workdir)){
  print(paste("change R work directory to ",opt$workdir,sep=""))
  setwd(opt$workdir)
 }
project_dir=getwd()


#read sample information matrix
#if more than one sample, should give the list csv used in cellranger aggr
cellranger_csv=paste(project_dir,"02_Cellranger","all.csv",sep="/")
if(! is.null(opt$aggr_csv)){
  print(paste("set aggr_csv to ",opt$aggr_csv,sep=""))
  cellranger_csv=opt$aggr_csv
 }

if(exists("cellranger_csv")){
    sample_list=read.csv(cellranger_csv,stringsAsFactors =F)
}else{
    stop("no cellranger_csv file was found\n")
}

if(nrow(sample_list) < 2){
  stop("this pipline need at least 2 sample\n")
}


#import DE contrast information
if(! is.null(opt$contrast)){
  print(paste("get compare contrast file ",opt$contrast,sep=""))
  sample_contrast=read.table(opt$contrast)
 }else{
#sample_contrast=switch(EXPR=as.character(nrow(sample_list)),"2"=combn(sample_list[,1],2), t(combn(sample_list[,1],2)))
 sample_contrast=t(combn(sample_list[,1],2))
 }

##################################
#library(SPOTlight)
library(Seurat)
library(dplyr)
library(SeuratData)
library(biomaRt)
library(RColorBrewer)
library(rlang)
library(ggplot2)
library(scales)
library(patchwork)
library(glmGamPoi)
library(reshape2)
library(future)
library(data.table)
library(harmony)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
###############################

if(opt$mt_regress){
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208,vars.to.regress ="percent.mt", verbose = opt$verbose)
    }
}else{
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose)   
    }
}
# 调色板
blue_green_red <- scale_color_gradient2(low = "blue", mid = "green", high = "red")
byr <- CustomPalette(low = '#0000FF', mid = "#FFFF00",high = "#FF0000", k = 7000)
bgr <- CustomPalette(low = '#0000FF', mid = "#00FF00",high = "#FF0000", k = 7000)
#br
#rblack
#whitered
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
getPalette(18)

#########基础空转分析#############

# 读取空转数据

spatial.ob=list()
#VariableFeatures=c()
for(x in 1:nrow(sample_list)){
    datadir=gsub("filtered_feature_bc_matrix.h5","",sample_list[x,2])
    spatial.ob[[x]]<-Load10X_Spatial(data.dir = datadir)
    #spatial.ob <- Seurat::SCTransform(spatial.ob, assay ='Spatial',verbose = opt$verbose)
    spatial.ob[[x]]@meta.data$orig.ident <as.character(sample_list[x,1])
    spatial.ob[[x]]$sample=as.character(sample_list[x,1])     #添加在meta.data slot里面
    spatial.ob[[x]]$orig.ident=as.character(sample_list[x,1]) 
    spatial.ob[[x]]$group=as.character(sample_list[x,4])
    names(spatial.ob[[x]]@images)=as.character(sample_list[x,1])

    spatial.ob[[x]][["percent.mt"]] <- PercentageFeatureSet(spatial.ob[[x]], pattern = opt$mt)
   # saveRDS(spatial.ob[[x]],paste(project_dir,"/","02_Cellranger","/",sample,".rds",sep=""))
}

#spatial.ob=merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)],add.cell.ids=sample_list[,1])
#names(spatial.ob@images)=as.character(sample_list[,1])
#spatial.ob$sample=as.character(sample_list[,1])
#VariableFeatures(spatial.ob) <-VariableFeatures
#cell_num_raw=table(spatial.ob$sample)
#spatial.ob[["percent.mt"]] <- PercentageFeatureSet(spatial.ob, pattern = opt$mt)
names(spatial.ob)=as.character(sample_list[,1])

# 3 QC与SCTransform去除批次效应
spatial.ob <- lapply(X =spatial.ob, FUN =SCTransform2)
if( opt$RMbE ){
    #combined <- SplitObject(spatial.ob, split.by = "sample")

    features <- SelectIntegrationFeatures(object.list = spatial.ob, nfeatures = 3000)
    combined <- PrepSCTIntegration(object.list = spatial.ob, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", anchor.features = features)
    combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT") 
    # Choose the features to use when integrating multiple datasets. 
    # This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. 
    # It returns the top scoring features by this ranking.
    VariableFeatures(combined.sct)=rownames(combined.sct@assays$integrated@scale.data)
    DefaultAssay(combined.sct) <- "integrated"
}else{

    combined.sct=merge(spatial.ob[[1]],spatial.ob[2:length(spatial.ob)],add.cell.ids=sample_list[,1])
    VariableFeatures(combined.sct)=rownames(combined.sct@assays$SCT@data)
    DefaultAssay(combined.sct) <- "SCT"
}

j=0
for (i in 1:length(combined.sct@images)){
    j=j+1
   if( nrow(combined.sct@images[[j]]@coordinates)==0){
       combined.sct@images[[j]]=NULL
       j=j-1
   }
}
names(combined.sct@images)=as.character(sample_list[,1])

seurat_qc_dir=paste(project_dir,"03_QC_SCT",sep="/")
if(!file.exists(seurat_qc_dir)){
    dir.create(seurat_qc_dir,recursive=TRUE)
}

plot2=plot4=list()

for ( x in 1:nrow(sample_list) ){

    plot2[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x,1],features = "nCount_Spatial") #+ theme(legend.position = "right")
    pic1=wrap_plots(plot2,ncol=2)

    plot4[[x]] <- SpatialFeaturePlot(combined.sct, images=sample_list[x,1],features = "nCount_SCT") #+ theme(legend.position = "right")
    pic2=wrap_plots(plot4,ncol=2)
}

plot1 <- VlnPlot(combined.sct, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), ncol = 3,group.by="sample", pt.size = 0.01)
plot3 <- VlnPlot(combined.sct, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), ncol = 3,group.by="sample", pt.size = 0.01)
## 图1：标准化前的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图
ggsave(plot1, width=12, height=7, filename=paste(seurat_qc_dir,"pre-QC-VlnPlot.pdf",sep="/"))
ggsave(plot3, width=12, height=7, filename=paste(seurat_qc_dir,"post-QC-VlnPlot.pdf",sep="/"))
## 图2：标准化后的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图
nrow1=ceiling(length(plot2)/2);nrow2=ceiling(length(plot4)/2)

ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.pdf",sep="/"))
ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.png",sep="/"))
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.pdf",sep="/"))
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.png",sep="/"))

fwrite(as.data.frame(combined.sct[["Spatial"]]@counts),paste(seurat_qc_dir,"cells_GeneCounts.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["SCT"]]@data),paste(seurat_qc_dir,"cells_SCT_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["integrated"]]@data),paste(seurat_qc_dir,"cells_integrated_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")
if(!file.exists(rds_dir)){
    dir.create(rds_dir,recursive=TRUE)
}

combined.sct <- RunPCA(combined.sct,npcs=100, verbose =  opt$verbose)
combined.sct <- RunHarmony(combined.sct, group.by.vars="orig.ident", assay.use="integrated")
## 选择后续分析使用的维度
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_1,"ElbowPlot.pdf",sep="/"),  dpi=1200)
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_1,"ElbowPlot.png",sep="/"),dpi=1200)
## 构建KNN图 
combined.sct <- FindNeighbors(combined.sct,reduction = "harmony", dims = 1:ndim)

saveRDS(spatial.ob,paste(rds_dir,"separated.sct.rds",sep="/"))

saveRDS(combined.sct,paste(rds_dir,"combined.sct.rds",sep="/"))
vctrs, cli, spatstat.data, spatstat.geom, caTools, cpp11, BiocManager, tiff, xfun, broom, Rhdf5lib, GenomeInfoDb, HDF5Array, sparseMatrixStats, IRanges, MatrixGenerics, beachmat, DelayedMatrixStats, DelayedArray, BiocParallel, S4Vectors, SingleCellExperiment, limma, spatstat.core, scuttle, edgeR, DropletUtils

# 4 降维、聚类、可视化
## 整合的分群分析
for (res in c(0.2,0.4,0.6,0.7,0.8,0.9,1)) {
    seurat_cluster_dir=paste(project_dir,paste0("04_Cluster_",res),sep="/")
    if(!file.exists(seurat_cluster_dir)){
        dir.create(seurat_cluster_dir,recursive=TRUE)
    }

    seurat_cluster_dir_1=paste(seurat_cluster_dir,"Integrated",sep="/")
    if(!file.exists(seurat_cluster_dir_1)){
        dir.create(seurat_cluster_dir_1,recursive=TRUE)
    }

    ## 降维

    ## 分群
    combined.sct <- FindClusters(combined.sct, resolution=res, verbose = opt$verbose)
    cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(combined.sct),"sample"=combined.sct$sample))),sample~cluster)
    fwrite(cluster_summary,paste(seurat_cluster_dir_1,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

    ## 分群结果展示UMAP VS TSNE
    combined.sct <- RunUMAP(combined.sct, reduction = "harmony", dims = 1:ndim)
    combined.sct <- RunTSNE(combined.sct, reduction = "harmony", dims = 1:ndim)
    plot1 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
    plot2 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
    plot3 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
    plot4 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
    pic3=wrap_plots(plot1, plot2,plot3,plot4,ncol=2)
    ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.pdf",sep="/"),dpi=1200)
    ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.png",sep="/"), dpi=1200)

    plot5 <- DimPlot(combined.sct, reduction = "umap",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
    plot6 <- DimPlot(combined.sct, reduction = "tsne",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
    nrow=ceiling(length(sample_list[,1])/2)
    ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.pdf",sep="/"),dpi=1200)
    ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.png",sep="/"),dpi=300)
    ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.pdf",sep="/"), dpi=1200)
    ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.png",sep="/"), dpi=300)

    ## 空间染色片上的分群结果
    plot7=list()
    width_add=0
    for ( x in 1:nrow(sample_list) ){
        plot7[[x]] <- SpatialDimPlot(combined.sct,image=sample_list[x,1], label = TRUE, label.size = 5,repel=TRUE)+ 

        labs(title=as.character(sample_list[x,1]))+
        theme(plot.title=element_text(color = "black", size = 20,hjust=0.5,vjust=0.5),
            legend.text=element_text(face="plain",size=12),
            legend.title = element_text(size=12,face = "bold",hjust=0.5,vjust=0.5),
            legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm') )

        plot7[[x]] <- plot7[[x]] + guides(fill= guide_legend(nrow =min(15,length(unique( plot7[[x]]$data$ident))), title = "Clusters"))
        names(plot7)[x]=sample_list[x,1]
        width_add=max(width_add,ceiling(length(unique( plot7[[x]]$data$ident))/15)/2)
    }
    pic4=wrap_plots(plot7,ncol=2)
    nrow=ceiling(length(sample_list[,1])/2)
    ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.pdf",sep="/"),  dpi=1200)
    ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.png",sep="/"),  dpi=300)


    
    ##  查找cluster biomarkers特征基因，获取其表达数据
    DefaultAssay(combined.sct) <- "SCT"
    combined.sct <- ScaleData(combined.sct, verbose =  opt$verbose)  ##scale.data将会被替换
    cluster_markers_all <- FindAllMarkers(object = combined.sct, 
                                                assay = "SCT",
                                                slot="data",
                                                verbose = opt$verbose, 
                                                only.pos = TRUE, 
                                                logfc.threshold = 0.25,
                                                min.pct = 0.1,
                                                test.use = "wilcox")   ###还有别的方法
    cluster_markers_all <- cluster_markers_all %>%  filter(p_val_adj <= 0.05 )
    fwrite(cluster_markers_all,paste(seurat_cluster_dir_1,"all_markers.xls",sep="/"),row.names=TRUE,col.names=TRUE,sep="\t")

    all_top10_markers=cluster_markers_all %>% top_n(n = 10, wt = avg_log2FC)

    plot8=VlnPlot(combined.sct, features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)
    ggsave(plot8,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7)
    ggsave(plot8,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7)

    plot9=FeaturePlot(combined.sct, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
    ggsave(plot9,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
    ggsave(plot9,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)

    plot10=DotPlot(combined.sct, features = rev(unique(all_top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
    ggsave(plot10,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15)
    ggsave(plot10,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15)

    if(!file.exists(paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",sep="/"))){
    dir.create(paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",sep="/"),recursive=TRUE)
    }

    for ( y in 1:length(all_top10_markers$gene)){
        plot11 <- SpatialFeaturePlot(combined.sct,slot="data",features=all_top10_markers$gene[y],combine=F,alpha=c(0.1,1),stroke = 0.5)
        pic5=wrap_plots(plot11,ncol=2,widths=3,heights=3)
        nrow=ceiling(length(sample_list[,1])/2)
        ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200)
        ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300)
    }

    DefaultAssay(combined.sct) <- "SCT"
    cluster_top10_markers=cluster_markers_all %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

    plot11=DoHeatmap(combined.sct, slot="scale.data", features = rev(unique(cluster_top10_markers$gene)))

    ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.pdf",sep="/"),width = 20,height = 14)
    ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.png",sep="/"),width = 20,height = 14)

    for( clust_num in unique(Idents(combined.sct))){
        seurat_cluster_dir_tmp=paste(seurat_cluster_dir_1,paste("cluster",clust_num,sep="_"),sep="/")
        if(!file.exists(seurat_cluster_dir_tmp)){
        dir.create(seurat_cluster_dir_tmp,recursive=TRUE)
        }
        cluster_markers=subset(cluster_markers_all,cluster==clust_num)
        C_num=nrow(cluster_markers)
        if(C_num==0){
            system(paste0("echo 'NO marker was found for this cluster!' > ",seurat_cluster_dir_tmp,paste("/cluster",clust_num,"markers.xls",sep="_")))
            next()
        }
        fwrite(cluster_markers,paste(seurat_cluster_dir_tmp,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)

        top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)


        plot12=VlnPlot(combined.sct, features = top10_markers$gene,pt.size = 0.1 ,ncol=5)

        ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.pdf",sep="/"),width =20,height = 7)
        ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.png",sep="/"),width =20,height = 7)

        plot13=FeaturePlot(combined.sct, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
        ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
        ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)

        plot14=DotPlot(combined.sct, features = rev(unique(top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()

        ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.pdf",sep="/"),width = 15,height =15)
        ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.png",sep="/"),width = 15,height =15)

        for ( y in 1:length(top10_markers$gene)){
            plot15 <- SpatialFeaturePlot(combined.sct,slot="data",features=top10_markers$gene[y],alpha=c(0.1,1),combine=FALSE,stroke = 0.1)
            
            pic6=wrap_plots(plot15,ncol=2,widths=3,heights=3)
            nrow=ceiling(length(sample_list[,1])/2)
            ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200)
            ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300)
        }
    }

    ## 分样本的分群分析

    plot7=list()
    width_add=0
    seurat_cluster_dir_2=paste(seurat_cluster_dir,"separated",sep="/")
    if(!file.exists(seurat_cluster_dir_2)){
        dir.create(seurat_cluster_dir_2,recursive=TRUE)
    }

    for(y in 1:nrow(sample_list)){
        seurat_cluster_dir_2_tmp=paste(seurat_cluster_dir_2,sample_list[y,1],sep="/")
        if(!file.exists(seurat_cluster_dir_2_tmp)){
            dir.create(seurat_cluster_dir_2_tmp,recursive=TRUE)
        }
        ## 降维
        spatial.ob[[y]]  <- RunPCA(spatial.ob[[y]],npcs=100, verbose =  opt$verbose)
        #spatial.ob[[y]] <- RunHarmony(spatial.ob[[y]], group.by.vars="orig.ident", assay.use="integrated")
        ## 选择后续分析使用的维度
        plot=ElbowPlot(spatial.ob[[y]],ndim=100,reduction="pca")
        ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
        ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_2_tmp,"ElbowPlot.pdf",sep="/"),  dpi=1200)
        ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_2_tmp,"ElbowPlot.png",sep="/"),dpi=1200)
        ## 构建KNN图 
        spatial.ob[[y]] <- FindNeighbors(spatial.ob[[y]],reduction = "pca", dims = 1:ndim)
        ## 分群
        spatial.ob[[y]] <- FindClusters(spatial.ob[[y]], resolution=res, verbose = opt$verbose)
        cluster_summary=dcast(as.data.frame(table(data.frame("cluster"=Idents(spatial.ob[[y]]),"sample"=spatial.ob[[y]]$sample))),sample~cluster)
        fwrite(cluster_summary,paste(seurat_cluster_dir_2_tmp,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

        ## 分群结果展示UMAP VS TSNE
        spatial.ob[[y]] <- RunUMAP(spatial.ob[[y]], reduction = "pca", dims = 1:ndim)
        spatial.ob[[y]] <- RunTSNE(spatial.ob[[y]], reduction = "pca", dims = 1:ndim)
        plot1 <- DimPlot(spatial.ob[[y]], reduction = "umap",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
        plot2 <- DimPlot(spatial.ob[[y]], reduction = "tsne",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
        #plot3 <- DimPlot(spatial.ob[[y]], reduction = "umap",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
        #plot4 <- DimPlot(spatial.ob[[y]], reduction = "tsne",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
        pic3=wrap_plots(plot1, plot2,ncol=2)
        ggsave(pic3, width=16, height=8,  filename=paste(seurat_cluster_dir_2_tmp,"UMAP_VS_tSNE.pdf",sep="/"),dpi=1200)
        ggsave(pic3, width=16, height=8,  filename=paste(seurat_cluster_dir_2_tmp,"UMAP_VS_tSNE.png",sep="/"), dpi=1200)

        #plot5 <- DimPlot(spatial.ob[[y]], reduction = "umap",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
        #plot6 <- DimPlot(spatial.ob[[y]], reduction = "tsne",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
        #nrow=ceiling(length(sample_list[,1])/2)
        #ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"UMAP_eachsample.pdf",sep="/"),dpi=1200)
        #ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"UMAP_eachsample.png",sep="/"),dpi=300)
        #ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"tSNE_eachsample.pdf",sep="/"), dpi=1200)
        #ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"tSNE_eachsample.png",sep="/"), dpi=300)

        ## 空间染色片上的分群结果 ##结果在循环结束后最后画出
        plot7[[y]] <- SpatialDimPlot(spatial.ob[[y]],image=sample_list[y,1], label = TRUE, label.size = 5,repel=TRUE)+ 

        labs(title=as.character(sample_list[y,1]))+
        theme(plot.title=element_text(color = "black", size = 20,hjust=0.5,vjust=0.5),
            legend.text=element_text(face="plain",size=12),
            legend.title = element_text(size=12,face = "bold",hjust=0.5,vjust=0.5),
            legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm') )

        plot7[[y]] <- plot7[[y]] + guides(fill= guide_legend(nrow =min(15,length(unique( plot7[[y]]$data$ident))), title = "Clusters"))
        names(plot7)[y]=sample_list[y,1]
        width_add=max(width_add,ceiling(length(unique( plot7[[y]]$data$ident))/15)/2)
    

        ##  查找cluster biomarkers特征基因，获取其表达数据
        DefaultAssay(spatial.ob[[y]] ) <- "SCT"
        spatial.ob[[y]]  <- ScaleData(spatial.ob[[y]] , verbose =  opt$verbose)  ##scale.data将会被替换
        cluster_markers_all <- FindAllMarkers(object = spatial.ob[[y]] , 
                                                    assay = "SCT",
                                                    slot="data",
                                                    verbose = opt$verbose, 
                                                    only.pos = TRUE, 
                                                    logfc.threshold = 0.25,
                                                    min.pct = 0.1,         ###不可太小，否则像一些特殊细胞的marker就找不到了。
                                                    test.use = "wilcox")   ###还有别的方法
        cluster_markers_all <- cluster_markers_all %>%  filter(p_val_adj <= 0.05 )
        fwrite(cluster_markers_all,paste(seurat_cluster_dir_2_tmp,"all_markers.xls",sep="/"),row.names=TRUE,col.names=TRUE,sep="\t")

        all_top10_markers=cluster_markers_all %>% top_n(n = 10, wt = avg_log2FC)

        plot8=VlnPlot(spatial.ob[[y]] , features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)
        ggsave(plot8,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7)
        ggsave(plot8,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7)

        plot9=FeaturePlot(spatial.ob[[y]] , features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
        ggsave(plot9,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
        ggsave(plot9,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)

        plot10=DotPlot(spatial.ob[[y]] , features = rev(unique(all_top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
        ggsave(plot10,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15)
        ggsave(plot10,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15)

        if(!file.exists(paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",sep="/"))){
            dir.create(paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",sep="/"),recursive=TRUE)
        }

        for ( z in 1:length(all_top10_markers$gene)){
            plot11 <- SpatialFeaturePlot(spatial.ob[[y]] ,slot="data",features=all_top10_markers$gene[z],combine=F,alpha=c(0.1,1),stroke = 0.5)
            pic5=wrap_plots(plot11,ncol=2,widths=3,heights=3)
            nrow=ceiling(length(sample_list[,1])/2)
            ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[z],"pdf",sep="."),sep="/"),  dpi=1200)
            ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[z],"png",sep="."),sep="/"),  dpi=300)
        }

        DefaultAssay(spatial.ob[[y]] ) <- "SCT"
        cluster_top10_markers=cluster_markers_all %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

        plot11=DoHeatmap(spatial.ob[[y]] , slot="scale.data", features = rev(unique(cluster_top10_markers$gene)))

        ggsave(plot11,filename=paste(seurat_cluster_dir_2_tmp,"top10_marker_echo_cluster_heatmap.pdf",sep="/"),width = 20,height = 14)
        ggsave(plot11,filename=paste(seurat_cluster_dir_2_tmp,"top10_marker_echo_cluster_heatmap.png",sep="/"),width = 20,height = 14)

        for( clust_num in unique(Idents(spatial.ob[[y]] ))){
            seurat_cluster_dir_tmp=paste(seurat_cluster_dir_2_tmp,paste("cluster",clust_num,sep="_"),sep="/")
            if(!file.exists(seurat_cluster_dir_tmp)){
                dir.create(seurat_cluster_dir_tmp,recursive=TRUE)
            }
            cluster_markers=subset(cluster_markers_all,cluster==clust_num)
            C_num=nrow(cluster_markers)
            if(C_num==0){
                system(paste0("echo 'NO marker was found for this cluster!' > ",seurat_cluster_dir_tmp,paste("/cluster",clust_num,"markers.xls",sep="_")))
                next()
            }
            fwrite(cluster_markers,paste(seurat_cluster_dir_tmp,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=T)

            top10_markers=cluster_markers %>%  top_n(n = 10, wt = avg_log2FC)


            plot12=VlnPlot(spatial.ob[[y]] , features = top10_markers$gene,pt.size = 0.1 ,ncol=5)

            ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.pdf",sep="/"),width =20,height = 7)
            ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.png",sep="/"),width =20,height = 7)

            plot13=FeaturePlot(spatial.ob[[y]] , features = top10_markers$gene, min.cutoff = "q9",ncol=5)
            ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7)
            ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7)

            plot14=DotPlot(spatial.ob[[y]] , features = rev(unique(top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()

            ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.pdf",sep="/"),width = 15,height =15)
            ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.png",sep="/"),width = 15,height =15)

            for ( z in 1:length(top10_markers$gene)){
                plot15 <- SpatialFeaturePlot(spatial.ob[[y]] ,slot="data",features=top10_markers$gene[z],alpha=c(0.1,1),combine=FALSE,stroke = 0.1)
                
                pic6=wrap_plots(plot15,ncol=2,widths=3,heights=3)
                nrow=ceiling(length(sample_list[,1])/2)
                ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[z],"pdf",sep="."),sep="/"),  dpi=1200)
                ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[z],"png",sep="."),sep="/"),  dpi=300)
            }
        }

    }

    pic4=wrap_plots(plot7,ncol=2)
    nrow=ceiling(length(sample_list[,1])/2)
    ggsave(pic4, width=13.5+width_add, height=6*nrow, filename=paste(seurat_cluster_dir_2,"SpatialDimPlot.pdf",sep="/"),  dpi=1200)
    ggsave(pic4, width=13.5+width_add, height=6*nrow, filename=paste(seurat_cluster_dir_2,"SpatialDimPlot.png",sep="/"),  dpi=300)

}

#################################################################################################################################
#################################################################################################################################
#5 find different gene between sample for each cluster
seurat_diff_cluster_dir=paste(project_dir,"05_DiffAnalysis_perCluster",sep="/")
if(!file.exists(seurat_diff_cluster_dir)){
   dir.create(seurat_diff_cluster_dir,recursive=TRUE)
}

combined.sct$cluster_sample <- paste(combined.sct$sample,paste("cluster",Idents(combined.sct),sep="") , sep = "_")
combined.sct$celltype <- Idents(combined.sct)
Idents(combined.sct) <- "cluster_sample"

sample_cluster_avg=AverageExpression(combined.sct,"SCT")$SCT
fwrite(data.frame(genes=rownames(sample_cluster_avg),sample_cluster_avg),paste(seurat_diff_cluster_dir,"avgExpression_per_sample_per_cluster.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

cluster_sample=as.vector(unique(Idents(combined.sct)))
nrow=ceiling(length(sample_list[,1])/2)
for(x in 1:nrow(sample_contrast)){

	compare=paste(sample_contrast[x,1],sample_contrast[x,2],sep="_vs_")
	compare_dir=paste(seurat_diff_cluster_dir,compare,sep="/")
	if(!file.exists(compare_dir)){
	   dir.create(compare_dir,recursive=TRUE)
	}
	for(sub_cluster in as.character(unique(combined.sct$celltype))){
		cp1=paste(sample_contrast[x,2],paste("cluster",sub_cluster,sep=""),sep="_")
		cp2=paste(sample_contrast[x,1],paste("cluster",sub_cluster,sep=""),sep="_")
		if(!cp1 %in% cluster_sample || ! cp2 %in% cluster_sample){
			next;
		}

		cluster_compare_dir=paste(compare_dir,paste("cluster",sub_cluster,sep=""),sep="/")
		if(!file.exists(cluster_compare_dir)){
		   dir.create(cluster_compare_dir,recursive=TRUE)
		}
		ncell1=CellsByIdentities(combined.sct,cp1)
		ncell2=CellsByIdentities(combined.sct,cp2)
		if (length(ncell1[[1]])<3 || length(ncell2[[1]])<3){
            system(paste0("echo 'too few cells of cluster ",sub_cluster," in sample ",sample_contrast[x,1], " or sample ", sample_contrast[x,2],".' >  ",paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep=""))
			next;
		}

		diff.cluster=FindMarkers(combined.sct, ident.1 = cp1, ident.2 = cp2, verbose = opt$verbose) %>%  filter(p_val_adj < 0.05 )
        if(nrow(diff.cluster)==0){
            system(paste0("echo 'no diff gene was found for this cluster' >  ",paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep=""))
        }else{
            fwrite(diff.cluster,paste(cluster_compare_dir,"diff_gene.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=T)
            sub_compare_cluster1=paste(sample_contrast[x,1],paste("cluster",unique(combined.sct$celltype),sep=""),sep="_")
            sub_compare_cluster2=paste(sample_contrast[x,2],paste("cluster",unique(combined.sct$celltype),sep=""),sep="_")
            diff.cluster$gene=row.names(diff.cluster)

            #top3_diff_gene=diff.cluster$gene[1:3]
            top10_diff.cluster=diff.cluster  %>% top_n(n = 10, wt = avg_log2FC)
            top10_diff_gene=top10_diff.cluster$gene

            ncol=length(top10_diff_gene)
            plots <- VlnPlot(combined.sct, features = top10_diff_gene, split.by = "sample", group.by = "celltype",pt.size = 0.01, combine = FALSE,idents = c(sub_compare_cluster1,sub_compare_cluster2))
            pic7=wrap_plots(plots = plots, ncol = 1)
            ggsave(pic7,filename=paste(cluster_compare_dir,"top10_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 17)
            ggsave(pic7,filename=paste(cluster_compare_dir,"top10_diffgene_exp_vilion.png",sep="/"),width = 8,height = 17)

            plot16=FeaturePlot(combined.sct, features = top10_diff_gene, split.by = "sample", max.cutoff = 3, pt.size = 0.01,  cols = c("grey", "red"))
            ggsave(plot16,filename=paste(cluster_compare_dir,"top10_diffgene_umap.pdf",sep="/"),width = 4*nrow,height = 2.4*ncol)
            ggsave(plot16,filename=paste(cluster_compare_dir,"top10_diffgene_umap.png",sep="/"),width = 4*nrow,height = 2.4*ncol)
		}

	}
}




