
suppressPackageStartupMessages(library(optparse))
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
        make_option(c("-l", "--resolution"),type="double", help="resolution value for FindClusters, default %default",default="1.0"),
        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),
        make_option(c("-d", "--rds"), help="rds for input, default %default",default="suerat_Rdata/separated.sct.Rdata"),
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
#suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(harmony))
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
#getPalette(18)

#########基础空转分析#############

# 读取空转数据

load(opt$rds)

# 4 降维、聚类、可视化
## 整合的分群分析

res=opt$resolution
seurat_cluster_dir=paste(project_dir,paste0("04_Cluster_",res),sep="/")
if(!file.exists(seurat_cluster_dir)){
    dir.create(seurat_cluster_dir,recursive=TRUE)
}

seurat_cluster_dir_2=paste(seurat_cluster_dir,"separated",sep="/")
if(!file.exists(seurat_cluster_dir_2)){
    dir.create(seurat_cluster_dir_2,recursive=TRUE)
}
## 分样本的分群分析

plot7=list()
width_add=0

for(y in 1:nrow(sample_list)){
    seurat_cluster_dir_2_tmp=paste(seurat_cluster_dir_2,sample_list[y,1],sep="/")
    sample=sample_list[y,1]
    if(!file.exists(seurat_cluster_dir_2_tmp)){
        dir.create(seurat_cluster_dir_2_tmp,recursive=TRUE)
    }
    ## 降维
    spatial.ob[[sample]]  <- RunPCA(spatial.ob[[sample]],npcs=100, verbose =  opt$verbose)
    #spatial.ob[[sample]] <- RunHarmony(spatial.ob[[sample]], group.by.vars="orig.ident", assay.use="integrated")
    ## 选择后续分析使用的维度
    plot=ElbowPlot(spatial.ob[[sample]],ndim=100,reduction="pca")
    ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
    ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_2_tmp,"ElbowPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
    ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir_2_tmp,"ElbowPlot.png",sep="/"),dpi=1200,limitsize = FALSE)
    ## 构建KNN图 
    spatial.ob[[sample]] <- FindNeighbors(spatial.ob[[sample]],reduction = "pca", dims = 1:ndim, verbose =  opt$verbose)
    ## 分群
    spatial.ob[[sample]] <- FindClusters(spatial.ob[[sample]], resolution=res, verbose = opt$verbose)
    cluster_summary=reshape2::dcast(as.data.frame(table(data.frame("cluster"=Idents(spatial.ob[[sample]]),"sample"=spatial.ob[[sample]]$sample))),sample~cluster)
    fwrite(cluster_summary,paste(seurat_cluster_dir_2_tmp,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

    ## 分群结果展示UMAP VS TSNE
    spatial.ob[[sample]] <- RunUMAP(spatial.ob[[sample]], reduction = "pca", dims = 1:ndim, verbose =  opt$verbose)
    spatial.ob[[sample]] <- RunTSNE(spatial.ob[[sample]], reduction = "pca", dims = 1:ndim, verbose =  opt$verbose)
    plot1 <- DimPlot(spatial.ob[[sample]], reduction = "umap",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
    plot2 <- DimPlot(spatial.ob[[sample]], reduction = "tsne",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
    #plot3 <- DimPlot(spatial.ob[[sample]], reduction = "umap",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
    #plot4 <- DimPlot(spatial.ob[[sample]], reduction = "tsne",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
    pic3=wrap_plots(plot1, plot2,ncol=2)
    ggsave(pic3, width=16, height=8,  filename=paste(seurat_cluster_dir_2_tmp,"UMAP_VS_tSNE.pdf",sep="/"),dpi=1200,limitsize = FALSE)
    ggsave(pic3, width=16, height=8,  filename=paste(seurat_cluster_dir_2_tmp,"UMAP_VS_tSNE.png",sep="/"), dpi=1200,limitsize = FALSE)

    #plot5 <- DimPlot(spatial.ob[[sample]], reduction = "umap",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
    #plot6 <- DimPlot(spatial.ob[[sample]], reduction = "tsne",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
    #nrow=ceiling(length(sample_list[,1])/2)
    #ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"UMAP_eachsample.pdf",sep="/"),dpi=1200,limitsize = FALSE)
    #ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"UMAP_eachsample.png",sep="/"),dpi=300,limitsize = FALSE)
    #ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"tSNE_eachsample.pdf",sep="/"), dpi=1200,limitsize = FALSE)
    #ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_2,"tSNE_eachsample.png",sep="/"), dpi=300,limitsize = FALSE)

    ## 空间染色片上的分群结果 ##结果在循环结束后最后画出
    plot7[[y]] <- SpatialDimPlot(spatial.ob[[sample]],images=sample_list[y,1], label = TRUE, label.size = 5,repel=TRUE)+ 

    labs(title=as.character(sample_list[y,1]))+
    theme(plot.title=element_text(color = "black", size = 20,hjust=0.5,vjust=0.5),
        legend.text=element_text(face="plain",size=12),
        legend.title = element_text(size=12,face = "bold",hjust=0.5,vjust=0.5),
        legend.key.width=unit(0.8,'cm'),legend.key.height=unit(0.8,'cm') )

    plot7[[y]] <- plot7[[y]] + guides(fill= guide_legend(nrow =min(15,length(unique( plot7[[y]]$data$ident))), title = "Clusters"))
    names(plot7)[y]=sample_list[y,1]
    width_add=max(width_add,ceiling(length(unique( plot7[[y]]$data$ident))/15)/2)


    ##  查找cluster biomarkers特征基因，获取其表达数据
    DefaultAssay(spatial.ob[[sample]] ) <- "SCT"
    spatial.ob[[sample]] <- PrepSCTFindMarkers(spatial.ob[[sample]])
    spatial.ob[[sample]]  <- ScaleData(spatial.ob[[sample]] , verbose =  opt$verbose)  ##scale.data将会被替换
    cluster_markers_all <- FindAllMarkers(object = spatial.ob[[sample]] , 
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

    plot8=VlnPlot(spatial.ob[[sample]] , features = all_top10_markers$gene,pt.size = 0.1 ,ncol=5)
    ggsave(plot8,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
    ggsave(plot8,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

    plot9=FeaturePlot(spatial.ob[[sample]] , features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
    ggsave(plot9,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
    ggsave(plot9,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

    plot10=DotPlot(spatial.ob[[sample]] , features = rev(unique(all_top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
    ggsave(plot10,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15,limitsize = FALSE)
    ggsave(plot10,filename=paste(seurat_cluster_dir_2_tmp,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15,limitsize = FALSE)

    if(!file.exists(paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",sep="/"))){
        dir.create(paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",sep="/"),recursive=TRUE)
    }

    for ( z in 1:length(all_top10_markers$gene)){
        plot11 <- SpatialFeaturePlot(spatial.ob[[sample]] ,slot="data",features=all_top10_markers$gene[z],combine=F,alpha=c(0.1,1),stroke = 0.5)
        pic5=wrap_plots(plot11,ncol=2,widths=3,heights=3)
        nrow=ceiling(length(sample_list[,1])/2)
        ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[z],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
        ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_2_tmp,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[z],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
    }

    DefaultAssay(spatial.ob[[sample]] ) <- "SCT"
    cluster_top10_markers=cluster_markers_all %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

    plot11=DoHeatmap(spatial.ob[[sample]] , slot="scale.data", features = rev(unique(cluster_top10_markers$gene)))

    ggsave(plot11,filename=paste(seurat_cluster_dir_2_tmp,"top10_marker_echo_cluster_heatmap.pdf",sep="/"),width = 20,height = 14,limitsize = FALSE)
    ggsave(plot11,filename=paste(seurat_cluster_dir_2_tmp,"top10_marker_echo_cluster_heatmap.png",sep="/"),width = 20,height = 14,limitsize = FALSE)

    for( clust_num in unique(Idents(spatial.ob[[sample]] ))){
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


        plot12=VlnPlot(spatial.ob[[sample]] , features = top10_markers$gene,pt.size = 0.1 ,ncol=5)

        ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.pdf",sep="/"),width =20,height = 7,limitsize = FALSE)
        ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.png",sep="/"),width =20,height = 7,limitsize = FALSE)

        plot13=FeaturePlot(spatial.ob[[sample]] , features = top10_markers$gene, min.cutoff = "q9",ncol=5)
        ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
        ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

        plot14=DotPlot(spatial.ob[[sample]] , features = rev(unique(top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()

        ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.pdf",sep="/"),width = 15,height =15,limitsize = FALSE)
        ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.png",sep="/"),width = 15,height =15,limitsize = FALSE)

        for ( z in 1:length(top10_markers$gene)){
            plot15 <- SpatialFeaturePlot(spatial.ob[[sample]] ,slot="data",features=top10_markers$gene[z],alpha=c(0.1,1),combine=FALSE,stroke = 0.1)
            
            pic6=wrap_plots(plot15,ncol=2,widths=3,heights=3)
            nrow=ceiling(length(sample_list[,1])/2)
            ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[z],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
            ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[z],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
        }
    }

}

pic4=wrap_plots(plot7,ncol=2)
nrow=ceiling(length(sample_list[,1])/2)
ggsave(pic4, width=13.5+width_add, height=6*nrow, filename=paste(seurat_cluster_dir_2,"SpatialDimPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(pic4, width=13.5+width_add, height=6*nrow, filename=paste(seurat_cluster_dir_2,"SpatialDimPlot.png",sep="/"),  dpi=300,limitsize = FALSE)



