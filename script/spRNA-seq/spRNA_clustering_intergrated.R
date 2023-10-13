
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
        make_option(c("-d", "--rds"), help="rds for input, default %default",default="suerat_Rdata/combined.sct.rds"),
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

combined.sct=readRDS(opt$rds)

# 4 降维、聚类、可视化
## 整合的分群分析
res=opt$resolution
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
cluster_summary=reshape2::dcast(as.data.frame(table(data.frame("cluster"=Idents(combined.sct),"sample"=combined.sct$sample))),sample~cluster)
fwrite(cluster_summary,paste(seurat_cluster_dir_1,"cluster_summary.xls",sep="/"),sep="\t",quote=F,row.names=F,col.names=T)

## 分群结果展示UMAP VS TSNE
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
combined.sct <- RunUMAP(combined.sct, reduction = "harmony", dims = 1:ndim, verbose = opt$verbose)
combined.sct <- RunTSNE(combined.sct, reduction = "harmony", dims = 1:ndim, verbose = opt$verbose)
plot1 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
plot2 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=0.1)
plot3 <- DimPlot(combined.sct, reduction = "umap",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
plot4 <- DimPlot(combined.sct, reduction = "tsne",shuffle=TRUE, group.by = "sample",label = FALSE,seed=19900208,pt.size=0.1)
pic3=wrap_plots(plot1, plot2,plot3,plot4,ncol=2)
ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.pdf",sep="/"),dpi=1200,limitsize = FALSE)
ggsave(pic3, width=16, height=15,  filename=paste(seurat_cluster_dir_1,"UMAP_VS_tSNE.png",sep="/"), dpi=300,limitsize = FALSE)

plot5 <- DimPlot(combined.sct, reduction = "umap",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
plot6 <- DimPlot(combined.sct, reduction = "tsne",split.by = "sample",shuffle=TRUE, label = TRUE,seed=19900208,pt.size=1,ncol=2)
nrow=ceiling(length(sample_list[,1])/2)
ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.pdf",sep="/"),dpi=1200,limitsize = FALSE)
ggsave(plot5, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"UMAP_eachsample.png",sep="/"),dpi=300,limitsize = FALSE)
ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.pdf",sep="/"), dpi=1200,limitsize = FALSE)
ggsave(plot6, width=16, height=7.5*nrow,  filename=paste(seurat_cluster_dir_1,"tSNE_eachsample.png",sep="/"), dpi=300,limitsize = FALSE)

## 空间染色片上的分群结果
plot7=list()
width_add=0
for ( x in 1:nrow(sample_list) ){
    plot7[[x]] <- SpatialDimPlot(combined.sct,images=sample_list[x,1], label = TRUE, label.size = 5,repel=TRUE)+ 

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
ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(pic4, width=13.5+width_add, height=6*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialDimPlot.png",sep="/"),  dpi=300,limitsize = FALSE)



##  查找cluster biomarkers特征基因，获取其表达数据
DefaultAssay(combined.sct) <- "SCT"
combined.sct <- PrepSCTFindMarkers(combined.sct)
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
ggsave(plot8,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_vilion.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
ggsave(plot8,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_vilion.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

plot9=FeaturePlot(combined.sct, features = all_top10_markers$gene, min.cutoff = "q9",ncol=5)
ggsave(plot9,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
ggsave(plot9,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

plot10=DotPlot(combined.sct, features = rev(unique(all_top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()
ggsave(plot10,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_exp_pct.pdf",sep="/"),width = 15,height = 15,limitsize = FALSE)
ggsave(plot10,filename=paste(seurat_cluster_dir_1,"allmarkers_top10_exp_pct.png",sep="/"),width = 15,height = 15,limitsize = FALSE)

if(!file.exists(paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",sep="/"))){
dir.create(paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",sep="/"),recursive=TRUE)
}

for ( y in 1:length(all_top10_markers$gene)){
    plot11 <- SpatialFeaturePlot(combined.sct,slot="data",features=all_top10_markers$gene[y],combine=F,alpha=c(0.1,1),stroke = 0.5)
    pic5=wrap_plots(plot11,ncol=2,widths=3,heights=3)
    nrow=ceiling(length(sample_list[,1])/2)
    ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
    ggsave(pic5, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_1,"SpatialFeaturePlot_all_top10",paste("all_top10",all_top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
}

DefaultAssay(combined.sct) <- "SCT"
cluster_top10_markers=cluster_markers_all %>% group_by(cluster)%>% top_n(n = 10, wt = avg_log2FC)

plot11=DoHeatmap(combined.sct, slot="scale.data", features = rev(unique(cluster_top10_markers$gene)))

ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.pdf",sep="/"),width = 20,height = 14,limitsize = FALSE)
ggsave(plot11,filename=paste(seurat_cluster_dir_1,"top10_marker_echo_cluster_heatmap.png",sep="/"),width = 20,height = 14,limitsize = FALSE)

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

    ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.pdf",sep="/"),width =20,height = 7,limitsize = FALSE)
    ggsave(plot12,filename=paste(seurat_cluster_dir_tmp,"top10_vilion.png",sep="/"),width =20,height = 7,limitsize = FALSE)

    plot13=FeaturePlot(combined.sct, features = top10_markers$gene, min.cutoff = "q9",ncol=5)
    ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.pdf",sep="/"),width = 20,height = 7,limitsize = FALSE)
    ggsave(plot13,filename=paste(seurat_cluster_dir_tmp,"top10_cell_exp_distribution.png",sep="/"),width = 20,height = 7,limitsize = FALSE)

    plot14=DotPlot(combined.sct, features = rev(unique(top10_markers$gene)), cols = rainbow(nrow(sample_list)), dot.scale = 8,split.by = "sample") + RotatedAxis()

    ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.pdf",sep="/"),width = 15,height =15,limitsize = FALSE)
    ggsave(plot14,filename=paste(seurat_cluster_dir_tmp,"top10_exp_pct.png",sep="/"),width = 15,height =15,limitsize = FALSE)

    for ( y in 1:length(top10_markers$gene)){
        plot15 <- SpatialFeaturePlot(combined.sct,slot="data",features=top10_markers$gene[y],alpha=c(0.1,1),combine=FALSE,stroke = 0.1)
        
        pic6=wrap_plots(plot15,ncol=2,widths=3,heights=3)
        nrow=ceiling(length(sample_list[,1])/2)
        ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"pdf",sep="."),sep="/"),  dpi=1200,limitsize = FALSE)
        ggsave(pic6, width=6, height=3*nrow,  filename=paste(seurat_cluster_dir_tmp,paste("top10",top10_markers$gene[y],"png",sep="."),sep="/"),  dpi=300,limitsize = FALSE)
    }
}



#################################################################################################################################
#################################################################################################################################
#5 find different gene between sample for each cluster
seurat_diff_cluster_dir=paste(project_dir,paste0("05_DiffAnalysis_perCluster_",res),sep="/")
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
            ggsave(pic7,filename=paste(cluster_compare_dir,"top10_diffgene_exp_vilion.pdf",sep="/"),width = 8,height = 17,limitsize = FALSE)
            ggsave(pic7,filename=paste(cluster_compare_dir,"top10_diffgene_exp_vilion.png",sep="/"),width = 8,height = 17,limitsize = FALSE)

            plot16=FeaturePlot(combined.sct, features = top10_diff_gene, split.by = "sample", max.cutoff = 3, pt.size = 0.01,  cols = c("grey", "red"))
            ggsave(plot16,filename=paste(cluster_compare_dir,"top10_diffgene_umap.pdf",sep="/"),width = 4*nrow,height = 2.4*ncol,limitsize = FALSE)
            ggsave(plot16,filename=paste(cluster_compare_dir,"top10_diffgene_umap.png",sep="/"),width = 4*nrow,height = 2.4*ncol,limitsize = FALSE)
		}

	}
}




