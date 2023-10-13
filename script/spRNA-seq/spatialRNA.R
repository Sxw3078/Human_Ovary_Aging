
suppressPackageStartupMessages(library(optparse))
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
	    make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-b", "--RMbE"), help="remove batch effect by function IntegrateData , default %default",default=TRUE),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT-"),
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
#suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
#suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
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
options(future.globals.maxSize = 100000 * 1024^2)
options(future.seed=TRUE)
###############################

if(opt$mt_regress){
    SCTransform2=function(x){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208,vars.to.regress ="percent.mt", verbose = opt$verbose, vst.flavor = "v2")
    }
}else{
    SCTransform2=function(x){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose,vst.flavor = "v2")
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

spatial.ob=list()
#VariableFeatures=c()
for(x in 1:nrow(sample_list)){
    datadir=gsub("filtered_feature_bc_matrix.h5","",sample_list[x,2])
    spatial.ob[[x]]<-suppressWarnings(Load10X_Spatial(data.dir = datadir))
    #spatial.ob <- Seurat::SCTransform(spatial.ob, assay ='Spatial',verbose = opt$verbose)
  #  spatial.ob[[x]]@meta.data$orig.ident<-as.character(sample_list[x,1])
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

# 按照分组进行SCT
if(2==2){
    spatial.group=list()
    for (group in unique(sample_list[,"group"])){
        group.samples=sample_list[which(sample_list[,"group"]==group),"sample"]
        spatial.group[[group]]=merge(spatial.ob[[group.samples[1]]],spatial.ob[group.samples[-1]],add.cell.ids= group.samples) %>% SCTransform(.,method= "glmGamPoi", assay = "Spatial",new.assay.name = "SCT1", return.only.var.genes = FALSE, seed.use=19900208,vars.to.regress ="percent.mt", verbose = opt$verbose, vst.flavor = "v2")
    }
    spatial.combined=merge(spatial.group[[1]],spatial.group[2:length(spatial.group)])
    combined.sct= SCTransform(spatial.combined,method= "glmGamPoi", assay = "SCT1", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose,vst.flavor = "v2")
}else{
    spatial.combined=merge(spatial.ob[[1]],spatial.ob[2:length(spatial.ob)])
    combined.sct= SCTransform(spatial.combined,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose,vst.flavor = "v2")
}

if( opt$RMbE ){

    spatial.ob <- SplitObject(combined.sct, split.by = "sample")

    features <- SelectIntegrationFeatures(object.list = spatial.ob, nfeatures = 3000,verbose = opt$verbose)
    combined <- PrepSCTIntegration(object.list = spatial.ob, anchor.features = features,verbose = opt$verbose)
    anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", anchor.features = features,verbose = opt$verbose)
    combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",verbose = opt$verbose) 
    # Choose the features to use when integrating multiple datasets. 
    # This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. 
    # It returns the top scoring features by this ranking.
    VariableFeatures(combined.sct)=rownames(combined.sct@assays$integrated@scale.data)
    DefaultAssay(combined.sct) <- "integrated"
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
ggsave(plot1, width=12, height=7, filename=paste(seurat_qc_dir,"pre-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
ggsave(plot3, width=12, height=7, filename=paste(seurat_qc_dir,"post-QC-VlnPlot.pdf",sep="/"),limitsize = FALSE)
## 图2：标准化后的空间表达分布图  molecular counts across spots before normalization   原始数据空间表达分布图

nrow1=ceiling(length(plot2)/2);
nrow2=ceiling(length(plot4)/2)

ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic1, width=10, height=5*nrow1, filename=paste(seurat_qc_dir,"pre-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.pdf",sep="/"),limitsize = FALSE)
ggsave(pic2, width=10, height=5*nrow2, filename=paste(seurat_qc_dir,"post-QC-SpatialFeaturePlot.png",sep="/"),limitsize = FALSE)

fwrite(as.data.frame(combined.sct[["Spatial"]]@counts),paste(seurat_qc_dir,"cells_GeneCounts.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["SCT"]]@data),paste(seurat_qc_dir,"cells_SCT_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

fwrite(as.data.frame(combined.sct[["integrated"]]@data),paste(seurat_qc_dir,"cells_integrated_expression.xls.gzip",sep="/"),sep="\t",quote=F,row.names=T,col.names=T,compress="gzip")

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")
if(!file.exists(rds_dir)){
    dir.create(rds_dir,recursive=TRUE)
}
seurat_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
if(!file.exists(seurat_cluster_dir)){
    dir.create(seurat_cluster_dir,recursive=TRUE)
}

combined.sct <- RunPCA(combined.sct,npcs=100, verbose =  opt$verbose)
combined.sct <- RunHarmony(combined.sct, group.by.vars="orig.ident", assay.use="integrated", verbose =  opt$verbose)
## 选择后续分析使用的维度
plot=ElbowPlot(combined.sct,ndim=100,reduction="harmony")
ndim=min(which(plot$data$stdev<(0.05*(max(plot$data$stdev)-min(plot$data$stdev))+min(plot$data$stdev))))
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.pdf",sep="/"),  dpi=1200,limitsize = FALSE)
ggsave(plot, width=16, height=10, units = "cm", filename=paste(seurat_cluster_dir,"ElbowPlot.png",sep="/"),dpi=300,limitsize = FALSE)
## 构建KNN图 
combined.sct <- FindNeighbors(combined.sct,reduction = "harmony", dims = 1:ndim, verbose =  opt$verbose)

save(spatial.ob,file=paste(rds_dir,"separated.sct.Rdata",sep="/"))

saveRDS(combined.sct,paste(rds_dir,"combined.sct.rds",sep="/"))









if(1==2){


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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
combined.sct=readRDS("suerat_Rdata/combined.sct.rds")
#combined.sct2=readRDS("../2021_03_19/suerat_Rdata/combined.sct.rds")
df=combined.sct@meta.data

samples=c("Y1_1", "M1_1", "O1_1", "O1_2", "Y2_1", "M2_1", "O2_1", "Y3_1", "Y3_2", "Y3_3", "M3_1", "M3_2", "M3_3", "O3_1", "O3_2")
df1=df[which(str_sub(df$orig.ident,2,2)==1),]
df2=df[which(str_sub(df$orig.ident,2,2)==2),]
df3=df[which(str_sub(df$orig.ident,2,2)==3),]

meangenes=mediangenes=meanreads=medianreads=nspot=c()
for (sample in samples){
    dff=df[which(df$orig.ident==sample),]
    dfs=df[which(str_sub(df$orig.ident,2,2)==str_sub(sample,2,2)),]

    plot1=ggplot(data=dff,mapping=aes(x=nCount_Spatial,y=nFeature_Spatial))+
        geom_point(color="red",fill="red")+
        xlab("Reads_per_Spot")+
        ylab("Geness_per_Spot")+
        theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(size=1.8))+
        scale_x_continuous(expand=expand_scale(mult = c(0.02, 0.05)),limits = c(0.0, max(dfs$nCount_Spatial)))+
        scale_y_continuous(expand=expand_scale(mult = c(0.02, 0.05)),,limits = c(0.0, max(dfs$nFeature_Spatial)))

    plot2=ggplot(data=dff,mapping=aes(x=nCount_Spatial))+
        geom_histogram(color="white",fill="white", bins = 40 )+ 
        theme(panel.grid.minor=element_blank(),
                panel.grid.major.x=element_line(size=1.8), 
                panel.grid.major.y=element_blank(),
                axis.text=element_blank(),
                axis.title.y=element_blank(),
                axis.ticks=element_blank()
                )+ xlab("")+
        scale_x_continuous(expand=expand_scale(mult = c(0.02, 0.05)),limits = c(0.0, max(dfs$nCount_Spatial)))+
        scale_y_continuous(expand=expand_scale(mult = c(0, 0)))


    plot3=ggplot(data=dff,mapping=aes(x=nFeature_Spatial))+
        geom_histogram(color="white",fill="white", bins = 40)+ 
        theme(panel.grid.minor=element_blank(),
                panel.grid.major.y=element_line(size=1.8),
                panel.grid.major.x=element_blank(),
                axis.text=element_blank(),
                axis.title.x=element_blank(),
                axis.ticks=element_blank()
            )+ xlab("")+
        coord_flip()+
        scale_x_continuous(expand=expand_scale(mult = c(0.02, 0.05)),limits = c(0.0, max(dfs$nFeature_Spatial)))+
        scale_y_continuous(expand=expand_scale(mult = c(0, 0)))
    
    layout <- "
    BBBBBBB#
    AAAAAAAC
    AAAAAAAC
    AAAAAAAC
    AAAAAAAC
    AAAAAAAC
    AAAAAAAC
    AAAAAAAC
    "
    pic = wrap_plots(A = plot1, B = plot2, C = plot3, design = layout)

    ggsave(pic, width=5, height=5, filename=paste0(sample,"_basic_stats.pdf"))
    meangenes=c(meangenes,mean(dff$nFeature_Spatial))
    mediangenes=c(mediangenes,median(dff$nFeature_Spatial))
    meanreads=c(meanreads,mean(dff$nCount_Spatial))
    medianreads=c(medianreads,median(dff$nCount_Spatial))
    nspot=c(nspot,length(dff$nCount_Spatial))
}




basic_stat=data.frame("sample"=samples,"mean_genes"=meangenes,"median_genes"=mediangenes,"mean_reads"=meanreads,"median_reads"=medianreads,"nspot"=nspot)
write.table(basic_stat,"basic_stats.tsv",sep="\t")






#combined.sct=readRDS("suerat_Rdata/combined.sct.rds")
#combined.sct2=readRDS("../2021_03_19/suerat_Rdata/combined.sct.rds")

for (sample in samples){

    sub_dataset=subset(combined.sct, subset=orig.ident == eval(sample))


    plot <- VlnPlot(sub_dataset, "percent.mt", group.by="sample", pt.size = 0.01, adjust = 1, y.max=25,cols="white")+
                xlab("")+ ylim(-0.1,25)+
                ggtitle(sample)+
                theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank() ,legend.position = 'none'     ) 

   

    ggsave(plot, width=2, height=5, filename=paste0(sample,"_mt.pdf"))

}




##基本信息表达
for (sample in samples){
    spatial.ob=subset(combined.sct, subset=orig.ident == eval(sample))
    plot0=SpatialPlot(
              object = spatial.ob,
              features = "nFeature_Spatial",
              #alpha = c(0.1, 1),
              crop=TRUE,
              pt.size.factor = 0,
              stroke=0,
              images=sample,
              image.alpha=1,
              min.cutoff=0.001
              )  +theme(legend.position="none")

    plot1=SpatialPlot(
              object = spatial.ob,
              features = "nFeature_Spatial",
              #alpha = c(0.1, 1),
              crop=TRUE,
              pt.size.factor = 1.8,
              stroke=0,
              images=sample,
              image.alpha=0,
              min.cutoff=0.001
              )  +theme(legend.position="right")+
              guides(fill=guide_colourbar(title="nGene",barwidth=0.8,barheight=5,label.theme=element_text(size=12),title.theme=element_text(size=14)))

    plot2=SpatialPlot(
              object = spatial.ob,
              features = "nCount_Spatial",
              #alpha = c(0.1, 1),
              crop=TRUE,
              pt.size.factor = 1.8,
              stroke=0,
              images=sample,
              image.alpha=0,
              min.cutoff=0.001
              )   +theme(legend.position="right")+
              guides(fill=guide_colourbar(title="nUMI",barwidth=0.8,barheight=5,label.theme=element_text(size=12),title.theme=element_text(size=14)))
   
    plot3=plot0+plot2+plot1
 
    ggsave(plot3, width=15, height=5,dpi=300,units ="in",bg="white", filename=paste0("03_QC_SCT/SpatialFeaturePlot_",sample,"_gene_and_UMI.png"),limitsize = FALSE)

}


  ##基因表达
  suppressPackageStartupMessages(library(doParallel))

for (sample in samples){
  spatial.ob=subset(combined.sct, subset=orig.ident == eval(sample))
  
  DefaultAssay(spatial.ob) <- "SCT"
   for (i in c("ACTA2","AMH","APOA1","CCL5","CDKN1A","CDKN1B","CDKN2A","CDKN2B","DCN","FIGLA","GSTA1",
     "HSD17B1","IFI30","IL1B","IL6","IL7R","KLRB1","MUSTN1","NKG7","STAG3","STAG","TM4SF1","TUBB8","TYROBP","VWF","ZP3")){
    if( ! i %in% spatial.ob@assays$SCT@data@Dimnames[[1]]){print("not in");next()}
    plot2=SpatialFeaturePlot(
              object = spatial.ob,
              slot="data",
              features = i,
              #alpha = c(0.1, 1),
              crop=TRUE,
              image=sample,
              pt.size.factor = 1.2,
              stroke=0.1,
              min.cutoff=0.001
              )  
    plot3=plot2+ 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
    if(!file.exists(paste0("genes_exp/",sample))){
        dir.create(paste0("genes_exp/",sample))
    }
    #ggsave(plot2, width=5, height=5,dpi=900,units ="in",bg="white", filename=paste0("genes_exp/SpatialFeaturePlot_",i,".color1.pdf"),limitsize = FALSE)
    ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename=paste0("genes_exp/",sample,"/SpatialFeaturePlot_",i,".png"),limitsize = FALSE)
    }
}


}



if(1==2){

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
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(stringr))
    combined.sct=readRDS("suerat_Rdata/combined.sct.rds")
    #combined.sct2=readRDS("../2021_03_19/suerat_Rdata/combined.sct.rds")
    df=combined.sct@meta.data

   

        
    spatial.ob=list()
    #VariableFeatures=c()
    for(x in 1:nrow(sample_list)){
        datadir=gsub("filtered_feature_bc_matrix.h5","",sample_list[x,2])
        spatial.ob[[x]]<-suppressWarnings(Load10X_Spatial(data.dir = datadir))
        #spatial.ob <- Seurat::SCTransform(spatial.ob, assay ='Spatial',verbose = opt$verbose)
        spatial.ob[[x]]@meta.data$orig.ident<-as.character(sample_list[x,1])
        spatial.ob[[x]]$sample=as.character(sample_list[x,1])     #添加在meta.data slot里面
        spatial.ob[[x]]$orig.ident=as.character(sample_list[x,1]) 
        spatial.ob[[x]]$group=as.character(sample_list[x,4])
        names(spatial.ob[[x]]@images)=as.character(sample_list[x,1])

        spatial.ob[[x]][["percent.mt"]] <- PercentageFeatureSet(spatial.ob[[x]], pattern = opt$mt)
    # saveRDS(spatial.ob[[x]],paste(project_dir,"/","02_Cellranger","/",sample,".rds",sep=""))
    }



    names(spatial.ob)=as.character(sample_list[,1])

    # 3 QC与SCTransform去除批次效应
    spatial.ob <- lapply(X =spatial.ob, FUN =SCTransform2)
    if( opt$RMbE ){
        #combined <- SplitObject(spatial.ob, split.by = "sample")

        features <- SelectIntegrationFeatures(object.list = spatial.ob, nfeatures = 3000,verbose = opt$verbose)
        combined <- PrepSCTIntegration(object.list = spatial.ob, anchor.features = features,verbose = opt$verbose)
        anchors <- FindIntegrationAnchors(object.list = combined, normalization.method = "SCT", anchor.features = features,verbose = opt$verbose)
        combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT",verbose = opt$verbose) 
        # Choose the features to use when integrating multiple datasets. 
        # This function ranks features by the number of datasets they are deemed variable in, breaking ties by the median variable feature rank across datasets. 
        # It returns the top scoring features by this ranking.
        VariableFeatures(combined.sct)=rownames(combined.sct@assays$integrated@scale.data)
        DefaultAssay(combined.sct) <- "integrated"
    }else{

        merged.rna=merge(spatial.ob[[1]],spatial.ob[2:length(spatial.ob)],add.cell.ids=sample_list[,1])
        #VariableFeatures(merged.rna)=rownames(combined.sct@assays$SCT@data)
        #DefaultAssay(merged.sct) <- "SCT"
        merged.sct <- SCTransform2(merged.rna)
    }

    j=0
    for (i in 1:length(merged.sct@images)){
        j=j+1
    if( nrow(merged.sct@images[[j]]@coordinates)==0){
        merged.sct@images[[j]]=NULL
        j=j-1
    }
    }
    names(merged.sct@images)=as.character(sample_list[,1])
    
    ##基因表达
    spatial.ob2=list()
    for (sample in samples){
        spatial.ob2[[sample]]=subset(merged.sct, subset=orig.ident == eval(sample))
        DefaultAssay(spatial.ob2[[sample]]) <- "SCT"
    }

    spatial.ob3=list()
    for (group in c("Y","M","O")){
        spatial.ob3[[group]]=subset(merged.sct, subset=group == eval(group))
        DefaultAssay(spatial.ob3[[group]]) <- "SCT"
    }
    #combined.sct=merge(spatial.ob2[[1]],spatial.ob2[2:length(spatial.ob2)],add.cell.ids=samples)
    #VariableFeatures(combined.sct)=rownames(combined.sct@assays$SCT@data)
    #DefaultAssay(combined.sct) <- "SCT"
    

    plot3=VlnPlot(merged.sct, features = c("FOXP1"),group.by="group",assay="SCT",slot="data")

    ggsave(plot3, width=10, height=10,dpi=600,units ="in",bg="white", filename=paste0("genes_exp/VlnPlot_",i,".pdf"),limitsize = FALSE)
        
    for (sample in samples){
    
        for (i in c("FOXP1")){
            DefaultAssay(spatial.ob2[[sample]]) <- "SCT"

            plot2=SpatialPlot(
                    object =spatial.ob2[[sample]],
                    slot="scale.data",
                    features = i,
                    alpha = c(1, 1),
                    crop=TRUE,
                    image.alpha=0,
                   # pt.size.factor = 1.2,
                    stroke=0.1,
                    min.cutoff=0.001,
                   # cols=15,
                    images=sample
                    )  +
                theme(plot.title = element_text(color = "black", size = 10,hjust = 0.5), 
                    axis.text.x=element_text(color = "black",size = 10),
                    axis.text.y=element_text(color = "black", size = 10),
                    axis.title.x=element_text(color = "black",size = 12),
                    axis.title.y=element_text(color = "black", size = 12),
                    axis.ticks.length = unit(.1, "cm"),
                    axis.ticks=element_line(color="black"),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    legend.justification = c(0.5,0.5),
                    legend.position="right",
                    legend.text=element_text(face="plain",size=12),
                    legend.title = element_text(face="plain",size=12),
                    legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
                    #plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
                    strip.text = element_text(face = "bold"),
                    panel.background = element_rect(fill = "black")
                )
            # plot3=plot2+ 
            #         scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            #         scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
            if(!file.exists(paste0("genes_exp/",sample))){
                dir.create(paste0("genes_exp/",sample))
            }
            #ggsave(plot2, width=5, height=5,dpi=900,units ="in",bg="white", filename=paste0("genes_exp/SpatialFeaturePlot_",i,".color1.pdf"),limitsize = FALSE)
            ggsave(plot2, width=10, height=10,dpi=600,units ="in",bg="white", filename=paste0("genes_exp/",sample,"/SpatialFeaturePlot_",i,".pdf"),limitsize = FALSE)
        }
    }

}
