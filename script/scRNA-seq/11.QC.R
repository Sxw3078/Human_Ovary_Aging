library(Seurat)
library(grid)
library(ggplot2)
library(paletteer)
library(scales)
library(patchwork)
library(RColorBrewer)
source("/PERSONALBIO/work/singlecell/s00/software/script/1.source/color/color.R")
colors = colorls$"NPG"

single.ob=readRDS("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/ruanchao.wumeng.YMO.rds")

##########UMI和Genes的分布图##########
data=single.ob@meta.data
note=summary(data)
note_2=strsplit(note[4,3],":",fixed=T)
note_3=strsplit(note[4,4],":",fixed=T)
data$sample=gsub("Y18_HW","Y18",data$sample)
data$sample=gsub("Y22_LXJ","Y22",data$sample)
data$sample=gsub("Y28_XZR","Y28",data$sample)
data$sample=gsub("Y38_LHM","Y37",data$sample)
data$sample=gsub("Y38_LWY","Y38",data$sample)
data$sample=gsub("Y40_LRX","Y39",data$sample)
data$sample=gsub("Y48_FXY","Y47",data$sample)
data$sample=gsub("Y48_JJY","Y48",data$sample)
data$sample=gsub("Y48_WWH","Y49",data$sample)

data$nCount_RNA=data$nCount_RNA/1000
pic<-ggplot()+
     geom_histogram(data=data,aes(nCount_RNA,..count..),color="#0099B3",fill="#0099B3",bins=52)+
     geom_density(data=data,aes(nCount_RNA,0.96*..count..,color="red"),size=1)+
     theme_bw()+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border=element_blank())+
     theme(axis.line=element_line(colour="black"))+
     labs(x="the number of UMIs detected in each cell(10^3)",y="cell count")+
     annotate("text",x=30,y=6000,label=paste0("Mean=",as.numeric(note_2[[1]][2])),col="black")+
     scale_x_continuous(limits = c(0,50))+
     theme(legend.title=element_blank())+theme(legend.position = 'none')
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","the_number_of_UMIs_detected_in_each_cell.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","the_number_of_UMIs_detected_in_each_cell.png",sep="/"),pic)

data$nFeature_RNA=data$nFeature_RNA/1000
pic<-ggplot()+
     geom_histogram(data=data,aes(nFeature_RNA,..count..),color="#0099B3",fill="#0099B3",bins=100)+
     geom_density(data=data,aes(nFeature_RNA,0.05*..count..,color="red"),size=1)+
     theme_bw()+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border=element_blank())+
     theme(axis.line=element_line(colour="black"))+
     labs(x="the number of genes detected in each cell(10^3)",y="cell count")+
     annotate("text",x=3.5,y=1800,label=paste0("Mean=",as.numeric(note_3[[1]][2])),col="black")+
     scale_x_continuous(limits = c(0,5))+
     theme(legend.title=element_blank())+theme(legend.position = 'none')
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","the_number_of_genes_detected_in_each_cell.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","the_number_of_genes_detected_in_each_cell.png",sep="/"),pic)

##########UMI和Genes的分组分布图##########
data=single.ob@meta.data
data$sample=gsub("Y18_HW","Y18",data$sample)
data$sample=gsub("Y22_LXJ","Y22",data$sample)
data$sample=gsub("Y28_XZR","Y28",data$sample)
data$sample=gsub("Y38_LHM","Y37",data$sample)
data$sample=gsub("Y38_LWY","Y38",data$sample)
data$sample=gsub("Y40_LRX","Y39",data$sample)
data$sample=gsub("Y48_FXY","Y47",data$sample)
data$sample=gsub("Y48_JJY","Y48",data$sample)
data$sample=gsub("Y48_WWH","Y49",data$sample)
data=within(data, Sample<-factor(Sample,levels=c("Y","M","O")))
with(data,levels(Sample))

pic<-ggplot(data=data,aes(nCount_RNA))+
     geom_density(fill="#5FC6FF",color="#5FC6FF")+
     theme_bw()+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border=element_blank())+
     theme(axis.line=element_line(colour="black"))+
     labs(x="UMI counts",y="cells")+
     scale_x_continuous(limits = c(0,30000))+
     facet_wrap(~ Sample)+
     theme(strip.text=element_text(colour='black',size=rel(1.2)),strip.background=element_rect(fill='white',colour='white',size=rel(0)))+
     theme(legend.title=element_blank())+theme(legend.position = 'none')
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","UMI_counts.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","UMI_counts.png",sep="/"),pic)

pic<-ggplot(data=data,aes(nFeature_RNA))+
     geom_density(fill="#5FC6FF",color="#5FC6FF")+
     theme_bw()+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(panel.border=element_blank())+
     theme(axis.line=element_line(colour="black"))+
     labs(x="Expressed Genes",y="cells")+
     scale_x_continuous(limits = c(0,5000))+
     facet_wrap(~ Sample)+
     theme(strip.text=element_text(colour='black',size=rel(1.2)),strip.background=element_rect(fill='white',colour='white',size=rel(0)))+
     theme(legend.title=element_blank())+theme(legend.position = 'none')
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","Expressed_Genes.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","Expressed_Genes.png",sep="/"),pic)

##########Percentage of mitochondrial genes##########
data=single.ob@meta.data
data$sample=gsub("Y18_HW","Y18",data$sample)
data$sample=gsub("Y22_LXJ","Y22",data$sample)
data$sample=gsub("Y28_XZR","Y28",data$sample)
data$sample=gsub("Y38_LHM","Y37",data$sample)
data$sample=gsub("Y38_LWY","Y38",data$sample)
data$sample=gsub("Y40_LRX","Y39",data$sample)
data$sample=gsub("Y48_FXY","Y47",data$sample)
data$sample=gsub("Y48_JJY","Y48",data$sample)
data$sample=gsub("Y48_WWH","Y49",data$sample)
pic<-ggplot(data=data,aes(x=sample,y=percent.mt,col=sample))+
     geom_violin(trim = FALSE)+
     geom_boxplot(width = 0.2)+
     labs(x="Sample",y="Percentage of mitochondrial genes")+
     theme_bw()+
     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
     theme(legend.title=element_blank())+theme(legend.position = 'none')+
     theme(panel.border=element_blank())+
     theme(axis.line=element_line(size=0.5,colour="black"))
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","Percentage_of_mitochondrial_genes.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","Percentage_of_mitochondrial_genes.png",sep="/"),pic)

##########UMI和Genes的UMAP图##########
pic <- FeaturePlot(single.ob,features=c("nCount_RNA"),raster=FALSE,order=T)+
       theme(legend.position = c(0.8,0.3))+
       ggtitle("nUMI")+
       theme(plot.title = element_text(hjust = 0.5))
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","nUMI.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","nUMI.png",sep="/"),pic)

pic <- FeaturePlot(single.ob,features=c("nFeature_RNA"),raster=FALSE,order=T)+
       theme(legend.position = c(0.8,0.3))+
       ggtitle("nGene")+
       theme(plot.title = element_text(hjust = 0.5))
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","nGene.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/QC","nGene.png",sep="/"),pic)
