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
Idents(single.ob) <- single.ob@meta.data[,"celltype"]
DefaultAssay(single.ob) <- "RNA"

##########全部celltype的点图##########
list_genes =c("CD2","KLRB1","IL7R","GZMA","NKG7","CCL5","IFI30",
              "LYZ","TYROBP","CLDN5","VWF","TM4SF1","RGS5",
              "MUSTN1","ACTA2","COL1A2","STAR","DCN","TUBB8",
              "ZP3","FIGLA","HSD17B1","AMH","GSTA1")

pic<-DotPlot(single.ob,
             features=list_genes,
             cols=c("grey", "red"),
             cluster.idents = T)+ 
  coord_flip()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5))+
  labs(x=NULL,y=NULL) + 
  guides(size=guide_legend(order=1,title="Percent Expressed"),colour=guide_colorbar(order=2,title="Average Expression"))+
  scale_color_gradientn(colours=c('#330066','#336699','#66CC66','#FFCC33'))
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","All_celltype_dotplot.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","All_celltype_dotplot.png",sep="/"),pic)

##########FeaturePlot##########
pic <- FeaturePlot(single.ob,features=c("AMH","ZP3","DCN","ACTA2","GSTA1","TUBB8","STAR","MUSTN1",
                                        "TM4SF1","TYROBP","CCL5","IL7R","VWF","IFI30","NKG7","KLRB1"),raster=FALSE,order=T)+
       theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","Merge_featureplot.pdf",sep="/"),pic,width=10.5,height=6.5)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","Merge_featureplot.png",sep="/"),pic,width=10.5,height=6.5)

##########Theca & stroma的点图##########
single.ob=readRDS("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Theca_A_stroma_subcls.rds")
Idents(single.ob) <- single.ob@meta.data[,"celltype"]
DefaultAssay(single.ob) <- "RNA"
list_genes =c("CD74","HLA-DRB1","SRGN","HLA-DRA","CCL4","FBLN1",
              "NR4A2","CXCL2","IGFBP6","CNN1","ACTG2","TAGLN","ACTA2",
              "CCN2","COL1A2","COL3A1","COL1A1","CYB5A","STAR")
pic<-DotPlot(single.ob,
             features=list_genes,
             cols=c("grey", "blue"),
             cluster.idents = T)+ 
  coord_flip()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5))+
  labs(x="Features",y=NULL) + 
  guides(size=guide_legend(order=2,title="Percent Expressed"),colour=guide_colorbar(order=1,title="Average Expression"))
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","TandS_dotplot.pdf",sep="/"),pic)
ggsave(paste("/PERSONALBIO/work/singlecell/s04/Analysis/Plot/Pictures/Expression","TandS_dotplot.png",sep="/"),pic)