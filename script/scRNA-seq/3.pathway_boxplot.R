library(dplyr)
library(Seurat)
library(UCell)
library(ggplot2)
library(clusterProfiler)
library(ggpubr)
rds = "./multi.rds"
geneset = "pathway.gmt"
outdir = "./"
plot_idents = "group"
needcelltype <- "ALL"
source("color.R")
clustcol<- colorls[["NPG"]]

sc <- readRDS(rds)
if (needcelltype == "ALL") {
sc <- sc
} else {
sc <- subset(sc,celltype %in% needcelltype)
}

x <- readLines(geneset)
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
gene_gmt <- lapply(res, "[", -c(1:2))

sc <- AddModuleScore(object = sc,
  features = gene_gmt,
  name = names(res)
  )
pathname <- names(res)
colnames(sc@meta.data)[length(colnames(sc@meta.data))] <- pathname
saveRDS(sc,"./AddmoduleScore.rds")
data<- FetchData(sc,vars = c("seurat_clusters","group",pathname))
P<- ggplot(data, aes(x=group,y = data[,pathname])) +
      theme_bw()+
      theme(panel.grid = element_blank(),axis.text.x=element_text(angle=1,hjust = 1,vjust=0.5),plot.title = element_text(hjust = 0.5))+
      labs(x=NULL,y="Score",title = newname)+
      geom_boxplot(outlier.size=1,position=position_dodge(0),aes(color = factor(group)))+
      NoLegend()
print(P)
ggsave(paste0(outdir_this,"/",pathname,"_","Score.box_celltype.png"),height = 8,width = 8)
ggsave(paste0(outdir_this,"/",pathname,"_","Score.box_celltype.pdf"),height = 8,width = 8)
sig <- compare_means(Cellular_senescence_signaling_pathway~group, data=data)
write.table(sig,"sig.xls",sep = "\t",col.names =NA)
