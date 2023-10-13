library(Seurat)
library(tidyverse)
library(clusterProfiler)
sc <- readRDS("sub_Theca_A_stroma.rds")
Idents(sc) <- "celltype"
plotgene <- ""
for (celltype in c("Theca_A_sroma_1","Theca_A_sroma_2","Theca_A_sroma_3","Theca_A_sroma_4","Theca_A_sroma_5")){
gene <-  read.table(paste0("./","cluster_",celltype,"_markers.xls"),sep = "\t",header =T,row.names = 1)
diff.cluster=arrange(gene, desc(avg_logFC))
geneList <- subset(diff.cluster,p_val < 0.05 & avg_logFC > 0.25)
geneList <- geneList$gene
en <- enrichGO(gene       = geneList,
                 OrgDb         = "org.Hs.eg.db",
                 keyType       = 'SYMBOL',
                 ont           = "ALL",
                 minGSSize     = 1 ,
                 maxGSSize     = 100000 ,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 1,
                 qvalueCutoff  = 1)
                 
getBGnumber <- function(ratio,split="/"){
  list<-strsplit(ratio, split = split)[[1]]
  return(as.numeric(list[1]))
}
if(!is.null(en)){
res=en@result
res$Total=apply(res,1,function(x){getBGnumber(x[5])})

df=data.frame( Category=res$ONTOLOGY,
                   GO=res$ID,
                   Term=res$Description,
                   List=res$Count,
                   Total=res$Total,
                   Pvalue=res$pvalue,
                   adjustPvalue=res$p.adjust,
                   Gene=res$geneID
    )
    df<-df[order(df$Pvalue),]
    write.table(df,file=paste0("./",celltype,"GO_enrichment.xls"),col.names=T,row.names=F,quote=F,sep='\t')
}
plotgene <- c(plotgene,gene$gene)
}


p <- DoHeatmap(sc,
          features = plotgene[-1],
          group.by = "subcls",
          assay = 'RNA') + scale_fill_gradientn(colors = c("blue","white","firebrick3"))
ggsave("heatmap.png",p)
ggsave("heatmap.pdf",p)