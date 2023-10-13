library(clusterProfiler)
library(org.Hs.eg.db)

diff.cluster <- read.table("/markers.xls",sep = "\t",header=T,row.names=1)
diff.cluster=arrange(diff.cluster, desc(avg_logFC))
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
    write.table(df,file="./GO_enrichment.xls",col.names=T,row.names=F,quote=F,sep='\t')
}