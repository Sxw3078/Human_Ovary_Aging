library(Seurat)
library(ggplot2)
sc <- readRDS("/PERSONALBIO/work/singlecell/s05/Analysis/T202009028/code/test/Theca_A_stroma_subcls.rds")
Idents(sc) <- "group"
DefaultAssay(sc) <- "RNA"
sc_all <- sc
sc <- subset(sc,subcls %in% "Theca & stroma 1")
p <- VlnPlot(sc, 
        features = c("CDKN1A"),
        pt.size = 0) +theme_bw()
        
        


gene_test <- function(SerautObj, 
                           genes.use, 
                           group.by=NULL, 
                           assay = "RNA", 
                           comp = NULL, 
                           alpha_start = .05, 
                           Bonferroni = T,
                           only_postive =F) {
  p_val.out <- c()
  stat.out <- c()
  condition.out <- c()
  gene.out <- c()
  if (only_postive == F){
    for (gene in genes.use){
      group1_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[1],])
      group1_exp = SerautObj@assays[[assay]]@data[gene, group1_cellname] 

      group2_cellname = rownames(SerautObj@meta.data[SerautObj@meta.data[[group.by]] == comp[2],])
      group2_exp = SerautObj@assays[[assay]]@data[gene, group2_cellname]
      t_out = t.test(group1_exp, group2_exp)
      cond = paste(comp[1], comp[2], sep = "_")
      condition.out <- c(condition.out, cond)
      stat.out <- c(stat.out, t_out[["statistic"]])
      p_val.out <- c(p_val.out, t_out[["p.value"]])
      gene.out <- c(gene.out, gene)
    }
  }
    new_alpha = alpha_start/(2*length(genes.use))
    cat(paste("\n", "P-value for significance: p <", new_alpha, "\n"))
    sig_out = p_val.out < new_alpha
    dfOUT<- data.frame(gene=gene.out, condition = condition.out, p_val = p_val.out, statistic = stat.out, significant = sig_out)

    dfOUT$sig = ifelse(dfOUT$p_val > 0.05, "ns",
                       ifelse(dfOUT$p_val > 0.01, '*',
                              ifelse(dfOUT$p_val > 0.001, "**", "****")))

  return(dfOUT)
}
A <- gene_test(sc, 
                    genes.use = c("CDKN1A"),
                    group.by = 'group', 
                    comp = c("Y", "M"))

B <- gene_test(sc, 
                    genes.use = c("CDKN1A"),
                    group.by = 'group', 
                    comp = c("Y", "O"))
C <- gene_test(sc, 
                    genes.use = c("CDKN1A"),
                    group.by = 'group', 
                    comp = c("M", "O"))                  

