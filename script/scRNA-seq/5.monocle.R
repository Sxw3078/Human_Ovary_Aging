bb <- "sub.rds"
idents = "celltype"
print('monocle2')
library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)
Idents(bb) = idents
DefaultAssay(bb) ='RNA'
data <- as(as.matrix(bb@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = bb@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocds <- newCellDataSet(data,
								  phenoData = pd,
								  featureData = fd,
								  lowerDetectionLimit = 0.5,
								  expressionFamily = negbinomial.size())
print("format data done , filter select genes ")
#pData(monocds)$Cluster<-as.factor(pData(monocds)$celltype) 
pData(monocds)['Cluster']=bb@active.ident	
monocds <- estimateSizeFactors(monocds)
monocds <- estimateDispersions(monocds)
monocds <- detectGenes(monocds, min_expr = 0.1)
print(head(fData(monocds)))
expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 10)) # nolint
monocds <- monocds[expressed_genes, ]
disp_table <- dispersionTable(monocds)
express_genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
dir=paste0('monocle_output')
trajectory_cluster_dir=paste('./',dir,sep="")
monocds <- setOrderingFilter(monocds, express_genes)
monocds <- reduceDimension(
monocds,
max_components = 2,
method = "DDRTree")
monocds <- orderCells(monocds)
saveRDS(monocds,file=paste(trajectory_cluster_dir,"monocle.rds",sep="/"))
A = monocds@featureData@data %>% filter(use_for_ordering == 'TRUE')

plot_cell_trajectory(monocds,color_by="seurat_clusters", size=1,show_backbone=TRUE)+ facet_wrap("~seurat_clusters", nrow = 2)+ scale_color_manual(values = colour)
ggsave(paste0(trajectory_cluster_dir,'/','seurat_clusters.pdf'),width = 12,height = 8) 
ggsave(paste0(trajectory_cluster_dir,'/','seurat_clusters.png'),width = 12,height = 8, dpi = 300,bg = 'white')
p1<-plot_cell_trajectory(monocds, color_by = "Cluster")+ scale_color_manual(values = colour)
ggsave(p1,file=paste(trajectory_cluster_dir,"monocle_Cluster.pdf",sep="/"))
ggsave(p1,file=paste(trajectory_cluster_dir,"monocle_Cluster.png",sep="/"))
p2<-plot_cell_trajectory(monocds, color_by = "State")+ scale_color_manual(values = colour)
ggsave(p2,file=paste(trajectory_cluster_dir,"monocle_State.pdf",sep="/"))
ggsave(p2,file=paste(trajectory_cluster_dir,"monocle_State.png",sep="/"))
p3<-plot_cell_trajectory(monocds, color_by = "Pseudotime")
ggsave(p3,file=paste(trajectory_cluster_dir,"monocle_Pseudotime.pdf",sep="/"))
ggsave(p3,file=paste(trajectory_cluster_dir,"monocle_Pseudotime.png",sep="/"))
p<-plot_cell_trajectory(monocds, color_by = "sample")+ scale_color_manual(values = colour)
ggsave(p,file=paste(trajectory_cluster_dir,"monocle_sample.pdf",sep="/"))
ggsave(p,file=paste(trajectory_cluster_dir,"monocle_sample.png",sep="/"))
pn1<-plot_cell_trajectory(monocds, color_by = "State")+facet_wrap(~State, nrow = 1)+ scale_color_manual(values = colour)
ggsave(pn1,file=paste(trajectory_cluster_dir,"monocle_State_split.pdf",sep="/"))
ggsave(pn1,file=paste(trajectory_cluster_dir,"monocle_State_split.png",sep="/"))
pn2<-plot_cell_trajectory(monocds, color_by = "Cluster")+facet_wrap(~Cluster, nrow = nrow)+ scale_color_manual(values = colour)
ggsave(pn2,file=paste(trajectory_cluster_dir,"monocle_Cluster_split.pdf",sep="/"),width=24,height=8)
ggsave(pn2,file=paste(trajectory_cluster_dir,"monocle_Cluster_split.png",sep="/"),width=24,height=8)
pn4<-plot_cell_trajectory(monocds, color_by = "group")+facet_wrap(~group, nrow = nrow)+ scale_color_manual(values = colour)
ggsave(pn4,file=paste(trajectory_cluster_dir,"monocle_group_split.pdf",sep="/"),width=4)
ggsave(pn4,file=paste(trajectory_cluster_dir,"monocle_group_split.png",sep="/"),width=4)




gbm_cds = readRDS(paste(trajectory_cluster_dir,"monocle.rds",sep="/"))
BEAM_res <- BEAM(gbm_cds, branch_point = 1, cores = 2)  ## branch_point: alternative
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
monocle_outdir = './'
write.table(BEAM_res, file = paste(monocle_outdir,'gene_related_to_branch.txt',sep='/'), row.names = T, quote = F)
heatmap_gene<-row.names(BEAM_res)[order(BEAM_res$qval)][1:50]
pdf(file = paste(monocle_outdir,"branch_dependent_gene_heatmap.pdf",sep='/'), width = 9, height = 5)
plot_genes_branched_heatmap(gbm_cds[heatmap_gene,],
                                          branch_point = 1,
                                          num_clusters = 3,
                                          cores = 2,
                                          use_gene_short_name = T,
                                          show_rownames = T)
dev.off()
branched_genes <- row.names(BEAM_res)[order(BEAM_res$qval)][1:12]
for (i in c(3,6,9,12)){
fit = try(plot_genes_branched_pseudotime(gbm_cds[branched_genes[i-2:i],],
                       branch_point = 1,
                       color_by = "Cluster",
                       ncol = 1) +ggtitle(sample)+
                       theme(plot.title = element_text(hjust = 0.5),legend.position = "right"))