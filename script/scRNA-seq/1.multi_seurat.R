#multi all sample
library(dplyr)
library(Seurat)
library(cowplot)
library(ggplot2)
library(reshape2)
library(DoubletFinder)
library(harmony)
library(dittoSeq)
library(scales)
project_dir=getwd()
mt="^MT-"
datadir=paste0("UMIcounts_HGSOC.tsv")
single.data<- read.table(datadir,sep = '\t',header = T ,row.names = 1 ,check.names = F)
single.ob<-CreateSeuratObject(counts = single.data, project = "single", min.cells = 3, min.features = 200)
anno = read.table('HGSOC.tsv',sep = '\t',row.names = 1,header = T,check.names = F)
single.ob = single.ob[,rownames(anno)]
single.ob$sample=anno$sample
single.ob$treatment_phase = anno$treatment_phase
single.ob$Seurat_clusters = anno$cell_subtype
single.ob$cluster = anno$cell_type 
seurat_exp_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
ifnb.list <- SplitObject(single.ob, split.by = "sample")
#integrat sample object
immune.combined<-merge(ifnb.list[[1]],ifnb.list[2:length(ifnb.list)])%>%NormalizeData() %>% FindVariableFeatures()
immune.combined<-ScaleData(immune.combined,feature=rownames(immune.combined)) %>% RunPCA(verbose = FALSE) %>% RunHarmony( group.by.vars = "sample")
immune.combined<-RunUMAP(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined<-RunTSNE(immune.combined, reduction = "harmony", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "harmony", dims = 1:20) %>% FindClusters()
immune.combined$seurat_clusters <- anno$cell_subtype
Idents(immune.combined) <- anno$cell_subtype
seurat_exp_cluster_dir=paste(project_dir,"04_Cluster",sep="/")
ElbowPlot(immune.combined)
cluster_number<-table(immune.combined$seurat_clusters)%>% reshape2::melt()
colors=hue_pal()(nrow(cluster_number))
cluster_number$Cluster <- as.factor(cluster_number$Var1)
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p2)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(seurat_exp_cluster_dir,"cluster_umap.png",sep="/"),width = 14,height = 7)
#find marker
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- ScaleData(immune.combined, feature=rownames(immune.combined),verbose = FALSE)
saveRDS(immune.combined,"use.rds")
markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(seurat_exp_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

for( clust_num in  unique(Idents(immune.combined))){
        cluster_dir=paste(seurat_exp_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
        cluster_dir_enrich=paste(cluster_dir,"enrichment",sep="/")
        if(!file.exists(cluster_dir_enrich)){dir.create(cluster_dir_enrich)}
        cluster_markers=subset(markers,cluster==clust_num)
        rownames(cluster_markers)<-cluster_markers$gene
        if(nrow(cluster_markers)>1){
        genelist=cluster_markers$gene
        write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
        }
}





#addanno




seurat_anno <- readRDS("use.rds")
seurat_anno <- read.table("anno.txt",sep="\t",header=F,col.names = c("cluster","anno"))
seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
seurat_anno
new.cluster.ids <-as.vector(seurat_anno$anno)
names(new.cluster.ids) <-seurat_anno$cluster

print(seurat_anno$cluster)
print(unique(seurat_obj@active.ident))
seurat_obj <- subset(seurat_obj,idents=as.vector(seurat_anno$cluster))
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltype <- Idents(seurat_obj)
anno_cluster_dir <- "./anno"
p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
plot_grid(p2)
ggsave(paste(anno_cluster_dir,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(anno_cluster_dir,"cluster_umap.png",sep="/"),width = 14,height = 7)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ScaleData(seurat_obj, feature=rownames(seurat_obj),verbose = FALSE)
saveRDS(seurat_obj,"anno.rds")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(anno_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

for( clust_num in unique(Idents(seurat_obj))){
        cluster_dir=paste(anno_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
        cluster_markers=subset(markers,cluster==clust_num)
        rownames(cluster_markers)<-cluster_markers$gene
        if(nrow(cluster_markers)>1){
        genelist=cluster_markers$gene
        write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
        }
}








#subcluster
seurat_obj <- readRDS("anno.rds")
cluster <- seurat_obj$celltype[1]
print("subset data")
seurat_obj <-subset(seurat_obj,celltype %in% cluster)
ifnb.list <- SplitObject(seurat_obj, split.by = "sample")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
sample_list <- unique(seurat_obj@meta.data$sample)
seurat_obj <- ScaleData(seurat_obj,feature=rownames(seurat_obj), verbose = FALSE)
seurat_obj <- FindVariableFeatures(object = seurat_obj,selection.method = 'vst', nfeatures = 2000)
seurat_obj <- RunPCA(seurat_obj,  features = VariableFeatures(object = seurat_obj) ,verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:20)
seurat_obj <- FindClusters(seurat_obj,resolution = resol)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
saveRDS(seurat_obj,"sub.rds")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sub_cluster_dir <- "./sub"
write.table(markers,paste(sub_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

for( clust_num in  unique(Idents(seurat_obj))){
        cluster_dir=paste(sub_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
        cluster_markers=subset(markers,cluster==clust_num)
        rownames(cluster_markers)<-cluster_markers$gene
        if(nrow(cluster_markers)>1){
        genelist=cluster_markers$gene
        write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
        }
}




#subcluster_addanno
seurat_anno <- readRDS("sub.rds")
seurat_anno <- read.table("subanno.txt",sep="\t",header=F,col.names = c("cluster","anno"))
seurat_anno <-seurat_anno %>% as_tibble() %>% 
  separate_rows(cluster, sep = ",")
seurat_anno
new.cluster.ids <-as.vector(seurat_anno$anno)
names(new.cluster.ids) <-seurat_anno$cluster

print(seurat_anno$cluster)
print(unique(seurat_obj@active.ident))
seurat_obj <- subset(seurat_obj,idents=as.vector(seurat_anno$cluster))
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$celltype <- Idents(seurat_obj)
anno_cluster_dir <- "./sub_anno"
p2 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
plot_grid(p2)
ggsave(paste(anno_cluster_dir,"cluster_umap.pdf",sep="/"),width = 14,height = 7)
ggsave(paste(anno_cluster_dir,"cluster_umap.png",sep="/"),width = 14,height = 7)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- ScaleData(seurat_obj, feature=rownames(seurat_obj),verbose = FALSE)
saveRDS(seurat_obj,"sub_anno.rds")
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,paste(anno_cluster_dir,"all_markers.xls",sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)

for( clust_num in unique(Idents(seurat_obj))){
        cluster_dir=paste(anno_cluster_dir,paste("cluster",clust_num,sep="_"),sep="/")
        cluster_markers=subset(markers,cluster==clust_num)
        rownames(cluster_markers)<-cluster_markers$gene
        if(nrow(cluster_markers)>1){
        genelist=cluster_markers$gene
        write.table(cluster_markers,paste(cluster_dir,paste("cluster",clust_num,"markers.xls",sep="_"),sep="/"),sep="\t",quote=F,row.names=T,col.names=NA)
        }
}



