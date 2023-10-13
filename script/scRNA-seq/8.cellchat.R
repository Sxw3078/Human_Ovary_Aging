####zhaoyunfei
#####https://www.jianshu.com/p/f7e56cba6121

suppressMessages({
library(CellChat)
library(ggalluvial)
library(ggplot2)
library(patchwork)
library(ggalluvial)
options(stringsAsFactors = FALSE)
library(ComplexHeatmap)
library(igraph)
library(Seurat)
library(reticulate)
library(cowplot)
})
use_python("/PERSONALBIO/work/singlecell/s04/conda/envs/scRNA/bin/python")


rds = "/PERSONALBIO/work/singlecell/s02/Analysis/scirpt-prepare/Y.rds"
cluster = NULL
outdir = "/PERSONALBIO/work/singlecell/s02/Analysis/scirpt-prepare/Y"
number_Patterns = 5
sample = "Y"
species = "Human"
Search = "Secreted Signaling"
pvalue = 0.05
raw = T
#source('/TJPROJ6/SC/personal_dir/zhaoyunfei/CellChat/function.R')

Seurat.obj = readRDS(rds)
fit = try(Seurat.obj@assays$RNA@data)  ####输入均一化的矩阵
if('try-error' %in% class(fit)){
data.input = Seurat.obj@assays$Spatial@data}else{
data.input = Seurat.obj@assays$RNA@data}

if (!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}

if (!is.null(cluster)){
anno <- read.csv(cluster,header=T,check.names=F)
if (length(anno$Cluster) != length(colnames(Seurat.obj))){print("warning ~~~ ,The length of celltype file is not match the sc data,subset will be acrry out ~~~")}
Seurat.obj = Seurat.obj[,anno$Barcode]
data.input = data.input[,anno$Barcode]
identity = data.frame(group = anno$Cluster, row.names = anno$Barcode)
}else{
identity = data.frame(group = Seurat.obj$celltype, row.names = names(Seurat.obj$celltype))
}

cellchat <- createCellChat(object = data.input,meta = identity,group.by = "group")
cellchat <- addMeta(cellchat, meta = identity, meta.name = "Cluster")
cellchat <- setIdent(cellchat, ident.use = "Cluster")
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

if (species == 'Human'){CellChatDB <- CellChatDB.human}else if (species == 'Mouse'){CellChatDB <- CellChatDB.mouse}else{CellChatDB <- CellChatDB.zebrafish}
pdf(paste(outdir,paste('DatabaseCategory',species,'pdf',sep = '.'),sep = '/'),width = 12,height = 7)
showDatabaseCategory(CellChatDB) ##Show the structure of the database
dev.off()
####dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = Search) # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object ,cellchat@DB包含了配受体库的信息
###Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # Subset the expression data of signaling genes in CellChatDB,截取后的矩阵存放在cellchat@data.signaling
#future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat,thresh.p = pvalue)###挑选高变基因,从截取后的矩阵进行计算，存放在cellchat@var.features
cellchat <- identifyOverExpressedInteractions(cellchat)  ####挑选特异的配受体对，存放在cellchat@LR

##A diffusion process is used to smooth genes’ expression values based on their neighbors’ defined in a high-confidence experimentally validated protein-protein network.
if (species == 'Human'){cellchat <- projectData(cellchat, PPI.human)}else if (species == 'Mouse'){cellchat <- projectData(cellchat, PPI.mouse)}else{cellchat <- projectData(cellchat, PPI.human)}
###Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat,raw.use = FALSE,trim = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10) 
###Extract the inferred cellular communication network as a data frame
#输出受配体对
df.net <- subsetCommunication(cellchat,slot.name = "net")
write.table(df.net,paste0(outdir,sprintf("/%s.Communication.lr.xls",sample)),quote=F,row.names=F,col.names=T,sep='\t')

###Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat,thresh = pvalue)  ###存放在cellchat@netP
###Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat , thresh = pvalue)  ##3Calculate the aggregated network by counting the number of links or summarizing the communication probability
#####Visualization and systems analysis of cell-cell communication network
if (!dir.exists(paste(outdir,'circle',sep = '/'))){dir.create(paste(outdir,'circle',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'circle',sep = '/'))

groupSize <- as.numeric(table(cellchat@idents))
pdf(paste(sample,'singaling.network.circle.count.pdf',sep='.'),width = 9,height = 9)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",arrow.size = 0.5,arrow.width = 2)
dev.off()
pdf(paste(sample,'singaling.network.circle.weights.pdf',sep='.'),width = 9,height = 9)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

#mat <- cellchat@net$weight
#for (i in 1:nrow(mat)) {
#  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#  mat2[i, ] <- mat[i, ]
#  pdf(paste0(paste(outdir,'circle/',sep = '/'),sample,".communcation.network.",rownames(mat)[i],".pdf"))
#  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = 25, title.name = rownames(mat)[i])
#  dev.off()
#}

if (!dir.exists(paste(outdir,'pathway',sep = '/'))){dir.create(paste(outdir,'pathway',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'pathway',sep = '/'))

n <- round(length(levels(cellchat@idents))/2)
vertex.receiver = seq(1,n)

if (length(cellchat@netP$pathways) < 15){num = length(cellchat@netP$pathways)}else(num = 15)

for (way in cellchat@netP$pathways[1:num]){
pdf(paste(sample,way,'singaling.pathway.network.pdf',sep='.'),width = 9,height = 9)
netVisual_aggregate(cellchat, signaling = way,  vertex.receiver = vertex.receiver, vertex.size = groupSize,thresh = 0.05)
dev.off()
pdf(paste(sample,way,'singaling.pathway.circle.pdf',sep='.'),width = 9,height = 9)
netVisual_aggregate(cellchat, signaling = way, layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7,thresh = pvalue)
dev.off()
pdf(paste(sample,way,'singaling.pathway.chord.pdf',sep='.'),width = 9,height = 9)
netVisual_aggregate(cellchat, signaling = way, layout = "chord", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7,thresh = pvalue)
dev.off()
#Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
fit<-try(netAnalysis_contribution(cellchat, signaling = way,thresh = pvalue),silent=TRUE)
if(!'try-error' %in% class(fit)){
pdf(paste(sample,way,'ligand.receptor.contribution.pdf',sep='.'),width = 8,height = 8)
print(netAnalysis_contribution(cellchat, signaling = way,thresh = pvalue))
dev.off()
}

fit<-try(netVisual_heatmap(cellchat, signaling = way, color.heatmap = "Reds"),silent=TRUE)
if(!'try-error' %in% class(fit)){
pdf(paste(sample,way,'ligand.receptor.heatmap.pdf',sep='.'),width = 10,height = 8)
print(netVisual_heatmap(cellchat, signaling = way, color.heatmap = "Reds"))
dev.off()
}
}


###a single ligand-receptor pair
if (!dir.exists(paste(outdir,'lr_pair',sep = '/'))){dir.create(paste(outdir,'lr_pair',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'lr_pair',sep = '/'))
for (way in cellchat@netP$pathways[1:num]){
print("hahah")
print(way)
pairLR <- extractEnrichedLR(cellchat, signaling = way, geneLR.return = FALSE,thresh = pvalue)

if (dim(pairLR)[1] > 5){pairLR = pairLR[1:5,]}

dir.create(way)
setwd(paste0('./',way))
print("hahahah")
print(pairLR)
for (lr in pairLR){
pdf(paste(sample,lr,'singaling.pathway.network.pdf',sep='.'),width = 10,height = 9)
netVisual_individual(cellchat, signaling = way,  pairLR.use = lr, vertex.receiver = vertex.receiver,thresh = pvalue)
dev.off()

pdf(paste(sample,lr,'singaling.pathway.circle.pdf',sep='.'),width = 9,height = 9)
netVisual_individual(cellchat, signaling = way, pairLR.use = lr, layout = "circle",thresh = pvalue)
dev.off()

pdf(paste(sample,lr,'singaling.pathway.chord.pdf',sep='.'),width = 9,height = 9)
p = netVisual_individual(cellchat, signaling = way, pairLR.use = lr, layout = "chord",thresh = pvalue)
print(p)
dev.off()
}
setwd('../')
}

if (!dir.exists(paste(outdir,'cluster',sep = '/'))){dir.create(paste(outdir,'cluster',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'cluster',sep = '/'))
for (cell in levels(cellchat@idents)){
dir.create(paste(outdir,'cluster',paste0('cluster_',cell),sep = '/'),recursive = TRUE)
setwd(paste(outdir,'cluster',paste0('cluster_',cell),sep = '/'))
pdf(paste(sample,'cluster',cell,'sender.lr.bubble.pdf',sep='.'),width = 9,height = 16)
print(netVisual_bubble(cellchat, sources.use = cell, targets.use = levels(cellchat@idents),thresh = pvalue))
dev.off()

for (way in cellchat@netP$pathways[1:num]){
pdf(paste(sample,'cluster',cell,way,'sender.lr.bubble.pdf',sep='.'),width = 9,height = 6)
print(netVisual_bubble(cellchat, sources.use = cell, targets.use = levels(cellchat@idents), signaling = way,thresh = pvalue))
dev.off()
}

pdf(paste(sample,'cluster',cell,'receiver.lr.bubble.pdf',sep='.'),width = 8,height = 10)
print(netVisual_bubble(cellchat, sources.use = levels(cellchat@idents), targets.use = cell,thresh = pvalue))
dev.off()

for (way in cellchat@netP$pathways){
pdf(paste(sample,'cluster',cell,way,'receiver.lr.bubble.pdf',sep='.'),width = 9,height = 7)
print(netVisual_bubble(cellchat, sources.use = levels(cellchat@idents), targets.use = cell, signaling = way,thresh = pvalue))
dev.off()

}
}

if (!dir.exists(paste(outdir,'dot_violin',sep = '/'))){dir.create(paste(outdir,'dot_violin',sep = '/'),recursive = TRUE)}

setwd(paste(outdir,'dot_violin',sep = '/'))


for (way in cellchat@netP$pathways[1:num]){

pdf(paste(sample,way,'violin.pdf',sep = '.'))

print(plotGeneExpression(cellchat, signaling = way,type = 'violin'))

dev.off()

fit = try(plotGeneExpression(cellchat, signaling = way,type = 'dot'))

if (! 'try-error' %in% class(fit)){

pdf(paste(sample,way,'dot.pdf',sep = '.'),width = 8,height = 6)

print(plotGeneExpression(cellchat, signaling = way,type = 'dot'))

dev.off()
}}

###Compute and visualize the network centrality scores
if (!dir.exists(paste(outdir,'centrality_scores_pathway',sep = '/'))){dir.create(paste(outdir,'centrality_scores',sep = '/'),recursive = TRUE)}

setwd(paste(outdir,'centrality_scores',sep = '/'))

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

for (way in cellchat@netP$pathways[1:num]){
pdf(paste(sample,way,'heatmap.pdf',sep='.'), width = 9, height = 6)
netAnalysis_signalingRole_network(cellchat, signaling = way)
dev.off()

pdf(paste(sample,way,'strength.dot.pdf',sep='.'), width = 9, height = 7)
print(netAnalysis_signalingRole_scatter(cellchat, signaling = way))
dev.off()

}

pdf(paste(sample,'heatmap.pathway.outgoing.pdf',sep='.'), width = 9, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
dev.off()

pdf(paste(sample,'heatmap.pathway.incoming.pdf',sep='.'), width = 9, height = 12)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
dev.off()

saveRDS(cellchat, file=paste0(outdir,'/',sample,'.cellchat.rds'))

if (!dir.exists(paste(outdir,'Pattern',sep = '/'))){dir.create(paste(outdir,'Pattern',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'Pattern',sep = '/'))

cellchat <- identifyCommunicationPatterns(cellchat, pattern = 'outgoing', k = number_Patterns)
###Identify global communication patterns and major signals for specific cell groups
pdf(paste(sample,'outgoing','singaling.pattern.pdf',sep='.'),width=12,height = 9)
identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = number_Patterns)
dev.off()

pdf(paste(sample,'outgoing','singaling.pathway.river.pdf',sep='.'),width=12,height = 9)
netAnalysis_river(cellchat, pattern = 'outgoing')
dev.off()

pdf(paste(sample,'outgoing','singaling.pathway.dot.pdf',sep='.'),width=12,height = 9)
netAnalysis_dot(cellchat, pattern = 'outgoing')
dev.off()

cellchat <- identifyCommunicationPatterns(cellchat, pattern = 'incoming', k = number_Patterns)

pdf(paste(sample,'incoming','singaling.pattern.pdf',sep='.'),width=12,height = 9)
identifyCommunicationPatterns(cellchat, pattern = 'incoming', k = number_Patterns)
dev.off()

pdf(paste(sample,'incoming','singaling.pathway.river.pdf',sep='.'),width=12,height = 9)
netAnalysis_river(cellchat, pattern = 'incoming')
dev.off()
pdf(paste(sample,'incoming','singaling.pathway.dot.pdf',sep='.'),width=9,height = 8)
netAnalysis_dot(cellchat, pattern = 'incoming')
dev.off()


###Identify signaling groups based on their functional similarity

if (!dir.exists(paste(outdir,'similarity',sep = '/'))){dir.create(paste(outdir,'similarity',sep = '/'),recursive = TRUE)}
setwd(paste(outdir,'similarity',sep = '/'))

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Visualization in 2D-space
pdf(paste(sample,'pathway.group.functional.similarity.pdf',sep='.'),width=9,height = 7)
netVisual_embedding(cellchat, type = "functional",pathway.remove.show = F, label.size = 3.5)
dev.off()

pdf(paste(sample,'pathway.group.functional.similarity.grid.pdf',sep='.'),width=9,height = 7)
netVisual_embeddingZoomIn(cellchat, type = "functional",nCol = 2)
dev.off()

cellchat <- computeNetSimilarity(cellchat, type = "structural", thresh = 0.25)
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
pdf(paste(sample,'pathway.group.structural.similarity.grid.pdf',sep='.'),width=9,height = 7)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()

pdf(paste(sample,'pathway.group.structural.similarity.grid.pdf',sep='.'),width=9,height = 7)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
dev.off()





