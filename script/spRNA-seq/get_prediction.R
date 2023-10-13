suppressPackageStartupMessages(library(Seurat))

sample_list=read.csv("../02_Cellranger/all.csv",stringsAsFactors =F)


for(sample in sample_list$sample ){
    sc_data=readRDS(paste0("../suerat_Rdata/",sample,"_spatial_cell_locations.mapping.ref.rds"))
    cell_list=c("Endothelial","Theca & stroma","Granulosa","Smooth muscle","NK","monocyte","oocyte","T")
    cell_list_tmp=cell_list[which(cell_list %in% colnames(sc_data@meta.data))]
    predicttion=sc_data@meta.data[,c(cell_list_tmp,"max")]
    predicttion$"predicted.id"=apply(predicttion[,cell_list_tmp],1,function(x){cell_list[which.max(x)]})
    predicttion=predicttion[,c("predicted.id",cell_list_tmp,"max")]
    write.table(predicttion,paste0(sample,"_label_transfer_bc.tsv"),sep="\t")
}