################################输入参数
################################

suppressPackageStartupMessages(library(optparse))
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-a", "--aggr"), help="cellranger aggr output directory, default 02_Cellranger/aggr"),
	      make_option(c("-s", "--aggr_csv"), help="csv file list sample information,also used when run cellrange aggr, default 02_Cellranger/all.csv"),
        make_option(c("-c", "--contrast"), help="compare group used when different expresion between sample/group, default pairwise "),
        make_option(c("-b", "--RMbE"), help="remove batch effect by function IntegrateData , default %default",default=TRUE),
        make_option(c("-m", "--mt"), help="mt gene pattern , default %default",default="^MT"),
        make_option(c("-r", "--mt_regress"), help="regree the mt gene in the SCTransform, default %default",default=FALSE),
        make_option(c("-f", "--mt_cutoff"),type="integer", help="mt percent filter cutfoff , default %default",default="10"),
       # make_option(c("-l", "--resolution"),type="double", help="resolution value for FindClusters, default %default",default="1.0"),
        make_option(c("-p", "--npc"),type="integer", help="number of dim used in pca  , default %default",default=30),
        make_option(c("-d", "--spatial_data"), help="rds for spatial data with celltypes, default %default",default="suerat_Rdata/M1_1_spatial_cell_locations.mapping.ref.rds"),
        make_option(c("-z", "--singcell"), help="rds for singcell data, default %default",default="../2021_03_19_added2/suerat_Rdata/sc_combined.subcluster2.rds"),  
        make_option(c("-t", "--celltype"), help="which celltype or celltypes to be subclustered, default %default",default="Granulosa"),     
        make_option(c("-x", "--celltype2"), help="which celltype or celltypes to be drawed, default %default",default="oocyte"),           
        make_option(c("--colors"), help="color set for scatterpie, choose from NPG,NEJM,JAMA,Lancet,AAAS,report,DEF, default %default",default="Lancet"),        
        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default=FALSE)

)
opt_parser=OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

################################
################################

if(! is.null(opt$workdir)){
  print(paste("change R work directory to ",opt$workdir,sep=""))
  setwd(opt$workdir)
 }
project_dir=getwd()

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")
if(!file.exists(rds_dir)){
    dir.create(rds_dir,recursive=TRUE)
}
#read sample information matrix
#if more than one sample, should give the list csv used in cellranger aggr
cellranger_csv=paste(project_dir,"02_Cellranger","all.csv",sep="/")
if(! is.null(opt$aggr_csv)){
  print(paste("set aggr_csv to ",opt$aggr_csv,sep=""))
  cellranger_csv=opt$aggr_csv
 }

if(exists("cellranger_csv")){
    sample_list=read.csv(cellranger_csv,stringsAsFactors =F)
}else{
    stop("no cellranger_csv file was found\n")
}

if(nrow(sample_list) < 2){
  stop("this pipline need at least 2 sample\n")
}


#import DE contrast information
if(! is.null(opt$contrast)){
  print(paste("get compare contrast file ",opt$contrast,sep=""))
  sample_contrast=read.table(opt$contrast)
 }else{
#sample_contrast=switch(EXPR=as.character(nrow(sample_list)),"2"=combn(sample_list[,1],2), t(combn(sample_list[,1],2)))
 sample_contrast=t(combn(sample_list[,1],2))
}

################################载入R包
################################

suppressPackageStartupMessages(library(SPOTlight))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressWarnings(suppressPackageStartupMessages(library(SeuratData)))
suppressMessages(suppressPackageStartupMessages(library(biomaRt)))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(igraph))
options(future.globals.maxSize = 10000 * 1024^2)

###############################载入自定义功能
###############################开始

source("scatterpie_plot.r")
source("spatial_scatterpie.r")


if(opt$mt_regress){
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208,vars.to.regress ="percent.mt", verbose = opt$verbose)
    }
}else{
    SCTransform2=function(x,y){
        SCTransform(x,method= "glmGamPoi", assay = "Spatial", return.only.var.genes = FALSE, seed.use=19900208, verbose = opt$verbose)
    }
}

# 调色板
blue_green_red <- scale_color_gradient2(low = "blue", mid = "green", high = "red")
byr <- CustomPalette(low = '#0000FF', mid = "#FFFF00",high = "#FF0000", k = 7000)
bgr <- CustomPalette(low = '#0000FF', mid = "#00FF00",high = "#FF0000", k = 7000)
#br
#rblack
#whitered
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
color.ls <- list(
    'NPG' = c(
        "#E54B34", "#4CBAD4", "#009F86", "#3B5387", "#F29A7F", "#8491B3", "#91D1C1",
		"#DC0000", "#7E6047", "#CCCCCC", "#BC8B83", "#33ADAD", "#347988", "#9F7685",
		"#C1969A", "#8BB0BB", "#CE8662", "#B04929", "#A59487", "#E3907E", "#D46F5B",
		"#41B4C1", "#278C87", "#726486", "#DA988C", "#88A0B7", "#B9AC91", "#C63517",
		"#927A66", "#DBAEA4", "#97A4AB", "#21A69A", "#3A6688", "#C98882", "#A593A7",
		"#8EC0BE", "#D85935", "#985738", "#B9AFA9", "#E67059", "#E5BFB9", "#B2CED4",
		"#779F99", "#747A87", "#F2DCD5", "#A7ABB3", "#C1D1CD", "#DCA5A5", "#7E7770",
		"#CCCCCC"
    ),
    'NEJM' = c(
		"#BB3B28", "#0072B4", "#E08626", "#1F854D", "#7876B1", "#6E99AC", "#DC91FF",
		"#ED4C97", "#905B6E", "#A07C74", "#928A3C", "#587F7F", "#7487AF", "#A897D5",
		"#E970C9", "#D6435F", "#A86463", "#7C7E95", "#BA9554", "#52826E", "#858BB0",
		"#99A2C1", "#E39AE4", "#E26E95", "#747291", "#C2916E", "#6E8856", "#758298",
		"#8198AE", "#CCAAEA", "#EC82BF", "#C96265", "#BB7B72", "#5A93B4", "#E0B383",
		"#528569", "#9494B1", "#8DA3AC", "#EEC8FF", "#ED9DC2", "#90767F", "#A08E8A",
		"#928E67", "#6C7F7F", "#929BAF", "#BFB6D5", "#E9ACD9", "#D68D9B", "#A88686",
		"#898A95"
    ),
    'JAMA' = c(
		"#374D54", "#DF8E44", "#00A0D4", "#B24645", "#79AE97", "#6A6599", "#7F796B",
		"#8E6E4F", "#A7988E", "#91778A", "#A07F6C", "#748998", "#776E82", "#5C6360",
		"#655E52", "#C6946A", "#6F8CAF", "#AB6558", "#779C98", "#71698D", "#6E6D65",
		"#B67E4A", "#7B9DB1", "#A56167", "#919781", "#6F7798", "#7C7376", "#49585A",
		"#465154", "#DFB792", "#6ABAD4", "#B27C7B", "#93AEA3", "#827F99", "#7F7C75",
		"#8E7E6F", "#A7A09B", "#91848E", "#A09086", "#869198", "#7D7882", "#606362",
		"#65625C", "#C6AD98", "#8F9EAF", "#AB8882", "#8A9C9A", "#7F7B8D", "#6E6E6A",
		"#B69A80"
    ),
    'Lancet' = c(
		"#00468B", "#EC0000", "#41B43F", "#0099B3", "#66426F", "#FDAE91", "#AC002A",
		"#ACB6B6", "#9E3B4F", "#B2811E", "#40A67D", "#6B7EA9", "#C98599", "#D7675A",
		"#B86E6B", "#6A7BA1", "#6E4D6D", "#D1793E", "#5EAD73", "#6694AE", "#AE80A1",
		"#EAA392", "#B46264", "#949CAB", "#C65356", "#8D9D4B", "#4D9F9B", "#8A7CA4",
		"#E3ABA9", "#C26161", "#B69C9A", "#596D96", "#46688B", "#EC7676", "#7BB47A",
		"#5AA6B3", "#997F9F", "#FDD6C7", "#AC566B", "#B1B6B6", "#9E6D76", "#B29A68",
		"#73A692", "#8A93A9", "#C9A7B1", "#D79F99", "#B89392", "#868EA1", "#6E5E6E",
		"#D1A588"
    )[-8],
    'AAAS' = c(
		"#3A4892", "#ED0000", "#008B45", "#631879", "#00817F", "#BA0020", "#5F549A",
		"#A10055", "#7F807F", "#A93B52", "#9F6D25", "#545A62", "#50557D", "#8A5C4E",
		"#9A3C5D", "#883A77", "#96546A", "#626489", "#7E4372", "#C85113", "#407354",
		"#5C3B7B", "#657167", "#AC293F", "#774989", "#9D3860", "#727285", "#CC2C31",
		"#6E7E35", "#5E3F6D", "#3C6C7E", "#A54137", "#83497B", "#962866", "#8D6B75",
		"#51568D", "#666D92", "#ED7777", "#468B68", "#6E4979", "#418180", "#BA5D6D",
		"#7D779A", "#A1517B", "#808080", "#A9727D", "#9F8662", "#5B5E62", "#67697D",
		"#8A736C"
    ),
    'report' = c(
		"#BF5A17", "#F0017F", "#386CB0", "#FDBF85", "#BEADD3", "#7FC97F", "#FA7F72",
		"#666666", "#B4B1B1", "#D8434F", "#A95597", "#AA949D", "#E1B6AD", "#A2BCAA",
		"#C7A978", "#B1746C", "#8C8A8A", "#C28665", "#CC5036", "#CE3F8A", "#7B7FA7",
		"#EFBB9A", "#B1B5BF", "#A7BA7B", "#D67A6F", "#787878", "#BE9B8A", "#E42F67",
		"#7D63A3", "#D4A992", "#D0B2C1", "#92C394", "#E29675", "#8C6D69", "#A09D9D",
		"#C27140", "#BF8C6B", "#F079B7", "#748EB0", "#FDDEC1", "#C9C0D3", "#A4C9A4",
		"#FABCB6", "#666666", "#B4B3B3", "#D88D93", "#A97FA0", "#AA9FA3", "#E1CCC7",
		"#AFBCB3"
    ),
    'DEF' = c(
		"#FF6600", "#FFFF66", "#009966", "#FF6666", "#666600", "#CCFFCC", "#669933",
		"#339966", "#FFB637", "#9ACB68", "#A78A65", "#B26B39", "#9BAF69", "#99CA7E",
		"#52984E", "#B0893E", "#FFAC57", "#D3E587", "#7A9371", "#D88671", "#83894D",
		"#BEE4B4", "#6C9857", "#849262", "#FFE47A", "#77B27A", "#D49281", "#8C723C",
		"#BDD6A8", "#8BB16F", "#5A986A", "#D99353", "#FFB380", "#FFFFB3", "#4D9980",
		"#FFB3B3", "#666633", "#E6FFE6", "#809966", "#669980", "#FFDB9B", "#B3CB9A",
		"#A79986", "#B28F76", "#A5AF8C", "#B2CAA4", "#759873", "#B09C77", "#FFD6AB",
		"#DCE5B6"
    )
)


color.continuous.ls <- list(
    'BrBG' = c('#8C510A', '#F5F1E7', '#01665E'),
    'PiYG' = c('#C51B7D', '#F8F0F4', '#4D9221'),
    'RdBu' = c('#B2182B', '#F8EFE9', '#2166AC'),
    'bl2rd' = c('#0000FF', '#0CE2F2', '#FF0000'),
    'gnbu' = c('#084081', '#6AC1C8', '#F7FCF0'),
    'matlablike' = c('#0000AA', '#A1FFDB', '#AA0000'),
    'Gwr' = c('#0000FF', '#ffffff', '#FF0000'),
    'Bwr' = c('#3399CC', '#ffffff', '#FF6666')
)

#getPalette(18)

###############################载入自定义功能
###############################结束

###############################载入单细胞RNA数据
###############################

sc_data_raw=readRDS(opt$singcell)

if(length(unique(sc_data_raw@meta.data$sample))>1){
    sc_data.list <- SplitObject(sc_data_raw, split.by = "sample")
}else{
    sc_data.list=sc_data_raw
}

sc_sample_list=unique(sc_data_raw@meta.data$sample)



cols2=color.continuous.ls[['RdBu']]

###############################载入空间RNA数据
###############################

spatial.ct=readRDS(opt$spatial_data)
if(is.null(spatial.ct)){
    {writeLines("spatial.ct was not exist in the .rds file specified by the \"-d\" or \"--spatial_data\" parameter, please check it manually."); q("no")}
}
sampleID=spatial.ct$sample[1]

sampleID2=sample_list[which(sample_list$sample==sampleID),5]
if(sampleID2 %in% sc_sample_list){
  sc_data=sc_data.list[[sampleID2]]
}else{
  #sc_data=merge(sc_data.list[[opt$samples]])
  sc_data=sc_data_raw
}

rm(sc_data.list)
###############################循环进行spotlight
###############################开始

spotlight_dir=paste(project_dir,"07_mapping_ref_sub",sep="/")
if(!file.exists(spotlight_dir)){
    dir.create(spotlight_dir,recursive=TRUE)
}

sc4subcluster=cluster_markers_all=spotlight=decon_mtrx=nmf_mod=cell_types_all=decon_mtrx_sub=decon_mtrx_new=list()

print(paste0("Now sample ",sampleID," is processing!\n"))

spotlight_dir_sample=paste(spotlight_dir,sampleID,sep="/")
if(!file.exists(spotlight_dir_sample)){
    dir.create(spotlight_dir_sample,recursive=TRUE)
}

# if(file.exists(paste(rds_dir,paste0("spotlight.",sampleID,".rds"),sep="/"))){
   
#   spotlight[[sampleID]]=readRDS(paste(rds_dir,paste0("spotlight.",sampleID,".rds"),sep="/"))
 
# }else{



cells=unique(sc_data@meta.data$celltype)
cells_other=cells[cells!=opt$celltype]
sc_tmp1 = subset(sc_data, subset=celltype %in% cells_other)
sc_tmp1$new.ident=sc_tmp1$celltype
sc_tmp2 = subset(sc_data, subset=celltype == opt$celltype)
sc_tmp2$new.ident=sc_tmp2$celltype3

sc4subcluster = merge(sc_tmp1,sc_tmp2)
mapper=unique(sc4subcluster$new.ident)
#sc4subcluster$new.ident=paste0("celltype.",as.numeric(factor(sc4subcluster$new.ident,levels=mapper)))

Idents(sc4subcluster)="new.ident"

sc4subcluster=SCTransform(sc4subcluster, verbose = opt$verbose) %>%
RunPCA(., verbose = opt$verbose) %>%
RunUMAP(., dims = 1:30, verbose = opt$verbose)   # 在前面一起标准化

spatial.ct@meta.data$contained=(spatial.ct@meta.data[,opt$celltype]>0)
spatial_tmp1=subset(spatial.ct, subset=contained == "TRUE")
spatial_tmp2=subset(spatial.ct, subset=contained != "TRUE")
cells_cmbn_list=list()
for (j in 1:length(which(colSums(spatial_tmp1@meta.data[,cells_other])>0))){

  cells_cmbn_list[[j]]=c(0,1)

}
cells_cmbn=expand.grid( cells_cmbn_list)==0
colnames(cells_cmbn)=names(which(colSums(spatial_tmp1@meta.data[,cells_other])>0))

curr.df=list()
df=c()
vv=1
spatial=spatial_tmp1 %>%     
        SCTransform( ., assay = "Spatial", verbose = FALSE) %>%
        RunPCA(verbose = FALSE)
for (v in 1:nrow(cells_cmbn)){
  spots_list=c()
  for (h in 1:ncol(cells_cmbn)){
    if(h!=1){
      spots_list=intersect(spots_list, which((spatial_tmp1@meta.data[,names(which(colSums(spatial_tmp1@meta.data[,cells_other])>0))]==0)[,h]==cells_cmbn[v,h]))
    }else{
      spots_list=which((spatial_tmp1@meta.data[,names(which(colSums(spatial_tmp1@meta.data[,cells_other])>0))]==0)[,h]==cells_cmbn[v,h])
    }
  }

  if(length(spots_list)==0){print(paste0("v = ",v," skiped"));next()}

  #curr.spatial=subset(spatial_tmp1,cells=unique(c(1,2,3,4,5,6,7,8,9,spots_list)))   ###增加几个spot，否则会报错
  # 空转数据SCT转化
  # curr.spatial=subset(spatial_tmp1,cells=unique(c(1:10,spots_list))) %>%     ###增加几个spot，否则会报错
  #              SCTransform( ., assay = "Spatial", verbose = FALSE) %>%
  #              RunPCA(verbose = FALSE)
  curr.spatial=spatial
  curr.cell=colnames(cells_cmbn)[cells_cmbn[v,]]
  curr.sc=subset( sc4subcluster, subset=celltype %in% c(curr.cell,opt$celltype)) #%>%
  # SCTransform(., verbose = opt$verbose) %>%
  # RunPCA(., verbose = opt$verbose) %>%
  # RunUMAP(., dims = 1:30, verbose = opt$verbose)   # 在前面一起标准化
  Idents(curr.sc)="new.ident"
  #mapper=unique(factor(curr.sc@meta.data$new.ident))
  curr.sc$new.ident=factor(curr.sc$new.ident,levels=mapper)
  # curr.cluster_markers <- FindAllMarkers(object = curr.sc, 
  #                                     assay = "SCT",
  #                                     slot = "data",
  #                                     verbose = opt$verbose, 
  #                                     only.pos = TRUE,             # 只选择在cluster上调的基因，only.pos=T；
  #                                     logfc.threshold = 0.5,      # 从SCT的data数据中获取差异基因，为了获得所有可能的差异基因，不需要太高的倍数变化；
  #                                     min.pct = 0.5)               # 选择在同类细胞中半数细胞都表达的基因，min.pct=0.5；如果效果不好，则继续降低
 # 获取每个点spot的预测分数
  anchors <- FindTransferAnchors(reference = curr.sc, query = curr.spatial, normalization.method = "SCT")
  #sc4spotlight[[sample_list[i,5]]]@meta.data$celltype3=paste0("celltype.",as.numeric(sc4spotlight[[sample_list[i,5]]]@meta.data$celltype))
  predictions.assay <- TransferData(anchorset = anchors, refdata = curr.sc$new.ident, prediction.assay = TRUE,weight.reduction = curr.spatial[["pca"]], dims = 1:30)
  rownames(predictions.assay@data)=gsub("-","_",rownames(predictions.assay@data))
  curr.spatial[["predictions"]] <- predictions.assay
  mapper2 = mapper[which(mapper %in% rownames(predictions.assay@data ))]
  # 绘图     
  DefaultAssay(curr.spatial) <- "predictions"


  curr.ratio=predictions.assay@data[intersect(rownames(predictions.assay@data),setdiff(mapper,colnames(cells_cmbn))),colnames(spatial_tmp1)[spots_list]]   ### 注意，这里又去掉了增加的几个spot
  #nmf_mod <- curr.spotlight[[1]]

  if( length(spots_list)==1 ){
    curr.df[[v]]=curr.ratio/sum(curr.ratio)*curr.spatial@meta.data[colnames(spatial_tmp1)[spots_list],opt$celltype]
    curr.df[[v]] <-curr.df[[v]] %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
    colnames(curr.df[[v]])=c("barcodes",colnames(spatial_tmp1)[spots_list])
  }else{
    curr.df[[v]]=t(t(apply(curr.ratio[,colnames(spatial_tmp1)[spots_list]],2,function(x){x/sum(x)}))*curr.spatial@meta.data[colnames(spatial_tmp1)[spots_list],opt$celltype])
    curr.df[[v]] <-curr.df[[v]] %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
    #colnames( decon_df)[which(colnames( decon_df)!="barcodes")]=rownames( curr.spatial@assays$predictions@data)
  }

  if(vv==1){
    df=curr.df[[v]] 
    vv=vv+1
  }else{
    df=df %>% data.frame() %>% 
    dplyr::full_join(curr.df[[v]], by = "barcodes") 
  }
 
  print(paste0("v = ",v," finished"))
}
df=df  %>% tibble::column_to_rownames("barcodes") %>% t()
rownames(df)=gsub(".","-",rownames(df),fixed=T)

##将信息添加到空间数据中
df_tmp=as.data.frame(df)
key=0
if(nchar(opt$celltype2)>0 & (opt$celltype2 %in% colnames(spatial_tmp1@meta.data))){
  df_tmp$others=as.numeric(1-spatial_tmp1@meta.data[,opt$celltype2]-rowSums(df))
  key=1
}else{
  df_tmp$others=as.numeric(1-rowSums(df,na.rm=T))
}

df_tmp <- df_tmp %>%
data.frame() %>%
tibble::rownames_to_column("barcodes")

colnames(df_tmp)=c("barcodes",colnames(df),"others")


spatial_tmp1@meta.data <- spatial_tmp1@meta.data %>%
tibble::rownames_to_column("barcodes") %>%
dplyr::left_join(df_tmp, by = "barcodes") %>%
tibble::column_to_rownames("barcodes")

spatial_tmp2@meta.data[,colnames(df)] = 0
if(nchar(opt$celltype2)>0 & (opt$celltype2 %in% colnames(spatial_tmp2@meta.data))){
  spatial_tmp2@meta.data$others=1-spatial_tmp2@meta.data[,opt$celltype2]
  key=1
}else{
  spatial_tmp2@meta.data$others=1
}

spatial.sub=merge(spatial_tmp1,spatial_tmp2) 
spatial.sub@meta.data = spatial.sub@meta.data %>% tibble::rownames_to_column("barcodes")

spatial.ct@meta.data = spatial.ct@meta.data %>%
tibble::rownames_to_column("barcodes") %>%
dplyr::left_join(spatial.sub@meta.data[,c(colnames(df),"others","barcodes")], by = "barcodes") %>%
tibble::column_to_rownames("barcodes")



## 所有细胞类型的空间展示

if(key==1){
  cell_types_interest=c(sort(colnames(df)),opt$celltype2)
  cols=c(color.ls[[opt$colors]][1:(ncol(df)+1)],"grey")
  names(cols)=c(cell_types_interest,"others")
  cols2=c(color.ls[[opt$colors]][1:(ncol(df)+1)])
  names(cols2)=c(cell_types_interest)
}else{
  cell_types_interest=sort(colnames(df))
  cols=c(color.ls[[opt$colors]][1:ncol(df)],"grey")
  names(cols)=c(cell_types_interest,"others")
  cols2=color.ls[[opt$colors]][1:ncol(df)]
  names(cols2)=cell_types_interest
}




# plot5=spatial_scatterpie(se_obj = spatial.ct,
#                             cell_types_all = c(cell_types_interest,"others"),
#                              cell_types_interest =cell_types_interest,
#                             #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
#                             img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
#                           #   img_alpha = 0.2,
#                             pie_scale = 0.39) +
#                             scale_color_manual(values=cols, limits = names(cols))+
#                             scale_fill_manual(values=cols,limits = names(cols))

# ggsave(plot5, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sampleID,"_",opt$celltype,".png"),sep="/"),limitsize = FALSE)

# plot4=spatial_scatterpie(se_obj = spatial.ct,
#                             cell_types_all = cell_types_interest,
#                              cell_types_interest = cell_types_interest,
#                             #img_path = paste(project_dir,"02_Cellranger/photos",sampleID,sep="/"),
#                             img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
#                           #   img_alpha = 0.2,
#                             pie_scale = 0.39) +
#                             scale_color_manual(values=cols,limits = names(cols2))+
#                             scale_fill_manual(values=cols,limits = names(cols2))

# ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie2_",sampleID,"_",opt$celltype,".png"),sep="/"),limitsize = FALSE)


# for(j in 1:ncol(df)){
#   cell=colnames(df)[j]
#   if(sum(spatial.ct@meta.data[,cell],na.rm=T)==0) {next()}
#   plot3=SpatialFeaturePlot(
#         object = spatial.ct,
#         image = sampleID,
#         features = cell,
#         #alpha = c(0.1, 1)
#         pt.size.factor = 1.2,
#         stroke=0.1,
#         min.cutoff=0.001
#         ) 
#   plot4=plot3+
#         scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
#         scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
#   ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color1.png"),sep="/"),limitsize = FALSE)
#   ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color2.png"),sep="/"),limitsize = FALSE)
# }

# saveRDS(spatial.ct,paste(rds_dir,paste0(sampleID,"_",opt$celltype,"_spatial_subcell_locations.mapping.ref.rds"),sep="/"))
# system(paste0("touch ",spotlight_dir_sample,"/spotlight_for_subclusters",sampleID,"_",gsub("&","and",gsub(" ","_",opt$celltype)),".finished"))


##亚分群合并
if(2==2){
spatial.ct2=spatial.ct
mum1=1
if("Granulosa_1" %in% colnames(spatial.ct2@meta.data) & "Granulosa_5" %in% colnames(spatial.ct2@meta.data)){
  spatial.ct$Granulosa_1 =spatial.ct2$Granulosa_1 + spatial.ct2$Granulosa_5
  num1=2
}
if("Granulosa_1" %in% colnames(spatial.ct2@meta.data) & (!"Granulosa_5" %in% colnames(spatial.ct2@meta.data))){
  spatial.ct$Granulosa_1 =spatial.ct2$Granulosa_1 
  num1=2
}
if((! "Granulosa_1" %in% colnames(spatial.ct2@meta.data)) & "Granulosa_5" %in% colnames(spatial.ct2@meta.data)){
  spatial.ct$Granulosa_1 = spatial.ct2$Granulosa_5
  num1=2
}
if("Granulosa_2" %in% colnames(spatial.ct2@meta.data) & "Granulosa_3" %in% colnames(spatial.ct2@meta.data) &  "Granulosa_4" %in% colnames(spatial.ct2@meta.data)){
  spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_2 + spatial.ct2$Granulosa_3 + spatial.ct2$Granulosa_4
  num=num1+1
} 
if("Granulosa_2" %in% colnames(spatial.ct2@meta.data) & "Granulosa_3" %in% colnames(spatial.ct2@meta.data) & (! "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_2 + spatial.ct2$Granulosa_3 
  num=num1+1
} 
if("Granulosa_2" %in% colnames(spatial.ct2@meta.data) & (! "Granulosa_3" %in% colnames(spatial.ct2@meta.data)) & (! "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_2  
  num=num1+1
} 
if("Granulosa_2" %in% colnames(spatial.ct2@meta.data) & (! "Granulosa_3" %in% colnames(spatial.ct2@meta.data)) & ( "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_2  + spatial.ct2$Granulosa_4 
  num=num1+1
} 
if((! "Granulosa_2" %in% colnames(spatial.ct2@meta.data)) & (! "Granulosa_3" %in% colnames(spatial.ct2@meta.data)) & ( "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_4 
  num=num1+1
} 
if((! "Granulosa_2" %in% colnames(spatial.ct2@meta.data)) & ( "Granulosa_3" %in% colnames(spatial.ct2@meta.data)) & (! "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_3
  num=num1+1
} 
if((! "Granulosa_2" %in% colnames(spatial.ct2@meta.data)) & ( "Granulosa_3" %in% colnames(spatial.ct2@meta.data)) & ( "Granulosa_4" %in% colnames(spatial.ct2@meta.data))){
   spatial.ct$Granulosa_2 = spatial.ct2$Granulosa_3 + spatial.ct2$Granulosa_4 
  num=num1+1
} 
#spatial.ct$Granulosa_3 = spatial.ct2$Granulosa_4 
if("Granulosa_3" %in% colnames(spatial.ct2@meta.data)) {
  spatial.ct$Granulosa_3=NULL
}
if("Granulosa_4" %in% colnames(spatial.ct2@meta.data)) {
  spatial.ct$Granulosa_4=NULL
}
if("Granulosa_5" %in% colnames(spatial.ct2@meta.data)) {
  spatial.ct$Granulosa_5=NULL
}
cell_types_interest=cell_types_interest[c(1:num,length(cell_types_interest))]
cols=c(color.ls[[opt$colors]][1:length(cell_types_interest)],"grey")
names(cols)=c(cell_types_interest,"others")
cols2=c(color.ls[[opt$colors]][1:length(cell_types_interest)])
names(cols2)=c(cell_types_interest)

spotlight_dir_sample=paste0(spotlight_dir_sample,"_new2")
if(!file.exists(spotlight_dir_sample)){
    dir.create(spotlight_dir_sample,recursive=TRUE)
}

if( strsplit(sample_list[which(sample_list[,1]==sampleID),6],"[.]")[[1]][2] %in% c("jpg","jpeg")){
  plot5=spatial_scatterpie2(se_obj = spatial.ct,
                              cell_types_all = c(cell_types_interest,"others"),
                              cell_types_interest =cell_types_interest,
                              #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                              img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                            #   img_alpha = 0.2,
                            slice=sampleID,
                              pie_scale = 0.39) +
                              scale_color_manual(values=cols, limits = names(cols))+
                              scale_fill_manual(values=cols,limits = names(cols))
}else{
  plot5=spatial_scatterpie(se_obj = spatial.ct,
                              cell_types_all = c(cell_types_interest,"others"),
                              cell_types_interest =cell_types_interest,
                              #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                              img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                            #   img_alpha = 0.2,
                            slice=sampleID,
                              pie_scale = 0.39) +
                              scale_color_manual(values=cols, limits = names(cols))+
                              scale_fill_manual(values=cols,limits = names(cols))
}

ggsave(plot5, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sampleID,"_",opt$celltype,".png"),sep="/"),limitsize = FALSE)


 
# plot5=scatterpie_plot(se_obj = spatial.ct,
#                             cell_types_all = c(cell_types_interest,"others"),
#                              cell_types_interest =cell_types_interest,
#                             #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
#                             #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
#                             slice="M1_1",
#                           #   img_alpha = 0.2,
#                             pie_scale = 0.39) +
#                             scale_color_manual(values=cols, limits = names(cols))+
#                             scale_fill_manual(values=cols,limits = names(cols))

# ggsave(plot5, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sampleID,"_",opt$celltype,".png"),sep="/"),limitsize = FALSE)

if( strsplit(sample_list[which(sample_list[,1]==sampleID),6],"[.]")[[1]][2] %in% c("jpg","jpeg")){
plot4=spatial_scatterpie2(se_obj = spatial.ct,
                            cell_types_all = cell_types_interest,
                             cell_types_interest = cell_types_interest,
                            #img_path = paste(project_dir,"02_Cellranger/photos",sampleID,sep="/"),
                            img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                            slice=sampleID,
                          #   img_alpha = 0.2,
                            pie_scale = 0.39) +
                            scale_color_manual(values=cols,limits = names(cols2))+
                            scale_fill_manual(values=cols,limits = names(cols2))
}else{

plot4=spatial_scatterpie(se_obj = spatial.ct,
                            cell_types_all = cell_types_interest,
                             cell_types_interest = cell_types_interest,
                            #img_path = paste(project_dir,"02_Cellranger/photos",sampleID,sep="/"),
                            img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                            slice=sampleID,
                          #   img_alpha = 0.2,
                            pie_scale = 0.39) +
                            scale_color_manual(values=cols,limits = names(cols2))+
                            scale_fill_manual(values=cols,limits = names(cols2))
}

ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie2_",sampleID,"_",opt$celltype,".png"),sep="/"),limitsize = FALSE)


for(j in 1:length(which(cell_types_interest!=opt$celltype2))){
  cell=cell_types_interest[which(cell_types_interest!=opt$celltype2)][j]
  if(sum(spatial.ct@meta.data[,cell],na.rm=T)==0) {next()}
  plot3=SpatialFeaturePlot(
        object = spatial.ct,
        image = sampleID,
        features = cell,
        #alpha = c(0.1, 1)
        pt.size.factor = 1.2,
        stroke=0.1,
        min.cutoff=0.001
        ) 
  plot4=plot3+
        scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
        scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
  ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color1.png"),sep="/"),limitsize = FALSE)
  ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color2.png"),sep="/"),limitsize = FALSE)
}


}
###
if(1==2){
spatial.ct2=spatial.ct

spatial.ct$venous_SMCs=NULL
spatial.ct$arterial_SMCs=NULL

spatial.ct$arterial_SMCs= spatial.ct2$"Smooth muscle_0" +spatial.ct2$"Smooth muscle_1" +spatial.ct2$"Smooth muscle_2"
spatial.ct$venous_SMCs= spatial.ct2$"Smooth muscle_3" +spatial.ct2$"Smooth muscle_4" +spatial.ct2$"Smooth muscle_5" +spatial.ct2$"Smooth muscle_6"

spatial.ct$others=1-spatial.ct$arterial_SMCs-spatial.ct$venous_SMCs
#spatial.ct$Granulosa_3 = spatial.ct2$Granulosa_4 
#spatial.ct$Granulosa_4=NULL
#spatial.ct$Granulosa_5=NULL
cell_types_interest=c("arterial_SMCs","venous_SMCs",opt$celltype2)
cols=c(color.ls[[opt$colors]][1:length(cell_types_interest)],"grey")
names(cols)=c(cell_types_interest,"others")
cols2=c(color.ls[[opt$colors]][1:length(cell_types_interest)])
names(cols2)=c(cell_types_interest)

spotlight_dir_sample=paste0(spotlight_dir_sample,"_new")
if(!file.exists(spotlight_dir_sample)){
    dir.create(spotlight_dir_sample,recursive=TRUE)
}


plot5=spatial_scatterpie(se_obj = spatial.ct,
                            cell_types_all = c(cell_types_interest,"others"),
                             cell_types_interest =cell_types_interest,
                            #img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                            img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                          #   img_alpha = 0.2,
                            pie_scale = 0.39) +
                            scale_color_manual(values=cols, limits = names(cols))+
                            scale_fill_manual(values=cols,limits = names(cols))

ggsave(plot5, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sampleID,"_",opt$celltype,"_rename.png"),sep="/"),limitsize = FALSE)

plot4=spatial_scatterpie(se_obj = spatial.ct,
                            cell_types_all = cell_types_interest,
                             cell_types_interest = cell_types_interest,
                            #img_path = paste(project_dir,"02_Cellranger/photos",sampleID,sep="/"),
                            img_path = paste(project_dir,"02_Cellranger/photos",sample_list[which(sample_list[,1]==sampleID),6],sep="/"),
                          #   img_alpha = 0.2,
                            pie_scale = 0.39) +
                            scale_color_manual(values=cols,limits = names(cols2))+
                            scale_fill_manual(values=cols,limits = names(cols2))

ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie2_",sampleID,"_",opt$celltype,"_rename.png"),sep="/"),limitsize = FALSE)


for(j in 1:length(which(cell_types_interest!=opt$celltype2))){
  cell=cell_types_interest[which(cell_types_interest!=opt$celltype2)][j]
  if(sum(spatial.ct@meta.data[,cell],na.rm=T)==0) {next()}
  plot3=SpatialFeaturePlot(
        object = spatial.ct,
        image = sampleID,
        features = cell,
        #alpha = c(0.1, 1)
        pt.size.factor = 1.2,
        stroke=0.1,
        min.cutoff=0.001
        ) 
  plot4=plot3+
        scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
        scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
  ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color1.png"),sep="/"),limitsize = FALSE)
  ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sampleID,"_",cell,".color2.png"),sep="/"),limitsize = FALSE)
}


}