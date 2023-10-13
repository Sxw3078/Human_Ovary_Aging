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
        make_option(c("-d", "--spatial_sep"), help="rds for spatial data of separated samples, default %default",default="suerat_Rdata/separated.sct.Rdata"),
        make_option(c("-d", "--spatial_int"), help="rds for spatial data of intergrated samples, default %default",default="suerat_Rdata/separated.sct.rds"),
        make_option(c("-z", "--singcell"), help="rds for singcell data, default %default",default="suerat_Rdata/seuratObject_celltype.rds"),      
        make_option(c("--colors"), help="color set for scatterpie, choose from NPG,NEJM,JAMA,Lancet,AAAS,report,DEF, custom, default %default",default="custom"),        
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

if(file.exists(cellranger_csv)){
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

scatterpie_plot <- function(se_obj,
                            cell_types_all,
                            slice = NULL,
                            scatterpie_alpha = 1,
                            cell_types_interest = NULL,
                            pie_scale = 0.4) {

  # Check variables
  if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
  if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is.numeric(scatterpie_alpha)) stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale)) stop("ERROR: pie_scale must be numeric between 0 and 1!")

  # Loading libraries
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

  metadata_ds <- data.frame(se_obj@meta.data)

  colnames(metadata_ds) <- colnames(se_obj@meta.data)

  if (is.null(cell_types_interest)) {
    cell_types_interest <- cell_types_all
  }

  # If not all cell types are in the cell types of interest we only want to keep those spots which have at least one of the cell types of interest
  if (!all(cell_types_all %in% cell_types_interest)) {

    metadata_ds <- metadata_ds %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::mutate(rsum = rowSums(.[, cell_types_interest, drop = FALSE])) %>%
      dplyr::filter(rsum != 0) %>%
      dplyr::select("barcodeID") %>%
      dplyr::left_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),by = "barcodeID") %>%
      tibble::column_to_rownames("barcodeID")
  }

  ## If slice is not selected set it to the first element in the list of slices
  if (is.null(slice) | (!is.null(slice) && !slice %in% names(se_obj@images))) {
    slice <- names(se_obj@images)[1]
    print(sprintf("Using slice %s", slice))
  }

  ## Preprocess data
  spatial_coord <- data.frame(se_obj@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::mutate(imagerow_scaled = imagerow * se_obj@images[[slice]]@scale.factors$lowres,
                  imagecol_scaled = imagecol * se_obj@images[[slice]]@scale.factors$lowres) %>%
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"), by = "barcodeID")

  # Plot the scatterplot
  scatterpie_plt <- suppressMessages(ggplot() +
                     scatterpie::geom_scatterpie(data = spatial_coord,
                                                 ggplot2::aes(
                                                   x = imagecol_scaled,
                                                   y = imagerow_scaled),
                                                 cols = cell_types_all,
                                                 color = NA,
                                                 alpha = scatterpie_alpha,
                                                 pie_scale = pie_scale) +
                     ggplot2::scale_y_reverse() +
                     ggplot2::theme_void())

  return(scatterpie_plt)

}

spatial_scatterpie <- function(se_obj,
                               cell_types_all,
                               img_path,
                               cell_types_interest = NULL,
                               slice = NULL,
                               scatterpie_alpha = 1,
                               pie_scale = 1) {

  # Check variables
  if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
  if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
  if (!is.character(img_path)) stop("ERROR: must be a character string!")
  if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
  if (!is.numeric(scatterpie_alpha)) stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
  if (!is.numeric(pie_scale)) stop("ERROR: pie_scale must be numeric between 0 and 1!")

  # Loading libraries
  suppressMessages(require(ggplot2))
  suppressMessages(require(cowplot))
  suppressMessages(require(imager))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))
  suppressMessages(require(png))
  suppressMessages(require(jpeg))
  suppressMessages(require(grid))

  metadata_ds <- data.frame(se_obj@meta.data)

  colnames(metadata_ds) <- colnames(se_obj@meta.data)

  if (is.null(cell_types_interest)) {
    cell_types_interest <- cell_types_all
  }

  # If not all cell types are in the cell types of interest we only want to keep those spots which have at least one of the cell types of interest
  if (!all(cell_types_all %in% cell_types_interest)) {

    metadata_ds <- metadata_ds %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest,
                                           drop = FALSE])) %>%
      dplyr::filter(rsum != 0) %>%
      dplyr::select("barcodeID") %>%
      dplyr::left_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                       by = "barcodeID") %>%
      tibble::column_to_rownames("barcodeID")
  }

  ## If slice is not selected set it to the first element in the list of slices
  if (is.null(slice) | (!is.null(slice) && !slice %in% names(se_obj@images))) {
    slice <- names(se_obj@images)[1]
    warning(sprintf("Using slice %s", slice))
  }

  ## Preprocess data
  spatial_coord <- data.frame(se_obj@images[[slice]]@coordinates) %>%
    tibble::rownames_to_column("barcodeID") %>%
    dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                      by = "barcodeID")

  ### Load histological image into R
  #### Extract file format, JPEG or PNG
  img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))

  if(img_frmt %in% c(".jpg", "jpeg")) {
    img <- jpeg::readJPEG(img_path)
  } else if (img_frmt == ".png") {
    img <- png::readPNG(img_path)
  }

   # Convert image to grob object
  img_grob <- grid::rasterGrob(img,
                               interpolate = FALSE,
                               width = grid::unit(1, "npc"),
                               height = grid::unit(1, "npc"))

  #img_grob$raster = as.raster(matrix(alpha(img_grob$raster,0.5),nrow=nrow(img_grob$raster),byrow=T)) #没有用

  ## Plot spatial scatterpie plot
  scatterpie_plt <- suppressMessages(
    ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grob = img_grob,
        xmin = 0,
        xmax = ncol(img),
        ymin = 0,
        ymax = -nrow(img)) +
      scatterpie::geom_scatterpie(
        data = spatial_coord,
        ggplot2::aes(x = imagecol,
                     y = imagerow),
                     cols = cell_types_all,
                     color = NA,
                     alpha = scatterpie_alpha,
                     pie_scale = pie_scale) +
      ggplot2::scale_y_reverse() +
      ggplot2::ylim(nrow(img), 0) +
      ggplot2::xlim(0, ncol(img)) +
      cowplot::theme_half_open(11, rel_small = 1) +
      ggplot2::theme_void() +
      ggplot2::coord_fixed(ratio = 1,
                           xlim = NULL,
                           ylim = NULL,
                           expand = TRUE,
                           clip = "on"))
  return(scatterpie_plt)
}

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
		"#00468B", "#EC0000", "#41B43F", "#0099B3", "#925E9F", "#FDAE91", "#AC002A",
		"#ACB6B6", "#9E3B4F", "#B2811E", "#40A67D", "#6B7EA9", "#C98599", "#D7675A",
		"#B86E6B", "#6A7BA1", "#6E4D6D", "#D1793E", "#5EAD73", "#6694AE", "#AE80A1",
		"#EAA392", "#B46264", "#949CAB", "#C65356", "#8D9D4B", "#4D9F9B", "#8A7CA4",
		"#E3ABA9", "#C26161", "#B69C9A", "#596D96", "#46688B", "#EC7676", "#7BB47A",
		"#5AA6B3", "#997F9F", "#FDD6C7", "#AC566B", "#B1B6B6", "#9E6D76", "#B29A68",
		"#73A692", "#8A93A9", "#C9A7B1", "#D79F99", "#B89392", "#868EA1", "#6E5E6E",
		"#D1A588"
    )[-5],
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
    ),
    'custom'=c(
    "#DC7B1E","#209EBB","#19679A","#C54733","#FDAE91","#548F40","#ACB6B6","#F0E8B3"
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


###############################载入自定义功能
###############################结束

###############################载入单细胞RNA数据
###############################

sc_data=readRDS(opt$singcell)

if(length(unique(sc_data@meta.data$sample))>1){
    sc_data.list <- SplitObject(sc_data, split.by = "sample")
}else{
    sc_data.list=sc_data
}

sc_sample_list=unique(sc_data@meta.data$sample)

cols=color.ls[[opt$colors]][1:length(unique(sc_data$celltype))]
names(cols)=as.character(unique(sc_data$celltype))

cols2=color.continuous.ls[['RdBu']]

###############################载入空间RNA数据
###############################

load(opt$spatial_sep)
if(is.null(spatial.ob)){
    {writeLines("spatial.ob was not exist in the .Rdata space specified by the \"-d\" or \"--spatial_sep\" parameter, please check it manually."); q("no")}
}

###############################循环进行spotlight
###############################开始

spotlight_dir=paste(project_dir,"06_mapping_ref",sep="/")
if(!file.exists(spotlight_dir)){
    dir.create(spotlight_dir,recursive=TRUE)
}

sc4spotlight=cluster_markers_all=spotlight=decon_mtrx=nmf_mod=cell_types_all=decon_mtrx_sub=decon_mtrx_new=list()

for(i in 1:nrow(sample_list)){

  print(paste0("Now sample ",sample_list[i,1]," is processing!\n"))

  spotlight_dir_sample=paste(spotlight_dir,sample_list[i,1],sep="/")
  if(!file.exists(spotlight_dir_sample)){
      dir.create(spotlight_dir_sample,recursive=TRUE)
  }
  
  if(file.exists(paste(spotlight_dir_sample,"mapping_ref.finished",sep="/"))){
      next()
  }

    if(sample_list[i,5] %in% sc_sample_list){
      # 单细胞数据SCT转化
      sc4spotlight[[sample_list[i,5]]] = SCTransform(sc_data.list[[sample_list[i,5]]],ncells=3000, verbose = opt$verbose) %>%
      RunPCA(., verbose = opt$verbose) %>%
      RunUMAP(., dims = 1:30, verbose = opt$verbose)
      mapper=levels(sc4spotlight[[sample_list[i,5]]]@meta.data$celltype)
      # 空转数据SCT转化
      spatial.ob[[sample_list[i,1]]] <- SCTransform( spatial.ob[[sample_list[i,1]]], assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
      # 获取每个点spot的预测分数
      anchors <- FindTransferAnchors(reference = sc4spotlight[[sample_list[i,5]]] , query = spatial.ob[[sample_list[i,1]]], normalization.method = "SCT",reduction="cca")
      #sc4spotlight[[sample_list[i,5]]]@meta.data$celltype3=paste0("celltype.",as.numeric(sc4spotlight[[sample_list[i,5]]]@meta.data$celltype))
      predictions.assay <- TransferData(anchorset = anchors, refdata = sc4spotlight[[sample_list[i,5]]]$celltype, prediction.assay = TRUE,weight.reduction = spatial.ob[[sample_list[i,1]]][["pca"]], dims = 1:30)
      spatial.ob[[sample_list[i,1]]][["predictions"]] <- predictions.assay
      mapper = mapper[which(mapper %in% rownames(predictions.assay@data ))]
      # 绘图     
      DefaultAssay(spatial.ob[[sample_list[i,1]]]) <- "predictions"
      plot1 =SpatialFeaturePlot(spatial.ob[[sample_list[i,1]]], features = mapper, pt.size.factor = 1.6, ncol = 2, crop = TRUE)
      nrow=ceiling(length(mapper)/2)
      ggsave(plot1, width=13.5, height=6*nrow,dpi=300, filename=paste(spotlight_dir_sample,paste0("SpatialFeaturePlot_all_cells_",sample_list[i,1],".png"),sep="/"),limitsize = FALSE)
 
      # 查找空间分布限制的细胞类型（一般没有必要）
      #spatial.ob[[sample_list[i,1]]] <- FindSpatiallyVariableFeatures(spatial.ob[[sample_list[i,1]]], assay = "predictions", selection.method = "markvariogram",  features = rownames(spatial.ob[[sample_list[i,1]]]), r.metric = 5, slot = "data")
      #top.clusters <- head(SpatiallyVariableFeatures(spatial.ob[[sample_list[i,1]]]), 4)
      #SpatialPlot(object = spatial.ob[[sample_list[i,1]]], features = top.clusters, ncol = 2)
      #ggsave(plot2, width=5, height=5,dpi=300, filename=paste(spotlight_dir_sample,paste0("test2",sample_list[i,1],".png"),sep="/"),limitsize = FALSE)


    }else{
      if(is.null(cluster_markers_all[["all"]] )){
          if(file.exists(paste(rds_dir,"sc_combined.sct.rds",sep="/"))){
            combined.sct=readRDS(paste(rds_dir,"sc_combined.sct.rds",sep="/"))
          }else{
            combined=merge(sc_data.list[[1]],sc_data.list[2:length(sc_data.list)],add.cell.ids=names(sc_data.list))
            # combined=merge(sc_data.list[[3]],sc_data.list[[4]],add.cell.ids=names(sc_data.list)[3:4])  ## 自定制
            VariableFeatures(combined)=rownames(combined@assays$RNA@data)
            combined.sct = SCTransform(combined,ncells=3000, verbose = opt$verbose) %>%
            RunPCA(., verbose = opt$verbose) %>%
            RunUMAP(., dims = 1:30, verbose = opt$verbose)
            
            # saveRDS(combined.sct,paste(rds_dir,"sc_combined.Y38_LHM.Y38_LWY.sct.rds",sep="/"))
            saveRDS(spotlight[[sample_list[i,1]]],paste(rds_dir,paste0("mapping_ref.",sample_list[i,1],".rds"),sep="/"))
          }
       }
                                          
   
      mapper=as.character(unique(factor(combined.sct@meta.data$celltype)))
      # 空转数据SCT转化
      spatial.ob[[sample_list[i,1]]] <- SCTransform( spatial.ob[[sample_list[i,1]]], assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
      # 获取每个点spot的预测分数
      anchors <- FindTransferAnchors(reference = combined.sct , query = spatial.ob[[sample_list[i,1]]], normalization.method = "SCT")
      #combined.sct@meta.data$celltype3=paste0("celltype.",as.numeric( combined.sct@meta.data$celltype))
      predictions.assay <- TransferData(anchorset = anchors, refdata = combined.sct$celltype, prediction.assay = TRUE,weight.reduction = spatial.ob[[sample_list[i,1]]][["pca"]], dims = 1:30)
      spatial.ob[[sample_list[i,1]]][["predictions"]] <- predictions.assay
       mapper = mapper[which(mapper %in% rownames(predictions.assay@data ))]
      # 绘图     
      DefaultAssay(spatial.ob[[sample_list[i,1]]]) <- "predictions"
      plot1 =SpatialFeaturePlot(spatial.ob[[sample_list[i,1]]], features = mapper, pt.size.factor = 1.6, ncol = 2, crop = TRUE)
      nrow=ceiling(length(mapper)/2)
      ggsave(plot1, width=12.5, height=7*nrow,dpi=300, filename=paste(spotlight_dir_sample,paste0("SpatialFeaturePlot_all_cells_",sample_list[i,1],".png"),sep="/"),limitsize = FALSE)
    }
 
    decon_df<-spatial.ob[[sample_list[i,1]]]@assays$predictions@data %>%
    t() %>% data.frame() %>%
    tibble::rownames_to_column("barcodes")
    colnames( decon_df)[which(colnames( decon_df)!="barcodes")]=rownames( spatial.ob[[sample_list[i,1]]]@assays$predictions@data)


    spatial.ob[[sample_list[i,1]]]@meta.data <-  spatial.ob[[sample_list[i,1]]]@meta.data %>% 
    #  dplyr::select(!(colnames( decon_mtrx_sub[[sample_list[i,1]]],))) %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")        

    ## 所有细胞类型的空间展示

    # plot2=spatial_scatterpie(se_obj = spatial.ob[[sample_list[i,1]]],
    #                             cell_types_all = mapper,
    #                             img_path = paste(project_dir,"02_Cellranger/photos",sample_list[i,6],sep="/"),
    #                             #   img_alpha = 0.2,
    #                             pie_scale = 0.39) +
    #                             scale_color_manual(values=cols)+
    #                             scale_fill_manual(values=cols)

    plot2 <- plotSpatialScatterpie(
        x = spatial.ob[[sample_list[i,1]]],
        y = spatial.ob[[sample_list[i,1]]]@meta.data[mapper],
        cell_types = mapper,
        img = paste(project_dir,"02_Cellranger/photos",sample_list[i,6],sep="/"),
        scatterpie_alpha = 1,
        pie_scale = 0.4,
        degrees = -90,
        # Pivot the image on its x axis
        axis = "h") +
        scale_fill_manual(
            values = cols,
            breaks = mapper)

    ggsave(plot2, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sample_list[i,1],".png"),sep="/"),limitsize = FALSE)
    #ggsave(plot5, width=5, height=5,dpi=300,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sample_list[i,1],".pdf"),sep="/"),limitsize = FALSE)

    for(j in 1:length( mapper)){
        cell=mapper[j]
        if(sum(spatial.ob[[sample_list[i,1]]]@meta.data[[cell]])==0) {next()}
        plot3=SpatialFeaturePlot(
                object = spatial.ob[[sample_list[i,1]]],
                features = cell,
                #alpha = c(0.1, 1),
                crop=TRUE,
                pt.size.factor = 1.2,
                stroke=0.1,
                min.cutoff=0.001
                ) 
        plot4=plot3 + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))

        plot5= spatial_scatterpie(se_obj = spatial.ob[[sample_list[i,1]]],
                                cell_types_all =mapper,
                                img_path = paste(project_dir,"02_Cellranger/photos",sample_list[i,6],sep="/"),
                                cell_types_interest = cell,
                                #  img_alpha = 0.2,
                                pie_scale = 0.39) +
                                scale_color_manual(values=cols)+
                                scale_fill_manual(values=cols)
        # plot5 <- plotSpatialScatterpie(
        #     x = spatial.ob[[sample_list[i,1]]],
        #     y = spatial.ob[[sample_list[i,1]]]@meta.data[mapper],
        #     cell_types = cell,
        #     img = paste(project_dir,"02_Cellranger/photos",sample_list[i,6],sep="/"),
        #     scatterpie_alpha = 1,
        #     pie_scale = 0.4,
        #     degrees = -90,
        #     # Pivot the image on its x axis
        #     axis = "h") +
        #     scale_fill_manual(
        #         values = cols,
        #         breaks = mapper)

        ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sample_list[i,1],"_",cell,".color1.png"),sep="/"),limitsize = FALSE)
        ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sample_list[i,1],"_",cell,".color2.png"),sep="/"),limitsize = FALSE)
        ggsave(plot5, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sample_list[i,1],"_",cell,".png"),sep="/"),limitsize = FALSE)
        #ggsave(plot3, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatialDimPlot_",sample_list[i,1],"_",cell,".pdf"),sep="/"),limitsize = FALSE)
        #ggsave(plot4, width=5, height=5,dpi=600,units ="in",bg="white", filename= paste(spotlight_dir_sample,paste0("spatial_scatterpie_",sample_list[i,1],"_",cell,".pdf"),sep="/"),limitsize = FALSE)
    }

    saveRDS(spatial.ob[[sample_list[i,1]]],paste(rds_dir,paste0(sample_list[i,1],"_spatial_cell_locations.mapping.ref.rds"),sep="/"))

}

###############################循环进行mapping_ref
###############################结束


