################################输入参数
################################

suppressPackageStartupMessages(library(optparse))
options(bitmapType='cairo')
option_list <- list(
        make_option(c("-w", "--workdir"), help="script work directory ,defualt is run directory " ),
        make_option(c("-s", "--sampleID"), help="which sample, default %default",default="Y48"),        
        make_option(c("-g", "--geneset"), help="the geneset file, default %default",default="09_pathway_mapping/Pathway_GeneName.xls"), 
        make_option(c("-v", "--verbose"), help="Shown more processing information , default %default",default=FALSE),
        make_option(c("-p", "--ptsize"), help="size of the points , default %default",default=0.85),   # Y18:0.85 Y28:0.91  Y38:1  Y48:0.85
        make_option(c("-y", "--ymove"), help="size of the points , default %default",default=0)       # Y18:0    Y28:18   Y38:-32  Y48：-25
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

result_dir=paste0(project_dir,"/09_pathway_mapping/",opt$sampleID)
if(!file.exists(result_dir)){
    dir.create(result_dir,recursive=TRUE)
}


hlpr_image_add_on3=function (object, display_image, of_sample) 
{
    if (base::isTRUE(display_image)) {
        sample_image <- getImage(object, of_sample)
        if ("Image" %in% base::class(sample_image)) {
            image_raster <- grDevices::as.raster(x = sample_image)
            img_info <- image_raster %>% magick::image_read() %>% 
                magick::image_info()
            st_image <- image_raster %>% magick::image_read() 
            image_add_on <- ggplot2::annotation_raster(raster = st_image, 
                xmin = 0, ymin = 0, xmax = img_info$width, ymax = img_info$height)
        }
        else {
            base::warning(glue::glue("Content of slot 'image' for sample '{of_sample}' must be of class 'Image' not of class '{base::class(sample_image)}'."))
            image_add_on <- list()
        }
    }
    else {
        image_add_on <- list()
    }
    base::return(image_add_on)
}


plotSurface2=function (object, color_by = NULL, method_gs = NULL, normalize = NULL, 
    smooth = NULL, smooth_span = NULL, pt_alpha = NULL, pt_clr = NULL, 
    pt_clrp = NULL, pt_clrsp = NULL, pt_size = NULL, clrp_adjust = NULL, 
    display_image = NULL, display_title = NULL, complete = NULL, 
    verbose = NULL, of_sample = NA, ...) 
{
    SPATA2:::hlpr_assign_arguments(object)
    check_pt(pt_size, pt_alpha, pt_clrsp)
    SPATA2:::check_display(display_title, display_image)
    of_sample <- check_sample(object = object, of_sample = of_sample, 
        desired_length = 1)
    if (!base::is.null(color_by)) {
        color_to <- SPATA2:::check_color_to(color_to = color_by, all_genes = getGenes(object, 
            of_sample = of_sample), all_gene_sets = getGeneSets(object), 
            all_features = getFeatureNames(object, of_sample = of_sample))
    }
    else {
        color_to <- list(color = pt_clr)
    }
    coords_df <- getCoordsDf(object, of_sample = of_sample)

   # coords_df$x= max(coords_df$x)+min(coords_df$x)-coords_df$x 
    coords_df$y= 1.01*(max(coords_df$y)+min(coords_df$y))-coords_df$y + opt$ymove

    plot_list <- SPATA2:::hlpr_scatterplot(object = object, spata_df = coords_df, 
        color_to = color_to, pt_size = pt_size, pt_alpha = pt_alpha, 
        pt_clrp = pt_clrp, pt_clrsp = pt_clrsp, method_gs = method_gs, 
        normalize = normalize, smooth = smooth, smooth_span = smooth_span, 
        verbose = verbose, complete = complete, clrp.adjust = c(subs.by.segm = "lightgrey", 
            clrp_adjust), display_title = display_title, ...)
    fig=ggplot2::ggplot(data = plot_list$data, mapping = ggplot2::aes(x = x, 
        y = y)) + hlpr_image_add_on3(object, display_image, of_sample) + 
        plot_list$add_on + ggplot2::coord_equal() + ggplot2::theme_void()
    if(!display_image){
                 fig + theme_bw()+
                 theme(axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank(),panel.grid=element_blank())
    }else{
      fig
    }

}

library(SPATA2)
library(magrittr)
library(ggplot2)
library(patchwork)


geneset=read.table(opt$geneset,sep='\t')

spata_obj <-
    initiateSpataObject_10X(directory_10X=paste0("02_Cellranger/",opt$sampleID),sample_name=opt$sampleID)

# for (i in 1:nrow(geneset)){
#     spata_obj=addGeneSet(object=spata_obj,class_name="custom",gs_name=geneset[i,1],genes=strsplit(geneset[i,2],split=";")[[1]])
# }

# create a spata-object from scratch

for (i in 1:nrow(geneset)){
    spata_obj=addGeneSet(object=spata_obj,class_name="custom",gs_name=geneset[i,1],genes=strsplit(geneset[i,2],split=";")[[1]],overwrite=TRUE)

# plot gene-set expression 
    p2 = plotSurface2(object = spata_obj,
                    of_sample = opt$sampleID,
                    color_by = paste0("custom_",geneset[i,1]),
                    pt_size = opt$ptsize,
                    display_image =TRUE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = FALSE
                    ) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
                    

    # plot gene-set expression (spatially smoothed by 0.1)
    p3 <- plotSurface2(object = spata_obj,
                        of_sample =  opt$sampleID,
                    color_by = paste0("custom_",geneset[i,1]),
                    pt_size = opt$ptsize,
                    display_image =FALSE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = TRUE,
                    smooth_span = 0.1) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
                    

    # plot gene-set expression (spatially smoothed by 0.2)
    p4 <- plotSurface2(object = spata_obj,
                        of_sample =  opt$sampleID,
                    color_by = paste0("custom_",geneset[i,1]),
                    pt_size = opt$ptsize,
                    display_image =FALSE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = TRUE,
                    smooth_span = 0.2) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))

                    

    # combine with patchwork 

    ggsave(p2, width=10, height=10, units = "cm", filename=paste0(result_dir,"/custom_",geneset[i,1],".png"),  dpi=900,limitsize = FALSE)

    ggsave(p3, width=10, height=10, units = "cm", filename=paste0(result_dir,"/custom_",geneset[i,1],"_smooth01.png"),  dpi=900,limitsize = FALSE)

    ggsave(p4, width=10, height=10, units = "cm", filename=paste0(result_dir,"/custom_",geneset[i,1],"_smooth02.png"),  dpi=900,limitsize = FALSE)
} 

for ( geneset2 in c("RCTM_DNA_REPAIR","HM_DNA_REPAIR" ,"BP.GO_DNA_REPAIR")){

  p2 = plotSurface2(object = spata_obj,
                    of_sample = opt$sampleID,
                    color_by = geneset2,
                    pt_size = opt$ptsize,
                    display_image =TRUE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = FALSE
                    ) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
                    

    # plot gene-set expression (spatially smoothed by 0.1)
    p3 <- plotSurface2(object = spata_obj,
                        of_sample =  opt$sampleID,
                    color_by =geneset2,
                    pt_size = opt$ptsize,
                    display_image =FALSE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = TRUE,
                    smooth_span = 0.1) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))
                    

    # plot gene-set expression (spatially smoothed by 0.2)
    p4 <- plotSurface2(object = spata_obj,
                        of_sample =  opt$sampleID,
                    color_by = geneset2,
                    pt_size = opt$ptsize,
                    display_image =FALSE,
                    display_title =FALSE,
                    pt_clrsp = "Terrain 2",
                    smooth = TRUE,
                    smooth_span = 0.2) + 
            scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
            scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))

                    

    # combine with patchwork 

    ggsave(p2, width=10, height=10, units = "cm", filename=paste0(result_dir,"/",geneset2,".png"),  dpi=900,limitsize = FALSE)

    ggsave(p3, width=10, height=10, units = "cm", filename=paste0(result_dir,"/",geneset2,"_smooth01.png"),  dpi=900,limitsize = FALSE)

    ggsave(p4, width=10, height=10, units = "cm", filename=paste0(result_dir,"/",geneset2,"_smooth02.png"),  dpi=900,limitsize = FALSE)
}