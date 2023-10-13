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


project_dir=getwd()

rds_dir=paste(project_dir,"suerat_Rdata",sep="/")

#read sample information matrix
#if more than one sample, should give the list csv used in cellranger aggr
cellranger_csv=paste(project_dir,"02_Cellranger","all.csv",sep="/")

if(file.exists(cellranger_csv)){
    sample_list=read.csv(cellranger_csv,stringsAsFactors =F)
}else{
    stop("no cellranger_csv file was found\n")
}

if(nrow(sample_list) < 2){
  stop("this pipline need at least 2 sample\n")
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

for(i in 1:nrow(sample_list)){
    spatial.ob <- readRDS(paste(rds_dir,paste0(sample_list[i,1],"_spatial_cell_locations.mapping.ref.rds"),sep="/"))
    mapper=c("Theca & stroma","Endothelial","monocyte","Smooth muscle","NK","T","oocyte","Granulosa")
    mapper=mapper[which(mapper %in%  names(spatial.ob@meta.data))]
    cols=c("#DC7B1E","#209EBB","#19679A","#C54733","#FDAE91","#548F40","#ACB6B6","#F0E8B3")
    meta_add=fread(paste(rds_dir,paste0(sample_list[i,1],"_meta_df.tsv"),sep="/"))
    spatial.ob@meta.data <- spatial.ob@meta.data %>% tibble::rownames_to_column("V1") %>%
      dplyr::left_join(meta_add, by = "V1") %>%  tibble::column_to_rownames("V1") 

    if (length(which(spatial.ob@meta.data$oocyte>0.0))>0){
        oocytes <- rownames(spatial.ob@meta.data)[which(spatial.ob@meta.data$oocyte>0.05)]
        neighbours <- stringr::str_split(spatial.ob@meta.data[which(spatial.ob@meta.data$oocyte>0.05),"neighbour_bcs"],",",simplify=T)%>% as.vector()
        
        bc_remain=c(neighbours,oocytes) %>%unique()
        spatial.sub=subset(spatial.ob,cells=bc_remain)

        plot1 <- SpatialFeaturePlot(spatial.sub, features = "CCL5_CCR5_scores",pt.size.factor = 1.2, crop=F) + theme(legend.position = "left")
        plot2 <- SpatialFeaturePlot(spatial.sub, features = "MDK_LRP1_scores",pt.size.factor = 1.2, crop=F) + theme(legend.position = "left")
        plot3 <- SpatialFeaturePlot(spatial.sub, features = "PTN_NCL_scores",pt.size.factor = 1.2, crop=F) + theme(legend.position = "left")
        

        plot4 <- plotSpatialScatterpie(
          x = spatial.sub,
          y = spatial.sub@meta.data[mapper],
          cell_types = mapper,
          img = paste(project_dir,"02_Cellranger/photos",sample_list[i,6],sep="/"),
          scatterpie_alpha = 1,
          pie_scale = 0.6,
          degrees = -90,
          # Pivot the image on its x axis
          axis = "h") +
          scale_fill_manual(
              values = cols,
              breaks = mapper)

        P1 <- wrap_plots(plot1, plot4)
        P2 <- wrap_plots(plot2, plot4)
        P3 <- wrap_plots(plot3, plot4)
        ggsave(P1, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_CCL5_CCR5_scores",".png"),sep="/"),limitsize = FALSE)
        ggsave(P2, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_MDK_LRP1_scores",".png"),sep="/"),limitsize = FALSE)
        ggsave(P3, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_PTN_NCL_scores",".png"),sep="/"),limitsize = FALSE)
      
        ggsave(P1, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_CCL5_CCR5_scores",".pdf"),sep="/"),limitsize = FALSE)
        ggsave(P2, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_MDK_LRP1_scores",".pdf"),sep="/"),limitsize = FALSE)
        ggsave(P3, width=10, height=5,dpi=600,units ="in",bg="white", filename= paste("08_stlearn",sample_list[i,1],paste0(sample_list[i,1],"_PTN_NCL_scores",".pdf"),sep="/"),limitsize = FALSE)
    }
   

} 

combined.sct=readRDS("suerat_Rdata/combined.sct.rds")
#combined.sct2=readRDS("../2021_03_19/suerat_Rdata/combined.sct.rds")
df=combined.sct@meta.data
samples=c("Y1_1", "Y2_1", "Y3_1", "Y3_2", "Y3_3", "M1_1", "M2_1", "M3_1", "M3_2", "M3_3", "O1_1", "O1_2", "O2_1", "O3_1", "O3_2")
j=0
for (i in 1:length(combined.sct@images)){
  j=j+1
  if( nrow(combined.sct@images[[j]]@coordinates)==0){
      combined.sct@images[[j]]=NULL
      j=j-1
  }
}
names(combined.sct@images)=as.character(sample_list[,1])

##基因表达
spatial.ob2=list()
for (sample in samples){
    spatial.ob2[[sample]]=subset(combined.sct, subset=orig.ident == eval(sample))
    DefaultAssay(spatial.ob2[[sample]]) <- "SCT"
}

for (sample in sample_list[,1]){
    DefaultAssay(spatial.ob2[[sample]]) <- "SCT"


    j=0
    for (i in 1:length(spatial.ob2[[sample]]@images)){
      j=j+1
      if( nrow(spatial.ob2[[sample]]@images[[j]]@coordinates)==0){
          spatial.ob2[[sample]]@images[[j]]=NULL
          j=j-1
      }
    }
    names(spatial.ob2[[sample]]@images)=as.character(sample)

    plot2=SpatialPlot(
            object =spatial.ob2[[sample]],
            slot="data",
            features = "LMNA",
            alpha = c(1, 1),
            crop=TRUE,
            image.alpha=0,
            # pt.size.factor = 1.2,
            stroke=0.1,
            min.cutoff=0.001,
            # cols=15,
            images=sample,

            )  +
        theme(plot.title = element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.ticks.length = unit(.1, "cm"),
            axis.ticks=element_line(color="black"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.justification = c(0.5,0.5),
            legend.position="right",
            legend.text=element_text(face="plain",size=10),
            legend.title = element_text(face="plain",size=12),
            legend.key.width=unit(0.5,'cm'),legend.key.height=unit(0.5,'cm'),
            #plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
            strip.text = element_text(face = "bold"),
            panel.background = element_rect(fill = "black")
        )+scale_fill_distiller(palette = "Spectral",limits = c(0, 4))

    # plot3=plot2+ 
    #         scale_color_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))+
    #         scale_fill_gradientn(colours=c("#371C47","#3B528B","#21908C","#5DC863","#FDE529"))

    #ggsave(plot2, width=5, height=5,dpi=900,units ="in",bg="white", filename=paste0("genes_exp/SpatialFeaturePlot_",i,".color1.pdf"),limitsize = FALSE)
    ggsave(plot2, width=6, height=5, dpi=600,units ="in",bg="white", filename=paste0("genes_exp/LMNA/",sample,"_SpatialFeaturePlot_LMNA.pdf"),limitsize = FALSE)

}

