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
    dplyr::mutate(
      imagerow_scaled =
        imagerow * se_obj@images[[slice]]@scale.factors$lowres,
      imagecol_scaled =
        imagecol * se_obj@images[[slice]]@scale.factors$lowres) %>%
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
        ggplot2::aes(x = imagecol_scaled,
                     y = imagerow_scaled),
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


spatial_scatterpie2 <- function(se_obj,
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