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
