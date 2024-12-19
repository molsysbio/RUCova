#' Calculated mean of normalised Iridium isotopes
#'
#' @param data SingleCellExperiment object with counts in linear scale 
#' @param dna_channels Vector specifying the names of the DNA channels
#' @param q Quantile for normalisation.
#' @return The SingleCellExperiment object with an extra column "mean_BC" in the counts assay.#'
#' @export
#' @import dplyr
#' @import ComplexHeatmap
#' @import circlize
#' @import Matrix
#' @import grid

calc_mean_DNA <- function(data,dna_channels, q) {
  
  # Check if the input is a SingleCellExperiment
  if (inherits(data, "SingleCellExperiment")){
    sce <- data
    dna_data <- assay(sce)[dna_channels, , drop = FALSE]
  } else {
    # Handle input as a data frame or matrix
    dna_data <- t(data[,dna_channels])
    
  }
  # Perform calculations
    asinh_input <- asinh(dna_data)
    quantiles <- apply(asinh_input, 1, quantile, probs = q, names = FALSE)
    scaling_factors <- quantiles[[1]] / quantiles
    scaled_data <- asinh_input * scaling_factors
    mean_DNA <- sinh(apply(scaled_data, 2, mean))
  
  # Add the mean_DNA to the appropriate structure
    if(inherits(data, "SingleCellExperiment")){
      # Add the mean_DNA as a new row in the "counts" assay
      sce <- SingleCellExperiment(
        assays = list(counts = rbind(sce@assays@data$counts,mean_DNA)),    # Assay data
        rowData = DataFrame(measured_channels = c(rowData(sce)$measured_channels,"mean_DNA")),# Metadata for rows
        colData = colData(sce)# Metadata for columns
      )
      
      return(sce)
    } else {
    data$mean_DNA <- mean_DNA
    
    return(data)
  }
}



#' Calculated mean of normalised highest BC per cell
#'
#' @param data SingleCellExperiment object with counts in linear scale 
#' @param bc_channels Vector specifying the names of the BC channels
#' @param q Quantile for normalisation.
#' @param n_bc number of barcoding isotopes per cell.
#' @return The SingleCellExperiment object with an extra column "mean_BC" in the counts assay.
#' @export
#'

calc_mean_BC <- function(data,bc_channels, n_bc, q) {
  
  # Check if the input is a SingleCellExperiment
  if (inherits(data, "SingleCellExperiment")){
    sce <- data
    bc_data <- assay(sce)[bc_channels, , drop = FALSE]
  } else {
    # Handle input as a data frame or matrix
    bc_data <- t(data[,bc_channels])
    
  }
  # Perform calculations
  asinh_input <- asinh(bc_data)
  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = FALSE)
  scaling_factors <- quantiles[[1]] / quantiles
  scaled_data <- asinh_input * scaling_factors
  mean_BC <- scaled_data |>
    apply(2, sort.int, decreasing = TRUE, method = "shell") |>
    #Rfast::colSort(descending = TRUE) |>
    head(n_bc) |>
    colMeans() |>
    sinh()
  
  # Add the mean_BC to the appropriate structure
  if(inherits(data, "SingleCellExperiment")){
    # Add the mean_BC as a new row in the "counts" assay
    sce <- SingleCellExperiment(
      assays = list(counts = rbind(sce@assays@data$counts,mean_BC)),    # Assay data
      rowData = DataFrame(measured_channels = c(rowData(sce)$measured_channels,"mean_BC")),# Metadata for rows
      colData = colData(sce)# Metadata for columns
    )
    
    return(sce)
  } else {
    data$mean_BC <- mean_BC

    return(data)
  }
}

#' Plot pearson correlation coefficients between markers on a double triangular heatmap (lower triangle: before RUCova, upper triangle: after RUCova).
#' @param lower Matrix with pearson correlation coefficient between markers eg.: before RUCova, to be plotted in the lower triangle.
#' @param upper  Matrix with pearson correlation coefficient between markers eg.: after RUCova, to be plotted in the upper triangle.
#' @return #A heatmap
#' @export
#'
#'
heatmap_compare_corr <- function(lower, upper){

tmp <-  as.matrix(Matrix::tril(as.matrix(lower)) + Matrix::triu(as.matrix(upper)))

diag(tmp) <- NA
hm <- tmp |> 
  ComplexHeatmap::Heatmap(col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                          name = "Pearson corr.coef",
                          na_col = "grey",
                          rect_gp = gpar(col = "black", lwd = 0),
                          cluster_rows = FALSE,
                          cluster_columns = FALSE,
                          row_names_side = "left",
                          column_names_side = "top",
                          row_title_rot = 0,
                          row_title_side = "right")

draw(hm)

}

