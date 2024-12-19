#' Calculated mean of normalised Iridium isotopes
#'
#' @param sce SingleCellExperiment object with counts in linear scale 
#' @param dna_channels Vector specifying the names of the DNA channels
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples 
#'
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
#' @param sce SingleCellExperiment object with counts in linear scale 
#' @param dna_channels Vector specifying the names of the DNA channels
#' @param q Quantile for normalisation.
#' @param n_bc number of barcoding isotopes per cell.
#' @return A vector.
#' @examples
#' data <- RUCova::HNSCC_data
#' data |> dplyr::mutate(mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di, n_bc = 4, q = 0.95))
#' 
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
#' @examples 
#' data <- RUCova::HNSCC_data
#' data <- data |> dplyr::mutate(cell_id = 1:n(),
#'                        mean_DNA = RUCova::calc_mean_DNA(DNA_191Ir, DNA_193Ir, q = 0.95),
#'                        mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di,
#'                        Dead_cells_194Pt,Dead_cells_198Pt, n_bc = 4, q = 0.95))
#' m <- c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2","pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8","pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")
#' out <- RUCova::rucova(data, markers = m, SUCs = c("mean_DNA", "mean_BC", "total_ERK", "pan_Akt"), apply_asinh_SUCs = TRUE, col_name_sample = "line", 
#' center_SUCs = "across_samples", model = "interaction", keep_offset = TRUE)
#' data_reg <- out$data_reg
#' lower <-  data |> filter(line == "Cal33") |> 
#' dplyr::mutate_at(vars(m,x), asinh) |> 
#' dplyr::select(m,x) |> 
#' cor(method= "pearson")
#' upper <-  data_reg |> filter(line == "Cal33") |> #regressed data set
#' dplyr::mutate_at(vars(m,x), asinh) |> 
#' dplyr::select(m,x) |> 
#' cor(method= "pearson")
#' heatmap_compare_corr(lower,upper)
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

