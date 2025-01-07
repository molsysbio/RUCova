#' Calculated mean of normalised Iridium isotopes
#' @param sce A SingleCellExperiment object with markers and SUCs in linear scale stored in the assay "name_assay". Asinh transformation is applied within the function.
#' @param name_assay A string specifying the name of the assay including the DNA channels in linear scale.
#' @param dna_channels Vector specifying the names of the DNA channels
#' @param q Quantile for normalisation.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import dplyr
#' @return The SingleCellExperiment object with an extra column "mean_BC" in the corresponding assay.
#' @examples 
#' sce <- RUCova::sce
#' dna_channels <- c("DNA_191Ir", "DNA_193Ir")
#' sce <- RUCova::calc_mean_DNA(sce, name_assay = "counts", dna_channels, q = 0.95)
#' @export

calc_mean_DNA <- function(sce, name_assay = "counts", dna_channels, q) {
  
  # Check if the input is a SingleCellExperiment
  if (inherits(sce, "SingleCellExperiment")){
    dna_data <- assay(sce,name_assay)[dna_channels, , drop = FALSE]
  } else {
    # Handle input as a data frame or matrix
    data <- sce
    dna_data <- t(data[,dna_channels])
    
  }
  # Perform calculations
  asinh_input <- asinh(dna_data)
  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = FALSE)
  scaling_factors <- quantiles[[1]] / quantiles
  scaled_data <- asinh_input * scaling_factors
  mean_DNA <- sinh(apply(scaled_data, 2, mean))
  
  # Add the mean_DNA to the appropriate structure
  if(inherits(sce, "SingleCellExperiment")){
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
#' @param sce A SingleCellExperiment object with markers and SUCs in linear scale stored in the assay "name_assay". Asinh transformation is applied within the function.
#' @param name_assay A string specifying the name of the assay including the BC channels in linear scale. Default is "counts".
#' @param bc_channels Vector specifying the names of the BC channels
#' @param q Quantile for normalisation. Default is 0.95.
#' @param n_bc number of barcoding isotopes per cell. n_bc = 3 for the Fluidigm kit.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import dplyr
#' @return The SingleCellExperiment object with an extra column "mean_BC" in the corresponding assay.
#' @examples 
#' sce <- RUCova::sce
#' bc_channels <- c(c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di"),
#'   c("Dead_cells_194Pt", "Dead_cells_198Pt")
#' )
#' sce <- RUCova::calc_mean_BC(sce, name_assay = "counts", bc_channels, n_bc = 4, q = 0.95)
#' @export
#'

calc_mean_BC <- function(sce, name_assay = "counts", bc_channels, n_bc, q = 0.95) {
  
  # Check if the input is a SingleCellExperiment
  if (inherits(sce, "SingleCellExperiment")){
    bc_data <- assay(sce)[bc_channels, , drop = FALSE]
  } else {
    # Handle input as a data frame or matrix
    data <- sce
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
  if(inherits(sce, "SingleCellExperiment")){
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

#' Plot pearson correlation coefficients between markers on a double triangular heatmap (lower triangle: before RUCova, upper triangle: after RUCova). If RUCova has not been applied, the output is a symmetric heatmap.
#' @param sce A SingleCellExperiment object with markers and SUCs in linear scale stored in the assay "name_assay". Asinh transformation is applied within the function.
#' @param name_assay_before A string specifying the name of the assay before RUCova (with original counts in linear scale).
#' @param name_assay_after A string specifying the name of the assay before RUCova (with original counts in linear scale).
#' @param name_reduced_dim A string specifying the name of the dimensionality reduction data stored under ``reducedDim()``. If "PCA", then PCs will be included in the heatmao. 
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @import grid
#' @import circlize
#' @import ComplexHeatmap
#' @import tidyverse
#' @import tidyr
#' @return #A heatmap with pearson correlation coefficients.
#' @examples
#' sce <- RUCova::sce
#' bc_channels <- c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di",
#' "Dead_cells_194Pt", "Dead_cells_198Pt")
#' sce <- RUCova::calc_mean_BC(sce, name_assay = "counts", bc_channels, n_bc = 4, q = 0.95)
#' dna_channels <- c("DNA_191Ir", "DNA_193Ir")
#' sce <- RUCova::calc_mean_DNA(sce, name_assay = "counts", dna_channels, q = 0.95)
#' # Markers:
#' m <- c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38",
#' "pChk2","pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1",
#' "pSmad1.8","pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")
#' # SUCs::
#' x <- c("total_ERK", "pan_Akt", "mean_DNA", "mean_BC")
#' sce <- RUCova::rucova(sce = sce, name_assay_before = "counts", markers = m, SUCs = x, 
#' apply_asinh_SUCs = TRUE,  model = "interaction", center_SUCs = "across_samples", 
#' col_name_sample = "line", name_assay_after = "counts_interaction")
#' heatmap_compare_corr(sce, name_assay_before = "counts", name_assay_after = "counts_interaction")
#' @export
#'
#'
heatmap_compare_corr <- function(sce, name_assay_before = "counts", name_assay_after = NULL, name_reduced_dim = NULL){
  
  #### before: no models pars needed
  data_before <- t(assay(sce,name_assay_before)) |>  as.tibble()
  

  ## if it is an initial evaluation with no assay after RUCova:
  if(is.null(name_assay_after) || name_assay_before == name_assay_after){
    data_before <- data_before |> mutate_all(asinh)
    
    ## add PCA
    if (!is.null(name_reduced_dim)){
      data_before <- data_before |> cbind(reducedDim(sce, type = name_reduced_dim))
    }
    
    corr_before <-  data_before |>  cor(method= "pearson")
    corr_after <-   corr_before
    
  } else { #if RUCova was already applied:
    
    data_before <- data_before |>  cbind(colData(sce)) 
    data_after <-  t(assay(sce,name_assay_after)) |>   as.tibble() |> cbind(colData(sce)) 

    ## add PCA
    if (!is.null(name_reduced_dim)){
      data_before <- data_before |> cbind(reducedDim(sce, type = name_reduced_dim))
      data_after <- data_after |> cbind(reducedDim(sce, type = name_reduced_dim))
    }
    
    #transform the data as done for the model:
    
    markers <- sce@metadata[[paste0("model_", name_assay_after)]]$markers
    SUCs <-  sce@metadata[[paste0("model_", name_assay_after)]]$SUCs
    sample <- sce@metadata[[paste0("model_", name_assay_after)]]$col_name_sample
    
    #asinh markers and apply_asinh_SUCs?
    
    if(sce@metadata[[paste0("model_", name_assay_after)]]$apply_asinh_SUCs){ 
      data_before <- data_before |> mutate_at(vars(markers,SUCs), asinh) 
      data_after <- data_after |> mutate_at(vars(markers,SUCs), asinh)
    } else {
      data_before <- data_before |> mutate_at(vars(markers), asinh) 
      data_after <- data_after |> mutate_at(vars(markers), asinh)
    }
    
    #center_SUCs

    if(sce@metadata[[paste0("model_", name_assay_after)]]$center_SUCs == "per_sample"){
      data_before <- data_before |> group_by(!!sym(sample)) |> mutate(across(all_of(SUCs), ~ .x - mean(.x))) |> ungroup()
      data_after <- data_after |>  group_by(!!sym(sample)) |> mutate(across(all_of(SUCs), ~ .x - mean(.x))) |> ungroup()
    } else {
      data_before <- data_before |>  mutate(across(all_of(SUCs), ~ .x - mean(.x)))
      data_after <- data_after |> mutate(across(all_of(SUCs), ~ .x - mean(.x)))
    }
    
    if(!is.null(name_reduced_dim)){
      DRs <- colnames(reducedDim(sce, type = name_reduced_dim))
      corr_before <-  data_before |> select(markers,SUCs, DRs) |>  cor(method= "pearson")
      corr_after <-   data_after |> select(markers,SUCs, DRs) |> cor(method= "pearson")
    } else {
      corr_before <-  data_before |> select(markers,SUCs) |>  cor(method= "pearson")
      corr_after <-   data_after |> select(markers,SUCs) |> cor(method= "pearson")
    }
    
    
  } 
  
  
  tmp <-  as.matrix(Matrix::tril(as.matrix(corr_before)) + Matrix::triu(as.matrix(corr_after)))
  
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
                            row_title_side = "left",
                            column_title = "after RUCova",
                            row_title = "before RUCova")
  
  draw(hm)
  
}
