# RUC function
#library(dplyr)
#library(fastDummies)
#library(tidyr)
#library(stringr)
#library(tibble)


#' Remove unwanted covariance
#'
#' @param sce A SingleCellExperiment object with markers and SUCs in linear scale stored in the assay "name_assay_before". Asinh transformation is applied within the function.
#' @param name_assay_before A string specifying the name of the assay before RUCova (with original counts in linear scale).
#' @param markers Vector of marker names to normalise, y (in linear scale).
#' @param SUCs Vector of surrogates of unwanted covariance to use for normalisation, x (in linear scale).
#' @param name_reduced_dim string specifying the name of the dimensionality reduction result in the SingleCellExperiment sce.
#' @param apply_asinh_SUCs Apply (TRUE) or not (FALSE) asinh transformation to the SUCs. TRUE if SUCs are the measured surrogates, FALSE if SUCs are PCs.
#' @param model A character: "simple", "offset" or "interaction" defining the model.
#' @param col_name_sample A character indicating the column name in "data" defining each sample.
#' @param center_SUCs A character "across_samples" or "per_sample" defining how to center the SUCs in zero.
#' @param keep_offset Keep (TRUE) or not (FALSE) the offset intercept between samples.+
#' @param name_assay_after A string specifying the name of the assay after RUCova (with regressed counts in linear scale).
#' @return The input SingleCellExperiment object with an additional assay (name_assay_after) and a list in the metadata containing all the model details. 
#' @import dplyr
#' @import fastDummies
#' @import tidyr
#' @import stringr
#' @import tibble
#' @import SingleCellExperiment
#' @import magrittr
#' @examples
#' library(SingleCellExperiment)
#' sce <- RUCova::sce
#' bc_channels <- c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di", 
#' "Dead_cells_194Pt", "Dead_cells_198Pt")
#' sce <- RUCova::calc_mean_BC(sce, name_assay = "counts", bc_channels, n_bc = 4, q = 0.95)
#' dna_channels <- c("DNA_191Ir", "DNA_193Ir")
#' sce <- RUCova::calc_mean_DNA(sce, name_assay = "counts", dna_channels, q = 0.95)
#' # Markers:
#' m <- c("pH3","IdU","Cyclin_D1","Cyclin_B1", "Ki.67","pRb","pH2A.X","p.p53","p.p38","pChk2",
#' "pCDC25c","cCasp3","cPARP","pAkt","pAkt_T308","pMEK1.2","pERK1.2","pS6","p4e.BP1","pSmad1.8",
#' "pSmad2.3","pNFkB","IkBa", "CXCL1","Lamin_B1", "pStat1","pStat3", "YAP","NICD")
#' # SUCs::
#' x <- c("total_ERK", "pan_Akt", "mean_DNA", "mean_BC")
#' sce <- RUCova::rucova(sce = sce, name_assay_before = "counts", markers = m, SUCs = x, 
#' apply_asinh_SUCs = TRUE,  model = "interaction", center_SUCs = "across_samples", 
#' col_name_sample = "line", name_assay_after = "counts_interaction")
#' @export
rucova <- function(sce, name_assay_before = "counts",  markers, SUCs = c("mean_DNA", "mean_BC", "total_ERK", "pan_Akt"), name_reduced_dim = "PCA", apply_asinh_SUCs = TRUE, model = "interaction", col_name_sample = "line",
                               center_SUCs = "across_samples", keep_offset = TRUE, name_assay_after = "counts_rucova") {

 
   if (missing(sce) == TRUE || !inherits(sce, "SingleCellExperiment")){
     stop("Please provide a SingleCellExperiment class")
     
   }
    data <- t(assay(sce,name_assay_before)) |> cbind(colData(sce)) |> as.data.frame()
    
    if (grepl("PC",SUCs)[1]){# if model is based on PCs, add this info to the data
      data <- data |> cbind(reducedDim(sce, name_reduced_dim))
      } 
  
    # Type of model ------------------------------------
  if (model == "offset" || model == "interaction" || center_SUCs == "per_sample") {
    if (missing(col_name_sample) == TRUE){
      stop("Please specify argument `col_name_sample`")

    }
          if(apply_asinh_SUCs == TRUE) {
            dt <- data |>
              dplyr::rename(sample = col_name_sample) |>
              mutate(across(all_of(c(markers, SUCs)), asinh))
          } else {
            dt <- data |>
              dplyr::rename(sample = col_name_sample) |>
              mutate(across(all_of(markers), asinh))
          }


     # dt$sample <- droplevels(dt$sample)

  } else {
          if(apply_asinh_SUCs == TRUE) {

          dt <- data |>
            mutate(across(all_of(c(markers, SUCs)), asinh))

          } else {

            dt <- data |>
              mutate(across(all_of(markers), asinh))
          }
  }

  if (center_SUCs == "per_sample") {
    dt <- dt |>
      group_by(sample) |>
      mutate(across(all_of(SUCs), ~ .x - mean(.x))) |>
      ungroup()
  } else if (center_SUCs == "across_samples") {
    dt <- dt |>
      mutate(across(all_of(SUCs), ~ .x - mean(.x)))
  } else {
    stop("Please specify argument 'center_SUCs'")
  }

  # Add dummy variables if necessary  ------------------------------------
  if (model == "interaction" || model == "offset") {
    dt <- dummy_cols(.data = dt,
                     select_columns = "sample",
                     remove_first_dummy = TRUE)

    dummy_sample_var <- colnames(dt)[grepl("sample_", colnames(dt))]
    dummy_values <- dt |>
      group_by(sample) |>
      summarise(across(all_of(dummy_sample_var), max))
  }

  #dt <- data |> mutate(cell_id = 1:n())
  
  # Model function and coefficients  ------------------------------------
  if (model == "interaction") {
    slope_dummy <- levels(interaction(SUCs, dummy_sample_var, sep = " : "))
    n_coeff <- 1 + length(dummy_sample_var) + length(SUCs) + (length(SUCs) * length(dummy_sample_var))
  } else if (model == "offset") {
    slope_dummy <- NULL
    n_coeff <- 1 + length(dummy_sample_var) + length(SUCs)
  } else {
    dummy_sample_var <- NULL
    slope_dummy <- NULL
    n_coeff <- 1 + length(SUCs)
  }

  # Regression ------------------------------------
  fits <- lapply(markers, function(marker_to_fit) {
    print(paste0("Fitting ", marker_to_fit))
    formula <- reformulate(termlabels = c(dummy_sample_var, SUCs, slope_dummy),
                           response = marker_to_fit)
    lm(formula = formula, data = dt)
    }) |>
    setNames(markers)

  
  model_coefficients.new <- sapply(fits, coef) |> t() |>
    as_tibble(rownames = "marker")

  model_residuals.new <- sapply(fits, resid) |>
    as_tibble() |>
    mutate(cell_id = dt$cell_id, .before = 1)

  adjr2.new <- vapply(fits, function(fit) {
    SSres <- sum(residuals(fit)^2)
    SStot <- sum((fit$model[[1]] - mean(fit$model[[1]]))^2)
    fit_rsquared <- 1 - (SSres / SStot)
    n <- nrow(fit$model)
    p <- length(coef(fit)) - 1
    1 - (1 - fit_rsquared) * (n - 1) / (n - p - 1)
  }, FUN.VALUE = numeric(1)) |>
    tibble::enframe(name = "marker", value = "adj_r_squared")

  # For output
  model_formula <- reformulate(
    termlabels = c(dummy_sample_var, SUCs, slope_dummy),
    response = "y")

  # Regressed values ------------------------------------
  new_values <- model_residuals.new

  for (m in (as.vector(markers))) {
    new_values[m] <-
      as.numeric(as.vector(unlist(model_residuals.new[m]))) + # residuals
      as.numeric(model_coefficients.new[model_coefficients.new$marker == m, 2])   # intercept for all

    if (keep_offset == TRUE) {
      for (i in dummy_sample_var) {
        new_values[m] <- pull(new_values,m) +
          as.numeric(model_coefficients.new[model_coefficients.new$marker == m, i]) * pull(dt, i) # offset_model
      }
    }
  }

  data_reg <-
    data |>
    select(-all_of(markers)) |>
    left_join(new_values, by = "cell_id") |>
    select(colnames(data)) |> # same order
    mutate(across(all_of(markers), sinh))

  # Effective coefficients ------------------------------------
  eff_coefficients <- model_coefficients.new
  if (model == "interaction" || model == "offset") {
    baseline_sample <-  as.character(dummy_values[[1]][1])
  }

  eff_coefficients <- eff_coefficients |>
      pivot_longer(names_to = "coef_key", -marker) |>
      mutate(surrogate = ifelse(str_detect(coef_key, paste(SUCs, collapse = "|")),
                                str_extract(coef_key, paste(SUCs, collapse = "|")), as.logical(FALSE)))

    if (model == "interaction") {
      eff_coefficients <- eff_coefficients |>
        mutate(sample = str_remove(coef_key, "sample_"),
               sample = str_remove(sample, paste0(":", surrogate)),
               sample = ifelse(sample %in% str_remove(dummy_sample_var, "sample_"), sample, baseline_sample)) |>
        group_by(surrogate, marker) |>
        mutate(eff_value = ifelse(sample == baseline_sample,
                                  value,
                                  value + value[sample == baseline_sample])) |>
        ungroup()
    } else if (model == "offset") {
      eff_coefficients <- eff_coefficients |>
        mutate(sample = str_remove(coef_key, "sample_"),
               sample = ifelse(sample %in% str_remove(dummy_sample_var, "sample_") &
                                 surrogate == "FALSE", sample, baseline_sample),
               sample = ifelse(surrogate %in% SUCs, "all", sample)) |>
        group_by(surrogate, marker) |>
        mutate(eff_value = ifelse(sample == baseline_sample | surrogate %in% SUCs,
                                  value,
                                  value + value[sample == baseline_sample])) |>
        ungroup()
    } else {
      eff_coefficients <- eff_coefficients |>
        mutate(eff_value = value)
    }

  # Standardized slopes (effect size) ------------------------------------
    if (model == "simple" || model == "offset") { # 1 slope across all cells
      sd_values <- dt |>
        summarise(across(all_of(c(markers, SUCs)), sd)) |>
        ungroup() |>
        pivot_longer(names_to = "marker", values_to = "sd_y", markers) |>
        pivot_longer(names_to = "surrogate", values_to = "sd_x", SUCs)

      stand_slopes <- eff_coefficients |>
        filter(surrogate != FALSE) |>
        left_join(sd_values, by = c("marker", "surrogate")) |>
        mutate(stand_value = eff_value * sd_x / sd_y)
    } else { # 1 slope per sample
      sd_values <- dt |>
        group_by(sample) |>
        summarise_at(vars(markers, SUCs), sd) |>
        ungroup() |>
        pivot_longer(names_to = "marker", values_to = "sd_y", markers) |>
        pivot_longer(names_to = "surrogate", values_to = "sd_x", SUCs)

      stand_slopes <- eff_coefficients |>
        filter(surrogate != FALSE) |>
        left_join(sd_values, by = c("marker", "surrogate", "sample")) |>
        mutate(stand_value = eff_value * sd_x / sd_y)
    }

    
      assay(sce, name_assay_after) <- data_reg |> select(rownames(sce)) |> t()
  
      out_ruc <- list(name_assay_before,markers, SUCs, name_reduced_dim, apply_asinh_SUCs, model,col_name_sample,center_SUCs, keep_offset, name_assay_after, model_formula, model_coefficients.new, eff_coefficients, model_residuals.new, adjr2.new, stand_slopes)
      names(out_ruc) <- c("name_assay_before", "markers", "SUCs", "name_reduced_dim","apply_asinh_SUCs", "model", "col_name_sample", "center_SUCs", "keep_offset", "name_assay_after","model_formula", "model_coefficients","eff_coefficients", "model_residuals", "adjr2", "stand_slopes")
      
      metadata(sce)[[paste0("model_",name_assay_after)]] <- out_ruc
      
  return(sce)
}
  