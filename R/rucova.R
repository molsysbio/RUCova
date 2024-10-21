# RUC function
library(dplyr)
library(fastDummies)
library(tidyr)
library(stringr)
library(tibble)


#' Remove unwanted covariance
#'
#' @param data A tibble with markers and SUCs in linear scale. Asinh transformation is applied within the function.
#' @param markers Vector of marker names to normalise, y (in linear scale).
#' @param SUCs Vector of surrogates of unwanted covariance to use for normalisation, x (in linear scale).
#' @param apply_asinh_SUCs Apply (TRUE) or not (FALSE) asinh transformation to the SUCs. TRUE if SUCs are the measured surrogates, FALSE if SUCs are PCs.
#' @param model A character: "simple", "offset" or "interaction" defining the model.
#' @param col_name_sample A character indicating the column name in "data" defining each sample.
#' @param center_SUCs A character "across_samples" or "per_sample" defining how to center the SUCs in zero.
#' @param keep_offset Keep (TRUE) or not (FALSE) the offset intercept between samples.
#' @return Normalised tibble with marker and surrogate values in linear scale (as the input).
#' @examples 
#' dontrun{
#' RUCova::rucova(data, markers_to_norm,SUCs = c("PC1","PC2","PC3","PC4"), apply_asinh_SUCs = FALSE, col_name_sample = "line", 
#' center_SUCs = "across_samples", model = "interaction", keep_offset = TRUE)
#' RUCova::rucova(data, markers_to_norm,SUCs = c("mean_DNA", "mean_BC", "total_ERK", "pan_Akt"), apply_asinh_SUCs = TRUE, col_name_sample = "line", 
#' center_SUCs = "across_samples", model = "interaction", keep_offset = TRUE)
#' }
#' @import dplyr
#' @import fastDummies
#' @import tidyr
#' @import stringr
#' @import tibble
#' @export
rucova <- function(data, markers, SUCs, apply_asinh_SUCs, model = "interaction", col_name_sample = "line",
                               center_SUCs = "across_samples", keep_offset = TRUE) {

  # model = c("simple","offset","interaction"), interaction = slope+offset
  # keep_offset = TRUE or FALSE

  # Type of model ------------------------------------
  if (model == "offset" || model == "interaction" || center_SUCs == "per_sample") {
    if (missing(col_name_sample)){
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

    if (isTRUE(keep_offset)) {
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

  # Output ------------------------------------
  
  out_ruc <- list(data_reg, markers, SUCs, apply_asinh_SUCs, model,col_name_sample,center_SUCs, keep_offset, col_name_sample, model_formula, model_coefficients.new,eff_coefficients, model_residuals.new, adjr2.new, stand_slopes)
  names(out_ruc) <- c("data_reg", "markers", "SUCs", "apply_asinh_SUCs", "model", "col_name_sample", "center_SUCs", "keep_offset", "model_formula", "model_coefficients","eff_coefficients", "model_residuals", "adjr2", "stand_slopes")
  return(out_ruc)
}
