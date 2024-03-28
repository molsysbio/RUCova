# RUC function
library(dplyr)
library(fastDummies)
library(tidyr)
library(stringr)
library(tibble)


#' Regress surrogates of unwanted covariance
#'
#' @param data A tibble.
#' @param markers Vector of marker values to normalise.
#' @param surrogates Vector of surrogates to use for normalisation.
#' @return Normalised tibble.
#' @examples
#'
#' @import dplyr
#' @import fastDummies
#' @import tidyr
#' @import stringr
#' @import tibble
#' @export
rucova <- function(data, markers, surrogates, apply_asinh_SUC, col_name_sample = NULL,
                               center_surr = "per_sample", model = "simple", keep_offset = TRUE) {

  # model = c("simple","offset","interaction"), interaction = slope+offset
  # keep_offset = TRUE or FALSE

  # Type of model ------------------------------------
  if (model == "offset" || model == "interaction" || center_surr == "per_sample") {
    if (missing(col_name_sample)){
      stop("Please specify argument `col_name_sample`")

    }

          if(apply_asinh_SUC == TRUE) {
            dt <- data |>
              dplyr::rename(sample = col_name_sample) |>
              mutate(across(all_of(c(markers, surrogates)), asinh))
          } else {
            dt <- data |>
              dplyr::rename(sample = col_name_sample) |>
              mutate(across(all_of(markers), asinh))
          }


     # dt$sample <- droplevels(dt$sample)

  } else {
          if(apply_asinh_SUC == TRUE) {

          dt <- data |>
            mutate(across(all_of(c(markers, surrogates)), asinh))

          } else {

            dt <- data |>
              mutate(across(all_of(markers), asinh))
          }
  }

  if (center_surr == "per_sample") {
    dt <- dt |>
      group_by(sample) |>
      mutate(across(all_of(surrogates), ~ .x - mean(.x))) |>
      ungroup()
  } else if (center_surr == "across_samples") {
    dt <- dt |>
      mutate(across(all_of(surrogates), ~ .x - mean(.x)))
  } else {
    stop("Please specify argument 'center_surr'")
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

  # Model function and coefficients  ------------------------------------
  if (model == "interaction") {
    slope_dummy <- levels(interaction(surrogates, dummy_sample_var, sep = " : "))
    n_coeff <- 1 + length(dummy_sample_var) + length(surrogates) + (length(surrogates) * length(dummy_sample_var))
  } else if (model == "offset") {
    slope_dummy <- NULL
    n_coeff <- 1 + length(dummy_sample_var) + length(surrogates)
  } else {
    dummy_sample_var <- NULL
    slope_dummy <- NULL
    n_coeff <- 1 + length(surrogates)
  }

  # Regression ------------------------------------
  fits <- lapply(markers, function(marker_to_fit) {
    print(paste0("Fitting ", marker_to_fit))
    formula <- reformulate(termlabels = c(dummy_sample_var, surrogates, slope_dummy),
                           response = marker_to_fit)
    lm(formula = formula, data = dt)
    }) |>
    setNames(markers)

  model_coefficients.new <- sapply(fits, coef) |> t() |>
    as_tibble(rownames = "marker")

  model_residuals.new <- sapply(fits, resid) |>
    as_tibble() |>
    mutate(id = dt$id, .before = 1)

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
    termlabels = c(dummy_sample_var, surrogates, slope_dummy),
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
    left_join(new_values, by = "id") |>
    select(colnames(data)) |> # same order
    mutate(across(all_of(markers), sinh))

  # Effective coefficients ------------------------------------
  eff_coefficients <- model_coefficients.new
  if (model == "interaction" || model == "offset") {
    baseline_sample <-  as.character(dummy_values[[1]][1])
  }

  eff_coefficients <- eff_coefficients |>
      pivot_longer(names_to = "coef_key", -marker) |>
      mutate(surrogate = ifelse(str_detect(coef_key, paste(surrogates, collapse = "|")),
                                str_extract(coef_key, paste(surrogates, collapse = "|")), as.logical(FALSE)))

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
               sample = ifelse(surrogate %in% surrogates, "all", sample)) |>
        group_by(surrogate, marker) |>
        mutate(eff_value = ifelse(sample == baseline_sample | surrogate %in% surrogates,
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
        summarise(across(all_of(c(markers, surrogates)), sd)) |>
        ungroup() |>
        pivot_longer(names_to = "marker", values_to = "sd_y", markers) |>
        pivot_longer(names_to = "surrogate", values_to = "sd_x", surrogates)

      stand_slopes <- eff_coefficients |>
        filter(surrogate != FALSE) |>
        left_join(sd_values, by = c("marker", "surrogate")) |>
        mutate(stand_value = eff_value * sd_x / sd_y)
    } else { # 1 slope per sample
      sd_values <- dt |>
        group_by(sample) |>
        summarise_at(vars(markers, surrogates), sd) |>
        ungroup() |>
        pivot_longer(names_to = "marker", values_to = "sd_y", markers) |>
        pivot_longer(names_to = "surrogate", values_to = "sd_x", surrogates)

      stand_slopes <- eff_coefficients |>
        filter(surrogate != FALSE) |>
        left_join(sd_values, by = c("marker", "surrogate", "sample")) |>
        mutate(stand_value = eff_value * sd_x / sd_y)
    }

  # Output ------------------------------------
  out_ruc <- list(data_reg, markers, surrogates, center_surr, model, keep_offset, col_name_sample, model_formula, model_coefficients.new, model_residuals.new, adjr2.new, eff_coefficients, stand_slopes)
  names(out_ruc) <- c("data_reg", "markers", "surrogates", "center_surr", "model", "keep_offset", "col_name_sample", "model_formula", "model_coefficients", "model_residuals", "adjr2", "eff_coefficients", "stand_slopes")
  return(out_ruc)
}
