# fair_zscore <- function(x,vars,group_var){
#
#   x <- x %>%
#     mutate(across(vars,asinh))
#
#   mean_sd <- x %>%
#     group_by_(group_var) %>%
#     summarise(across(vars,sd)) %>% ungroup() %>%
#     summarise(across(vars,mean,.names = '{col}_mean_sd'))
#
#
#   x %>%
#     bind_cols(mean_sd) %>%
#     group_by_(group_var) %>%
#     mutate(across(.cols = vars, .fns = ~ (.x - mean(.x))/mean_sd(.x,)))
#
#
# }

# calc_mean_ir <- function(data, iridium_channels, q){
#
#
#   perc <- data %>%
#     select(iridium_channels) %>%
#     mutate_all(asinh) %>%
#     pivot_longer(names_to = "iridium_channels", values_to = "values", everything()) %>%
#     group_by(iridium_channels) %>%
#     summarise(percentile = quantile(values, probs = q)) %>%
#     ungroup() %>%
#     mutate(sf = percentile[iridium_channels == iridium_channels[1]]/percentile) %>% #first is the reference
#     select(iridium_channels, sf)
#
#
#   mean_ir <- data %>%
#     select(id,iridium_channels) %>%
#     mutate_at(vars(iridium_channels), asinh) %>%
#     pivot_longer(names_to = "iridium_channels", values_to = "value", -id) %>%
#     left_join(perc, by = "iridium_channels") %>%
#     mutate(value = value*sf) %>%
#     group_by(id) %>%
#     summarize(mean_ir = mean(value)) %>%
#     ungroup() %>%
#     mutate(mean_ir = sinh(mean_ir))
#
#   return(mean_ir)
# }

#' Calculated mean of normalised Iridium isotopes
#'
#' @param ... Channels to average.
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples
#'
#' @export
#'
calc_mean_DNA <- function(..., q) {
  asinh_input <- asinh(rbind(...))

  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = F)
  scaling_factors <- quantiles[[1]] / quantiles

  scaled_data <- asinh_input * scaling_factors
  means <- apply(scaled_data, 2, mean)

  sinh(means)
}

#' Calculated mean of normalised highest BC per cell
#'
#' @param ... Channels to average.
#' @param q Quantile for normalisation.
#' @param n_bc number of barcoding isotopes per cell
#' @return A vector.
#' @examples
#'
#' @export
#'
calc_mean_BC <- function(..., n_bc, q) {
  asinh_input <- asinh(rbind(...))

  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = F)
  scaling_factors <- quantiles[[1]] / quantiles

  scaled_data <- asinh_input * scaling_factors

  scaled_data |>
    apply(2, sort.int, decreasing = TRUE, method = "shell") |>
    #Rfast::colSort(descending = TRUE) |>
    head(n_bc) |>
    colMeans() |>
    sinh()
}



