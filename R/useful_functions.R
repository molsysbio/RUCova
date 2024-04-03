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

#' Calculated scaled mean of Iridium isotopes
#'
#' @param ... Channels to average.
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples
#'
#' @export
#'
calc_mean_DNA <- function(..., p) {
  asinh_input <- asinh(rbind(...))

  quantiles <- apply(asinh_input, 1, quantile, probs = p, names = F)
  scaling_factors <- quantiles[[1]] / quantiles

  scaled_data <- asinh_input * scaling_factors
  means <- apply(scaled_data, 2, mean)

  sinh(means)
}

#' Calculated scaled mean of Iridium isotopes
#'
#' @param ... Channels to average.
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples
#'
#' @export
#'
calc_mean_highest_bc <- function(..., bc_count, q) {
  asinh_input <- asinh(rbind(...))

  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = F)
  scaling_factors <- quantiles[[1]] / quantiles

  scaled_data <- asinh_input * scaling_factors

  scaled_data |>
    apply(2, sort.int, decreasing = TRUE, method = "shell") |>
    #Rfast::colSort(descending = TRUE) |>
    head(bc_count) |>
    colMeans() |>
    sinh()
}

identify_bc_isotope <- function(..., bc_count, q) {
  asinh_input <- asinh(rbind(...))

  quantiles <- apply(asinh_input, 1, quantile, probs = q, names = F)
  scaling_factors <- quantiles[[1]] / quantiles

  scaled_data <- asinh_input * scaling_factors

  scaled_data |>
    apply(2, sort.int, decreasing = TRUE, method = "shell") |>
    #Rfast::colSort(descending = TRUE) |>
    head(bc_count)
}

calc_pc1_bc <- function(...) {
  asinh_input <- asinh(cbind(...))
  pca_res <- asinh_input |>
    as.data.frame() |>
    mutate_all(scale) |>
    prcomp()
  -pca_res$x[,1] |>
    sinh()
}


# hella fast, but different result

calc_mean_used_bc <- function(data, bc_channels, bc_col, n_used_bc, q){

  perc <- data %>%
    pivot_longer(all_of(bc_channels),
                 names_to = "bc_ch",
                 values_transform = asinh) %>%
    select(bc_ch, value) %>%
    group_by(bc_ch) %>%
    summarise(percentile = quantile(value, probs = q)) %>%
    ungroup() %>%
    mutate(sf = percentile[bc_ch == bc_ch[1]] / percentile) %>% #first is the reference
    select(bc_ch, sf)


  usedBC <- data %>%
    pivot_longer(all_of(bc_channels),
                 names_to  = "bc_ch",
                 values_transform = asinh) %>%
    select(bc_col, bc_ch, value) %>%
    left_join(perc, by = "bc_ch") %>% # percentile normalization
    mutate(value = value * sf) %>%
    select(-sf) %>%
    group_by_(bc_col, "bc_ch") %>%
    summarise(value = mean(value)) %>%
    group_by_(bc_col) %>%
    mutate(ranks = rank(-value)) %>%
    mutate(use = ifelse(ranks %in% c(1:n_used_bc),1,0)) %>%
    ungroup() %>%
    select(all_of(bc_col), bc_ch, use)

  mean_used_bc <- data %>%
    select(id, all_of(bc_col), bc_channels) %>%
    mutate_at(vars(bc_channels), asinh) %>%
    pivot_longer(names_to = "bc_ch", values_to  ="value", -c({{ bc_col }}, id)) %>%
    left_join(usedBC, by = c(as.character(bc_col), "bc_ch")) %>%
    left_join(perc, by = "bc_ch") %>% # percentile normalization
    mutate(value = value * sf) %>%
    select(-sf) %>%
    filter(use == 1) %>% # we only care about used barcodes
    group_by(id) %>%
    summarise(mean_used_bc = mean(value)) %>%
    ungroup() %>%
    mutate(mean_used_bc = sinh(mean_used_bc))


  return(mean_used_bc)

}
