
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(Matrix)


#' Calculated mean of normalised Iridium isotopes
#'
#' @param ... Channels to average in linear scale. Asinh transformation is applied within the function.
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples 
#' data |> mutate(mean_DNA = RUCova::calc_mean_DNA(DNA_191Ir, DNA_193Ir, q = 0.95)))
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
#' @param ... Channels to average in linear scale. Asinh transformation is applied within the function.
#' @param q Quantile for normalisation.
#' @param n_bc number of barcoding isotopes per cell.
#' @return A vector.
#' @examples
#' data |> mutate(mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di, n_bc = 4, q = 0.95))
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

#' Plot pearson correlation coefficients between markers on a double triangular heatmap (lower triangle: before RUCova, upper triangle: after RUCova).
#' @param corr_values_before Matrix with pearson correlation coefficient between markers before RUCova.
#' @param corr_values_after  Matrix with pearson correlation coefficient between markers after RUCova.
#' @return #A heatmap
#' @examples 
#' plot_corr_before_after(corr_values_before,corr_values_after)
#'
#' @export
#'
#'
plot_corr_before_after <- function(corr_values_before, corr_values_after){

tmp <-  as.matrix(tril(corr_values_before,-1) + 
                    triu(corr_values_after))

diag(tmp) <- NA
hm <- tmp|> 
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

