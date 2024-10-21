#' Calculated mean of normalised Iridium isotopes
#'
#' @param ... Channels to average in linear scale. Asinh transformation is applied within the function.
#' @param q Quantile for normalisation.
#' @return A vector.
#' @examples 
#' data <- RUCova::HNSCC_data
#' data |> dplyr::mutate(mean_DNA = RUCova::calc_mean_DNA(DNA_191Ir, DNA_193Ir, q = 0.95))
#'
#' @export
#' @import dplyr
#' @import ComplexHeatmap
#' @import circlize
#' @import Matrix
#' @import grid

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
#' data <- RUCova::HNSCC_data
#' data |> dplyr::mutate(mean_BC = RUCova::calc_mean_BC(Pd102Di, Pd104Di, Pd105Di, Pd106Di, Pd108Di, Pd110Di, n_bc = 4, q = 0.95))
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
#' @param lower Matrix with pearson correlation coefficient between markers eg.: before RUCova, to be plotted in the lower triangle.
#' @param upper  Matrix with pearson correlation coefficient between markers eg.: after RUCova, to be plotted in the upper triangle.
#' @return #A heatmap
#' @examples 
#' library(dplyr)
#' lower <-  data |> 
#' mutate_at(vars(m,x), asinh) |> 
#' select(m,x) |> 
#' cor(method= "pearson")
#' upper <-  data_reg |> #regressed data set
#' mutate_at(vars(m,x), asinh) |> 
#' select(m,x) |> 
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

