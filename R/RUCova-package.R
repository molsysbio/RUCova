#' RUCova: Removing Unwanted Covariance in Mass Cytometry Data
#'
#' The RUCova package provides tools for removing unwanted covariance in mass cytometry (CyTOF) data.
#' It is designed to help researchers preprocess and analyze CyTOF data by identifying and correcting
#' for technical and biological sources of unwanted variation.
#'
#' @section Key Features:
#' \itemize{
#'   \item Correction of technical and biological variation.
#'   \item Visualization tools for assessing data quality.
#' }
#'
#' @section Getting Started:
#' To get started with RUCova, load the package and explore the example dataset:
#' \preformatted{
#' library(RUCova)
#' data(HNSCC_data)  # Example dataset
#' head(HNSCC_data)
#' }
#'
#' For a detailed workflow, see the vignette:
#' \preformatted{
#' vignette("RUCova")
#' }
#'
#' @section Author(s):
#' \itemize{
#'   \item Rosario Astaburuaga-García (maintainer)
#' }
#'
#' @section References:
#' For more information on the methods used in this package, see:
#' \itemize{
#'   \item Astaburuaga-García, R. et al. (2023). "RUCova: A tool for removing unwanted covariance in CyTOF data." \emph{Bioinformatics}.
#' }
#'
#' @keywords internal
"_PACKAGE"