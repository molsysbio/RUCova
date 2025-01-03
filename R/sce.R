#' SingleCellExperiment Object with HNSCC Data Set
#'
#' A `SingleCellExperiment` object containing mass cytometry data of single-cell marker signals.
#' The data includes one or more assays (e.g., "counts") with signals in linear scale. 
#' Rows represent markers, and columns represent cells. The data is clean, 
#' excluding calibration beads, debris, doublets, and dead cells. Single cells are demultiplexed, 
#' which is important for adapting linear fits to samples. Metadata such as samples and treatment conditions 
#' are stored in the `colData`.
#'
#' @docType data
#' @name sce
#' @usage SummarizedExperiment::assay(sce, "counts")
#' @format A `SingleCellExperiment` object with the following components:
#' \describe{
#'   \item{assays}{One or more assays, such as \code{"counts"}, containing the marker signals.}
#'   \item{colData}{Column metadata, including cell annotations such as \code{cell_id}, \code{line}, and \code{dose}.}
#'   \item{rowData}{Row metadata, including marker annotations.}
#' }
#' @param sce A `SingleCellExperiment` object containing the HNSCC dataset.
#' @param counts The assay name for marker signal data (e.g., "counts").
#' @aliases sce assay
"sce"