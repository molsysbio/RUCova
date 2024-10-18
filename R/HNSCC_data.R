#' HNSCC data set
#'
#' Is a tibble containing mass cytometry data of single-cell marker signals (rows = cells, columns = markers and metadata) in linear scale. 
#' This data set should be clean, meaning you excluded beads, debris, doublets, dead cells, and single-cells are demultiplexed (important if you want to adapt the linear fits to the samples).
#' In this example we offer a mass cytometry data set consisting of 8 Head-and-Neck Squamous Cell Carcinoma (HNSCC) lines in irradiated (10 Gy) and control (0 Gy) conditions (Figure 2 and Figure 3 in the manuscript).
#'
#' @docType data
#' @usage data(HNSCC_data)
#' @format A data frame with 108649 rows and 59 variables:
#' \describe{
#'   \item{dose}{0Gy and 10Gy}
#'   \item{line}{"Cal27","Cal33","UPCISCC099","UPCISCC131","UTSCC16A","UDSCC2","UPCISCC154","VUSCC147"  }
#'   ...
#' }
"HNSCC_data"