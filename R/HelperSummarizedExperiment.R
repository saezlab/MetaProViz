#' Process SummarizedExperiment object to extract assay, sample, and feature metadata
#'
#' This helper function extracts the numeric data matrix, sample metadata, and feature metadata from a SummarizedExperiment object.
#' If the input is not a SummarizedExperiment, an error is thrown.
#'
#' @param se_obj A SummarizedExperiment object.
#'
#' @return A list with elements:
#'   \item{data}{Numeric matrix with metabolites as columns and samples as rows}
#'   \item{metadata_sample}{Data frame with sample metadata (samples as row names)}
#'   \item{metadata_feature}{Data frame with feature metadata (metabolites as row names)}
#'
#' @examples
#' Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
#' mat <- Intra[-c(49:58), -c(1:5, 97)]
#' colData <- cellular_meta %>% dplyr::filter(Metabolite %in% colnames(mat)) %>% tibble::column_to_rownames("Metabolite")%>% .[colnames(mat), , drop = FALSE]
#' rowData <- Intra[-c(49:58), c(1:3)]
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat), colData = colData, rowData = rowData)
#' result <- process_se(se)
#' str(result)
#'
process_se <- function(se_obj) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package is required.")
  }
  if (!inherits(se_obj, "SummarizedExperiment"))  {
    stop("Provided object is not a SummarizedExperiment.")
  }
  assay_data <- SummarizedExperiment::assay(se_obj)
  metadata_sample <- SummarizedExperiment::rowData(se_obj)
  metadata_feature <- SummarizedExperiment::colData(se_obj)
  return(list(
    data = assay_data,
    metadata_sample = as.data.frame(metadata_sample),
    metadata_feature = as.data.frame(metadata_feature)
  ))
}
