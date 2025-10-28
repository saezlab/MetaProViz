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
#' data(intracell_raw_se)
#' result <- process_se(intracell_raw_se)
#' str(result)
#'
process_se <- function(se_obj) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment package is required.")
    }
    if (!inherits(se_obj, "SummarizedExperiment"))  {
        stop("Provided object is not a SummarizedExperiment.")
    }
    assay_mat <- SummarizedExperiment::assay(se_obj)
    assay_data <- as.data.frame(t(assay_mat))                # samples as rows, metabolites as columns
    metadata_sample <- as.data.frame(SummarizedExperiment::colData(se_obj))   # sample metadata
    metadata_feature <- as.data.frame(SummarizedExperiment::rowData(se_obj))  # feature metadata

    return(list(
        data = assay_data,
        metadata_sample = metadata_sample,
        metadata_feature = metadata_feature
    ))
}
