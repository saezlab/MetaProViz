# !/usr/bin/env Rscript

#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2022-2025
#  Saez Lab, Heidelberg University
#
#  Authors: see the file `README.md`
#
#  Distributed under the GNU GPLv3 License.
#  See accompanying file `LICENSE` or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#


#' Process SummarizedExperiment objects
#'
#' This helper function extracts the numeric data matrix, sample metadata, and
#' feature metadata from a SummarizedExperiment object. If the input is not a
#' SummarizedExperiment, an error is thrown.
#'
#' @param se_obj A SummarizedExperiment object.
#'
#' @return A list with elements: \item{data}{Numeric matrix with metabolites as
#'     columns and samples as rows} \item{metadata_sample}{Data frame with
#'     sample metadata (samples as row names)} \item{metadata_feature}{Data
#'     frame with feature metadata (metabolites as row names)}
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @noRd
process_se <- function(se_obj) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment package is required.")
    }
    if (!inherits(se_obj, "SummarizedExperiment")) {
        stop("Provided object is not a SummarizedExperiment.")
    }
    assay_mat <- assay(se_obj)
    # samples as rows, metabolites as columns
    assay_data <- as.data.frame(t(assay_mat))
    metadata_sample <- as.data.frame(colData(se_obj))  # sample metadata
    metadata_feature <- as.data.frame(rowData(se_obj))  # feature metadata

    # treat "no columns" as empty
    if (ncol(metadata_feature) == 0) metadata_feature <- NULL
    if (ncol(metadata_sample) == 0) metadata_sample <- NULL

    #return
    return(list(
        data = assay_data,
        metadata_sample = metadata_sample,
        metadata_feature = metadata_feature
    ))
}
