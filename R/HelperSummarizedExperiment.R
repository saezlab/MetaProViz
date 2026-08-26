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

#' Prepare standalone helper input
#'
#' Standardizes standalone helper handling of SummarizedExperiment input by
#' unpacking assay, sample metadata, and feature metadata when needed.
#'
#' @param data Input assay data frame or SummarizedExperiment.
#' @param metadata_sample Optional sample metadata for data frame input.
#'
#' @return A list containing the original input, assay data, sample metadata,
#'     and feature metadata.
#' @noRd
prepare_se_input <- function(data, metadata_sample = NULL) {
    input_data <- data
    metadata_feature <- NULL

    if (inherits(data, "SummarizedExperiment")) {
        log_info('Processing input SummarizedExperiment object.')
        se_list <- process_se(data)
        data <- se_list$data
        metadata_sample <- se_list$metadata_sample
        metadata_feature <- se_list$metadata_feature
    }

    list(
        input_data = input_data,
        data = data,
        metadata_sample = metadata_sample,
        metadata_feature = metadata_feature
    )
}

#' Build SummarizedExperiment from standalone helper output
#'
#' Reconstructs a SummarizedExperiment from a samples-by-features data frame
#' while aligning sample and feature metadata to the retained dimensions.
#'
#' @param assay_df Data frame with samples as rows and metabolites as columns.
#' @param metadata_sample Optional sample metadata.
#' @param metadata_feature Optional feature metadata.
#'
#' @return A SummarizedExperiment object.
#' @noRd
build_se_from_df <- function(
        assay_df,
        metadata_sample = NULL,
        metadata_feature = NULL
) {
    assay_df <- as.data.frame(assay_df)

    if (!is.null(metadata_sample)) {
        metadata_sample <- as.data.frame(metadata_sample)
        metadata_sample <- metadata_sample[rownames(assay_df), , drop = FALSE]
    }

    if (!is.null(metadata_feature)) {
        metadata_feature <- as.data.frame(metadata_feature)
        metadata_feature <- metadata_feature[colnames(assay_df), , drop = FALSE]
    }

    se_args <- list(
        assays = list(data = as.matrix(t(assay_df)))
    )

    if (!is.null(metadata_sample)) {
        se_args$colData <- S4Vectors::DataFrame(metadata_sample)
    }

    if (!is.null(metadata_feature)) {
        se_args$rowData <- S4Vectors::DataFrame(metadata_feature)
    }

    do.call(SummarizedExperiment::SummarizedExperiment, se_args)
}
