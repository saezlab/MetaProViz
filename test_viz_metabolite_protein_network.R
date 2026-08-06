#!/usr/bin/env Rscript

required_packages <- c(
    "dplyr",
    "tidyr",
    "stringr",
    "igraph",
    "ggraph",
    "ggplot2",
    "logger",
    "rlang",
    "tibble"
)
# library(MetaProViz)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
    stop(
        "Missing required package(s): ",
        paste(missing_packages, collapse = ", "),
        ". Install them before running this script."
    )
}

source("R/VizMetaboliteProteinNetwork.R")

feature_metadata <- data.frame(
    Metabolite = c("Lactate", "Succinate", NA),
    hmdb = c("HMDB0000190", "1544", "HMDB0000254; HMDB0000331"),
    stringsAsFactors = FALSE
)

metalinks_subset <- data.frame(
    hmdb = c("HMDB0000190", "HMDB0001544", "HMDB0000254"),
    gene_symbol = c("HCAR1", "SLC13A3", "OXGR1"),
    protein_type = c("receptor", "transporter", "receptor"),
    protein_type_clean = c("Receptor", "Transporter", "Receptor"),
    transport_direction = c(NA, "in", NA),
    type = c("Ligand-Receptor", "Production-Degradation", "Ligand-Receptor"),
    mode_of_regulation = c("Activating", NA, "Activating"),
    interaction_family = c("Receptor", "Transporter", "Receptor"),
    source = c("Example", "Example", "Example"),
    stringsAsFactors = FALSE
)

message("Running named-metabolite test...")
res_named <- viz_metabolite_protein_network(
    feature_metadata = feature_metadata,
    metalinks_df = metalinks_subset,
    metabolite_col = "Metabolite",
    hmdb_col = "hmdb",
    save_plot = NULL,
    print_plot = FALSE
)

stopifnot(is.list(res_named))
stopifnot(all(c("plot", "edges", "nodes", "matched_features", "unmatched_features", "saved_files") %in% names(res_named)))
stopifnot(inherits(res_named$plot, "ggplot"))
stopifnot(nrow(res_named$edges) > 0)
stopifnot(nrow(res_named$nodes) > 0)
stopifnot("Lactate" %in% res_named$matched_features$metabolite)
stopifnot("Succinate" %in% res_named$matched_features$metabolite)
stopifnot(length(res_named$saved_files) == 0)

message("Running HMDB fallback-label test...")
res_hmdb <- viz_metabolite_protein_network(
    feature_metadata = feature_metadata["hmdb"],
    metalinks_df = metalinks_subset,
    hmdb_col = "hmdb",
    save_plot = NULL,
    print_plot = FALSE
)

stopifnot("HMDB0000190" %in% res_hmdb$matched_features$metabolite)
stopifnot("HMDB0001544" %in% res_hmdb$matched_features$metabolite)
stopifnot("HMDB0000254" %in% res_hmdb$matched_features$metabolite)
stopifnot("HMDB0000190" %in% res_hmdb$nodes$name)
stopifnot("HMDB0001544" %in% res_hmdb$nodes$name)
stopifnot("HMDB0000254" %in% res_hmdb$nodes$name)

message("Running no-match test...")
res_empty <- viz_metabolite_protein_network(
    feature_metadata = data.frame(hmdb = "HMDB9999999", stringsAsFactors = FALSE),
    metalinks_df = metalinks_subset,
    hmdb_col = "hmdb",
    save_plot = NULL,
    print_plot = FALSE
)

stopifnot(is.null(res_empty$plot))
stopifnot(nrow(res_empty$edges) == 0)
stopifnot(nrow(res_empty$nodes) == 0)
stopifnot(nrow(res_empty$unmatched_features) == 1)

message("All viz_metabolite_protein_network checks passed.")


