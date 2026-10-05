#!/usr/bin/env Rscript
#
#  This file is part of the `MetaProViz` R package
#
#  Copyright 2023-2026
#  Saez Lab, Heidelberg University, University Hospital Heidelberg,
#  European Bioinformatics Institute (EMBL-EBI), University of Cologne
#
#  Authors: see the file `README.md`
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file `LICENSE.md` or copy at
#      https://opensource.org/license/bsd-3-clause
#
#  Website: https://saezlab.github.io/MetaProViz
#  Git repo: https://github.com/saezlab/MetaProViz
#

#
# Shared helpers for the prior knowledge networks: matching measured
# features to prior knowledge, building nodes and edges, set similarity,
# node/edge aesthetics and graph layout.
#


##
## Matching features to prior knowledge
##

#' Match measured features to a prior knowledge table
#'
#' Splits multi-ID cells of `feature_metadata`, normalises the IDs on both
#' sides and joins the features to `input_pk`. The returned associations keep
#' all columns of `input_pk`, so callers can derive edge and node attributes
#' from them.
#'
#' @param feature_metadata Data frame with one row per measured feature.
#' @param input_pk Long prior knowledge table with one ID per row.
#' @param metadata_info Named vector with at least `InputID`, `PriorID` and
#'     `PriorTerm`; `InputLabel` is optional.
#' @param id_type `"HMDB"` normalises HMDB IDs before matching; any other
#'     value matches the trimmed IDs exactly.
#' @param id_sep Separator of multiple IDs in one cell of
#'     `feature_metadata[[InputID]]`.
#'
#' @return A list with
#'     \item{associations}{One row per feature-ID-PK row match. Internal key
#'     columns `.feature_row_id`, `.metabolite`, `.term` and `.matched_id`
#'     are added in front of the `input_pk` columns.}
#'     \item{matched_features}{Expanded feature IDs with a match.}
#'     \item{unmatched_features}{Expanded feature IDs without a match,
#'     including features without any valid ID.}
#'
#' @noRd
.pk_associations <- function(
    feature_metadata,
    input_pk,
    metadata_info,
    id_type,
    id_sep
) {
    input_id <- metadata_info[["InputID"]]
    prior_id <- metadata_info[["PriorID"]]
    prior_term <- metadata_info[["PriorTerm"]]
    input_label <- if ("InputLabel" %in% names(metadata_info)) {
        metadata_info[["InputLabel"]]
    } else {
        NULL
    }

    raw_ids <- as.character(feature_metadata[[input_id]])
    normalized_cells <- vapply(
        raw_ids,
        function(cell) {
            ids <- .normalize_pk_ids(.split_id_cell(cell, sep = id_sep), id_type)
            ids <- unique(ids[!is.na(ids)])
            if (length(ids) == 0L) NA_character_ else paste(ids, collapse = id_sep)
        },
        character(1),
        USE.NAMES = FALSE
    )
    labels <- if (is.null(input_label)) {
        rep(NA_character_, nrow(feature_metadata))
    } else {
        as.character(feature_metadata[[input_label]])
    }
    # Features without a label are named by their (normalised) IDs
    labels <- .coalesce_network_label(labels, normalized_cells)
    labels <- .coalesce_network_label(labels, raw_ids)

    # One row per feature and ID
    id_lists <- lapply(raw_ids, .split_id_cell, sep = id_sep)
    expanded <- dplyr::tibble(
        .feature_row_id = rep(seq_len(nrow(feature_metadata)), lengths(id_lists)),
        metabolite = rep(labels, lengths(id_lists)),
        id_input = unlist(id_lists, use.names = FALSE)
    )
    expanded$id_normalized <- .normalize_pk_ids(expanded$id_input, id_type)

    # Features without any ID are kept in `unmatched_features`
    no_id <- setdiff(seq_len(nrow(feature_metadata)), expanded$.feature_row_id)
    expanded <- dplyr::bind_rows(
        expanded,
        dplyr::tibble(
            .feature_row_id = no_id,
            metabolite = labels[no_id],
            id_input = raw_ids[no_id],
            id_normalized = NA_character_
        )
    ) |>
        dplyr::distinct() |>
        dplyr::arrange(.data$.feature_row_id)

    pk <- input_pk
    pk$.matched_id <- .normalize_pk_ids(pk[[prior_id]], id_type)
    pk$.term <- as.character(pk[[prior_term]])
    pk <- pk[!is.na(pk$.matched_id) & !is.na(pk$.term) & pk$.term != "", , drop = FALSE]

    associations <- expanded |>
        dplyr::filter(!is.na(.data$id_normalized)) |>
        dplyr::transmute(
            .feature_row_id = .data$.feature_row_id,
            .metabolite = .data$metabolite,
            .matched_id = .data$id_normalized
        ) |>
        dplyr::distinct() |>
        dplyr::inner_join(pk, by = ".matched_id", relationship = "many-to-many")
    associations <- dplyr::relocate(
        associations,
        dplyr::all_of(c(".feature_row_id", ".metabolite", ".term", ".matched_id"))
    )

    is_matched <- paste(expanded$.feature_row_id, expanded$id_normalized) %in%
        paste(associations$.feature_row_id, associations$.matched_id)

    list(
        associations = associations,
        matched_features = expanded[is_matched, , drop = FALSE],
        unmatched_features = expanded[!is_matched, , drop = FALSE]
    )
}

#' Normalise metabolite IDs for matching
#'
#' @param x Character vector of IDs.
#' @param id_type `"HMDB"` brings HMDB IDs into the 7-digit `HMDB0000000`
#'     form; any other value only trims whitespace.
#'
#' @return Character vector of the same length; invalid IDs become `NA`.
#'
#' @noRd
.normalize_pk_ids <- function(x, id_type) {
    if (identical(toupper(id_type), "HMDB")) {
        return(.normalize_hmdb_token(x))
    }
    x <- stringr::str_trim(as.character(x))
    x[toupper(x) %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_
    x
}

#' @noRd
.normalize_hmdb_token <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x[toupper(x) %in% c("", "NA", "N/A", "NULL", "NAN")] <- NA_character_

    out <- rep(NA_character_, length(x))
    keep <- !is.na(x)
    if (!any(keep)) {
        return(out)
    }

    x_keep <- toupper(x[keep])
    x_keep <- gsub("^HMDB[: _-]*", "", x_keep)
    x_keep <- gsub("^0+(?=[0-9]+$)", "", x_keep, perl = TRUE)

    valid_digits <- grepl("^[0-9]+$", x_keep)
    digits <- x_keep
    digits[!valid_digits] <- NA_character_

    out[keep] <- ifelse(
        !is.na(digits),
        sprintf("HMDB%07d", as.integer(digits)),
        NA_character_
    )

    out
}

#' @noRd
.split_id_cell <- function(x, sep) {
    if (length(x) == 0L || is.null(x) || is.na(x)) {
        return(character(0))
    }

    parts <- unlist(strsplit(as.character(x), split = sep, fixed = TRUE), use.names = FALSE)
    parts <- stringr::str_trim(parts)
    parts[!is.na(parts) & parts != ""]
}

#' @noRd
.coalesce_network_label <- function(label, fallback) {
    label <- as.character(label)
    fallback <- as.character(fallback)
    use_fallback <- is.na(label) | stringr::str_trim(label) == ""
    label[use_fallback] <- fallback[use_fallback]
    label
}


##
## Node and edge attributes
##

#' Collapse a column to one value per group
#'
#' Numeric columns are averaged, all other columns take their first
#' non-missing value.
#'
#' @noRd
.summarise_attribute <- function(x) {
    if (is.numeric(x)) {
        if (all(is.na(x))) {
            return(NA_real_)
        }
        return(mean(x, na.rm = TRUE))
    }
    .first_non_missing(x)
}

#' @noRd
.first_non_missing <- function(x, default = NA_character_) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    if (length(x) == 0L) {
        return(default)
    }
    x[[1]]
}

#' Attach one attribute column per node
#'
#' @param nodes Data frame with a `name` column.
#' @param source Data frame holding the attribute.
#' @param key Column in `source` matching `nodes$name`.
#' @param column Attribute column in `source`.
#' @param target Name of the new column in `nodes`.
#'
#' @noRd
.add_node_attribute <- function(nodes, source, key, column, target) {
    values <- source |>
        dplyr::filter(!is.na(.data[[key]])) |>
        dplyr::group_by(name = as.character(.data[[key]])) |>
        dplyr::summarise(
            !!target := .summarise_attribute(.data[[column]]),
            .groups = "drop"
        )
    dplyr::left_join(nodes, values, by = "name")
}

#' Shorten labels and hide labels of low-degree nodes
#'
#' @noRd
.format_network_labels <- function(labels, degrees, label_mode, label_max_chars, label_degree_min) {
    labels <- as.character(labels)
    too_long <- !is.na(labels) & nchar(labels, type = "width") > label_max_chars
    labels[too_long] <- paste0(substr(labels[too_long], 1L, label_max_chars - 3L), "...")

    if (identical(label_mode, "all")) {
        return(labels)
    }

    labels[is.na(degrees) | degrees < label_degree_min] <- NA_character_
    labels
}


##
## Set similarity
##

#' Pairwise similarity of sets
#'
#' Used by [cluster_pk()] for term-term similarity and by
#' [viz_shared_pk_network()] for metabolite-metabolite similarity, so both
#' report the same values.
#'
#' @param sets Named list of character vectors without duplicates.
#' @param method `"shared"` returns the number of shared elements,
#'     `"jaccard"` |A intersect B| / |A union B| and `"overlap_coefficient"`
#'     |A intersect B| / min(|A|, |B|). The diagonal of the two coefficients
#'     is 1.
#'
#' @return Symmetric numeric matrix with the set names as dimnames.
#'
#' @noRd
.set_similarity <- function(sets, method = c("jaccard", "overlap_coefficient", "shared")) {
    method <- match.arg(method)

    elements <- unique(unlist(sets, use.names = FALSE))
    incidence <- matrix(
        0,
        nrow = length(sets),
        ncol = length(elements),
        dimnames = list(names(sets), elements)
    )
    incidence[cbind(
        rep(seq_along(sets), lengths(sets)),
        match(unlist(sets, use.names = FALSE), elements)
    )] <- 1

    shared <- tcrossprod(incidence)
    if (method == "shared") {
        return(shared)
    }

    sizes <- diag(shared)
    denom <- if (method == "jaccard") {
        outer(sizes, sizes, "+") - shared
    } else {
        outer(sizes, sizes, pmin)
    }
    similarity <- shared / denom
    similarity[denom == 0] <- 0
    diag(similarity) <- 1
    similarity
}


##
## Plotting
##

#' Graph layout
#'
#' Random layouts such as the default Fruchterman-Reingold ("fr") are
#' reproducible with a seed, without changing the global seed.
#'
#' @noRd
.network_layout <- function(graph, seed = NULL, layout = "fr") {
    if (is.null(seed)) {
        return(ggraph::ggraph(graph, layout = layout))
    }
    withr::with_seed(seed, ggraph::ggraph(graph, layout = layout))
}

#' Fill scale for a node attribute
#'
#' Numeric attributes get a continuous scale, which is diverging around 0
#' when the values have both signs (e.g. Log2FC). Other attributes get a
#' discrete palette, whose legend keys use the node `shape`.
#'
#' @noRd
.network_fill_scale <- function(values, name, shape = 21, option = "D") {
    if (is.numeric(values)) {
        if (any(values < 0, na.rm = TRUE) && any(values > 0, na.rm = TRUE)) {
            limit <- max(abs(values), na.rm = TRUE)
            return(ggplot2::scale_fill_gradient2(
                name = name,
                low = "#2166ac",
                mid = "white",
                high = "#b2182b",
                midpoint = 0,
                limits = c(-limit, limit),
                na.value = "grey80"
            ))
        }
        return(ggplot2::scale_fill_viridis_c(name = name, option = option, na.value = "grey80"))
    }

    ggplot2::scale_fill_manual(
        name = name,
        values = .discrete_palette(values),
        na.value = "grey80",
        guide = ggplot2::guide_legend(override.aes = list(shape = shape, size = 4))
    )
}

#' Named colours for the non-missing values of a discrete attribute
#'
#' @noRd
.discrete_palette <- function(values) {
    levels <- sort(unique(stats::na.omit(as.character(values))))
    if (length(levels) == 0L) {
        # Only missing values: they are drawn with `na.value`
        return("grey80")
    }
    stats::setNames(.network_palette(length(levels)), levels)
}

#' @noRd
.network_palette <- function(n) {
    # Okabe-Ito colours; yellow last as it is hard to see on white
    okabe_ito <- c(
        "#E69F00", "#56B4E9", "#009E73", "#0072B2",
        "#D55E00", "#CC79A7", "#999999", "#F0E442"
    )
    if (n <= length(okabe_ito)) {
        return(okabe_ito[seq_len(n)])
    }
    grDevices::hcl.colors(n, "Dynamic")
}
