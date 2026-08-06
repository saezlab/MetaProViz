#!/usr/bin/env Rscript
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

#' Visualize metabolite-protein networks from MetalinksDB
#'
#' Plot a directed metabolite-protein interaction network based on the tables
#' returned by [metsigdb_metalinks()]. The function also accepts pre-filtered
#' subsets of the MetalinksDB interaction table and subsets of metabolite
#' metadata.
#'
#' HMDB identifiers in `feature_metadata` and `metalinks_df` are normalized
#' before matching. Metabolite names are encouraged but not required. If
#' `metabolite_col` is omitted, or if some metabolite names are missing, the
#' corresponding normalized HMDB identifiers are used as metabolite labels.
#'
#' @param feature_metadata Data frame containing at least an HMDB identifier
#'     column and, optionally, a metabolite-name column.
#' @param metalinks_df Data frame returned by [metsigdb_metalinks()] or any
#'     subset of it. The input must contain at least the columns `hmdb` and
#'     `gene_symbol`. If present, the columns `metabolite`, `protein_type`,
#'     `protein_type_clean`, `transport_direction`, `type`,
#'     `mode_of_regulation`, `interaction_family`, and `source` are used to
#'     enrich the network output.
#' @param metabolite_col \emph{Optional: } Single column name in
#'     `feature_metadata` containing metabolite labels. If `NULL`, normalized
#'     HMDB identifiers are used as metabolite labels. \strong{Default = NULL}
#' @param hmdb_col Single column name in `feature_metadata` containing HMDB
#'     identifiers. \strong{Default = "hmdb"}
#' @param save_plot \emph{Optional: } File type(s) used to save the network
#'     plot. Options are `svg`, `pdf`, `png`, or `NULL`.
#'     \strong{Default = NULL}
#' @param path \emph{Optional: } Output directory for saved plots.
#'     \strong{Default = NULL}
#' @param plot_name \emph{Optional: } Title of the plot and basename for saved
#'     plot files. \strong{Default = "metalinks_network"}
#' @param width \emph{Optional: } Plot width in inches for saved plots.
#'     \strong{Default = 10}
#' @param height \emph{Optional: } Plot height in inches for saved plots.
#'     \strong{Default = 8}
#' @param hmdb_sep \emph{Optional: } Separator used for multiple HMDB IDs in a
#'     single cell of `feature_metadata[[hmdb_col]]`. \strong{Default = ";"}
#' @param return_data \emph{Optional: } If `TRUE`, return a list containing the
#'     plot object and intermediate network tables. If `FALSE`, return the plot
#'     invisibly. \strong{Default = TRUE}
#' @param print_plot \emph{Optional: } If `TRUE`, print the plot.
#'     \strong{Default = TRUE}
#' @param full_labels \emph{Optional: } Backward-compatible shortcut for
#'     `label_mode = "all"`. If `TRUE`, all node labels are shown regardless of
#'     degree. \strong{Default = FALSE}
#' @param label_mode \emph{Optional: } Labeling mode for network nodes.
#'     `"reduced"` only labels nodes with degree at least `label_degree_min`;
#'     `"all"` labels all nodes. \strong{Default = "reduced"}
#' @param label_max_chars \emph{Optional: } Maximum label width before truncation
#'     with `...`. \strong{Default = 20}
#' @param label_repel \emph{Optional: } If `TRUE`, use repelled text labels.
#'     \strong{Default = TRUE}
#' @param label_degree_min \emph{Optional: } Minimum node degree required for
#'     labeling when `label_mode = "reduced"`. \strong{Default = 2}
#'
#' @return If `return_data = TRUE`, a list with the following elements:
#' \describe{
#'   \item{plot}{The `ggraph`/`ggplot2` network plot, or `NULL` if no network
#'   could be built.}
#'   \item{edges}{A data frame of network edges.}
#'   \item{nodes}{A data frame of network nodes.}
#'   \item{matched_features}{Expanded and normalized feature-to-HMDB mappings
#'   that matched the Metalinks table.}
#'   \item{unmatched_features}{Expanded and normalized feature-to-HMDB mappings
#'   without a Metalinks match.}
#'   \item{saved_files}{Character vector of written plot files.}
#' }
#' If `return_data = FALSE`, the plot is returned invisibly.
#'
#' @examples
#' \dontrun{
#' metalinks_transporters <- metsigdb_metalinks(
#'     cell_location = c("Extracellular"),
#'     tissue_location = c("Kidney", "All Tissues"),
#'     biospecimen_location = c("Blood", "Urine")
#' ) |>
#'     dplyr::filter(stringr::str_detect(tolower(protein_type_clean), "transporter")) |>
#'     dplyr::distinct(hmdb, .keep_all = TRUE) |>
#'     dplyr::slice_head(n = 10)
#'
#' feature_metadata <- metalinks_transporters |>
#'     dplyr::transmute(
#'         Metabolite = dplyr::coalesce(.data$metabolite, .data$hmdb),
#'         hmdb = .data$hmdb
#'     )
#'
#' viz_metabolite_protein_network(
#'     feature_metadata = feature_metadata,
#'     metalinks_df = metalinks_transporters,
#'     metabolite_col = "Metabolite",
#'     hmdb_col = "hmdb",
#'     save_plot = NULL,
#'     print_plot = FALSE
#' )
#' }
#'
#' @importFrom logger log_info
#' @export
viz_metabolite_protein_network <- function(
    feature_metadata,
    metalinks_df,
    metabolite_col = NULL,
    hmdb_col = "hmdb",
    save_plot = NULL,
    path = NULL,
    plot_name = "metalinks_network",
    width = 10,
    height = 8,
    hmdb_sep = ";",
    return_data = TRUE,
    print_plot = TRUE,
    full_labels = FALSE,
    label_mode = c("reduced", "all"),
    label_max_chars = 20,
    label_repel = TRUE,
    label_degree_min = 2
) {
    check_param_VizMetaboliteProteinNetwork(
        feature_metadata = feature_metadata,
        metalinks_df = metalinks_df,
        metabolite_col = metabolite_col,
        hmdb_col = hmdb_col,
        save_plot = save_plot,
        path = path,
        plot_name = plot_name,
        width = width,
        height = height,
        hmdb_sep = hmdb_sep,
        return_data = return_data,
        print_plot = print_plot,
        full_labels = full_labels,
        label_mode = label_mode,
        label_max_chars = label_max_chars,
        label_repel = label_repel,
        label_degree_min = label_degree_min
    )

    logger::log_info("viz_metabolite_protein_network: Metalinks metabolite-protein network")

    label_mode <- match.arg(label_mode)
    if (isTRUE(full_labels)) {
        label_mode <- "all"
    }

    save_plot <- .normalize_metalinks_network_save_plot(save_plot)
    hmdb_sym <- rlang::sym(hmdb_col)

    feature_input <- feature_metadata |>
        dplyr::mutate(
            .feature_row_id = seq_len(nrow(feature_metadata)),
            .metabolite_name = .extract_metabolite_labels(
                feature_metadata = feature_metadata,
                metabolite_col = metabolite_col
            ),
            .hmdb_raw = as.character(!!hmdb_sym)
        )

    matched_features <- .expand_feature_hmdb(
        feature_input = feature_input,
        hmdb_sep = hmdb_sep
    )
    metalinks_norm <- .normalize_metalinks_df(metalinks_df)

    unmatched_features <- feature_input |>
        dplyr::transmute(
            .feature_row_id = .data$.feature_row_id,
            metabolite = .coalesce_network_label(
                label = .data$.metabolite_name,
                fallback = vapply(
                    .data$.hmdb_raw,
                    .normalize_hmdb_cell,
                    character(1),
                    sep = hmdb_sep
                )
            ),
            hmdb_input = .data$.hmdb_raw,
            hmdb_normalized = vapply(
                .data$.hmdb_raw,
                .normalize_hmdb_cell,
                character(1),
                sep = hmdb_sep
            )
        )

    if (nrow(matched_features) == 0L) {
        warning("No valid HMDB IDs found in `feature_metadata` after normalization.")
        result <- list(
            plot = NULL,
            edges = dplyr::tibble(),
            nodes = dplyr::tibble(),
            matched_features = dplyr::tibble(),
            unmatched_features = unmatched_features,
            saved_files = character(0)
        )
        return(if (isTRUE(return_data)) result else invisible(NULL))
    }

    interactions <- matched_features |>
        dplyr::inner_join(
            metalinks_norm,
            by = c("hmdb_normalized" = "hmdb_normalized"),
            relationship = "many-to-many"
        ) |>
        dplyr::filter(!is.na(.data$gene_symbol), .data$gene_symbol != "") |>
        dplyr::distinct()

    unmatched_features <- matched_features |>
        dplyr::distinct(
            .data$.feature_row_id,
            .data$metabolite,
            .data$hmdb_input,
            .data$hmdb_normalized
        ) |>
        dplyr::anti_join(
            interactions |>
                dplyr::distinct(
                    .data$.feature_row_id,
                    .data$metabolite,
                    .data$hmdb_input,
                    .data$hmdb_normalized
                ),
            by = c(".feature_row_id", "metabolite", "hmdb_input", "hmdb_normalized")
        )

    if (nrow(interactions) == 0L) {
        warning("No Metalinks interactions matched the normalized HMDB IDs.")
        result <- list(
            plot = NULL,
            edges = dplyr::tibble(),
            nodes = dplyr::tibble(),
            matched_features = matched_features,
            unmatched_features = unmatched_features,
            saved_files = character(0)
        )
        return(if (isTRUE(return_data)) result else invisible(NULL))
    }

    edges <- .build_metalinks_network_edges(interactions)
    nodes <- .build_metalinks_network_nodes(
        edges = edges,
        interactions = interactions,
        label_mode = label_mode,
        label_max_chars = label_max_chars,
        label_degree_min = label_degree_min
    )

    if (nrow(edges) == 0L || nrow(nodes) == 0L) {
        warning("No valid network edges could be constructed from the matched interactions.")
        result <- list(
            plot = NULL,
            edges = edges,
            nodes = nodes,
            matched_features = matched_features,
            unmatched_features = unmatched_features,
            saved_files = character(0)
        )
        return(if (isTRUE(return_data)) result else invisible(NULL))
    }

    plot_obj <- .make_metalinks_network_plot(
        nodes = nodes,
        edges = edges,
        plot_name = plot_name,
        label_repel = label_repel
    )

    if (isTRUE(print_plot)) {
        print(plot_obj)
    }

    saved_files <- .save_metalinks_network_plot(
        plot = plot_obj,
        save_plot = save_plot,
        path = path,
        plot_name = plot_name,
        width = width,
        height = height
    )

    result <- list(
        plot = plot_obj,
        edges = edges,
        nodes = nodes,
        matched_features = matched_features,
        unmatched_features = unmatched_features,
        saved_files = saved_files
    )

    if (isTRUE(return_data)) {
        return(result)
    }

    invisible(plot_obj)
}

#' @noRd
.extract_metabolite_labels <- function(feature_metadata, metabolite_col) {
    if (is.null(metabolite_col)) {
        return(rep(NA_character_, nrow(feature_metadata)))
    }

    as.character(feature_metadata[[metabolite_col]])
}

#' @noRd
.normalize_metalinks_network_save_plot <- function(save_plot) {
    if (is.null(save_plot)) {
        return(NULL)
    }

    save_plot <- tolower(as.character(save_plot))
    save_plot <- unique(save_plot[save_plot %in% c("png", "svg", "pdf")])

    if (length(save_plot) == 0L) {
        stop(
            "`save_plot` must contain only supported format(s): ",
            paste(c("png", "svg", "pdf"), collapse = ", "),
            ", or NULL."
        )
    }

    save_plot
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
.normalize_hmdb_cell <- function(x, sep = ";") {
    tokens <- .split_id_cell(x, sep = sep)
    if (length(tokens) == 0L) {
        return(NA_character_)
    }

    normalized <- .normalize_hmdb_token(tokens)
    normalized <- unique(normalized[!is.na(normalized) & normalized != ""])

    if (length(normalized) == 0L) {
        return(NA_character_)
    }

    paste(normalized, collapse = sep)
}

#' @noRd
.split_id_cell <- function(x, sep = ";") {
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

#' @noRd
.expand_feature_hmdb <- function(feature_input, hmdb_sep) {
    feature_input |>
        dplyr::transmute(
            .feature_row_id = .data$.feature_row_id,
            .metabolite_name = as.character(.data$.metabolite_name),
            hmdb_input = as.character(.data$.hmdb_raw)
        ) |>
        tidyr::separate_rows(hmdb_input, sep = hmdb_sep) |>
        dplyr::mutate(
            hmdb_input = stringr::str_trim(.data$hmdb_input),
            hmdb_normalized = .normalize_hmdb_token(.data$hmdb_input),
            metabolite = .coalesce_network_label(
                label = .data$.metabolite_name,
                fallback = .data$hmdb_normalized
            )
        ) |>
        dplyr::filter(!is.na(.data$hmdb_normalized), .data$hmdb_normalized != "") |>
        dplyr::filter(!is.na(.data$metabolite), .data$metabolite != "") |>
        dplyr::select(
            .data$.feature_row_id,
            .data$metabolite,
            .data$hmdb_input,
            .data$hmdb_normalized
        ) |>
        dplyr::distinct()
}

#' @noRd
.get_col_or_na <- function(df, col) {
    if (col %in% colnames(df)) {
        return(df[[col]])
    }
    rep(NA_character_, nrow(df))
}

#' @noRd
.normalize_metalinks_df <- function(metalinks_df) {
    protein_type <- .get_col_or_na(metalinks_df, "protein_type")
    protein_type_clean <- .get_col_or_na(metalinks_df, "protein_type_clean")

    metalinks_df |>
        dplyr::transmute(
            hmdb_normalized = .normalize_hmdb_token(.data$hmdb),
            metabolite_pk = as.character(.get_col_or_na(metalinks_df, "metabolite")),
            gene_symbol = as.character(.data$gene_symbol),
            protein_type = gsub("\"", "", as.character(protein_type)),
            protein_type_clean = gsub("\"", "", as.character(protein_type_clean)),
            transport_direction = as.character(.get_col_or_na(metalinks_df, "transport_direction")),
            type = as.character(.get_col_or_na(metalinks_df, "type")),
            mode_of_regulation = as.character(.get_col_or_na(metalinks_df, "mode_of_regulation")),
            interaction_family = as.character(.get_col_or_na(metalinks_df, "interaction_family")),
            source = as.character(.get_col_or_na(metalinks_df, "source"))
        ) |>
        dplyr::filter(!is.na(.data$hmdb_normalized), .data$hmdb_normalized != "")
}

#' @noRd
.classify_interaction <- function(transport_direction, type) {
    dplyr::case_when(
        !is.na(transport_direction) & transport_direction == "in" ~ "Transport_In",
        !is.na(transport_direction) & transport_direction == "out" ~ "Transport_Out",
        !is.na(type) & type == "Ligand-Receptor" ~ "Ligand-Receptor",
        !is.na(type) & type == "Production-Degradation" ~ "Production-Degradation",
        TRUE ~ "Unknown"
    )
}

#' @noRd
.classify_regulation <- function(mode_of_regulation) {
    dplyr::case_when(
        !is.na(mode_of_regulation) & mode_of_regulation == "Activating" ~ "Activating",
        !is.na(mode_of_regulation) & mode_of_regulation == "Inhibiting" ~ "Inhibiting",
        TRUE ~ NA_character_
    )
}

#' @noRd
.build_metalinks_network_edges <- function(interactions) {
    interactions |>
        dplyr::mutate(
            interaction = .classify_interaction(
                transport_direction = .data$transport_direction,
                type = .data$type
            ),
            regulation = .classify_regulation(.data$mode_of_regulation),
            from = dplyr::if_else(
                .data$interaction %in% c("Transport_In", "Production-Degradation"),
                .data$gene_symbol,
                .data$metabolite
            ),
            to = dplyr::if_else(
                .data$interaction %in% c("Transport_Out", "Ligand-Receptor"),
                .data$gene_symbol,
                .data$metabolite
            )
        ) |>
        dplyr::filter(
            !is.na(.data$from),
            !is.na(.data$to),
            .data$from != "",
            .data$to != "",
            .data$from != .data$to
        ) |>
        dplyr::group_by(.data$from, .data$to, .data$interaction, .data$regulation) |>
        dplyr::summarise(
            n_matches = dplyr::n(),
            hmdb_ids = paste(unique(.data$hmdb_normalized), collapse = ";"),
            sources = paste(unique(stats::na.omit(.data$source)), collapse = ";"),
            interaction_family = paste(
                unique(stats::na.omit(.data$interaction_family)),
                collapse = ";"
            ),
            .groups = "drop"
        ) |>
        dplyr::distinct()
}

#' @noRd
.build_metalinks_network_nodes <- function(edges, interactions, label_mode, label_max_chars, label_degree_min) {
    node_names <- unique(c(edges$from, edges$to))

    metabolite_nodes <- interactions |>
        dplyr::distinct(name = .data$metabolite) |>
        dplyr::mutate(
            node_type = "Metabolite",
            node_group = "Metabolite"
        )

    protein_nodes <- interactions |>
        dplyr::filter(!is.na(.data$gene_symbol), .data$gene_symbol != "") |>
        dplyr::group_by(name = .data$gene_symbol) |>
        dplyr::summarise(
            node_type = "Protein",
            node_group = .first_non_missing(
                c(
                    .data$protein_type_clean,
                    .data$protein_type,
                    .data$interaction_family
                ),
                default = "Protein"
            ),
            .groups = "drop"
        )

    nodes <- dplyr::bind_rows(metabolite_nodes, protein_nodes) |>
        dplyr::filter(.data$name %in% node_names) |>
        dplyr::distinct(.data$name, .keep_all = TRUE)

    graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)
    nodes$degree <- igraph::degree(graph, mode = "all")
    nodes$label <- .format_network_labels(
        labels = nodes$name,
        degrees = nodes$degree,
        label_mode = label_mode,
        label_max_chars = label_max_chars,
        label_degree_min = label_degree_min
    )
    nodes
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

#' @noRd
.format_network_labels <- function(labels, degrees, label_mode, label_max_chars, label_degree_min) {
    labels <- as.character(labels)

    rendered <- vapply(
        labels,
        function(label) {
            if (is.na(label)) {
                return(NA_character_)
            }
            if (nchar(label, type = "width") <= label_max_chars) {
                return(label)
            }
            paste0(substr(label, 1L, label_max_chars - 3L), "...")
        },
        character(1)
    )

    if (identical(label_mode, "all")) {
        return(rendered)
    }

    rendered[is.na(degrees) | degrees < label_degree_min] <- NA_character_
    rendered
}

#' @noRd
.make_metalinks_network_plot <- function(nodes, edges, plot_name, label_repel) {
    graph <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

    edge_colors <- c(
        "Transport_In" = "#1f78b4",
        "Transport_Out" = "#e31a1c",
        "Production-Degradation" = "#33a02c",
        "Ligand-Receptor" = "#6a3d9a",
        "Unknown" = "grey50"
    )
    reg_linetypes <- c("Activating" = "solid", "Inhibiting" = "dashed")
    node_shapes <- c("Metabolite" = 21, "Protein" = 22)
    node_fills <- c("Metabolite" = "#fdb863", "Protein" = "#80b1d3")

    set.seed(123)
    plot_obj <- ggraph::ggraph(graph, layout = "fr") +
        ggraph::geom_edge_link(
            ggplot2::aes(color = interaction, linetype = regulation, width = n_matches),
            arrow = grid::arrow(length = grid::unit(2.5, "mm")),
            alpha = 0.7,
            end_cap = ggraph::circle(2.5, "mm")
        ) +
        ggraph::geom_node_point(
            ggplot2::aes(shape = node_type, fill = node_type, size = degree),
            color = "black"
        ) +
        ggplot2::scale_shape_manual(name = "Node Type", values = node_shapes) +
        ggplot2::scale_fill_manual(name = "Node Type", values = node_fills) +
        ggplot2::scale_size_continuous(name = "Degree", range = c(3.5, 10)) +
        ggraph::scale_edge_color_manual(
            name = "Interaction",
            values = edge_colors,
            na.value = "grey50"
        ) +
        ggraph::scale_edge_linetype_manual(
            name = "Regulation",
            values = reg_linetypes,
            na.value = "solid"
        ) +
        ggraph::scale_edge_width_continuous(
            name = "Matched links",
            range = c(0.4, 1.3)
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(
            title = plot_name,
            subtitle = "Metalinks metabolite-protein network"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
            legend.position = "right",
            legend.box = "vertical",
            legend.title = ggplot2::element_text(size = 9),
            legend.text = ggplot2::element_text(size = 8)
        ) +
        ggraph::geom_node_text(
            ggplot2::aes(label = label),
            repel = label_repel,
            size = 2.8,
            max.overlaps = Inf
        )

    plot_obj
}

#' @noRd
.save_metalinks_network_plot <- function(plot, save_plot, path, plot_name, width, height) {
    if (is.null(save_plot)) {
        return(character(0))
    }

    out_dir <- if (is.null(path)) "." else path
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    saved_files <- character(0)
    for (fmt in save_plot) {
        out_file <- file.path(out_dir, paste0(plot_name, ".", fmt))
        if (fmt == "svg") {
            .check_required_package("svglite")
            ggplot2::ggsave(
                filename = out_file,
                plot = plot,
                width = width,
                height = height,
                units = "in",
                limitsize = FALSE,
                device = svglite::svglite
            )
        } else if (fmt == "png") {
            ggplot2::ggsave(
                filename = out_file,
                plot = plot,
                width = width,
                height = height,
                units = "in",
                limitsize = FALSE
            )
        } else if (fmt == "pdf") {
            ggplot2::ggsave(
                filename = out_file,
                plot = plot,
                width = width,
                height = height,
                units = "in",
                limitsize = FALSE
            )
        }
        saved_files <- c(saved_files, out_file)
    }

    saved_files
}

#' @noRd
.check_required_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package `", pkg, "` is required but not installed.")
    }
}
