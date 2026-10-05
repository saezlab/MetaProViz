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
# Networks of measured metabolites and prior knowledge
#

#' Network of metabolites and their prior knowledge terms
#'
#' Plot measured metabolites together with the prior knowledge (PK) terms they
#' are linked to, e.g. proteins from [metsigdb_metalinks()], pathways from
#' [metsigdb_kegg()] or cancer types from [metsigdb_macdb()]. Every link
#' between a metabolite and a term in `input_pk` becomes an edge.
#'
#' The metabolites in `feature_metadata` are matched to `input_pk` by ID. A
#' cell may hold several IDs separated by `id_sep`. With `id_type = "HMDB"`,
#' HMDB IDs are brought into the same 7-digit form on both sides (e.g.
#' "HMDB00123" and "HMDB0000123" match); all other IDs are matched exactly.
#' Use [translate_id()] beforehand if the data and the PK use different ID
#' types.
#'
#' Node colours and sizes and edge colours, line types and widths can be
#' mapped to columns via `metadata_info`. Numeric columns get a continuous
#' scale (diverging around 0 if the values have both signs, e.g. Log2FC),
#' other columns a discrete one. If several rows map to the same node or edge,
#' numeric values are averaged and otherwise the first non-missing value is
#' used. Nodes without a size value are drawn at the smallest size.
#'
#' Edges are undirected unless `metadata_info` contains `EdgeDirection`. This
#' column must hold `"to_term"` (arrow from metabolite to term) or
#' `"to_metabolite"` (arrow from term to metabolite); edges with other values
#' stay undirected. [metsigdb_metalinks()] provides such a column,
#' `direction`, together with the interaction category `interaction`.
#'
#' @param feature_metadata Data frame with one row per measured feature,
#'     holding the ID column and optional label and attribute columns.
#' @param input_pk Prior knowledge in long format, i.e. one metabolite ID per
#'     row, such as the tables returned by the `metsigdb_*()` functions.
#' @param metadata_info Named character vector mapping roles to column names:
#'     \describe{
#'       \item{InputID}{Required. ID column in `feature_metadata`.}
#'       \item{PriorID}{Required. ID column in `input_pk`.}
#'       \item{PriorTerm}{Required. Column in `input_pk` whose values become
#'       the term nodes, e.g. "gene_symbol" or "term".}
#'       \item{InputLabel}{Metabolite label column in `feature_metadata`.
#'       Without it, or for missing labels, the IDs are used.}
#'       \item{MetaboliteColor, MetaboliteSize}{Columns in `feature_metadata`
#'       for the fill and size of metabolite nodes. Size must be numeric.}
#'       \item{TermColor, TermSize}{Columns in `term_metadata`, or else in
#'       `input_pk`, for the fill and size of term nodes. Size must be
#'       numeric.}
#'       \item{EdgeColor, EdgeLinetype, EdgeWidth}{Columns in `input_pk` for
#'       the colour, line type and width of edges. Width must be numeric;
#'       by default it is the number of `input_pk` rows behind an edge.}
#'       \item{EdgeDirection}{Column in `input_pk` holding the edge direction
#'       (see Details).}
#'     }
#' @param term_metadata \emph{Optional: } Data frame with one row per term,
#'     e.g. the result of [cluster_ora()], containing the `PriorTerm` column.
#'     Used for `TermColor` and `TermSize`. \strong{Default = NULL}
#' @param id_type \emph{Optional: } `"HMDB"` to normalise HMDB IDs before
#'     matching. Any other value (e.g. "KEGG", "PubChem") matches IDs exactly.
#'     \strong{Default = "HMDB"}
#' @param id_sep \emph{Optional: } Separator of multiple IDs in one cell of
#'     the `InputID` column. \strong{Default = ";"}
#' @param label_mode \emph{Optional: } `"reduced"` labels only nodes with a
#'     degree of at least `label_degree_min`; `"all"` labels all nodes.
#'     \strong{Default = "reduced"}
#' @param label_max_chars \emph{Optional: } Labels longer than this are
#'     shortened with "...". \strong{Default = 20}
#' @param label_degree_min \emph{Optional: } Minimum degree of labelled nodes
#'     if `label_mode = "reduced"`. \strong{Default = 2}
#' @param label_repel \emph{Optional: } If TRUE, labels are repelled from
#'     each other. \strong{Default = TRUE}
#' @param seed \emph{Optional: } Seed for random layouts such as "fr". With
#'     `NULL` the layout changes between calls. The global random seed is not
#'     changed. \strong{Default = NULL}
#' @param layout \emph{Optional: } Graph layout passed to [ggraph::ggraph()],
#'     e.g. "fr" (force-directed), "kk" or "stress". \strong{Default = "fr"}
#' @param plot_name \emph{Optional: } Plot title and name of the saved files.
#'     \strong{Default = "PK_Network"}
#' @param save_plot \emph{Optional: } File type of the saved plot: "svg",
#'     "pdf", "png" or NULL. \strong{Default = "svg"}
#' @param save_table \emph{Optional: } File type of the saved tables: "csv",
#'     "xlsx", "txt" or NULL. \strong{Default = NULL}
#' @param print_plot \emph{Optional: } If TRUE, the plot is printed.
#'     \strong{Default = TRUE}
#' @param path \emph{Optional: } Path to the folder the results are saved in.
#'     \strong{Default = NULL}
#' @param plot_width,plot_height,plot_unit \emph{Optional: } Size of the saved
#'     plot. \strong{Default = 25 x 20 cm}
#'
#' @return A list with
#'     \item{DF}{List of `edges` (one row per metabolite-term link with the
#'     mapped attributes and `n_pk_rows`), `nodes` (one row per node with
#'     `node_type`, `degree` and the mapped attributes), `matched_features`
#'     and `unmatched_features` (the split and normalised feature IDs with or
#'     without a match in `input_pk`).}
#'     \item{Plot}{List with the ggraph plot `pk_network`, or `NULL` if no
#'     feature matched `input_pk`.}
#'
#' @seealso [viz_shared_pk_network()] to connect metabolites by the terms
#'     they share; [cluster_pk()] to connect terms by the metabolites they
#'     share.
#'
#' @examples
#' # Biocrates amino acids and the transporters they interact with
#' data(biocrates_features)
#' amino_acids <- biocrates_features |>
#'     dplyr::filter(Class == "Aminoacids") |>
#'     dplyr::select(TrivialName, HMDB)
#' aa_hmdb <- trimws(unlist(strsplit(amino_acids$HMDB, ",")))
#'
#' transporters <- metsigdb_metalinks(
#'     hmdb_ids = aa_hmdb,
#'     save_table = NULL,
#'     exclude_metabolites = NULL
#' ) |>
#'     dplyr::filter(interaction_family == "Transporter-metabolite")
#'
#' network <- viz_pk_network(
#'     feature_metadata = amino_acids,
#'     input_pk = transporters,
#'     metadata_info = c(
#'         InputID = "HMDB",
#'         InputLabel = "TrivialName",
#'         PriorID = "hmdb",
#'         PriorTerm = "gene_symbol",
#'         EdgeColor = "interaction",
#'         EdgeDirection = "direction"
#'     ),
#'     id_sep = ",",
#'     seed = 1,
#'     save_plot = NULL
#' )
#'
#' @importFrom logger log_info
#' @export
viz_pk_network <- function(
    feature_metadata,
    input_pk,
    metadata_info,
    term_metadata = NULL,
    id_type = "HMDB",
    id_sep = ";",
    label_mode = c("reduced", "all"),
    label_max_chars = 20,
    label_degree_min = 2,
    label_repel = TRUE,
    seed = NULL,
    layout = "fr",
    plot_name = "PK_Network",
    save_plot = "svg",
    save_table = NULL,
    print_plot = TRUE,
    path = NULL,
    plot_width = 25,
    plot_height = 20,
    plot_unit = "cm"
) {
    metaproviz_init()
    label_mode <- match.arg(label_mode)

    check_param_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = metadata_info,
        term_metadata = term_metadata,
        id_type = id_type,
        id_sep = id_sep,
        label_max_chars = label_max_chars,
        label_degree_min = label_degree_min,
        label_repel = label_repel,
        seed = seed,
        plot_name = plot_name,
        save_plot = save_plot,
        save_table = save_table,
        print_plot = print_plot,
        path = path
    )

    log_info("viz_pk_network: Network of metabolites and prior knowledge terms")

    matches <- .pk_associations(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = metadata_info,
        id_type = id_type,
        id_sep = id_sep
    )

    edges <- dplyr::tibble()
    nodes <- dplyr::tibble()
    plot_obj <- NULL
    if (nrow(matches$associations) == 0L) {
        warning("No feature in `feature_metadata` matched `input_pk`.")
    } else {
        edges <- .pk_network_edges(matches$associations, metadata_info)
        nodes <- .pk_network_nodes(
            edges = edges,
            matches = matches,
            feature_metadata = feature_metadata,
            term_metadata = term_metadata,
            metadata_info = metadata_info
        )
        nodes$label <- .format_network_labels(
            labels = nodes$name,
            degrees = nodes$degree,
            label_mode = label_mode,
            label_max_chars = label_max_chars,
            label_degree_min = label_degree_min
        )
        plot_obj <- .make_pk_network_plot(
            nodes = nodes,
            edges = edges,
            metadata_info = metadata_info,
            plot_name = plot_name,
            label_repel = label_repel,
            seed = seed,
            layout = layout
        )
    }

    DF <- list(
        edges = dplyr::select(edges, -dplyr::any_of(c(".directed", ".to_metabolite"))),
        nodes = nodes,
        matched_features = matches$matched_features,
        unmatched_features = matches$unmatched_features
    )
    Plot <- list(pk_network = plot_obj)

    .save_pk_network(
        DF = DF,
        Plot = Plot,
        save_plot = save_plot,
        save_table = save_table,
        print_plot = print_plot,
        path = path,
        plot_name = plot_name,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_unit = plot_unit
    )

    invisible(list(DF = DF, Plot = Plot))
}


#' Network of metabolites connected by shared prior knowledge terms
#'
#' Connect measured metabolites that are linked to the same prior knowledge
#' (PK) terms, e.g. metabolites binding the same receptors in
#' [metsigdb_metalinks()], taking part in the same pathways in
#' [metsigdb_kegg()] or reported in the same cancer types in
#' [metsigdb_macdb()]. Node size and the number beneath each label show how
#' many terms a metabolite is linked to; edge width and labels show the
#' similarity of two metabolites.
#'
#' Metabolites are matched to `input_pk` as in [viz_pk_network()], so the
#' same `metadata_info` can be passed to both functions. Only `InputID`,
#' `InputLabel`, `PriorID`, `PriorTerm`, `MetaboliteColor` and
#' `MetaboliteSize` are used here; other entries are ignored.
#'
#' The raw number of shared terms favours metabolites with many terms. Use
#' `similarity = "jaccard"` to compare the term profiles of two metabolites
#' as a whole, or `similarity = "overlap_coefficient"` to see whether the
#' terms of one metabolite are mostly contained in those of the other. The
#' coefficients are calculated as in [cluster_pk()].
#'
#' To connect terms by the metabolites they share instead, use
#' [cluster_pk()]: with `input_format = "enrichment"` it takes an enrichment
#' result and only uses the measured metabolites behind each term.
#'
#' @inheritParams viz_pk_network
#' @param similarity \emph{Optional: } Edge weight: `"shared"` (number of
#'     shared terms), `"jaccard"` (shared terms divided by the terms of either
#'     metabolite) or `"overlap_coefficient"` (shared terms divided by the
#'     terms of the metabolite with fewer terms). \strong{Default = "shared"}
#' @param threshold \emph{Optional: } Minimum `similarity` for two
#'     metabolites to be connected, as in [cluster_pk()]. Metabolites without
#'     shared terms are never connected. For `similarity = "shared"` this is
#'     the minimum number of shared terms. \strong{Default = 0}
#' @param label_mode \emph{Optional: } `"all"` labels all metabolites;
#'     `"reduced"` labels only metabolites connected to at least
#'     `label_degree_min` other metabolites. \strong{Default = "all"}
#' @param label_degree_min \emph{Optional: } Minimum number of connected
#'     metabolites if `label_mode = "reduced"`. \strong{Default = 1}
#' @param plot_name \emph{Optional: } Plot title and name of the saved files.
#'     \strong{Default = "Shared_PK_Network"}
#' @param layout \emph{Optional: } Graph layout passed to [ggraph::ggraph()],
#'     e.g. "stress", "fr" (force-directed) or "kk". "stress" places
#'     unconnected metabolites next to the network instead of far away.
#'     \strong{Default = "stress"}
#'
#' @return A list with
#'     \item{DF}{List of `edges` (one row per connected metabolite pair with
#'     `shared`, `jaccard`, `overlap_coefficient`, the plotted `weight` and
#'     the `shared_terms`), `nodes` (one row per metabolite with `n_terms`
#'     and the mapped attributes), `associations` (the metabolite-term
#'     pairs the network is based on), `matched_features` and
#'     `unmatched_features`.}
#'     \item{Plot}{List with the ggraph plot `shared_pk_network`, or `NULL` if
#'     no feature matched `input_pk`.}
#'
#' @seealso [viz_pk_network()] to show the terms themselves.
#'
#' @examples
#' # Biocrates amino acids connected by the transporters they share
#' data(biocrates_features)
#' amino_acids <- biocrates_features |>
#'     dplyr::filter(Class == "Aminoacids") |>
#'     dplyr::select(TrivialName, HMDB)
#' aa_hmdb <- trimws(unlist(strsplit(amino_acids$HMDB, ",")))
#'
#' transporters <- metsigdb_metalinks(
#'     hmdb_ids = aa_hmdb,
#'     save_table = NULL,
#'     exclude_metabolites = NULL
#' ) |>
#'     dplyr::filter(interaction_family == "Transporter-metabolite")
#'
#' shared <- viz_shared_pk_network(
#'     feature_metadata = amino_acids,
#'     input_pk = transporters,
#'     metadata_info = c(
#'         InputID = "HMDB",
#'         InputLabel = "TrivialName",
#'         PriorID = "hmdb",
#'         PriorTerm = "gene_symbol"
#'     ),
#'     similarity = "jaccard",
#'     id_sep = ",",
#'     seed = 1,
#'     save_plot = NULL
#' )
#'
#' @importFrom logger log_info
#' @export
viz_shared_pk_network <- function(
    feature_metadata,
    input_pk,
    metadata_info,
    similarity = c("shared", "jaccard", "overlap_coefficient"),
    threshold = 0,
    id_type = "HMDB",
    id_sep = ";",
    label_mode = c("all", "reduced"),
    label_max_chars = 20,
    label_degree_min = 1,
    label_repel = TRUE,
    seed = NULL,
    layout = "stress",
    plot_name = "Shared_PK_Network",
    save_plot = "svg",
    save_table = NULL,
    print_plot = TRUE,
    path = NULL,
    plot_width = 25,
    plot_height = 20,
    plot_unit = "cm"
) {
    metaproviz_init()
    similarity <- match.arg(similarity)
    label_mode <- match.arg(label_mode)

    ignored <- intersect(
        names(metadata_info),
        c("TermColor", "TermSize", "EdgeColor", "EdgeLinetype", "EdgeWidth", "EdgeDirection")
    )
    if (length(ignored) > 0L) {
        log_info(
            "viz_shared_pk_network: Ignoring metadata_info entries: %s",
            paste(ignored, collapse = ", ")
        )
        metadata_info <- metadata_info[setdiff(names(metadata_info), ignored)]
    }

    check_param_pk_network(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = metadata_info,
        term_metadata = NULL,
        id_type = id_type,
        id_sep = id_sep,
        label_max_chars = label_max_chars,
        label_degree_min = label_degree_min,
        label_repel = label_repel,
        seed = seed,
        plot_name = plot_name,
        save_plot = save_plot,
        save_table = save_table,
        print_plot = print_plot,
        path = path
    )
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        is.na(threshold) || threshold < 0) {
        stop("`threshold` must be a single number greater than or equal to 0.")
    }

    log_info("viz_shared_pk_network: Network of metabolites sharing prior knowledge terms")

    matches <- .pk_associations(
        feature_metadata = feature_metadata,
        input_pk = input_pk,
        metadata_info = metadata_info,
        id_type = id_type,
        id_sep = id_sep
    )

    associations <- matches$associations |>
        dplyr::distinct(metabolite = .data$.metabolite, term = .data$.term)

    edges <- dplyr::tibble()
    nodes <- dplyr::tibble()
    plot_obj <- NULL
    if (nrow(associations) == 0L) {
        warning("No feature in `feature_metadata` matched `input_pk`.")
    } else {
        shared <- .shared_pk_network(
            associations = associations,
            similarity = similarity,
            threshold = threshold
        )
        edges <- shared$edges
        nodes <- .add_metabolite_attributes(
            nodes = shared$nodes,
            matched_features = matches$matched_features,
            feature_metadata = feature_metadata,
            metadata_info = metadata_info
        )
        nodes$label <- .format_network_labels(
            labels = nodes$name,
            degrees = nodes$degree,
            label_mode = label_mode,
            label_max_chars = label_max_chars,
            label_degree_min = label_degree_min
        )
        nodes$label <- ifelse(
            is.na(nodes$label),
            NA_character_,
            paste0(nodes$label, "\n(", nodes$n_terms, ")")
        )
        plot_obj <- .make_shared_pk_network_plot(
            nodes = nodes,
            edges = edges,
            metadata_info = metadata_info,
            similarity = similarity,
            plot_name = plot_name,
            label_repel = label_repel,
            seed = seed,
            layout = layout
        )
    }

    DF <- list(
        edges = edges,
        nodes = nodes,
        associations = associations,
        matched_features = matches$matched_features,
        unmatched_features = matches$unmatched_features
    )
    Plot <- list(shared_pk_network = plot_obj)

    .save_pk_network(
        DF = DF,
        Plot = Plot,
        save_plot = save_plot,
        save_table = save_table,
        print_plot = print_plot,
        path = path,
        plot_name = plot_name,
        plot_width = plot_width,
        plot_height = plot_height,
        plot_unit = plot_unit
    )

    invisible(list(DF = DF, Plot = Plot))
}


##
## Network tables
##

#' Edges of the metabolite-term network
#'
#' One edge per metabolite, term and combination of the discrete edge
#' attributes (colour, line type, direction). Numeric attributes are averaged.
#'
#' @noRd
.pk_network_edges <- function(associations, metadata_info) {
    edge_keys <- intersect(
        c("EdgeColor", "EdgeLinetype", "EdgeWidth", "EdgeDirection"),
        names(metadata_info)
    )
    edge_columns <- unique(unname(metadata_info[edge_keys]))
    group_columns <- edge_columns[!vapply(
        edge_columns,
        function(column) is.numeric(associations[[column]]),
        logical(1)
    )]
    numeric_columns <- setdiff(edge_columns, group_columns)

    edges <- associations |>
        dplyr::group_by(
            metabolite = .data$.metabolite,
            term = .data$.term,
            dplyr::across(dplyr::all_of(group_columns))
        ) |>
        dplyr::summarise(
            dplyr::across(dplyr::all_of(numeric_columns), .summarise_attribute),
            n_pk_rows = dplyr::n(),
            matched_ids = paste(unique(.data$.matched_id), collapse = ";"),
            .groups = "drop"
        )

    edges$.directed <- FALSE
    edges$.to_metabolite <- FALSE
    if ("EdgeDirection" %in% names(metadata_info)) {
        direction <- as.character(edges[[metadata_info[["EdgeDirection"]]]])
        edges$.directed <- !is.na(direction) & direction %in% c("to_term", "to_metabolite")
        edges$.to_metabolite <- !is.na(direction) & direction == "to_metabolite"
    }
    edges
}

#' Nodes of the metabolite-term network
#'
#' @noRd
.pk_network_nodes <- function(edges, matches, feature_metadata, term_metadata, metadata_info) {
    metabolite_nodes <- dplyr::tibble(name = unique(edges$metabolite), node_type = "Metabolite")
    metabolite_nodes <- .add_metabolite_attributes(
        nodes = metabolite_nodes,
        matched_features = matches$matched_features,
        feature_metadata = feature_metadata,
        metadata_info = metadata_info
    )

    term_nodes <- dplyr::tibble(name = unique(edges$term), node_type = "Term")
    for (key in intersect(c("TermColor", "TermSize"), names(metadata_info))) {
        column <- metadata_info[[key]]
        if (column %in% colnames(term_nodes)) {
            next
        }
        use_term_metadata <- !is.null(term_metadata) && column %in% colnames(term_metadata)
        term_nodes <- .add_node_attribute(
            nodes = term_nodes,
            source = if (use_term_metadata) term_metadata else matches$associations,
            key = if (use_term_metadata) metadata_info[["PriorTerm"]] else ".term",
            column = column,
            target = column
        )
    }

    nodes <- dplyr::bind_rows(metabolite_nodes, term_nodes)
    degree <- c(table(c(
        paste0("Metabolite\r", edges$metabolite),
        paste0("Term\r", edges$term)
    )))
    nodes$degree <- unname(degree[paste0(nodes$node_type, "\r", nodes$name)])
    dplyr::relocate(nodes, dplyr::all_of(c("name", "node_type", "degree")))
}

#' Add `MetaboliteColor` and `MetaboliteSize` columns to metabolite nodes
#'
#' @noRd
.add_metabolite_attributes <- function(nodes, matched_features, feature_metadata, metadata_info) {
    features <- matched_features |>
        dplyr::distinct(.data$.feature_row_id, .data$metabolite)
    for (key in intersect(c("MetaboliteColor", "MetaboliteSize"), names(metadata_info))) {
        column <- metadata_info[[key]]
        if (column %in% colnames(nodes)) {
            next
        }
        features[[column]] <- feature_metadata[[column]][features$.feature_row_id]
        nodes <- .add_node_attribute(
            nodes = nodes,
            source = features,
            key = "metabolite",
            column = column,
            target = column
        )
    }
    nodes
}

#' Metabolite-metabolite network of shared terms
#'
#' @noRd
.shared_pk_network <- function(associations, similarity, threshold) {
    sets <- split(associations$term, associations$metabolite)
    sets <- lapply(sets, unique)

    shared <- .set_similarity(sets, "shared")
    jaccard <- .set_similarity(sets, "jaccard")
    overlap <- .set_similarity(sets, "overlap_coefficient")

    weight <- list(shared = shared, jaccard = jaccard, overlap_coefficient = overlap)[[similarity]]
    pairs <- which(upper.tri(shared) & shared > 0 & weight >= threshold, arr.ind = TRUE)
    metabolites <- rownames(shared)
    edges <- dplyr::tibble(
        from = metabolites[pairs[, 1]],
        to = metabolites[pairs[, 2]],
        shared = as.integer(shared[pairs]),
        jaccard = jaccard[pairs],
        overlap_coefficient = overlap[pairs]
    )
    edges$weight <- edges[[similarity]]
    edges$shared_terms <- vapply(
        seq_len(nrow(edges)),
        function(i) paste(sort(intersect(sets[[edges$from[i]]], sets[[edges$to[i]]])), collapse = "; "),
        character(1)
    )

    nodes <- dplyr::tibble(
        name = metabolites,
        n_terms = unname(lengths(sets[metabolites]))
    )
    nodes$degree <- vapply(
        nodes$name,
        function(metabolite) sum(edges$from == metabolite | edges$to == metabolite),
        integer(1),
        USE.NAMES = FALSE
    )

    list(nodes = nodes, edges = edges)
}


##
## Plots
##

#' @noRd
.make_pk_network_plot <- function(
    nodes,
    edges,
    metadata_info,
    plot_name,
    label_repel,
    seed = NULL,
    layout = "fr"
) {
    metadata_info <- .drop_empty_mappings(
        metadata_info,
        tables = list(
            MetaboliteColor = nodes, MetaboliteSize = nodes,
            TermColor = nodes, TermSize = nodes,
            EdgeColor = edges, EdgeLinetype = edges, EdgeWidth = edges
        )
    )
    nodes$.id <- paste(nodes$node_type, nodes$name, sep = ":")
    metabolite_ids <- paste("Metabolite", edges$metabolite, sep = ":")
    term_ids <- paste("Term", edges$term, sep = ":")
    graph_edges <- edges
    graph_edges$from <- ifelse(edges$.to_metabolite, term_ids, metabolite_ids)
    graph_edges$to <- ifelse(edges$.to_metabolite, metabolite_ids, term_ids)
    graph_edges$.width <- if ("EdgeWidth" %in% names(metadata_info)) {
        edges[[metadata_info[["EdgeWidth"]]]]
    } else {
        edges$n_pk_rows
    }
    graph_edges <- dplyr::relocate(graph_edges, "from", "to")
    graph_nodes <- dplyr::relocate(nodes, ".id")

    graph <- igraph::graph_from_data_frame(graph_edges, directed = TRUE, vertices = graph_nodes)
    term_label <- if (metadata_info[["PriorTerm"]] == "term") {
        "Term"
    } else {
        paste0("Term (", metadata_info[["PriorTerm"]], ")")
    }

    plot_obj <- .network_layout(graph, seed, layout) +
        .pk_edge_layers(edges, metadata_info) +
        .pk_edge_scales(edges, metadata_info) +
        .pk_node_layers(nodes, metadata_info, term_label) +
        do.call(
            ggraph::geom_node_text,
            c(
                list(
                    mapping = ggplot2::aes(label = .data$label),
                    repel = label_repel,
                    size = 2.8,
                    na.rm = TRUE
                ),
                if (label_repel) list(max.overlaps = Inf)
            )
        ) +
        ggplot2::labs(title = plot_name) +
        .network_theme()

    plot_obj
}

#' Drop aesthetic mappings to columns without any non-missing value
#'
#' Such a column carries no information and breaks the plot legends.
#'
#' @param metadata_info Named character vector.
#' @param tables Named list: for each aesthetic key the table holding the
#'     mapped column.
#'
#' @noRd
.drop_empty_mappings <- function(metadata_info, tables) {
    for (key in intersect(names(tables), names(metadata_info))) {
        values <- tables[[key]][[metadata_info[[key]]]]
        if (all(is.na(values))) {
            warning(
                "Column `", metadata_info[[key]], "` (metadata_info[[\"", key,
                "\"]]) has only missing values for the plotted network and is not used."
            )
            metadata_info <- metadata_info[names(metadata_info) != key]
        }
    }
    metadata_info
}

#' Edge layers with arrows for directed and plain lines for undirected edges
#'
#' @noRd
.pk_edge_layers <- function(edges, metadata_info) {
    mapping <- list(edge_width = quote(.data$.width))
    if ("EdgeColor" %in% names(metadata_info)) {
        mapping$edge_colour <- rlang::expr(.data[[!!metadata_info[["EdgeColor"]]]])
    }
    if ("EdgeLinetype" %in% names(metadata_info)) {
        mapping$edge_linetype <- rlang::expr(.data[[!!metadata_info[["EdgeLinetype"]]]])
    }
    fixed <- if ("EdgeColor" %in% names(metadata_info)) list() else list(edge_colour = "grey50")

    layers <- list()
    if (any(edges$.directed)) {
        layers <- c(layers, list(do.call(
            ggraph::geom_edge_link,
            c(
                list(
                    mapping = do.call(ggplot2::aes, c(mapping, list(filter = quote(.data$.directed)))),
                    arrow = grid::arrow(length = grid::unit(2.5, "mm"), type = "closed"),
                    end_cap = ggraph::circle(3, "mm"),
                    edge_alpha = 0.8
                ),
                fixed
            )
        )))
    }
    if (any(!edges$.directed)) {
        layers <- c(layers, list(do.call(
            ggraph::geom_edge_link,
            c(
                list(
                    mapping = do.call(ggplot2::aes, c(mapping, list(filter = quote(!.data$.directed)))),
                    edge_alpha = 0.8
                ),
                fixed
            )
        )))
    }
    layers
}

#' @noRd
.pk_edge_scales <- function(edges, metadata_info) {
    width_name <- if ("EdgeWidth" %in% names(metadata_info)) {
        metadata_info[["EdgeWidth"]]
    } else {
        "PK rows per edge"
    }
    widths <- if ("EdgeWidth" %in% names(metadata_info)) edges[[metadata_info[["EdgeWidth"]]]] else edges$n_pk_rows
    scales <- list(ggraph::scale_edge_width_continuous(
        name = width_name,
        range = c(0.4, 1.6),
        # A legend with a single width carries no information
        guide = if (length(unique(stats::na.omit(widths))) > 1L) "legend" else "none"
    ))

    if ("EdgeColor" %in% names(metadata_info)) {
        column <- metadata_info[["EdgeColor"]]
        values <- edges[[column]]
        if (is.numeric(values)) {
            scales <- c(scales, list(ggraph::scale_edge_colour_viridis(name = column, na.value = "grey70")))
        } else {
            scales <- c(scales, list(ggraph::scale_edge_colour_manual(
                name = column,
                values = .discrete_palette(values),
                na.value = "grey70"
            )))
        }
    }
    if ("EdgeLinetype" %in% names(metadata_info)) {
        scales <- c(scales, list(ggraph::scale_edge_linetype_discrete(name = metadata_info[["EdgeLinetype"]])))
    }
    scales
}

#' Node layers: one per node type, each with its own fill and size scale
#'
#' @noRd
.pk_node_layers <- function(nodes, metadata_info, term_label) {
    default_fill <- c(Metabolite = "#fdb863", Term = "#80b1d3")
    color_keys <- c(Metabolite = "MetaboliteColor", Term = "TermColor")
    size_keys <- c(Metabolite = "MetaboliteSize", Term = "TermSize")
    shared_size <- !any(size_keys %in% names(metadata_info))

    layers <- list()
    for (type in c("Metabolite", "Term")) {
        type_nodes <- nodes[nodes$node_type == type, , drop = FALSE]
        mapping <- list(
            filter = rlang::expr(.data$node_type == !!type),
            shape = quote(.data$node_type)
        )
        fixed <- list(colour = "black")

        fill_column <- if (color_keys[[type]] %in% names(metadata_info)) metadata_info[[color_keys[[type]]]]
        size_column <- if (size_keys[[type]] %in% names(metadata_info)) metadata_info[[size_keys[[type]]]] else "degree"

        if (is.null(fill_column)) {
            fixed$fill <- default_fill[[type]]
        } else {
            mapping$fill <- rlang::expr(.data[[!!fill_column]])
        }
        mapping$size <- rlang::expr(.size_or_min(.data[[!!size_column]]))

        if (type == "Term" && !is.null(fill_column) && "MetaboliteColor" %in% names(metadata_info)) {
            layers <- c(layers, list(ggnewscale::new_scale_fill()))
        }
        if (type == "Term" && !shared_size) {
            layers <- c(layers, list(ggnewscale::new_scale("size")))
        }

        layers <- c(layers, list(do.call(
            ggraph::geom_node_point,
            c(list(mapping = do.call(ggplot2::aes, mapping)), fixed)
        )))
        if (!is.null(fill_column)) {
            layers <- c(layers, list(.network_fill_scale(
                values = type_nodes[[fill_column]],
                name = fill_column,
                shape = c(Metabolite = 21, Term = 22)[[type]],
                option = c(Metabolite = "D", Term = "C")[[type]]
            )))
        }
        if (!shared_size) {
            size_name <- if (size_column == "degree") paste(type, "degree") else size_column
            layers <- c(layers, list(ggplot2::scale_size_continuous(name = size_name, range = c(3, 9))))
        }
    }
    if (shared_size) {
        layers <- c(layers, list(ggplot2::scale_size_continuous(name = "Degree", range = c(3, 9))))
    }

    key_fill <- ifelse(color_keys %in% names(metadata_info), "white", default_fill)
    c(layers, list(
        ggplot2::scale_shape_manual(
            name = "Node type",
            values = c(Metabolite = 21, Term = 22),
            labels = c(Metabolite = "Metabolite", Term = term_label)
        ),
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(fill = unname(key_fill), size = 4)))
    ))
}

#' Replace missing sizes by the smallest size so no node is dropped
#'
#' @noRd
.size_or_min <- function(x) {
    if (all(is.na(x))) {
        return(rep(1, length(x)))
    }
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x
}

#' @noRd
.make_shared_pk_network_plot <- function(
    nodes,
    edges,
    metadata_info,
    similarity,
    plot_name,
    label_repel,
    seed = NULL,
    layout = "stress"
) {
    metadata_info <- .drop_empty_mappings(
        metadata_info,
        tables = list(MetaboliteColor = nodes, MetaboliteSize = nodes)
    )
    graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

    weight_name <- c(
        shared = "Shared terms",
        jaccard = "Jaccard index",
        overlap_coefficient = "Overlap coefficient"
    )[[similarity]]
    size_column <- if ("MetaboliteSize" %in% names(metadata_info)) metadata_info[["MetaboliteSize"]] else "n_terms"
    size_name <- if (size_column == "n_terms") "Terms per metabolite" else size_column
    fill_column <- if ("MetaboliteColor" %in% names(metadata_info)) metadata_info[["MetaboliteColor"]]

    node_mapping <- list(size = rlang::expr(.size_or_min(.data[[!!size_column]])))
    node_fixed <- list(shape = 21, colour = "black")
    if (is.null(fill_column)) {
        node_fixed$fill <- "#fdb863"
    } else {
        node_mapping$fill <- rlang::expr(.data[[!!fill_column]])
    }

    edge_layers <- list()
    if (nrow(edges) > 0L) {
        edge_layers <- list(
            ggraph::geom_edge_link(
                ggplot2::aes(
                    edge_width = .data$weight,
                    label = if (similarity == "shared") .data$weight else round(.data$weight, 2)
                ),
                edge_colour = "grey45",
                edge_alpha = 0.7,
                label_colour = "black",
                label_size = 3,
                angle_calc = "along",
                label_dodge = grid::unit(2, "mm"),
                check_overlap = TRUE
            ),
            ggraph::scale_edge_width_continuous(name = weight_name, range = c(0.5, 2.5))
        )
    }

    .network_layout(graph, seed, layout) +
        edge_layers +
        do.call(
            ggraph::geom_node_point,
            c(list(mapping = do.call(ggplot2::aes, node_mapping)), node_fixed)
        ) +
        (if (is.null(fill_column)) NULL else .network_fill_scale(nodes[[fill_column]], fill_column, shape = 21)) +
        ggplot2::scale_size_continuous(name = size_name, range = c(4, 12)) +
        ggraph::geom_node_text(
            ggplot2::aes(label = .data$label),
            repel = label_repel,
            size = 3,
            na.rm = TRUE
        ) +
        ggplot2::labs(
            title = plot_name,
            subtitle = paste0(
                "Node label: metabolite (terms); edge label: ",
                tolower(weight_name)
            )
        ) +
        .network_theme()
}

#' @noRd
.network_theme <- function() {
    list(
        ggplot2::theme_void(),
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", colour = NA),
            plot.title = ggplot2::element_text(hjust = 0.5),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
            legend.position = "right",
            legend.box = "vertical",
            legend.title = ggplot2::element_text(size = 9),
            legend.text = ggplot2::element_text(size = 8)
        )
    )
}

#' Print and save the network tables and plot
#'
#' @noRd
.save_pk_network <- function(
    DF,
    Plot,
    save_plot,
    save_table,
    print_plot,
    path,
    plot_name,
    plot_width,
    plot_height,
    plot_unit
) {
    Plot <- Plot[!vapply(Plot, is.null, logical(1))]

    if (isTRUE(print_plot)) {
        for (plot_obj in Plot) {
            print(plot_obj)
        }
    }

    if (is.null(save_plot) && is.null(save_table)) {
        return(invisible(NULL))
    }

    folder <- save_path(folder_name = "PKNetwork", path = path)
    save_res(
        inputlist_df = if (is.null(save_table)) NULL else DF[vapply(DF, nrow, integer(1)) > 0L],
        inputlist_plot = if (length(Plot) == 0L) NULL else Plot,
        save_table = save_table,
        save_plot = if (length(Plot) == 0L) NULL else save_plot,
        path = folder,
        file_name = plot_name,
        core = FALSE,
        print_plot = FALSE,
        plot_height = plot_height,
        plot_width = plot_width,
        plot_unit = plot_unit
    )
    invisible(NULL)
}
