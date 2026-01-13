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

#
# Graph visualization for clustered prior knowledge
#

#' Graph visualization for clustered terms
#'
#' @param similarity_matrix Square numeric matrix of term similarity values.
#'     Row/column names must match `clusters` names.
#' @param clusters Named vector of cluster labels for each term (e.g., "cluster1").
#' @param threshold Similarity threshold used to define edges. Values below the
#'     threshold are removed. Default = 0.5.
#' @param plot_name \emph{Optional: } String added to output files of the plot.
#'     Default = "ClusterGraph".
#' @param max_nodes \emph{Optional: } Maximum nodes for plotting. If set,
#'     keeps nodes from the largest component up to this limit (by degree).
#' @param min_degree \emph{Optional: } Minimum degree filter for graph plotting.
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'     Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting
#'     plots. \strong{Default = TRUE}
#' @param path {Optional:} String which is added to the resulting folder name.
#'     \strong{default: NULL}
#'
#' @return Graph plot as a ggplot object.
#'
#' @importFrom igraph graph_from_adjacency_matrix components degree induced_subgraph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggraph theme_graph
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 labs
#' @importFrom logger log_info log_trace
#' @noRd
viz_graph <- function(
    similarity_matrix,
    clusters,
    threshold = 0.5,
    plot_name = "ClusterGraph",
    max_nodes = NULL,
    min_degree = NULL,
    save_plot = "svg",
    print_plot = TRUE,
    path = NULL
) {
    # # ------------ Create log file ----------- ##
    # metaproviz_init()
    log_info("viz_graph: Graph plot visualization")

    # ---- Input checks ----------------------------------------------------
    if (!is.matrix(similarity_matrix)) {
        stop("similarity_matrix must be a numeric matrix.")
    }
    if (is.null(rownames(similarity_matrix)) || is.null(colnames(similarity_matrix))) {
        stop("similarity_matrix must have row and column names.")
    }
    if (!is.numeric(similarity_matrix)) {
        stop("similarity_matrix must be numeric.")
    }
    if (any(rownames(similarity_matrix) != colnames(similarity_matrix))) {
        stop("similarity_matrix must be square with matching row/column names.")
    }
    if (is.null(names(clusters))) {
        stop("clusters must be a named vector with term names.")
    }
    if (!all(rownames(similarity_matrix) %in% names(clusters))) {
        stop("All similarity_matrix terms must be present in clusters.")
    }
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
        stop("threshold must be a single numeric value.")
    }
    if (threshold < 0 || threshold > 1) {
        stop("threshold must be between 0 and 1.")
    }
    if (!is.null(max_nodes) && (!is.numeric(max_nodes) || max_nodes < 1)) {
        stop("max_nodes must be NULL or a positive integer.")
    }
    if (!is.null(min_degree) && (!is.numeric(min_degree) || min_degree < 0)) {
        stop("min_degree must be NULL or a non-negative integer.")
    }
    if (!is.logical(print_plot) || length(print_plot) != 1L) {
        stop("print_plot must be TRUE or FALSE.")
    }

    save_plot_options <- c("svg", "pdf", "png")
    if (!is.null(save_plot) && !(save_plot %in% save_plot_options)) {
        stop(
            "save_plot must be one of: ",
            paste(save_plot_options, collapse = ", "),
            ", or NULL."
        )
    }

    # ---- Create adjacency ------------------------------------------------
    adj <- similarity_matrix
    adj[adj < threshold] <- 0
    diag(adj) <- 0

    g <- graph_from_adjacency_matrix(
        adj,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
    )

    # Filter by degree if requested
    if (!is.null(min_degree)) {
        keep <- which(degree(g) >= min_degree)
        g <- induced_subgraph(g, vids = keep)
    }

    # If max_nodes is set, keep nodes from largest component up to the limit
    if (!is.null(max_nodes) && igraph::vcount(g) > max_nodes) {
        comps <- components(g)
        largest_comp <- which.max(comps$csize)
        comp_nodes <- which(comps$membership == largest_comp)
        comp_nodes <- comp_nodes[order(degree(g)[comp_nodes], decreasing = TRUE)]
        keep_nodes <- comp_nodes[seq_len(min(length(comp_nodes), max_nodes))]
        g <- induced_subgraph(g, vids = keep_nodes)
    }

    # Attach cluster labels to nodes for stable color mapping
    igraph::V(g)$cluster <- clusters[igraph::V(g)$name]

    graph_plot <-
        ggraph(g, layout = "fr") +
        geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
        geom_node_point(aes(color = cluster), size = 3) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +
        labs(title = plot_name) +
        theme_graph()

    # ---- Save and return -------------------------------------------------
    folder <- NULL
    if (!is.null(save_plot)) {
        folder <- save_path(
            folder_name = "Graph",
            path = path
        )
        log_trace(paste0("viz_graph results saved at ", folder))
    }

    save_res(
        inputlist_df = NULL,
        inputlist_plot = list(graph_plot = graph_plot),
        save_table = NULL,
        save_plot = save_plot,
        path = folder,
        file_name = plot_name,
        core = FALSE,
        print_plot = FALSE,
        plot_height = 3000,
        plot_width = 9000,
        plot_unit = "px"
    )

    if (print_plot) {
        print(graph_plot)
    }

    return(graph_plot)
}
