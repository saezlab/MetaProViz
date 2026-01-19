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
#' @param plot_threshold Similarity threshold used to define edges. Values below the
#'     threshold are removed. Default = 0.5.
#' @param plot_name \emph{Optional: } String added to output files of the plot.
#'     Default = "ClusterGraph".
#' @param max_nodes \emph{Optional: } Maximum nodes for plotting. If set,
#'     keeps nodes from the largest component up to this limit (by degree).
#' @param min_degree \emph{Optional: } Minimum degree filter for graph plotting.
#' @param node_sizes \emph{Optional: } Named numeric vector of node sizes,
#'     with names matching term IDs. Values are scaled for plotting.
#' @param show_density \emph{Optional: } If TRUE, add a hull background per
#'     cluster. Default = FALSE.
#' @param seed \emph{Optional: } Random seed for graph layout reproducibility.
#'     Default = NULL.
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'     Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting
#'     plots. \strong{Default = TRUE}
#' @param path {Optional:} String which is added to the resulting folder name.
#'     \strong{default: NULL}
#'
#' @return Graph plot as a ggplot object.
#' 
#' @examples
#' 
#' # Load example data
#' d <- metsigdb_kegg()
#' 
#' # Run clustering with graph plotting
#' r <- cluster_pk(
#'     d,
#'     metadata_info = c(
#'         metabolite_column = "MetaboliteID",
#'         pathway_column = "term"
#'     ),
#'     input_format = "long",
#'     similarity = "jaccard",
#'     threshold = 0.2,
#'     clust = "community",
#'     min = 2,
#'     plot_name = "GraphExample_long_format",
#'     save_plot = "png",
#'     min_degree = 1,
#'     seed = 123,
#'     show_density = FALSE,
#'     max_nodes = 1000
#' ) 
#' 
#' # run graph plotting separately from clustering output but without node_sizes
#' viz_graph(
#'     r$similarity_matrix,
#'     r$clusters,
#'     plot_threshold = 0.5,
#'     plot_name = "ClusterGraph",
#'     max_nodes = NULL,
#'     min_degree = NULL,
#'     node_sizes = NULL,
#'     show_density = TRUE,
#'     seed = 123,
#'     save_plot = NULL,
#'     print_plot = FALSE,
#'     path = NULL
#' )
#'
#' @importFrom igraph graph_from_adjacency_matrix components degree induced_subgraph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggraph theme_graph
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 labs
#' @importFrom grDevices chull
#' @importFrom logger log_info log_trace
#' @export
viz_graph <- function(
    similarity_matrix,
    clusters,
    plot_threshold = 0.5,
    plot_name = "ClusterGraph",
    max_nodes = NULL,
    min_degree = NULL,
    node_sizes = NULL,
    show_density = FALSE,
    seed = NULL,
    save_plot = "svg",
    print_plot = TRUE,
    path = NULL,
    plot_width = 3000,
    plot_height = 2000,
    plot_unit = "px"
) {
    # NSE vs. R CMD check workaround
    cluster <- name <- node_size <- weight <- x <- y <- .data <- NULL

    # # ------------ Create log file ----------- ##
    # metaproviz_init()
    log_info("viz_graph: Graph plot visualization")

    # ---- Input checks ----------------------------------------------------
    check_param_VizGraph(
        similarity_matrix = similarity_matrix,
        clusters = clusters,
        plot_threshold = plot_threshold,
        max_nodes = max_nodes,
        min_degree = min_degree,
        node_sizes = node_sizes,
        show_density = show_density,
        print_plot = print_plot,
        save_plot = save_plot,
        seed = seed
    )

    # ---- Create adjacency ------------------------------------------------
    adj <- similarity_matrix
    adj[adj < plot_threshold] <- 0
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

    if (!is.null(seed)) {
        set.seed(seed)
    }
    layout <- ggraph::create_layout(g, layout = "fr")

    node_size_vec <- NULL
    node_size_label <- NULL
    if (!is.null(node_sizes)) {
        node_size_vec <- node_sizes[igraph::V(g)$name]
        node_size_label <- attr(node_sizes, "label")
        if (!all(is.na(node_size_vec))) {
            med <- stats::median(node_size_vec, na.rm = TRUE)
            node_size_vec[is.na(node_size_vec)] <- med
        } else {
            node_size_vec <- NULL
        }
    }
    if (!is.null(node_size_vec)) {
        layout$node_size <- node_size_vec[layout$name]
    }

    hull_layer <- NULL
    if (isTRUE(show_density)) {
        hull_data <- layout[layout$cluster != "None", , drop = FALSE]
        hull_list <- lapply(
            split(hull_data, hull_data$cluster),
            function(df) {
                if (nrow(df) < 3) {
                    return(NULL)
                }
                idx <- chull(df$x, df$y)
                df <- df[idx, , drop = FALSE]
                rbind(df, df[1, , drop = FALSE])
            }
        )
        hull_df <- do.call(rbind, hull_list)
        if (!is.null(hull_df) && nrow(hull_df) > 0) {
            hull_layer <- ggplot2::geom_polygon(
                data = hull_df,
                aes(x = x, y = y, fill = cluster, group = cluster),
                alpha = 0.15,
                color = NA,
                show.legend = FALSE
            )
        }
    }

    cluster_levels <- sort(unique(layout$cluster))
    cluster_palette <- grDevices::hcl.colors(length(cluster_levels), "Dynamic")
    names(cluster_palette) <- cluster_levels
    if ("None" %in% cluster_levels) {
        cluster_palette["None"] <- "grey70"
    }

    if (is.null(node_size_vec)) {
        graph_plot <-
            ggraph(layout) +
            hull_layer +
            geom_edge_link(aes(edge_alpha = weight, edge_width = weight)) +
            ggraph::scale_edge_width(range = c(0.2, 1.2), name = "Edge weight (similarity)") +
            ggraph::scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
            geom_node_point(aes(color = cluster), size = 3) +
            ggplot2::scale_color_manual(values = cluster_palette, drop = FALSE) +
            ggplot2::scale_fill_manual(values = cluster_palette, drop = FALSE) +
            geom_node_text(aes(label = name), repel = TRUE, size = 3) +
            labs(title = plot_name) +
            theme_graph(base_family = "sans")
    } else {
        graph_plot <-
            ggraph(layout) +
            hull_layer +
            geom_edge_link(aes(edge_alpha = weight, edge_width = weight)) +
            ggraph::scale_edge_width(range = c(0.2, 1.2), name = "Edge weight (similarity)") +
            ggraph::scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
            geom_node_point(aes(color = cluster, size = node_size)) +
            ggplot2::scale_color_manual(values = cluster_palette, drop = FALSE) +
            ggplot2::scale_fill_manual(values = cluster_palette, drop = FALSE) +
            ggplot2::scale_size(
                range = c(1, 8),
                name = ifelse(is.null(node_size_label), "Node size", node_size_label)
            ) +
            geom_node_text(aes(label = name), repel = TRUE, size = 3) +
            labs(title = plot_name) +
            theme_graph(base_family = "sans")
    }

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
        plot_height = plot_height,
        plot_width = plot_width,
        plot_unit = plot_unit
    )

    if (print_plot) {
        print(graph_plot)
    }

    return(graph_plot)
}
