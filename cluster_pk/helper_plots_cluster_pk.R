##
## Plotting helper for cluster_pk results
##

#' Plot helper for cluster_pk results
#'
#' Generates a set of plots from a `cluster_pk` result:
#' similarity heatmap with dendrograms, cluster size bar plot, cluster size
#' treemap, similarity distribution with threshold line, graph view, and
#' dendrogram (if hierarchical).
#'
#' @param res Result list returned by `cluster_pk`.
#' @param threshold Similarity cutoff to show on plots (defaults to
#'     `res$threshold` if present).
#' @param hclust_method Linkage method if a dendrogram needs to be recomputed
#'     (defaults to `res$hclust_method` if present).
#' @param max_nodes Optional maximum nodes for the graph plot. If set, keeps
#'     nodes from the largest component up to this limit (by degree).
#' @param min_degree Optional minimum degree filter for graph plotting.
#'
#' @return A list with elements:
#'     \item{heatmap}{pheatmap object of the similarity matrix.}
#'     \item{cluster_bar}{ggplot bar plot of cluster sizes.}
#'     \item{cluster_treemap}{ggplot treemap of cluster sizes.}
#'     \item{sim_distribution}{ggplot histogram/density of similarities.}
#'     \item{graph}{ggraph plot of the thresholded graph.}
#'     \item{dendrogram}{ggplot dendrogram (hierarchical only), else NULL.}
#' @importFrom pheatmap pheatmap
#' @importFrom ggplot2 ggplot aes geom_col geom_histogram geom_density geom_text
#' @importFrom ggplot2 scale_y_continuous scale_fill_brewer coord_flip labs theme_bw
#' @importFrom ggplot2 theme geom_vline guides guide_legend
#' @importFrom treemapify geom_treemap geom_treemap_text
#' @importFrom dplyr mutate arrange desc
#' @importFrom igraph graph_from_adjacency_matrix components degree induced_subgraph
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#' @importFrom ggraph theme_graph
#' @importFrom stats as.dist hclust
#' @importFrom ggdendro dendro_data
#' @noRd
helper_plots_cluster_pk <- function(
        res,
        threshold = NULL,
        hclust_method = NULL,
        max_nodes = NULL,
        min_degree = NULL
    ) {

    # ---- Resolve parameters ---------------------------------------------
    if (is.null(threshold) && !is.null(res$threshold)) {
        threshold <- res$threshold
    }
    if (is.null(threshold)) {
        stop("threshold must be provided (or present in res$threshold).")
    }

    if (is.null(hclust_method) && !is.null(res$hclust_method)) {
        hclust_method <- res$hclust_method
    } 
    if (is.null(hclust_method)) {
        hclust_method <- "average"
    }

    sim <- res$similarity_matrix
    dist_mat <- res$distance_matrix
    cluster_labels <- res$clusters

    # ---- Similarity heatmap with dendrograms ----------------------------
    heatmap_plot <-
        pheatmap(
            sim,
            clustering_method = hclust_method,
            main = "Similarity heatmap",
            legend = TRUE,
            border_color = NA,
            show_rownames = FALSE,
            show_colnames = FALSE
        )

    # ---- Cluster size bar plot ------------------------------------------
    cs <- res$cluster_summary %>%
        arrange(desc(n_terms)) %>%
        mutate(cluster = factor(cluster, levels = cluster))

    cluster_bar <-
        ggplot(cs, aes(x = cluster, y = n_terms, fill = cluster)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        labs(
            title = "Cluster sizes",
            x = "Cluster",
            y = "Number of terms"
        ) +
        theme_bw()

    # ---- Cluster size treemap -------------------------------------------
    cluster_treemap <-
        ggplot(cs, aes(
            area = n_terms,
            fill = cluster,
            label = cluster
        )) +
        geom_treemap(show.legend = FALSE) +
        geom_treemap_text(
            colour = "black",
            place = "centre",
            grow = TRUE
        ) +
        labs(title = "Cluster size treemap")

    # ---- Similarity distribution ----------------------------------------
    # Upper triangle without diagonal
    upper_idx <- upper.tri(sim, diag = FALSE)
    sim_vals <- sim[upper_idx]
    sim_df <- data.frame(sim = sim_vals)

    sim_distribution <-
        ggplot(sim_df, aes(x = sim)) +
        geom_histogram(
            bins = 50,
            fill = "#3182bd",
            color = "white",
            alpha = 0.8
        ) +
        geom_vline(
            xintercept = threshold,
            color = "red",
            linetype = "dashed"
        ) +
        labs(
            title = "Similarity distribution",
            x = "Similarity",
            y = "Count"
        ) +
        theme_bw()

    # ---- Graph view -----------------------------------------------------
    # Build adjacency from thresholded similarity
    adj <- sim
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
    if (!is.null(max_nodes) && vcount(g) > max_nodes) {
        comps <- components(g)
        largest_comp <- which.max(comps$csize)
        comp_nodes <- which(comps$membership == largest_comp)
        # order by degree within the component
        comp_nodes <- comp_nodes[order(degree(g)[comp_nodes], decreasing = TRUE)]
        keep_nodes <- comp_nodes[seq_len(min(length(comp_nodes), max_nodes))]
        g <- induced_subgraph(g, vids = keep_nodes)
    }

    # Attach cluster labels to nodes for stable color mapping
    igraph::V(g)$cluster <- cluster_labels[igraph::V(g)$name]

    graph_plot <-
        ggraph(g, layout = "fr") +
        geom_edge_link(aes(edge_alpha = weight), show.legend = FALSE) +
        geom_node_point(aes(color = cluster), size = 3) +
        geom_node_text(aes(label = name), repel = TRUE, size = 3) +
        labs(title = "Graph view (thresholded similarity)") +
        theme_graph()

    # ---- Dendrogram -----------------------------------------------------
    if (!is.null(res$hclust)) {
        hc <- res$hclust
    } else {
        hc <- hclust(as.dist(dist_mat), method = hclust_method)
    }
    dendro <- ggdendro::dendro_data(hc, type = "rectangle")
    cut_height <- 1 - threshold
    dendrogram_plot <-
        ggplot() +
        geom_segment(
            data = dendro$segments,
            aes(x = x, y = y, xend = xend, yend = yend)
        ) +
        geom_segment(
            aes(
                x = min(dendro$labels$x),
                xend = max(dendro$labels$x),
                y = cut_height,
                yend = cut_height
            ),
            color = "red",
            linetype = "dashed"
        ) +
        theme_bw() +
        labs(title = "Dendrogram (hierarchical)")

    plots <- list(
        heatmap = heatmap_plot,
        cluster_bar = cluster_bar,
        cluster_treemap = cluster_treemap,
        sim_distribution = sim_distribution,
        graph = graph_plot,
        dendrogram = dendrogram_plot
    )

    return(plots)
}
