##
## Cluster prior knowledge by pathway overlap
##

#' Cluster terms in prior knowledge by set overlap
#'
#' @param data Long data frame with one ID per row, or enrichment-style table
#'     with a delimited metabolite list per term (see input_format).
#' @param metadata_info List with entries `metabolite_column` (metabolite ID
#'     column or delimited list column) and `pathway_column` (term column).
#'     Defaults to `c(metabolite_column = "MetaboliteID", pathway_column = "term")`.
#' @param similarity Similarity measure between term ID sets. Options:
#'     "jaccard" (default), "overlap_coefficient", or "correlation".
#'     Jaccard similarity is |A ∩ B| / |A ∪ B|. Overlap coefficient is
#'     |A ∩ B| / min(|A|, |B|). Jaccard is stricter for large sets, while
#'     overlap_coefficient is more permissive for nested sets.
#' @param correlation_method Correlation method when `similarity = "correlation"`.
#'     One of "pearson", "spearman", "kendall". Ignored otherwise.
#' @param input_format Input format of `data`. Use "long" for one ID per row
#'     (default) or "enrichment" for one term per row with a delimited ID list.
#'     The `metabolite_column` entry in metadata_info is interpreted accordingly.
#' @param delimiter Delimiter for metabolite ID lists when input_format =
#'     "enrichment". Ignored for input_format = "long". Default = "/".
#' @param threshold Similarity cutoff for keeping edges (applies to all
#'     clustering modes). Default = 0.5.
#' @param clust Clustering strategy: "components" (connected components on
#'     thresholded unweighted graph), "community" (Louvain on thresholded
#'     weighted graph), or "hierarchical" (hclust on distance = 1 - similarity).
#' @param hclust_method Linkage method for hierarchical clustering. One of
#'     "average" (default), "single", "complete", "ward.D", "ward.D2",
#'     "mcquitty", "median", "centroid". Used only when clust = "hierarchical".
#' @param min Minimum cluster size; smaller clusters are relabeled to "None".
#'     Default = 2.
#' @param plot_name \emph{Optional: } String added to output files of the plot.
#'     Default = "ClusterGraph".
#' @param max_nodes \emph{Optional: } Maximum nodes for plotting. If set,
#'     keeps nodes from the largest component up to this limit (by degree).
#'     Used only for the graph plot. Default = 10000.
#' @param min_degree \emph{Optional: } Minimum degree filter for graph plotting.
#'     Used only for the graph plot. Default = 1.
#' @param node_size_column \emph{Optional: } Numeric column name from `data`
#'     used to scale node sizes in the graph. Aggregated per term (mean) when
#'     multiple rows map to the same term. Default = NULL.
#' @param show_density \emph{Optional: } If TRUE, add a hull background per
#'     cluster to the graph. Default = FALSE.
#' @param seed \emph{Optional: } Random seed for graph layout reproducibility.
#'     Default = NULL.
#' @param save_plot \emph{Optional: } Select the file type of output plots.
#'     Options are svg, pdf, png or NULL. \strong{Default = "svg"}
#' @param print_plot \emph{Optional: } If TRUE prints an overview of resulting
#'     plots. \strong{Default = FALSE}
#' @param path {Optional:} String which is added to the resulting folder name.
#'     \strong{default: NULL}
#'
#' @return A list with:
#'     \item{data}{Input data with a `cluster` column added.}
#'     \item{cluster_summary}{Summary of cluster sizes and percentages.}
#'     \item{clusters}{Named vector of term -> cluster assignment.}
#'     \item{similarity_matrix}{Term-by-term similarity matrix.}
#'     \item{distance_matrix}{Term-by-term distance matrix (1 - similarity).}
#'     \item{graph_plot}{Graph plot returned by viz_graph.}
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
#'     show_density = TRUE,
#'     max_nodes = 1000
#' ) 
#' 
#' print(head(r$cluster_summary))
#' print(r$graph_plot)
#' 
#' ## add an example for an enrichment format result
#'
#' @importFrom dplyr group_by summarize ungroup mutate select left_join
#' @importFrom dplyr across n distinct filter tibble arrange
#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain
#' @importFrom stats cor as.dist hclust cutree
#' @importFrom rlang sym !! .data
#' @importFrom logger log_trace log_warn
#' @noRd
cluster_pk <- function(
    data,
    metadata_info = c(
        metabolite_column = "MetaboliteID",
        pathway_column = "term"
    ),
    similarity = c("jaccard", "overlap_coefficient", "correlation"),
    correlation_method = "pearson",
    input_format = c("long", "enrichment"),
    delimiter = "/",
    threshold = 0.5,
    clust = c("components", "community", "hierarchical"),
    hclust_method = "average",
    min = 2,
    plot_name = "ClusterGraph",
    max_nodes = 10000,
    min_degree = 1,
    node_size_column = NULL,
    show_density = FALSE,
    seed = NULL,
    save_plot = "png",
    print_plot = FALSE,
    path = NULL
) {

    # NSE vs. R CMD check workaround
    cluster <- .data <- NULL

    # ---- Input checks ----------------------------------------------------
    similarity <- match.arg(similarity)
    input_format <- match.arg(input_format)
    clust <- match.arg(clust)

    if (clust == "hierarchical") {
        hclust_method <- match.arg(
            hclust_method,
            c(
                "average", "single", "complete",
                "ward.D", "ward.D2",
                "mcquitty", "median", "centroid"
            )
        )
    }

    check_param_cluster_pk(
        data = data,
        metadata_info = metadata_info,
        input_format = input_format,
        delimiter = delimiter,
        threshold = threshold,
        min = min,
        node_size_column = node_size_column,
        show_density = show_density,
        seed = seed
    )

    id_col <- metadata_info[["metabolite_column"]]
    term_col <- metadata_info[["pathway_column"]]
    min <- as.integer(min)

    if (similarity == "correlation") {
        correlation_method <-
            match.arg(correlation_method, c("pearson", "spearman", "kendall"))
    }

    # ---- Build term -> ID list ------------------------------------------
    if (input_format == "long") {
        # Drop duplicated term-ID rows to avoid inflating overlaps
        data <- dplyr::distinct(data, !!sym(term_col), !!sym(id_col), .keep_all = TRUE)

        # Warn if terms have no IDs
        empty_terms <-
            data %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                n_ids = sum(!is.na(!!sym(id_col)) & !!sym(id_col) != ""),
                .groups = "drop"
            ) %>%
            dplyr::filter(n_ids == 0L)
        if (nrow(empty_terms) > 0L) {
            log_warn(
                "Terms with no IDs will be dropped: %s",
                paste(empty_terms[[term_col]], collapse = ", ")
            )
        }

        term_metabolites <-
            data %>%
            dplyr::filter(!is.na(!!sym(id_col)), !!sym(id_col) != "") %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                MetaboliteIDs = list(unique(!!sym(id_col))),
                .groups = "drop"
            )
    } else {
        # Enrichment input: split delimited metabolite IDs per term
        term_metabolites <-
            data %>%
            dplyr::mutate(
                MetaboliteIDs = strsplit(
                    as.character(.data[[id_col]]),
                    delimiter,
                    fixed = TRUE
                )
            ) %>%
            dplyr::mutate(
                MetaboliteIDs = lapply(
                    MetaboliteIDs,
                    function(x) {
                        x <- trimws(x)
                        x[x != ""]
                    }
                )
            ) %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                MetaboliteIDs = list(unique(unlist(MetaboliteIDs))),
                .groups = "drop"
            )

        empty_terms <- term_metabolites[lengths(term_metabolites$MetaboliteIDs) == 0L, , drop = FALSE]
        if (nrow(empty_terms) > 0L) {
            log_warn(
                "Terms with no IDs will be dropped: %s",
                paste(empty_terms[[term_col]], collapse = ", ")
            )
        }
        term_metabolites <- term_metabolites[lengths(term_metabolites$MetaboliteIDs) > 0L, , drop = FALSE]
    }

    terms <- term_metabolites[[term_col]]
    n_terms <- length(terms)
    if (n_terms == 0L) {
        stop("No terms with non-missing IDs after filtering.")
    }

    # ---- Similarity matrix ----------------------------------------------
    similarity_matrix <-
        matrix(
            1,
            nrow = n_terms,
            ncol = n_terms,
            dimnames = list(terms, terms)
        )

    if (similarity %in% c("jaccard", "overlap_coefficient")) {
        combs <- combn(seq_len(n_terms), 2L)
        for (j in seq_len(ncol(combs))) {
            i1 <- combs[1, j]
            i2 <- combs[2, j]
            set1 <- term_metabolites$MetaboliteIDs[[i1]]
            set2 <- term_metabolites$MetaboliteIDs[[i2]]
            inter <- length(intersect(set1, set2))
            if (similarity == "jaccard") {
                denom <- length(union(set1, set2))
            } else { # overlap coefficient
                denom <- min(length(set1), length(set2))
            }
            sim <- if (denom == 0) 0 else inter / denom
            t1 <- terms[i1]
            t2 <- terms[i2]
            similarity_matrix[t1, t2] <- sim
            similarity_matrix[t2, t1] <- sim
        }
    } else { # correlation
        metabolites <- unique(unlist(term_metabolites$MetaboliteIDs))
        binary_matrix <-
            matrix(
                0,
                nrow = n_terms,
                ncol = length(metabolites),
                dimnames = list(terms, metabolites)
            )
        for (i in seq_len(n_terms)) {
            ids <- term_metabolites$MetaboliteIDs[[i]]
            binary_matrix[i, colnames(binary_matrix) %in% ids] <- 1
        }
        corr <- cor(t(binary_matrix), method = correlation_method)
        corr[corr < 0] <- 0
        diag(corr) <- 1
        similarity_matrix <- corr
    }

    # Distance matrix
    distance_matrix <- 1 - similarity_matrix
    diag(distance_matrix) <- 0

    # Thresholded versions
    similarity_thr <- similarity_matrix
    similarity_thr[similarity_thr < threshold] <- 0
    diag(similarity_thr) <- 0

    distance_thr <- distance_matrix
    distance_thr[similarity_thr == 0] <- 1
    diag(distance_thr) <- 0

    # ---- Clustering ------------------------------------------------------
    clusters <- rep(NA_integer_, n_terms)
    names(clusters) <- terms

    if (clust == "components") {
        adj <- similarity_thr
        adj[adj > 0] <- 1
        g <- igraph::graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = NULL
        )
        mem <- igraph::components(g)$membership
        if (is.null(names(mem))) {
            names(mem) <- igraph::V(g)$name
        }
        clusters[names(mem)] <- mem
    } else if (clust == "community") {
        adj <- similarity_thr
        g <- igraph::graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = TRUE
        )
        mem <- igraph::cluster_louvain(g)$membership
        if (is.null(names(mem))) {
            names(mem) <- igraph::V(g)$name
        }
        clusters[names(mem)] <- mem
    } else if (clust == "hierarchical") {
        hc <- hclust(as.dist(distance_thr), method = hclust_method)
        cut_height <- 1 - threshold
        mem <- cutree(hc, h = cut_height)
        if (is.null(names(mem))) {
            names(mem) <- terms
        }
        clusters[names(mem)] <- mem
    }

    # ---- Apply minimum size filter --------------------------------------
    cluster_sizes <- table(clusters, useNA = "no")
    small <- names(cluster_sizes[cluster_sizes < min])
    clusters[as.character(clusters) %in% small] <- NA_integer_

    # Label clusters
    cluster_labels <- ifelse(
        is.na(clusters),
        "None",
        paste0("cluster", clusters)
    )
    names(cluster_labels) <- terms

    # Merge back to data
    term_metabolites$cluster <- cluster_labels[term_metabolites[[term_col]]]
    df <-
        data %>%
        dplyr::left_join(
            term_metabolites %>% dplyr::select(!!sym(term_col), cluster),
            by = term_col
        )

    # ---- Summary ---------------------------------------------------------
    cluster_summary <-
        term_metabolites %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(n_terms = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(pct_terms = 100 * n_terms / sum(n_terms))

    # ---- Node sizes ------------------------------------------------------
    node_sizes <- NULL
    if (!is.null(node_size_column)) {
        node_size_df <-
            data %>%
            dplyr::group_by(!!sym(term_col)) %>%
            dplyr::summarize(
                node_size = mean(.data[[node_size_column]], na.rm = TRUE),
                .groups = "drop"
            )
        node_size_df$node_size[is.nan(node_size_df$node_size)] <- NA_real_
        node_sizes <- node_size_df$node_size
        names(node_sizes) <- node_size_df[[term_col]]
        node_sizes <- node_sizes[terms]
        attr(node_sizes, "label") <- node_size_column
    }

    # ---- Graph plot ------------------------------------------------------
    graph_plot <- viz_graph(
        similarity_matrix = similarity_matrix,
        clusters = cluster_labels,
        threshold = threshold,
        plot_name = plot_name,
        max_nodes = max_nodes,
        min_degree = min_degree,
        node_sizes = node_sizes,
        show_density = show_density,
        seed = seed,
        save_plot = save_plot,
        print_plot = print_plot,
        path = path
    )

    out <- list(
        data = df,
        cluster_summary = cluster_summary,
        clusters = cluster_labels,
        similarity_matrix = similarity_matrix,
        distance_matrix = distance_matrix,
        graph_plot = graph_plot
    )

    log_trace(paste0(
        "cluster_pk completed with ",
        clust,
        " clustering; ",
        length(unique(cluster_labels)),
        " clusters (including None)."
    ))

    return(out)
}
