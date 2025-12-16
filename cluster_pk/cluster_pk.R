##
## Cluster prior knowledge by pathway overlap
##

#' Cluster terms in prior knowledge by set overlap
#'
#' @param data Long data frame with one ID per row; must contain columns for the
#'     target IDs and the grouping term.
#' @param metadata_info List with entries `InputID` (ID column) and
#'     `grouping_variable` (term column). Defaults to
#'     `c(InputID = "MetaboliteID", grouping_variable = "term")`.
#' @param similarity Similarity measure between term ID sets. Options:
#'     "jaccard" (default), "overlap" (overlap coefficient), or "correlation"
#'     (applied to the binary term-by-ID matrix using `matrix`).
#' @param matrix Correlation method when `similarity = "correlation"`.
#'     One of "pearson", "spearman", "kendall". Ignored otherwise.
#' @param threshold Similarity cutoff for keeping edges (applies to all
#'     clustering modes). Default = 0.7.
#' @param clust Clustering strategy: "components" (connected components on
#'     thresholded unweighted graph), "community" (Louvain on thresholded
#'     weighted graph), or "hierarchical" (hclust on distance = 1 - similarity).
#' @param hclust_method Linkage method for hierarchical c lustering. One of
#'     "average" (default), "single", "complete", "ward.D", "ward.D2",
#'     "mcquitty", "median", "centroid".
#' @param k Optional number of clusters for hierarchical clustering. If NULL,
#'     dendrogram is cut at height = 1 - threshold.
#' @param min Minimum cluster size; smaller clusters are relabeled to "None".
#'     Default = 2.
#' @param drop_negative If TRUE, negative correlations are set to zero before
#'     thresholding. Default = TRUE.
#' @return A list with:
#'     \item{data}{Input data with a `cluster` column added.}
#'     \item{similarity_matrix}{Term-by-term similarity matrix.}
#'     \item{distance_matrix}{Term-by-term distance matrix (1 - similarity).}
#'     \item{cluster_summary}{Summary of cluster sizes and percentages.}
#'     \item{clusters}{Named vector of term -> cluster assignment.}
#'
#' @importFrom dplyr group_by summarize ungroup mutate select left_join
#' @importFrom dplyr across n distinct filter tibble arrange
#' @importFrom igraph graph_from_adjacency_matrix components cluster_louvain
#' @importFrom stats cor as.dist hclust cutree
#' @importFrom rlang sym !!
#' @importFrom logger log_trace log_warn
#' 
#' @examples 
#' 
#' data <- metsigdb_kegg()
#' 
#' res <- cluster_pk(data, threshold = 0.2, clust = "components", min = 1)
#' table(res$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
#' res$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)
#' ## add plotting into example
#' 
#' @noRd
cluster_pk <- function(
    data,
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    similarity = c("jaccard", "overlap", "correlation"),
    matrix = "pearson",
    threshold = 0.7,
    clust = c("components", "community", "hierarchical"),
    hclust_method = "average",
    k = NULL,
    min = 2,
    drop_negative = TRUE
) {
    
    # NSE vs. R CMD check workaround
    cluster <- NULL
    
    # ---- Input checks ----------------------------------------------------
    ## move these to Helper functions probably
    similarity <- match.arg(similarity)
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
    
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.")
    }
    
    if (!all(c("InputID", "grouping_variable") %in% names(metadata_info))) {
        stop("metadata_info must contain InputID and grouping_variable.")
    }
    
    id_col <- metadata_info[["InputID"]]
    term_col <- metadata_info[["grouping_variable"]]
    
    if (!id_col %in% colnames(data)) {
        stop(sprintf("Column %s (InputID) not found in data.", id_col))
    }
    if (!term_col %in% colnames(data)) {
        stop(sprintf("Column %s (grouping_variable) not found in data.", term_col))
    }
    
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
        stop("`threshold` must be a single numeric value.")
    }
    if (threshold < 0 || threshold > 1) {
        stop("`threshold` must be between 0 and 1 (similarity scale).")
    }
    
    if (!is.null(k)) {
        if (!is.numeric(k) || length(k) != 1L || k < 1) {
            stop("`k` must be NULL or a single positive integer.")
        }
        k <- as.integer(k)
    }
    
    if (!is.numeric(min) || length(min) != 1L || min < 1) {
        stop("`min` must be a single integer >= 1.")
    }
    min <- as.integer(min)
    
    if (similarity == "correlation") {
        matrix <- match.arg(matrix, c("pearson", "spearman", "kendall"))
    }
    
    # Drop duplicated term-ID rows to avoid inflating overlaps
    data <- distinct(data, !!sym(term_col), !!sym(id_col), .keep_all = TRUE)
    
    # Warn if terms have no IDs
    empty_terms <-
        data %>%
        group_by(!!sym(term_col)) %>%
        summarize(
            n_ids = sum(!is.na(!!sym(id_col)) & !!sym(id_col) != ""),
            .groups = "drop"
        ) %>%
        filter(n_ids == 0L)
    if (nrow(empty_terms) > 0L) {
        log_warn(
            "Terms with no IDs will be dropped: %s",
            paste(empty_terms[[term_col]], collapse = ", ")
        )
    }
    
    # ---- Build term -> ID list ------------------------------------------
    term_metabolites <-
        data %>%
        filter(!is.na(!!sym(id_col)), !!sym(id_col) != "") %>%
        group_by(!!sym(term_col)) %>%
        summarize(
            MetaboliteIDs = list(unique(!!sym(id_col))),
            .groups = "drop"
        )
    
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
    
    if (similarity %in% c("jaccard", "overlap")) {
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
        corr <- cor(t(binary_matrix), method = matrix)
        if (drop_negative) {
            corr[corr < 0] <- 0
        }
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
        g <- graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = NULL
        )
        mem <- components(g)$membership
        clusters[names(mem)] <- mem
        hclust_obj <- NULL
    } else if (clust == "community") {
        adj <- similarity_thr
        g <- graph_from_adjacency_matrix(
            adj,
            mode = "undirected",
            weighted = TRUE
        )
        mem <- cluster_louvain(g)$membership
        clusters[names(mem)] <- mem
        hclust_obj <- NULL
    } else if (clust == "hierarchical") {
        hc <- hclust(as.dist(distance_thr), method = hclust_method)
        if (!is.null(k)) {
            mem <- cutree(hc, k = k)
        } else {
            cut_height <- 1 - threshold
            mem <- cutree(hc, h = cut_height)
        }
        clusters[names(mem)] <- mem
        hclust_obj <- hc
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
        left_join(
            term_metabolites %>% select(!!sym(term_col), cluster),
            by = term_col
        )
    
    # ---- Summary ---------------------------------------------------------
    cluster_summary <-
        term_metabolites %>%
        group_by(cluster) %>%
        summarize(n_terms = dplyr::n(), .groups = "drop") %>%
        mutate(pct_terms = 100 * n_terms / sum(n_terms))
    
    out <- list(
        data = df,
        similarity_matrix = similarity_matrix,
        distance_matrix = distance_matrix,
        cluster_summary = cluster_summary,
        clusters = cluster_labels,
        hclust = hclust_obj,
        hclust_method = hclust_method,
        threshold = threshold
    )
    
    log_trace(
        "cluster_pk completed with %s clustering; %i clusters (including None).",
        clust,
        length(unique(cluster_labels))
    )
    
    return(out)
}
