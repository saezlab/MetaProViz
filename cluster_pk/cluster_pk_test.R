
library(MetaProViz)
library(dplyr)
library(igraph)
library(logger)
library(rlang)
library(stats)

source("cluster_pk/cluster_pk.R")
source("cluster_pk/helper_plots_cluster_pk.R")


d <- metsigdb_kegg()
r <- cluster_pk(d)  # defaults; adjust threshold as needed


# 1) More permissive components: lower threshold, keep singletons
r1 <- cluster_pk(d, threshold = 0.2, clust = "components", min = 1)
table(r1$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r1$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

data_r1 <- r1$data
cluster_summary_r1 <- r1$cluster_summary
r1_sim_mat <- r1$similarity_matrix
r1_dist_mat <- r1$distance_matrix
r1_cluster_summary <- r1$cluster_summary
r1_clusters <- r1$clusters

## produce error due to unmatched hclust_method argument
# cluster_pk(d, threshold = 0.2, clust = "hierarchical", hclust_method = "gi")

# 2) Community detection on weighted graph, mid threshold
r2 <- cluster_pk(d, threshold = 0.3, clust = "community", min = 1)
table(r2$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r2$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 3) Hierarchical cut by height derived from threshold
r3 <- cluster_pk(d, threshold = 0.3, clust = "hierarchical", hclust_method = "average", k = NULL, min = 1)
table(r3$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r3$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 4) Hierarchical with fixed k (ignore threshold for cutting)
r4 <- cluster_pk(d, threshold = 0.3, clust = "hierarchical", hclust_method = "average", k = 10, min = 1)
table(r4$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r4$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 5) Correlation-based similarity (Spearman), community detection, drop negatives
r5 <- cluster_pk(d, similarity = "correlation", matrix = "spearman",
                 drop_negative = TRUE, threshold = 0.2, clust = "community", min = 1)
table(r5$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r5$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)



# 11) More permissive components: lower threshold, no singletons
r11 <- cluster_pk(d, threshold = 0.2, clust = "components", min = 2)
table(r11$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r11$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 12) Community detection on weighted graph, mid threshold
r12 <- cluster_pk(d, threshold = 0.3, clust = "community", min = 2)
table(r12$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r12$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 13) Hierarchical cut by height derived from threshold
r13 <- cluster_pk(d, threshold = 0.3, clust = "hierarchical", hclust_method = "average", k = NULL, min = 2)
table(r13$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r13$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 14) Hierarchical with fixed k (ignore threshold for cutting)
r14 <- cluster_pk(d, threshold = 0.3, clust = "hierarchical", hclust_method = "average", k = 10, min = 2)
table(r14$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r14$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# 15) Correlation-based similarity (Spearman), community detection, drop negatives
r15 <- cluster_pk(d, similarity = "correlation", matrix = "spearman",
                 drop_negative = TRUE, threshold = 0.2, clust = "community", min = 2)
table(r15$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r15$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)



# 
# Rebuilt cluster_pk to fix correctness, clarify behavior, and add options.
#     
#     What was wrong before
# 
# No input validation for required columns or parameters.
# Overlap vs distance naming was confusing; the “overlap_matrix” had 1 on the diagonal and was actually a distance.
# correlation_matrix was only defined in one branch, yet used elsewhere → errors.
# threshold was hard-coded (0.7), and term_clusters was unused.
# Graph clustering used a Gaussian kernel, making a fully connected graph → everything in one cluster.
# min was unused; no output return; no summaries.
#
#
#     What’s fixed and new
# 
# Added input checks for data frame, required columns (InputID, grouping_variable), parameter ranges, and duplication handling.
# Clear similarity/distance construction:
#     similarity options: "jaccard" (default), "overlap" (overlap coefficient), or "correlation" (pearson/spearman/kendall on binary term×ID matrix, with optional negative dropping).
# Distance = 1 − similarity; diagonals set correctly (sim=1, dist=0).
# Thresholding: edges with sim < threshold are zeroed; corresponding distances set to max (1).
# Clustering options (clust):
#     "components": connected components on thresholded unweighted graph.
# "community": Louvain on thresholded weighted graph.
# "hierarchical": hclust on distance; cut by height = 1 − threshold (default) or by k.
# min now applied: clusters smaller than min become “None”.
# Outputs: returns a list with clustered data (cluster column), similarity and distance matrices, cluster summary (counts and %), and term→cluster assignments; logging via log_trace/log_warn.
# Current parameters and behavior
# 
# similarity: choose Jaccard, overlap coefficient, or correlation (with matrix = pearson/spearman/kendall; drop_negative to zero negatives).
# threshold: similarity cutoff; controls sparsity and clustering across all methods.
# clust: components, community, or hierarchical (with hclust_method and optional k).
# min: minimum cluster size; smaller clusters labeled “None.”
# Defaults: similarity=jaccard, threshold=0.7, clust=components, min=2, drop_negative=TRUE.
#
#
#     What the function does now (flow)
# 
# Validate inputs; drop duplicate term–ID rows; warn on empty terms.
# Build term→ID lists.
# Compute term×term similarity (per chosen metric); derive distance.
# Threshold similarity; adjust distance accordingly.
# Cluster via selected method; cut dendrogram by height or k; relabel small clusters to “None.”
# Merge cluster labels back into the original data; produce cluster summary; return all objects.
# Empirical notes from metsigdb_kegg() tests
# 
# Jaccard overlaps are low; high thresholds (≥0.3) with community detection yield all “None.” Components/hierarchical with lower thresholds (≈0.2) and min=1 produce multiple clusters; min strongly affects “None.”
# Next ideas (if needed)
# 
# Add similarity distribution plots, heatmaps, dendrograms, and cluster size plots.
# Consider an “auto” threshold (e.g., percentile-based) or a helper to summarize similarity stats.
