# Simple test script for cluster_pk and helper_plots_cluster_pk

library(MetaProViz)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(treemapify)
library(ggraph)
library(igraph)
library(logger)

source("cluster_pk/cluster_pk.R")
source("cluster_pk/helper_plots_cluster_pk.R")

# Load example data
d <- metsigdb_kegg()

# Run clustering
r <- cluster_pk(d, threshold = 0.2, clust = "components", min = 2)

# Generate plots
plots <- helper_plots_cluster_pk(r)
plots2 <- helper_plots_cluster_pk(r, min_degree = 2)

# Print summaries
print(head(r$cluster_summary))

## create new dir in MetaProViz/cluster_pk/ to save plots
## if dir does not exist yet
output_dir <- "cluster_pk/plots_output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

png(filename = file.path(output_dir, "similarity_heatmap.png"), width = 10000, height = 4000, res = 300)
plots$heatmap
dev.off()

png(filename = file.path(output_dir, "cluster_bar.png"), width = 10000, height = 4000, res = 300)
plots$cluster_bar
dev.off()

png(filename = file.path(output_dir, "cluster_treemap.png"), width = 10000, height = 4000, res = 300)
plots$cluster_treemap
dev.off()

png(filename = file.path(output_dir, "similarity_distribution.png"), width = 10000, height = 4000, res = 300)
plots$sim_distribution
dev.off()

png(filename = file.path(output_dir, "graph.png"), width = 10000, height = 4000, res = 300)
plots$graph
dev.off()

png(filename = file.path(output_dir, "graph_min_degree_2.png"), width = 10000, height = 4000, res = 300)
plots2$graph
dev.off()


png(filename = file.path(output_dir, "dendrogram.png"), width = 10000, height = 4000, res = 300)
plots$dendrogram
dev.off()





