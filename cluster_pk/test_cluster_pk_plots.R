# Simple test script for cluster_pk with graph plotting

library(MetaProViz)
library(dplyr)
library(ggplot2)
library(logger)
library(igraph)
library(ggraph)
devtools::load_all()


source("R/VizGraph.R")
source("cluster_pk/cluster_pk.R")

# Load example data
d <- metsigdb_kegg()

# Run clustering with graph plotting
r <- cluster_pk(
    d,
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    similarity = "jaccard",
    threshold = 0.5,
    clust = "hierarchical",
    min = 2,
    debug = TRUE,
    save_plot = "png",
    min_degree = 1,
    max_nodes = 1000
) 

print(head(r$cluster_summary))
print(r$graph_plot)










