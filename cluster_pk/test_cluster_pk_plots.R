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
        metabolite_column = "MetaboliteID",
        pathway_column = "term"
    ),
    input_format = "long",
    similarity = "jaccard",
    threshold = 0.2,
    clust = "community",
    min = 2,
    plot_name = "Long",
    save_plot = "png",
    min_degree = 1,
    seed = 123,
    show_density = TRUE,
    max_nodes = 1000
) 

print(head(r$cluster_summary))
print(r$graph_plot)



source("R/VizGraph.R")
source("cluster_pk/cluster_pk.R")


## test for enrichment results

d_enrichment <- read.csv("D:/Masterarbeit_LOCAL/enrichment_results_test.csv", sep=";")
names(d_enrichment)[1] <- "ID"

# d_enrichment[1, "percentage_of_Pathway_detected"] <- 100

r_enrichment <- cluster_pk(
    d_enrichment,
    metadata_info = c(
        metabolite_column = "Metabolites_in_pathway",
        pathway_column = "ID"
    ),
    input_format = "enrichment",
    similarity = "jaccard",
    threshold = 0.5,
    clust = "community",
    min = 2,
    node_size_column = "percentage_of_Pathway_detected",
    save_plot = "png",
    plot_name = "Enrichment",
    print_plot = FALSE,
    min_degree = 1,
    seed = 42,
    show_density = TRUE,
    max_nodes = 1000
)


viz_graph(
    as.numeric(r_enrichment$similarity_matrix),
    r_enrichment$clusters,
    threshold = 0.5,
    plot_name = "ClusterGraph",
    max_nodes = NULL,
    min_degree = NULL,
    node_sizes = NULL,
    show_density = FALSE,
    seed = NULL,
    save_plot = "svg",
    print_plot = TRUE,
    path = NULL
) 

