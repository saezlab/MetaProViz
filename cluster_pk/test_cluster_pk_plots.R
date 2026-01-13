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
    threshold = 0.5,
    clust = "components",
    min = 2,
    save_plot = "png",
    min_degree = 1,
    max_nodes = 1000
) 

print(head(r$cluster_summary))
print(r$graph_plot)






## test for enrichment results

d_enrichment <- read.csv("D:/Masterarbeit_LOCAL/enrichment_results_test.csv", sep=";")

names(d_enrichment)[1] <- "ID"

r_enrichment <- cluster_pk(
    d_enrichment,
    metadata_info = c(
        metabolite_column = "Metabolites_in_pathway",
        pathway_column = "ID"
    ),
    input_format = "enrichment",
    similarity = "jaccard",
    threshold = 0.5,
    clust = "components",
    min = 2,
    save_plot = NULL,
    print_plot = TRUE,
    min_degree = 1,
    max_nodes = 1000
) 













# 
# 
# data(intracell_raw)
# Intra <- intracell_raw[-c(49:58), ] %>% tibble::column_to_rownames("Code")
# ResI <- dma(
#     data = Intra[, -c(1:3)],
#     metadata_sample = Intra[, c(1:3)],
#     metadata_info = c(
#         Conditions = "Conditions", Numerator = NULL, Denominator = "HK2"
#     )
# )
