# Basic cluster_pk checks across clustering modes

library(MetaProViz)
library(dplyr)
library(logger)

source("R/VizGraph.R")
source("cluster_pk/cluster_pk.R")

d <- metsigdb_kegg()

# Components
r1 <- cluster_pk(d, threshold = 0.2, clust = "components", min = 1, save_plot = NULL)
table(r1$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r1$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# Community
r2 <- cluster_pk(d, threshold = 0.2, clust = "community", min = 1, save_plot = NULL)
table(r2$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r2$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)

# Hierarchical
r3 <- cluster_pk(
    d,
    threshold = 0.2,
    clust = "hierarchical",
    hclust_method = "average",
    min = 1,
    save_plot = NULL
)
table(r3$data$cluster, useNA = "ifany") %>% sort(decreasing = TRUE) %>% head(10)
r3$cluster_summary %>% arrange(desc(n_terms)) %>% head(10)
