# Cluster terms in prior knowledge by set overlap

Cluster terms in prior knowledge by set overlap

## Usage

``` r
cluster_pk(
  data,
  metadata_info = c(metabolite_column = "MetaboliteID", pathway_column = "term"),
  similarity = c("jaccard", "overlap_coefficient", "correlation"),
  correlation_method = "pearson",
  input_format = c("long", "enrichment"),
  delimiter = "/",
  threshold = 0.5,
  plot_threshold = 0,
  clust = c("components", "community", "hierarchical"),
  hclust_method = "average",
  min = 2,
  plot_name = "ClusterGraph",
  max_nodes = 10000,
  min_degree = 1,
  node_size_column = NULL,
  show_density = FALSE,
  save_plot = "png",
  print_plot = FALSE,
  path = NULL
)
```

## Arguments

- data:

  Long data frame with one ID per row, or enrichment-style table with a
  delimited metabolite list per term (see input_format).

- metadata_info:

  List with entries `metabolite_column` (metabolite ID column or
  delimited list column) and `pathway_column` (term column). Defaults to
  `c(metabolite_column = "MetaboliteID", pathway_column = "term")`.

- similarity:

  Similarity measure between term ID sets. Options: "jaccard" (default),
  "overlap_coefficient", or "correlation". Jaccard similarity is \|A
  intersect B\| / \|A union B\|. Overlap coefficient is \|A intersect
  B\| / min(\|A\|, \|B\|). Jaccard is stricter for large sets, while
  overlap_coefficient is more permissive for nested sets.

- correlation_method:

  Correlation method when `similarity = "correlation"`. One of
  "pearson", "spearman", "kendall". Ignored otherwise.

- input_format:

  Input format of `data`. Use "long" for one ID per row (default) or
  "enrichment" for one term per row with a delimited ID list. The
  `metabolite_column` entry in metadata_info is interpreted accordingly.

- delimiter:

  Delimiter for metabolite ID lists when input_format = "enrichment".
  Ignored for input_format = "long". Default = "/".

- threshold:

  Similarity cutoff for keeping edges (applies to all clustering modes).
  Default = 0.5.

- plot_threshold:

  Similarity cutoff for plotting edges in viz_graph. Default = 0 (plot
  all edges with similarity \> 0).

- clust:

  Clustering strategy: "components" (connected components on thresholded
  unweighted graph), "community" (Louvain on thresholded weighted
  graph), or "hierarchical" (hclust on distance = 1 - similarity).

- hclust_method:

  Linkage method for hierarchical clustering. One of "average"
  (default), "single", "complete", "ward.D", "ward.D2", "mcquitty",
  "median", "centroid". Used only when clust = "hierarchical".

- min:

  Minimum cluster size; smaller clusters are relabeled to "None".
  Default = 2.

- plot_name:

  *Optional:* String added to output files of the plot. Default =
  "ClusterGraph".

- max_nodes:

  *Optional:* Maximum nodes for plotting. If set, keeps nodes from the
  largest component up to this limit (by degree). Used only for the
  graph plot. Default = 10000.

- min_degree:

  *Optional:* Minimum degree filter for graph plotting. Used only for
  the graph plot. Default = 1.

- node_size_column:

  *Optional:* Numeric column name from `data` used to scale node sizes
  in the graph. Aggregated per term (mean) when multiple rows map to the
  same term. Default = NULL.

- show_density:

  *Optional:* If TRUE, add a hull background per cluster to the graph.
  Default = FALSE.

- save_plot:

  *Optional:* Select the file type of output plots. Options are svg,
  pdf, png or NULL. **Default = "svg"**

- print_plot:

  *Optional:* If TRUE prints an overview of resulting plots. **Default =
  FALSE**

- path:

  Optional: String which is added to the resulting folder name.
  **default: NULL**

## Value

A list with:

- data:

  Input data with a `cluster` column added.

- cluster_summary:

  Summary of cluster sizes and percentages.

- clusters:

  Named vector of term -\> cluster assignment.

- similarity_matrix:

  Term-by-term similarity matrix.

- distance_matrix:

  Term-by-term distance matrix (1 - similarity).

- node_sizes:

  Named numeric vector of node sizes used in plotting (or NULL).

- graph_plot:

  Graph plot returned by viz_graph.

## Examples

``` r
# Create toy pathway data in long format
toy_pw <- data.frame(
    MetaboliteID = c("C1", "C2", "C3", "C1", "C2", "C4", "C3", "C4", "C5"),
    term = c("pA", "pA", "pA", "pB", "pB", "pB", "pC", "pC", "pC")
)

r <- cluster_pk(
    toy_pw,
    metadata_info = c(
        metabolite_column = "MetaboliteID",
        pathway_column = "term"
    ),
    input_format = "long",
    similarity = "jaccard",
    threshold = 0.1,
    clust = "community",
    min = 1,
    save_plot = NULL,
    print_plot = FALSE
)
```
