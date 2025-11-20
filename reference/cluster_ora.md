# Overrepresentation analysis by cluster

Uses enricher to run ORA on each of the metabolite cluster from any of
the MCA functions using a pathway list

## Usage

``` r
cluster_ora(
  data,
  metadata_info = c(ClusterColumn = "RG2_Significant", BackgroundColumn = "BG_method",
    PathwayTerm = "term", PathwayFeature = "Metabolite"),
  remove_background = TRUE,
  input_pathway,
  pathway_name = "",
  min_gssize = 10,
  max_gssize = 1000,
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  DF with metabolite names/metabolite IDs as row names. Metabolite
  names/IDs need to match the identifier type (e.g. HMDB IDs) in the
  input_pathway.

- metadata_info:

  *Optional:* Pass ColumnName of the column including the cluster names
  that ORA should be performed on (=ClusterColumn). BackgroundColumn
  passes the column name needed if remove_background=TRUE. Also pass
  ColumnName for input_pathway including term and feature names.
  (ClusterColumn= ColumnName data, BackgroundColumn = ColumnName data,
  PathwayTerm= ColumnName input_pathway, PathwayFeature= ColumnName
  input_pathway) **c(FeatureName="Metabolite",
  ClusterColumn="RG2_Significant", BackgroundColumn="BG_method",
  PathwayTerm= "term", PathwayFeature= "Metabolite")**

- remove_background:

  *Optional:* If TRUE, column BackgroundColumn name needs to be in
  metadata_info, which includes TRUE/FALSE for each metabolite to fall
  into background based on the chosen Background method for e.g.
  mca_2cond are removed from the universe. **default: TRUE**

- input_pathway:

  DF that must include column "term" with the pathway name, column
  "Feature" with the Metabolite name or ID and column "Description" with
  pathway description.

- pathway_name:

  *Optional:* Name of the pathway list used **default: ""**

- min_gssize:

  *Optional:* minimum group size in ORA **default: 10**

- max_gssize:

  *Optional:* maximum group size in ORA **default: 1000**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt" **default: "csv"**

- path:

  *Optional:* Path to the folder the results should be saved at.
  **default: NULL**

## Value

Saves results as individual .csv files.

## Examples

``` r
KEGG_Pathways <- metsigdb_kegg()
data(intracell_dma) # loads the object into your environment
DMAres <- intracell_dma %>%
    dplyr::filter(!is.na(KEGGCompound)) %>%
    tibble::column_to_rownames("KEGGCompound") %>%
    dplyr::select(-"Metabolite")
RES <- cluster_ora(
    data = DMAres,
    metadata_info = c(
        ClusterColumn = "Pathway",
        PathwayTerm = "term",
        PathwayFeature = "Metabolite"
    ),
    input_pathway = KEGG_Pathways,
    remove_background = FALSE
)
```
