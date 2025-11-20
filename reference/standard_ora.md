# Overrepresentation analysis of metabolite sets in pathways

Can be applied on the result of differential metabolite analysis (DMA),
requires a pathway list (from databases).

## Usage

``` r
standard_ora(
  data,
  metadata_info = c(pvalColumn = "p.adj", percentageColumn = "t.val", PathwayTerm =
    "term", PathwayFeature = "Metabolite"),
  cutoff_stat = 0.05,
  cutoff_percentage = 10,
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

  *Optional:* Pass ColumnName of the column including parameters to use
  for cutoff_stat and cutoff_percentage. Also pass ColumnName for
  input_pathway including term and feature names. (pvalColumn =
  ColumnName data, percentageColumn= ColumnName data, PathwayTerm=
  ColumnName input_pathway, PathwayFeature= ColumnName input_pathway)
  **c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term",
  PathwayFeature= "Metabolite")**

- cutoff_stat:

  *Optional:* p-adjusted value cutoff from ORA results. Must be a
  numeric value. **default: 0.05**

- cutoff_percentage:

  *Optional:* percentage cutoff of metabolites that should be considered
  for ORA. Selects top/Bottom % of selected percentage Column, usually
  t.val or Log2FC **default: 10**

- input_pathway:

  DF that must include column "term" with the pathway name, column
  "Metabolite" with the Metabolite name or ID and column "Description"
  with pathway description that will be depicted on the plots.

- pathway_name:

  *Optional:* Name of the input_pathway used **default: ""**

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
RES <- standard_ora(
    data = DMAres,
    input_pathway = KEGG_Pathways
)
```
