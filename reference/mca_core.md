# Metabolite clustering analysis for core experiments

Performs metabolite clustering analysis and computes clusters based on
regulatory rules between intracellular and culture media metabolomics in
core experiments.

## Usage

``` r
mca_core(
  data_intra,
  data_core,
  metadata_info_intra = c(ValueCol = "Log2FC", StatCol = "p.adj", cutoff_stat = 0.05,
    ValueCutoff = 1),
  metadata_info_core = c(DirectionCol = "core", ValueCol = "Log2(Distance)", StatCol =
    "p.adj", cutoff_stat = 0.05, ValueCutoff = 1),
  feature = "Metabolite",
  save_table = "csv",
  method_background = "Intra&core",
  path = NULL
)
```

## Arguments

- data_intra:

  DF for your data (results from e.g. dma) containing metabolites in
  rows with corresponding Log2FC and stat (p-value, p.adjusted) value
  columns.

- data_core:

  DF for your data (results from e.g. dma) containing metabolites in
  rows with corresponding Log2FC and stat (p-value, p.adjusted) value
  columns. Here we additionally require

- metadata_info_intra:

  *Optional:* Pass ColumnNames and Cutoffs for the intracellular
  metabolomics including the value column (e.g. Log2FC, Log2Diff, t.val,
  etc) and the stats column (e.g. p.adj, p.val). This must include:
  c(ValueCol=ColumnName_data_intra,StatCol=ColumnName_data_intra,
  cutoff_stat= NumericValue, ValueCutoff=NumericValue)
  **Default=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05,
  ValueCutoff=1)**

- metadata_info_core:

  *Optional:* Pass ColumnNames and Cutoffs for the consumption-release
  metabolomics including the direction column, the value column (e.g.
  Log2Diff, t.val, etc) and the stats column (e.g. p.adj, p.val). This
  must include: c(DirectionCol=
  ColumnName_data_core,ValueCol=ColumnName_data_core,StatCol=ColumnName_data_core,
  cutoff_stat= NumericValue,
  ValueCutoff=NumericValue)**Default=c(DirectionCol="core",
  ValueCol="Log2(Distance)",StatCol="p.adj", cutoff_stat= 0.05,
  ValueCutoff=1)**

- feature:

  *Optional:* Column name of Column including the Metabolite
  identifiers. This MUST BE THE SAME in each of your Input files.
  **Default="Metabolite"**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt" **default: "csv"**

- method_background:

  *Optional:* Background method \`Intra\|core, Intra&core, core, Intra
  or \* **Default="Intra&core"**

- path:

  *Optional:* Path to the folder the results should be saved at.
  **default: NULL**

## Value

List of two DFs: 1. summary of the cluster count and 2. the detailed
information of each metabolites in the clusters.

## Examples

``` r
data(intracell_dma)

# Create mock CoRe DMA results with required columns
core_dma <- data.frame(
    Metabolite = intracell_dma$Metabolite[1:50],
    `Log2(Distance)` = runif(50, -2, 2),
    p.adj = runif(50, 0, 0.1),
    core = sample(c("Consumption", "Release"), 50, replace = TRUE),
    check.names = FALSE
)

Res <- mca_core(
    data_intra = as.data.frame(intracell_dma),
    data_core = core_dma,
    save_table = NULL
)
```
