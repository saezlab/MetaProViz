# Merges the analytical replicates of an experiment

Merges the analytical replicates of an experiment

## Usage

``` r
replicate_sum(
  data,
  metadata_sample = NULL,
  metadata_info = c(Conditions = "Conditions", Biological_Replicates =
    "Biological_Replicates", Analytical_Replicates = "Analytical_Replicates"),
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  SummarizedExperiment (se) file including assay and colData. If se file
  is provided, metadata_sample is extracted from the colData of the se
  object. Alternatively provide a DF which contains unique sample
  identifiers as row names and metabolite numerical values in columns
  with metabolite identifiers as column names. Use NA for metabolites
  that were not detected.

- metadata_sample:

  *Optional:* Only required if you did not provide se file in parameter
  data. Provide DF which contains information about the samples, which
  will be combined with the input data based on the unique sample
  identifiers used as rownames. Must contain column with Conditions. If
  you do not have multiple conditions in your experiment assign all
  samples into the same condition. **Default = NULL**

- metadata_info:

  *Optional:* Named vector including the Conditions and Replicates
  information: c(Conditions="ColumnNameConditions",
  Biological_Replicates="ColumnName_metadata_sample",
  Analytical_Replicates="ColumnName_metadata_sample").**Default = NULL**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt", ot NULL **default: "csv"**

- path:

  *Optional:* Path to the folder the results should be saved at.
  **default: NULL**

## Value

DF with the merged analytical replicates

## Examples

``` r
data(intracell_raw_se)
Res <- replicate_sum(data = intracell_raw_se)

data(intracell_raw)
Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
Res <- replicate_sum(
    data = Intra[-c(49:58), -c(1:3)],
    metadata_sample = Intra[-c(49:58), c(1:3)]
)
```
