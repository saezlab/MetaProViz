# Meta prior-knowledge

Meta prior-knowledge

## Usage

``` r
meta_pk(
  data,
  metadata_sample,
  metadata_info = NULL,
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  SummarizedExperiment (se) file including assay and colData. If se file
  is provided, metadata_sample is extracted from the colData of the se
  object. Alternatively provide a DF with unique sample identifiers as
  row names and metabolite numerical values in columns with metabolite
  identifiers as column names. Use NA for metabolites that were not
  detected.

- metadata_sample:

  *Optional:* Only required if you did not provide se file in parameter
  data. Provide DF which contains metadata information about the
  samples, which will be combined with your input data based on the
  unique sample identifiers used as rownames. **Default = NULL**

- metadata_info:

  *Optional:* NULL or vector with column names that should be used, i.e.
  c("Age", "gender", "Tumour-stage"). **default: NULL**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  *Optional:* Path to the folder the results should be saved at.
  **default: NULL**

## Value

DF with prior knowledge based on patient metadata

## Examples

``` r
data(tissue_norm_se)
Res <- meta_pk(tissue_norm_se)

data(tissue_norm)
Tissue_Norm <- tissue_norm %>% tibble::column_to_rownames("Code")
Res <- meta_pk(
    data = Tissue_Norm[, -c(1:13)],
    metadata_sample = Tissue_Norm[, c(2, 4:5, 12:13)]
)
```
