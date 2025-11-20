# Create metabolite sets from existing genesets

Gene to metabolite translation is based on mappings in Recon-3D
(cosmosR).

## Usage

``` r
make_gene_metab_set(
  input_pk,
  metadata_info = c(Target = "gene"),
  pk_name = NULL,
  save_table = "csv",
  path = NULL
)
```

## Arguments

- input_pk:

  dataframe with two columns for source (=term) and Target (=gene), e.g.
  Hallmarks.

- metadata_info:

  *Optional:* Column name of Target in input_pk. **Default =
  c(Target="gene")**

- pk_name:

  *Optional:* Name of the prior knowledge resource. **default: NULL**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  Optional: String which is added to the resulting folder name
  **default: NULL**

## Value

List of two data frames: "GeneMetabSet" and "MetabSet".

## Examples

``` r
data(hallmarks)
make_gene_metab_set(hallmarks)
```
