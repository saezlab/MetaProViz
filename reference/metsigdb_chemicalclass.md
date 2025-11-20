# Metabolite chemical classes from RaMP DB

Metabolite chemical classes from RaMP DB

## Usage

``` r
metsigdb_chemicalclass(version = "2.5.4", save_table = "csv", path = NULL)
```

## Arguments

- version:

  *Optional:* Version of the RaMP database loaded from OmniPathR.
  **default: "2.5.4"**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  Optional: String which is added to the resulting folder name
  **default: NULL**

## Value

A data frame containing the Prior Knowledge.

## Examples

``` r
ChemicalClass <- metsigdb_chemicalclass()
```
