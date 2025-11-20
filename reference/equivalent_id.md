# Find additional potential IDs for "kegg", "pubchem", "chebi", "hmdb"

Find additional potential IDs for "kegg", "pubchem", "chebi", "hmdb"

## Usage

``` r
equivalent_id(
  data,
  metadata_info = c(InputID = "MetaboliteID"),
  from = "hmdb",
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  dataframe with at least one column with the detected metabolite IDs
  (one ID per row).

- metadata_info:

  *Optional:* Column name of metabolite IDs. **Default =
  list(InputID="MetaboliteID")**

- from:

  ID type that is present in your data. Choose between "kegg",
  "pubchem", "chebi", "hmdb". **Default = "hmdb"**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  Optional: Path to the folder the results should be saved at. **Default
  = NULL**

## Value

Input DF with additional column including potential additional IDs.

## Examples

``` r
data(cellular_meta)
DetectedIDs <- cellular_meta %>% tidyr::drop_na()
Res <- equivalent_id(
    data = DetectedIDs,
    metadata_info = c(InputID = "HMDB"),
    from = "hmdb"
)
#> Warning: The following IDs are duplicated and removed: HMDB0000725, HMDB0000267, HMDB0000755
#> chebi is used to find additional potential IDs for hmdb.
```
