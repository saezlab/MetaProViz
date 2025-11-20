# Translate IDs to/from KEGG, PubChem, Chebi, HMDB

Translate IDs to/from KEGG, PubChem, Chebi, HMDB

## Usage

``` r
translate_id(
  data,
  metadata_info = c(InputID = "MetaboliteID", grouping_variable = "term"),
  from = "kegg",
  to = c("pubchem", "chebi", "hmdb"),
  summary = FALSE,
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  dataframe with at least one column with the target (e.g. metabolite),
  you can add other columns such as source (e.g. term). Must be "long"
  DF, meaning one ID per row.

- metadata_info:

  *Optional:* Column name of Target in input_pk. **Default =
  list(InputID="MetaboliteID" , grouping_variable="term")**

- from:

  ID type that is present in your data. Choose between "kegg",
  "pubchem", "chebi", "hmdb". **Default = "kegg"**

- to:

  One or multiple ID types to which you want to translate your data.
  Choose between "kegg", "pubchem", "chebi", "hmdb". **Default =
  c("pubchem","chebi","hmdb")**

- summary:

  *Optional:* If TRUE a long summary tables are created. **Default =
  FALSE**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  Optional: Path to the folder the results should be saved at. **Default
  = NULL**

## Value

List with at least three DFs: 1) Original data and the new column of
translated ids spearated by comma. 2) Mapping information between
Original ID to Translated ID. 3) Mapping summary between Original ID to
Translated ID.

## Examples

``` r
if (FALSE) { # \dontrun{
KEGG_Pathways <- metsigdb_kegg()
Res <- translate_id(
    data = KEGG_Pathways,
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    from = c("kegg"),
    to = c("pubchem", "hmdb")
)
} # }
```
