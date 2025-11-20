# Create Mapping Ambiguities between two ID types

Create Mapping Ambiguities between two ID types

## Usage

``` r
mapping_ambiguity(
  data,
  from,
  to,
  grouping_variable = NULL,
  summary = FALSE,
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  Translated DF from translate_id reults or dataframe with at least one
  column with the target metabolite ID and another MetaboliteID type.
  One of the IDs can only have one ID per row, the other ID can be
  either

  by comma or a list. Optional: add other columns such as source (e.g.
  term).

- from:

  Column name of the secondary or translated metabolite identifier in
  data. Here can be multiple IDs per row either separated by comma " ,"
  or a list of IDs.

- to:

  Column name of original metabolite identifier in data. Here should
  only have one ID per row.

- grouping_variable:

  *Optional:* If NULL no groups are used. If TRUE provide column name in
  data containing the grouping_variable and features are grouped.
  **Default = NULL**

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

List with at least 4 DFs: 1-3) from-to-to: 1. MappingIssues, 2.
MappingIssues summary, 3. Long summary (If summary=TRUE) & 4-6)
to-to-from: 4. MappingIssues, 5. MappingIssues summary, 6. Long summary
(If summary=TRUE) & 7) Combined summary table (If summary=TRUE)

## Examples

``` r
if (FALSE) { # \dontrun{
KEGG_Pathways <- metsigdb_kegg()
InputDF <- translate_id(
    data = KEGG_Pathways,
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    from = c("kegg"),
    to = c("pubchem")
)[["TranslatedDF"]]
Res <- mapping_ambiguity(
    data = InputDF,
    from = "MetaboliteID",
    to = "pubchem",
    grouping_variable = "term",
    summary = TRUE
)
} # }
```
