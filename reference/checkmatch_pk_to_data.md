# Check and summarize relationship between prixor knowledge to measured

features

## Usage

``` r
checkmatch_pk_to_data(
  data,
  input_pk,
  metadata_info = c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term"),
  save_table = "csv",
  path = NULL
)
```

## Arguments

- data:

  dataframe with at least one column with the detected metabolite IDs
  (e.g. HMDB). If there are multiple IDs per detected peak, please
  separate them by comma ("," or ", " or chr list). If there is a main
  ID and additional IDs, please provide them in separate columns.

- input_pk:

  dataframe with at least one column with the metabolite ID (e.g. HMDB)
  that need to match data metabolite IDs "source" (e.g. term). If there
  are multiple IDs, as the original pathway IDs (e.g. KEGG) where
  translated (e.g. to HMDB), please separate them by comma ("," or ", "
  or chr list).

- metadata_info:

  Colum name of Metabolite IDs in data and input_pk as well as column
  name of grouping_variable in input_pk. **Default = c(InputID="HMDB",
  PriorID="HMDB", grouping_variable="term")**

- save_table:

  *Optional:* File types for the analysis results are: "csv", "xlsx",
  "txt". **Default = "csv"**

- path:

  Optional: Path to the folder the results should be saved at. **Default
  = NULL**

## Value

A `list` with three elements:

- `data_summary` — a data frame summarising matching results per input
  ID, including counts, conflicts, and recommended actions.

- `GroupingVariable_summary` — a detailed data frame showing matches
  grouped by the specified variable, with conflict annotations.

- `data_long` — a merged data frame of prior knowledge IDs and detected
  IDs in long format.

## Examples

``` r
if (FALSE) { # \dontrun{
data(cellular_meta)
DetectedIDs <- cellular_meta %>%
    dplyr::select("Metabolite", "HMDB") %>%
    tidyr::drop_na()
input_pathway <- translate_id(
    data = metsigdb_kegg(),
    metadata_info = c(
        InputID = "MetaboliteID",
        grouping_variable = "term"
    ),
    from = c("kegg"),
    to = c("hmdb")
)[["TranslatedDF"]] %>% tidyr::drop_na()
Res <- checkmatch_pk_to_data(
    data = DetectedIDs,
    input_pk = input_pathway,
    metadata_info = c(
        InputID = "HMDB",
        PriorID = "hmdb",
        grouping_variable = "term"
    )
)
} # }
```
