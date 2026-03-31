# Expand metabolite IDs by traversing RaMP ID mappings

Traverses pairwise RaMP mappings from
[`OmnipathR::ramp_id_mapping_table()`](https://saezlab.github.io/OmnipathR/reference/ramp_id_mapping_table.html)
across selected metabolite ID types until no new IDs are found.

## Usage

``` r
traverse_ids(
  data,
  id_types = c("HMDB", "KEGG", "CHEBI", "PUBCHEM"),
  delimiter = c(";", ","),
  save_table = "csv",
  path = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  Data frame with zero or more of the columns `HMDB`, `KEGG`, `CHEBI`,
  and `PUBCHEM`. Column names are matched case-insensitively against
  these exact names.

- id_types:

  Character vector of ID types to expand. Choose from `HMDB`, `KEGG`,
  `CHEBI`, and `PUBCHEM`.

- delimiter:

  Character string indicating whether multiple IDs within one cell are
  separated by semicolons or commas. Accepted values are `";"`, `","`,
  `"semicolon"`, or `"comma"`.

- save_table:

  *Optional:* File types for the analysis results are: `"csv"`,
  `"xlsx"`, `"txt"`. If `NULL`, no tables are saved. **Default = "csv"**

- path:

  *Optional:* Path to the folder the results should be saved at.
  **Default = NULL**

- verbose:

  Logical; if `TRUE`, prints pairwise mapping and edge construction
  diagnostics to the console.

## Value

Named list with three data frames:

- ExpandedDF:

  Input data with appended expanded ID columns and QC summary columns,
  including `all_seed_ids_compatible` (logical flag indicating whether
  all seed-ID pairs in each row are compatible).

- ID_pair_compatibility:

  Long-format table with one unique unordered seed-ID pair per input
  row. The first column `original_row_id` stores the original input row
  name. The table also includes `pair_compatible`, `compatibility_path`,
  and `all_seed_ids_compatible`.

- ID_Edges_prior_knowledge:

  Bidirectional ID edge table used for traversal and compatibility
  checks.

## Examples

``` r
input_df <- data.frame(
    name = c(
        "Acetone ; Propanal ; acetone",
        "Acetaldehyde oxime ; HMDB01122",
        "acetate",
        "Urea"
    ),
    all_ids = c(
        "HMDB01659 ; HMDB03366 ; C00207",
        "HMDB03656 ; HMDB01122",
        "C00033",
        "C00086"
     ),
     HMDB = c(
        "HMDB01659; HMDB03366",
        "HMDB03656; HMDB01122",
        NA,
        NA
     ),
     KEGG = c(
        "C00207",
        NA,
        "C00033",
        "C00086"
     ),
     CHEBI = NA,
     stringsAsFactors = FALSE
)

res <- traverse_ids(input_df)
#> Warning: Selected ID column 'PUBCHEM' not found in data and was created as NA.
#> Warning: Detected incompatible seed IDs in 2 row(s) (row_id: 1, 2). Not all seed IDs in these rows are mutually reachable through ID_Edges_prior_knowledge, suggesting they may map to different molecules. This can overexpand the ID space during traversal. Please manually remove or correct incompatible seed IDs and rerun traverse_ids() freshly on a clean input table.

df_translated <- res$ExpandedDF
head(df_translated)
#> # A tibble: 4 × 20
#>   name           all_ids HMDB  KEGG  CHEBI PUBCHEM all_seed_ids_compati…¹ row_id
#>   <chr>          <chr>   <chr> <chr> <lgl> <chr>   <lgl>                   <int>
#> 1 Acetone ; Pro… HMDB01… HMDB… C002… NA    NA      FALSE                       1
#> 2 Acetaldehyde … HMDB03… HMDB… NA    NA    NA      FALSE                       2
#> 3 acetate        C00033  NA    C000… NA    NA      TRUE                        3
#> 4 Urea           C00086  NA    C000… NA    NA      TRUE                        4
#> # ℹ abbreviated name: ¹​all_seed_ids_compatible
#> # ℹ 12 more variables: HMDB_translated <chr>, KEGG_translated <chr>,
#> #   CHEBI_translated <chr>, PUBCHEM_translated <chr>, n_seed_ids <int>,
#> #   n_HMDB_translated <int>, n_KEGG_translated <int>, n_CHEBI_translated <int>,
#> #   n_PUBCHEM_translated <int>, mapping_expanded <lgl>, ambiguous_seed <lgl>,
#> #   large_mapping <lgl>
```
