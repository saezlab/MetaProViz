# Check compatibility of seed ID pairs in input rows

Creates a long-format permutation table where each row represents one
unique unordered pair of seed IDs from the same input row, then flags
whether each pair is compatible via direct or secondary graph
connections.

## Usage

``` r
seed_id_compatibility_check(
  data,
  id_types = c("HMDB", "KEGG", "CHEBI", "PUBCHEM"),
  delimiter = c(";", ","),
  verbose = FALSE,
  edge_table = NULL,
  handle_partially_compatible = FALSE,
  handle_completely_incompatible = FALSE,
  completely_incompatible_priority = c("HMDB", "CHEBI", "PUBCHEM", "KEGG")
)
```

## Arguments

- data:

  Data frame with zero or more of the columns `HMDB`, `KEGG`, `CHEBI`,
  and `PUBCHEM`. Column names are matched case-insensitively against
  these exact names.

- id_types:

  Character vector of ID types to use. Choose from `HMDB`, `KEGG`,
  `CHEBI`, and `PUBCHEM`.

- delimiter:

  Character string indicating whether multiple IDs within one cell are
  separated by semicolons or commas. Accepted values are `";"`, `","`,
  `"semicolon"`, or `"comma"`.

- verbose:

  Logical; if `TRUE`, prints pairwise mapping and edge construction
  diagnostics to the console.

- edge_table:

  Optional precomputed bidirectional edge table with columns `id1`,
  `type1`, `id2`, `type2`. If `NULL`, the table is built internally.

- handle_partially_compatible:

  Logical; if `TRUE`, partially compatible features are cleaned by
  retaining only IDs from compatible pairs.

- handle_completely_incompatible:

  Logical; if `TRUE`, completely incompatible features are cleaned by
  retaining a single ID according to `completely_incompatible_priority`.

- completely_incompatible_priority:

  Character vector defining the namespace priority for resolving
  completely incompatible features. Supported values are `HMDB`, `KEGG`,
  `CHEBI`, and `PUBCHEM`. The default priority is
  `c("HMDB", "CHEBI", "PUBCHEM", "KEGG")`.

## Value

Named list with at least five data frames:

- ID_pair_compatibility:

  Long-format table with one unique unordered seed-ID pair per input
  row. The first column `original_row_id` stores the original input row
  name. The table also includes `pair_compatible`, `compatibility_path`
  (`direct`, `secondary`, `no_match`), and grouped
  `all_seed_ids_compatible`.

- data_with_compatibility:

  Original input data with appended `all_seed_ids_compatible` per input
  row (rows with fewer than two seed IDs are `TRUE`).

- feature_compatibility_summary:

  One row per input feature summarizing compatibility counts and
  assigning `fully_compatible`, `partially_compatible`, or
  `completely_incompatible`.

- data_after_handling:

  Feature-level table after optional automatic handling. If no handling
  is enabled, this matches the raw feature-level output aside from added
  summary columns.

- ID_pair_compatibility_after_handling:

  Pair-level compatibility table recomputed from `data_after_handling`.

If either handling option is enabled, the return object also includes
`handling_summary_text` and `handling_summary_metrics`.
