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
  edge_table = NULL
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

## Value

Named list with two data frames:

- ID_pair_compatibility:

  Long-format table with one unique unordered seed-ID pair per input
  row. The first column `original_row_id` stores the original input row
  name. The table also includes `pair_compatible`, `compatibility_path`
  (`direct`, `secondary`, `no_match`), and grouped
  `all_seed_ids_compatible`.

- data_with_compatibility:

  Original input data with appended `all_seed_ids_compatible` per input
  row (rows with fewer than two seed IDs are `TRUE`).
