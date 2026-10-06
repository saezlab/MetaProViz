# Process feature-metadata metabolite IDs in a staged workflow

Runs an ID-space QC and expansion workflow on feature metadata using the
exported MetaProViz ID helper functions. The workflow always starts with
ID exploration using
[`compare_pk()`](https://saezlab.github.io/MetaProViz/reference/compare_pk.md)
and
[`count_id()`](https://saezlab.github.io/MetaProViz/reference/count_id.md),
can optionally perform automatic seed-ID compatibility handling, and can
then run either
[`translate_id()`](https://saezlab.github.io/MetaProViz/reference/translate_id.md)
or
[`traverse_ids()`](https://saezlab.github.io/MetaProViz/reference/traverse_ids.md).

## Usage

``` r
id_processing(
  data,
  id_types = c("HMDB", "KEGG", "CHEBI", "PUBCHEM"),
  delimiter = c(";", ","),
  run_compatibility_check = FALSE,
  handle_partially_compatible = TRUE,
  handle_completely_incompatible = TRUE,
  completely_incompatible_priority = c("HMDB", "CHEBI", "PUBCHEM", "KEGG"),
  run_translation = FALSE,
  translation_from = NULL,
  translation_to = NULL,
  translation_summary = FALSE,
  run_traversal = FALSE,
  edge_table = NULL,
  compare_name_col = "Metabolite",
  save_plot = "svg",
  save_table = "csv",
  print_plot = TRUE,
  verbose = FALSE,
  path = NULL
)
```

## Arguments

- data:

  Data frame with feature metadata and zero or more of the columns
  `HMDB`, `KEGG`, `CHEBI`, and `PUBCHEM`. Column names are matched
  case-insensitively to these canonical names.

- id_types:

  Character vector of workflow-managed ID namespaces. Supported values
  are `HMDB`, `KEGG`, `CHEBI`, and `PUBCHEM`.

- delimiter:

  Character string indicating whether multiple IDs within one cell are
  separated by semicolons or commas. Accepted values are `";"`, `","`,
  `"semicolon"`, or `"comma"`.

- run_compatibility_check:

  Logical; if `TRUE`, run
  [`seed_id_compatibility_check()`](https://saezlab.github.io/MetaProViz/reference/seed_id_compatibility_check.md)
  before any expansion step. Default is `FALSE`.

- handle_partially_compatible:

  Logical; forwarded to
  [`seed_id_compatibility_check()`](https://saezlab.github.io/MetaProViz/reference/seed_id_compatibility_check.md).

- handle_completely_incompatible:

  Logical; forwarded to
  [`seed_id_compatibility_check()`](https://saezlab.github.io/MetaProViz/reference/seed_id_compatibility_check.md).

- completely_incompatible_priority:

  Character vector defining the namespace priority for resolving
  completely incompatible features.

- run_translation:

  Logical; if `TRUE`, run a translation step after the compatibility
  step. Mutually exclusive with `run_traversal`.

- translation_from:

  Source namespace for
  [`translate_id()`](https://saezlab.github.io/MetaProViz/reference/translate_id.md).
  Required when `run_translation = TRUE`.

- translation_to:

  One or more target namespaces for
  [`translate_id()`](https://saezlab.github.io/MetaProViz/reference/translate_id.md).
  Required when `run_translation = TRUE`.

- translation_summary:

  Logical; forwarded to
  [`translate_id()`](https://saezlab.github.io/MetaProViz/reference/translate_id.md).

- run_traversal:

  Logical; if `TRUE`, run
  [`traverse_ids()`](https://saezlab.github.io/MetaProViz/reference/traverse_ids.md)
  after the compatibility step. Mutually exclusive with
  `run_translation`. The seed-ID compatibility check inside
  [`traverse_ids()`](https://saezlab.github.io/MetaProViz/reference/traverse_ids.md)
  is always skipped here; use `run_compatibility_check` to check seed
  IDs before traversal. Default is `FALSE`.

- edge_table:

  Optional precomputed bidirectional edge table with columns `id1`,
  `type1`, `id2`, `type2`. If `NULL`, it is built internally when
  required for compatibility checking or traversal.

- compare_name_col:

  Optional feature-name column passed to
  [`compare_pk()`](https://saezlab.github.io/MetaProViz/reference/compare_pk.md)
  during stage-wise ID-space exploration. If missing from `data`, it is
  ignored automatically by
  [`compare_pk()`](https://saezlab.github.io/MetaProViz/reference/compare_pk.md).

- save_plot:

  Optional plot file type: `"svg"`, `"png"`, or `"pdf"`. If `NULL`,
  plots are not saved.

- save_table:

  Optional table file type: `"csv"`, `"xlsx"`, or `"txt"`. If `NULL`,
  tables are not saved.

- print_plot:

  Logical; whether saved plots should also be printed by `save_res()`.

- verbose:

  Logical; forwarded to compatibility and traversal helpers to control
  their detailed logging. `id_processing()` always prints its own
  workflow progress and result overview. Default is `FALSE`.

- path:

  Optional path where results should be saved.

## Value

Named list with three top-level entries:

- Data:

  Feature-metadata tables for each workflow stage, the aggregated
  ID-count summary, stage-wise QC tables, and any raw compatibility,
  translation, or traversal tables produced by enabled steps.

- Plot:

  Stage-wise
  [`compare_pk()`](https://saezlab.github.io/MetaProViz/reference/compare_pk.md)
  and
  [`count_id()`](https://saezlab.github.io/MetaProViz/reference/count_id.md)
  plots.

- Workflow:

  Chosen settings, steps run, collected stage messages, and a final
  result-overview text.

## Examples

``` r
data(tissue_meta)

qc_only <- id_processing(
    data = tissue_meta[seq_len(min(40, nrow(tissue_meta))), , drop = FALSE],
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE,
    verbose = FALSE
)
#> [id_processing] Starting workflow on 40 feature(s).
#> [id_processing] Selected namespaces: HMDB, KEGG, CHEBI, PUBCHEM
#> [id_processing] Steps enabled: compatibility=FALSE, translation=FALSE, traversal=FALSE
#> [id_processing] Defaults/assumptions: partial_auto=TRUE, complete_auto=TRUE, complete_priority=HMDB > CHEBI > PUBCHEM > KEGG
#> [id_processing] Quantifying ID space for stage 'input'.
#> [id_processing] Stage 'input': HMDB total_ids=14 no_id=26 single=14 multiple=0 | KEGG total_ids=9 no_id=31 single=9 multiple=0 | CHEBI total_ids=0 no_id=40 single=0 multiple=0 | PUBCHEM total_ids=22 no_id=18 single=22 multiple=0
#> [id_processing] Workflow steps run: initial_exploration. Returned Data tables: input. Returned Plot stages: input. Final stage 'input' contains 40 feature(s).
#> [id_processing] Suggestion: equivalent_id() is not part of this workflow yet. If you want additional ambiguity-aware within-namespace expansion, run it afterwards on the final feature metadata.
```
