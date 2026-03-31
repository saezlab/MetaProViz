# Retrieve WikiPathways metabolite mapping suitable for ORA.

Retrieves pathway to metabolite mappings from WikiPathways (via
`wikipathways_metabolites_sparql()`) via OmnipathR and returns a
long-format table with one metabolite identifier per row.

## Usage

``` r
metsigdb_wikipathways(species = "Homo sapiens", exclude_metabolites = "all")
```

## Arguments

- species:

  Character. Species name. Default is `"Homo sapiens"`.

- exclude_metabolites:

  Optional metabolite classes to exclude: NULL (exclude nothing), "all"
  (default), or any combination of c("ions", "small_molecules",
  "xenobiotics", "atoms").

## Value

A tibble in long format with columns `pathway_id`, `pathway_name`,
`pathway_url`, `n_metabolites_in_pathway`, and `metabolite_id`.
