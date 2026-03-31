# Retrieve MACDB metabolite-cancer associations.

Retrieves metabolite-cancer associations from MACDB via OmnipathR and
summarizes them to unique metabolite-cancer associations (by `term` and
`Metabolite_PubchemID`).

## Usage

``` r
metsigdb_macdb(exclude_metabolites = "all")
```

## Arguments

- exclude_metabolites:

  Optional metabolite classes to exclude: NULL (exclude nothing), "all"
  (default), or any combination of c("ions", "small_molecules",
  "xenobiotics", "atoms").

## Value

A data frame with one row per unique metabolite-cancer association,
collapsed metadata columns, and summary metrics (`evidence_count`,
`significance_count`, `association_score`).
