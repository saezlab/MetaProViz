# Retrieve Reactome metabolite sets suitable for ORA.

Queries the OmniPath resource through OmniPathR to obtail Reactome
pathway level metabolite sets.

## Usage

``` r
metsigdb_reactome(
  species = "Homo sapiens",
  pathway_ids = NULL,
  out_path = NULL,
  exclude_metabolites = "all"
)
```

## Arguments

- species:

  String. Optionally specify pathways to query from a species via full
  name or three letter code. Default = "Homo sapiens". NULL for all
  species.

- pathway_ids:

  String vector. Optionally provide pathway_ids to query. Default NULL
  to query all pathways.

- out_path:

  String. Optionally save results as csv into out_path. Default NULL.

- exclude_metabolites:

  Optional metabolite classes to exclude: NULL (exclude nothing), "all"
  (default), or any combination of c("ions", "small_molecules",
  "xenobiotics", "atoms").

## Value

A tibble in long format containing one row per metabolite for the
Reactome pathways.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- metsigdb_reactome()
head(df)
} # }
```
