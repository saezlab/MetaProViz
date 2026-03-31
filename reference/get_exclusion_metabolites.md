# Metabolites excluded from prior knowledge resources

Builds an exclusion table from hard-coded seed identifiers and
translates them across ID systems (HMDB, KEGG, CHEBI, PUBCHEM).

## Usage

``` r
get_exclusion_metabolites()
```

## Value

A data frame with columns `metabolite_id`, `class`, and `id_type`.

## Examples

``` r
if (FALSE) { # \dontrun{
get_exclusion_metabolites()
} # }
```
