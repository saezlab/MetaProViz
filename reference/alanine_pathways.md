# alanine_pathways

Manually curated table for the amino acid alanine toshowcase pathways
(wiki, reactome, etc.) and alanine IDs (chebi, hmdb, etc.) included in
those pathways

## Usage

``` r
alanine_pathways
```

## Format

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 204
rows and 5 columns.

## Examples

``` r
data(alanine_pathways)
head(alanine_pathways)
#> # A tibble: 6 × 5
#>   pathwayName                      pathwaySource pathwayId inputId    commonName
#>   <chr>                            <chr>         <chr>     <chr>      <chr>     
#> 1 Alanine and aspartate metabolism wiki          WP106     chebi:155… D-Alanine…
#> 2 Alanine and aspartate metabolism wiki          WP106     chebi:169… L-Alanine 
#> 3 Alanine and aspartate metabolism wiki          WP106     hmdb:HMDB… L-Alanine 
#> 4 Alanine and aspartate metabolism wiki          WP106     hmdb:HMDB… D-Alanine…
#> 5 Alanine and aspartate metabolism wiki          WP106     kegg:C000… L-Alanine 
#> 6 Alanine and aspartate metabolism wiki          WP106     kegg:C001… D-Alanine…
```
