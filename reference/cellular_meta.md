# cellular_meta

Metabolomics workbench project PR001418, study ST002226 and ST002224
measured metabolites were assigned HMDB and KEGG IDs as well as one main
metabolic pathway. metabolite trivial names. nitrogen supports renal
cancer progression , Nature Communications 2022,
[doi:10.1038/s41467-022-35036-4](https://doi.org/10.1038/s41467-022-35036-4)

## Usage

``` r
cellular_meta
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 199 rows and 5 columns.

## Examples

``` r
data(cellular_meta)
head(cellular_meta)
#> # A tibble: 6 × 5
#>   Metabolite                HMDB        KEGG.ID KEGGCompound             Pathway
#>   <chr>                     <chr>       <chr>   <chr>                    <chr>  
#> 1 N-acetylaspartate         HMDB0000812 C01042  N-Acetyl-L-aspartate     Alanin…
#> 2 argininosuccinate         HMDB0000052 C03406  N-(L-Arginino)succinate  Alanin…
#> 3 N-acetylaspartylglutamate HMDB0001067 C12270  N-Acetylaspartylglutama… Alanin…
#> 4 tyrosine                  HMDB0000158 C00082  L-Tyrosine               Amino …
#> 5 asparagine                HMDB0000168 C00152  L-Asparagine             Amino …
#> 6 glutamate                 HMDB0000148 C00025  L-Glutamate              Amino …
```
