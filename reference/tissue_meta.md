# Tissue_Metadata

In Hakimi et. al. metabolites were assigned to metabolite IDs, pathways,
platform, mass and other fetaure metainformation. row names being
metabolite trivial names.
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_meta
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 877 rows and 11 columns.

## Examples

``` r
data(tissue_meta)
head(tissue_meta)
#> # A tibble: 6 × 11
#>   Metabolite        SUPER_PATHWAY SUB_PATHWAY COMP_ID PLATFORM    RI  MASS CAS  
#>   <chr>             <chr>         <chr>         <dbl> <chr>    <dbl> <dbl> <chr>
#> 1 1,2-propanediol   Lipid         Ketone bod…   38002 GC/MS    1041   117  57-5…
#> 2 1,3-dihydroxyace… Carbohydrate  Glycolysis…   35963 GC/MS    1263   103  6214…
#> 3 1,5-anhydrogluci… Carbohydrate  Glycolysis…   20675 GC/MS    1789.  217  154-…
#> 4 10-heptadecenoat… Lipid         Long chain…   33971 LC/MS N… 5558   267. 2974…
#> 5 10-nonadecenoate… Lipid         Long chain…   33972 LC/MS N… 5775   295. 7303…
#> 6 13-HODE + 9-HODE  Lipid         Fatty acid…   37752 LC/MS N… 5247   295. NA   
#> # ℹ 3 more variables: PUBCHEM <dbl>, KEGG <chr>, Group_HMDB <chr>
```
