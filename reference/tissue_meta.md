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

An object of class `tbl_df` (inherits from `tbl`, `data.frame`) with 877
rows and 11 columns.

## Examples

``` r
data(tissue_meta)
head(tissue_meta)
#> # A tibble: 6 × 11
#>   Metabolite        CAS   SUPER_PATHWAY SUB_PATHWAY COMP_ID PLATFORM RI    MASS 
#>   <chr>             <chr> <chr>         <chr>       <chr>   <chr>    <chr> <chr>
#> 1 1,2-propanediol   57-5… Lipid         Ketone bod… 38002   GC/MS    1041  117  
#> 2 1,3-dihydroxyace… 6214… Carbohydrate  Glycolysis… 35963   GC/MS    1263  103  
#> 3 1,5-anhydrogluci… 154-… Carbohydrate  Glycolysis… 20675   GC/MS    1788… 217  
#> 4 1-arachidonoylgl… NA    Lipid         Lysolipid   33228   LC/MS P… 5554  544.…
#> 5 1-arachidonoylgl… NA    Lipid         Lysolipid   35186   LC/MS N… 5731  500.3
#> 6 1-arachidonoylgl… NA    Lipid         Lysolipid   34214   LC/MS N… 5479  619.4
#> # ℹ 3 more variables: PUBCHEM <chr>, KEGG <chr>, HMDB <chr>
```
