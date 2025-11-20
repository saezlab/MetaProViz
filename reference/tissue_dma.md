# tissue_dma

We performed differential metabolite analysis comparing ccRCC tissue
versus adjacent normal tissue using median normalised data from the
supplementary table 2 of Hakimi et. al.(="Tissue_Norm"). metabolite
values used as input with row names being metabolite trivial names.
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_dma
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 570 rows and 17 columns.

## Examples

``` r
data(tissue_dma)
head(tissue_dma)
#> # A tibble: 6 × 17
#>   Metabolite         Log2FC AveExpr  t.val    p.val    p.adj     B SUPER_PATHWAY
#>   <chr>               <dbl>   <dbl>  <dbl>    <dbl>    <dbl> <dbl> <chr>        
#> 1 1-arachidonoylgly… -0.433  0.0516  -2.64 8.70e- 3 1.22e- 2 -4.44 Lipid        
#> 2 1-arachidonoylgly… -1.26  -0.0204 -15.4  4.88e-39 6.04e-38 77.8  Lipid        
#> 3 1-arachidonoylgly… -0.968 -0.127  -11.5  2.95e-25 1.70e-24 46.2  Lipid        
#> 4 1-arachidonylglyc… -0.247 -1.24    -2.46 1.43e- 2 1.94e- 2 -4.89 Lipid        
#> 5 1-docosahexaenoyl… -0.915 -0.558   -6.52 3.27e-10 8.44e-10 11.9  Lipid        
#> 6 1-heptadecanoylgl… -1.62  -0.353   -9.86 7.51e-20 3.48e-19 33.8  Lipid        
#> # ℹ 9 more variables: SUB_PATHWAY <chr>, COMP_ID <dbl>, PLATFORM <chr>,
#> #   RI <dbl>, MASS <dbl>, CAS <chr>, PUBCHEM <dbl>, KEGG <chr>,
#> #   Group_HMDB <chr>
```
