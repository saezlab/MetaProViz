# tissue_dma_young

We performed differential metabolite analysis comparing ccRCC tissue
versus adjacent normal tissue of the patient's subset of young patient's
(age \<42 years) using median normalised data from the supplementary
table 2 of Hakimi et. al.(="Tissue_Norm"). metabolite values used as
input with row names being metabolite trivial names.
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_dma_young
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 570 rows and 17 columns.

## Examples

``` r
data(tissue_dma_young)
head(tissue_dma_young)
#> # A tibble: 6 × 17
#>   Metabolite          Log2FC AveExpr  t.val   p.val   p.adj      B SUPER_PATHWAY
#>   <chr>                <dbl>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl> <chr>        
#> 1 1-arachidonoylgly… -0.251   0.257  -0.332 0.745   0.902   -6.45  Lipid        
#> 2 1-arachidonoylgly… -1.39    0.170  -4.07  0.00102 0.00763 -0.753 Lipid        
#> 3 1-arachidonoylgly… -1.16    0.109  -3.24  0.00550 0.0257  -2.38  Lipid        
#> 4 1-arachidonylglyc… -0.0922 -1.31   -0.252 0.804   0.942   -6.48  Lipid        
#> 5 1-docosahexaenoyl… -1.29   -0.421  -1.80  0.0929  0.216   -4.99  Lipid        
#> 6 1-heptadecanoylgl… -1.93    0.0472 -2.68  0.0172  0.0631  -3.46  Lipid        
#> # ℹ 9 more variables: SUB_PATHWAY <chr>, COMP_ID <dbl>, PLATFORM <chr>,
#> #   RI <dbl>, MASS <dbl>, CAS <chr>, PUBCHEM <dbl>, KEGG <chr>,
#> #   Group_HMDB <chr>
```
