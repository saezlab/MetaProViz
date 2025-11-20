# tissue_dma_old

We performed differential metabolite analysis comparing ccRCC tissue
versus adjacent normal tissue of the patient's subset of old patient's
(age \> 58 years) using median normalised data from the supplementary
table 2 of Hakimi et. al.(="Tissue_Norm"). metabolite values used as
input with row names being metabolite trivial names.
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_dma_old
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 570 rows and 17 columns.

## Examples

``` r
data(tissue_dma_old)
head(tissue_dma_old)
#> # A tibble: 6 × 17
#>   Metabolite         Log2FC AveExpr  t.val    p.val    p.adj     B SUPER_PATHWAY
#>   <chr>               <dbl>   <dbl>  <dbl>    <dbl>    <dbl> <dbl> <chr>        
#> 1 1-arachidonoylgly… -0.663 -0.0529  -3.34 1.01e- 3 1.72e- 3 -2.23 Lipid        
#> 2 1-arachidonoylgly… -1.26  -0.0300 -12.4  7.73e-26 9.26e-25 47.9  Lipid        
#> 3 1-arachidonoylgly… -1.01  -0.165   -9.63 6.57e-18 4.02e-17 29.7  Lipid        
#> 4 1-arachidonylglyc… -0.256 -1.24    -2.03 4.36e- 2 5.79e- 2 -5.62 Lipid        
#> 5 1-docosahexaenoyl… -1.13  -0.617   -6.60 4.61e-10 1.50e- 9 11.9  Lipid        
#> 6 1-heptadecanoylgl… -1.83  -0.420   -9.22 9.04e-17 5.05e-16 27.1  Lipid        
#> # ℹ 9 more variables: SUB_PATHWAY <chr>, COMP_ID <dbl>, PLATFORM <chr>,
#> #   RI <dbl>, MASS <dbl>, CAS <chr>, PUBCHEM <dbl>, KEGG <chr>,
#> #   Group_HMDB <chr>
```
