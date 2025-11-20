# medium_raw

Metabolomics workbench project PR001418, study ST002226 where we
exported integrated raw peak values of intracellular metabolomics of HK2
and cccRCC cell lines 786-O, 786-M1A, 786-M2A, OS-RC-2, OS-LM1 and
RFX-631. a numeric column for each measured metabolite (raw data)
nitrogen supports renal cancer progression , Nature Communications 2022,
[doi:10.1038/s41467-022-35036-4](https://doi.org/10.1038/s41467-022-35036-4)

## Usage

``` r
medium_raw
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 44 rows and 77 columns.

## Examples

``` r
data(medium_raw)
head(medium_raw)
#> # A tibble: 6 × 77
#>   Code    Conditions Biological_Replicates GrowthFactor `valine-d8`
#>   <chr>   <chr>                      <dbl>        <dbl>       <dbl>
#> 1 MS51-01 blank                          1          NA   805641492.
#> 2 MS51-02 blank                          2          NA   774207449.
#> 3 MS51-03 blank                          3          NA   806230116.
#> 4 MS51-04 blank                          4          NA   799304084.
#> 5 MS51-05 blank                          5          NA   765139228.
#> 6 MS51-06 HK2                            1         249.  780552871.
#> # ℹ 72 more variables: `hipppuric acid-d5` <dbl>, `2-hydroxyglutarate` <dbl>,
#> #   `2-ketoglutarate` <dbl>, `3-Dehydro-L-threonate` <dbl>,
#> #   acetylcarnitine <dbl>, acetylcholine <dbl>, acetylornithine <dbl>,
#> #   aconitate <dbl>, alanine <dbl>, arginine <dbl>, asparagine <dbl>,
#> #   aspartate <dbl>, betaine <dbl>, `butyryl-carnitine` <dbl>,
#> #   `carbamoyl phosphate` <dbl>, carnitine <dbl>, citrulline <dbl>,
#> #   creatine <dbl>, creatinine <dbl>, cystine <dbl>, dihydroorotate <dbl>, …
```
