# intracell_raw

Metabolomics workbench project PR001418, study ST002224 where we
exported integrated raw peak values of intracellular metabolomics of HK2
and ccRCC cell lines 786-O, 786-M1A and 786-M2A. -`Conditions`:
Character vector indicating cell line identity - `Analytical_Replicate`:
Integer replicate number for analytical replicates
-`Biological_Replicate`: Integer replicate number for biological
replicates - Additional numeric columns (183 in total) containing raw
metabolite intensities nitrogen supports renal cancer progression,
Nature Communications 2022, DOI:10.1038/s41467-022-35036-4.

## Usage

``` r
intracell_raw
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 58 rows and 186 columns.

## Examples

``` r
data(intracell_raw)
head(intracell_raw)
#> # A tibble: 6 × 186
#>   Code    Conditions Analytical_Replicates Biological_Replicates `valine-d8`
#>   <chr>   <chr>                      <dbl>                 <dbl>       <dbl>
#> 1 MS55_01 HK2                            1                     1  1910140239
#> 2 MS55_02 HK2                            2                     1  2030901280
#> 3 MS55_03 HK2                            3                     1  2001950756
#> 4 MS55_04 HK2                            4                     1  1971520079
#> 5 MS55_05 786-O                          1                     1  2150817213
#> 6 MS55_06 786-O                          2                     1  2041935934
#> # ℹ 181 more variables: `hippuric acid-d5` <dbl>, `2/3-phosphoglycerate` <dbl>,
#> #   `2-aminoadipic acid` <dbl>, `2-hydroxyglutarate` <dbl>,
#> #   `2-ketoglutarate` <dbl>, `4-guanidinobutanoate` <dbl>,
#> #   `4-hydroxyphenyllactate` <dbl>, `5-aminolevulinic acid` <dbl>,
#> #   `5-aminovaleric acid*` <dbl>, `5-methylthioadenosine (MTA)` <dbl>,
#> #   `6-phosphogluconate` <dbl>, acetylcarnitine <dbl>, acetylcholine <dbl>,
#> #   `acetyl-CoA` <dbl>, adenosine <dbl>, ADP <dbl>, `ADP-ribose` <dbl>, …
```
