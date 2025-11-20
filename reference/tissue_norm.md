# tissue_norm

This is median normalised data from the supplementary table 2 of Hakimi
et al with metabolomic profiling on 138 matched clear cell renal cell
carcinoma (ccRCC)/normal tissue pairs. measured metabolite (normalised
data)
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_norm
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 276 rows and 584 columns.

## Examples

``` r
data(tissue_norm)
head(tissue_norm)
#> # A tibble: 6 × 584
#>   Code  CLIENT_IDENTIFIER TISSUE_TYPE MATCHING_INDEX GENDER RACE  AGE_AT_SURGERY
#>   <chr>             <dbl> <chr>                <dbl> <chr>  <chr>          <dbl>
#> 1 DIAG…                 3 TUMOR                    1 Male   White           74.8
#> 2 DIAG…                 4 NORMAL                   1 Male   White           74.8
#> 3 DIAG…                 5 TUMOR                    2 Male   White           77.2
#> 4 DIAG…                 6 NORMAL                   2 Male   White           77.2
#> 5 DIAG…                 7 TUMOR                    3 Female White           59.1
#> 6 DIAG…                 8 NORMAL                   3 Female White           59.1
#> # ℹ 577 more variables: PATHSTAGET_1 <chr>, `TYPE-STAGE` <chr>,
#> #   PATHSTAGEN_1 <chr>, PATHSTAGEM_1 <chr>, PATHGRADE1 <chr>, STAGE <chr>,
#> #   AGE <chr>, `1,2-propanediol` <dbl>, `1,3-dihydroxyacetone` <dbl>,
#> #   `1,5-anhydroglucitol (1,5-AG)` <dbl>, `10-heptadecenoate (17:1n7)` <dbl>,
#> #   `10-nonadecenoate (19:1n9)` <dbl>, `13-HODE + 9-HODE` <dbl>,
#> #   `15-methylpalmitate (isobar with 2-methylpalmitate)` <dbl>,
#> #   `16-hydroxypalmitate` <dbl>, …
```
