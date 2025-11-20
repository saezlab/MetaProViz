# mca_core_rules

Manually curated table defining the flow of information of the
Conusuption-Release and Intracellular metabolomics biological regulatory
clusters Regulatory labels from the different grouping methods. columns
(RG1-RG3)

## Usage

``` r
mca_core_rules
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 108 rows and 7 columns.

## Examples

``` r
data(mca_core_rules)
head(mca_core_rules)
#> # A tibble: 6 × 7
#>   Intra CoRe      Core_Direction RG1_All R2_Significant RG3_Change RG3_ChangeNum
#>   <chr> <chr>     <chr>          <chr>   <chr>          <chr>              <dbl>
#> 1 DOWN  DOWN      Released       Intra … Both_DOWN (Re… Both_DOWN…             1
#> 2 DOWN  Not Dete… Not Detected   Intra … None           None                   2
#> 3 DOWN  Not Sign… Released       Intra … None           None                   2
#> 4 DOWN  Signific… Released       Intra … Both_DOWN (Re… None                   2
#> 5 DOWN  Signific… Released       Intra … Opposite (Rel… None                   2
#> 6 DOWN  UP        Released       Intra … Opposite (Rel… Opposite …             3
```
