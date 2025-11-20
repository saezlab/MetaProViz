# mca_twocond_rules

Manually curated table defining the flow of information of the two
condition biological regulatory clusters Regulatory labels from the
different grouping methods. columns (RG1-RG3)

## Usage

``` r
mca_twocond_rules
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 36 rows and 5 columns.

## Examples

``` r
data(mca_twocond_rules)
head(mca_twocond_rules)
#> # A tibble: 6 × 5
#>   Cond1 Cond2                RG1_All       RG2_Significant RG3_SignificantChange
#>   <chr> <chr>                <chr>         <chr>           <chr>                
#> 1 DOWN  DOWN                 Cond1 DOWN +… Core_DOWN       Core_DOWN            
#> 2 DOWN  Not Detected         Cond1 DOWN +… Cond1_DOWN      Cond1_DOWN           
#> 3 DOWN  Not Significant      Cond1 DOWN +… Cond1_DOWN      Cond1_DOWN           
#> 4 DOWN  Significant Negative Cond1 DOWN +… Core_DOWN       Cond1_DOWN           
#> 5 DOWN  Significant Positive Cond1 DOWN +… Opposite        Cond1_DOWN           
#> 6 DOWN  UP                   Cond1 DOWN +… Opposite        Opposite             
```
