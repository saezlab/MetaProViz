# tissue_norm_se

This is median normalised data from the supplementary table 2 of Hakimi
et al with metabolomic profiling on 138 matched clear cell renal cell
carcinoma (ccRCC)/normal tissue pairs. coldata include patient metadata
(e.g. age, gender, sage, etc.)
[doi:10.1016/j.ccell.2015.12.004](https://doi.org/10.1016/j.ccell.2015.12.004)

## Usage

``` r
tissue_norm_se
```

## Format

An object of class `SummarizedExperiment` with 570 rows and 276 columns.

## Examples

``` r
data(tissue_norm_se)
head(tissue_norm_se)
#> class: SummarizedExperiment 
#> dim: 6 276 
#> metadata(0):
#> assays(1): data
#> rownames(6): 1,2-propanediol 1,3-dihydroxyacetone ... 10-nonadecenoate
#>   (19:1n9) 13-HODE + 9-HODE
#> rowData names(0):
#> colnames(276): DIAG-16076 DIAG-16077 ... DIAG-16354 DIAG-16355
#> colData names(5): TISSUE_TYPE GENDER RACE STAGE AGE
```
