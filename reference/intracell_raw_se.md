# intracell_raw_se

Metabolomics workbench project PR001418, study ST002224 where we
exported integrated raw peak values of intracellular metabolomics of HK2
and ccRCC cell lines 786-O, 786-M1A and 786-M2A converted into an se
object. -`Conditions`: Character vector indicating cell line identity -
`Analytical_Replicate`: Integer replicate number for analytical
replicates -`Biological_Replicate`: Integer replicate number for
biological replicates nitrogen supports renal cancer progression, Nature
Communications 2022, DOI:10.1038/s41467-022-35036-4.

## Usage

``` r
intracell_raw_se
```

## Format

An object of class `SummarizedExperiment` with 182 rows and 58 columns.

## Examples

``` r
data(intracell_raw_se)
head(intracell_raw_se)
#> class: SummarizedExperiment 
#> dim: 6 58 
#> metadata(0):
#> assays(1): counts
#> rownames(6): valine-d8 hippuric acid-d5 ... 2-hydroxyglutarate
#>   2-ketoglutarate
#> rowData names(0):
#> colnames(58): MS55_01 MS55_02 ... POOL8 POOL9
#> colData names(3): Conditions Analytical_Replicates
#>   Biological_Replicates
```
