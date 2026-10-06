# medium_raw_se

Metabolomics workbench project PR001418, study ST002226 where we
exported integrated raw peak values of intracellular metabolomics of HK2
and cccRCC cell lines 786-O, 786-M1A, 786-M2A, OS-RC-2, OS-LM1 and
RFX-631 converted into an se object. -`Conditions`: Character vector
indicating cell line identity -`Biological_Replicates`: Integer
replicate number for biological replicates -`GrowthFactor`: Character
vector indicating growth factor condition nitrogen supports renal cancer
progression, Nature Communications 2022,
[doi:10.1038/s41467-022-35036-4](https://doi.org/10.1038/s41467-022-35036-4)

## Usage

``` r
medium_raw_se
```

## Format

An object of class `SummarizedExperiment` with 73 rows and 44 columns.

## Examples

``` r
data(medium_raw_se)
head(medium_raw_se)
#> class: SummarizedExperiment 
#> dim: 6 44 
#> metadata(0):
#> assays(1): data
#> rownames(6): valine-d8 hipppuric acid-d5 ... 3-Dehydro-L-threonate
#>   acetylcarnitine
#> rowData names(4): HMDB KEGG.ID KEGGCompound Pathway
#> colnames(44): MS51-01 MS51-02 ... POOL5 POOL6
#> colData names(3): Conditions Biological_Replicates GrowthFactor
```
