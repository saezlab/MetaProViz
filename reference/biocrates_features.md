# biocrates_features

Biocrates kit feature information of the "MxP Quant 500 XL kit" that
covers more than 1,000 metabolites with biochemical class information
and the exported different metabolite IDs (HMDB, KEGG, etc.).
information (INCHI, Key, etc.).

## Usage

``` r
biocrates_features
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 1019 rows and 14 columns.

## Examples

``` r
data(biocrates_features)
head(biocrates_features)
#> # A tibble: 6 × 14
#>   TrivialName TrivialName_Prior2023 Class    CAS   CHEBI CID   HMDB  INCHI Key  
#>   <chr>       <chr>                 <chr>    <chr> <chr> <chr> <chr> <chr> <chr>
#> 1 1-Met-His   1-Met-His             Aminoac… 332-… 50599 92105 HMDB… InCh… BRMW…
#> 2 3-IAA       3-IAA                 Indoles… 87-5… 16411 802   HMDB… InCh… SEOV…
#> 3 3-IPA       3-IPA                 Indoles… 830-… 43580 3744  HMDB… InCh… GOLX…
#> 4 3-Met-His   3-Met-His             Aminoac… 368-… 27596 64969 HMDB… InCh… JDHI…
#> 5 5-AVA       5-AVA                 Aminoac… 660-… 15887 138   HMDB… InCh… JJMD…
#> 6 AABA        AABA                  Aminoac… 1492… 2879… 6657… HMDB… InCh… QWCK…
#> # ℹ 5 more variables: IUPAC <chr>, LIMID <chr>, MESH <chr>, Molecule <chr>,
#> #   SMILES <chr>
```
