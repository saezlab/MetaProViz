# intracell_dma

Metabolomics workbench project PR001418, study ST002224 where we
performed differential metabolite analysis comparing intracellular
metabolomics of 786-M1A versus HK2 cells. metabolite values used as
input with row names being metabolitetrivial names. nitrogen supports
renal cancer progression , Nature Communications 2022,
[doi:10.1038/s41467-022-35036-4](https://doi.org/10.1038/s41467-022-35036-4)

## Usage

``` r
intracell_dma
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 179 rows and 14 columns.

## Examples

``` r
data(intracell_dma)
head(intracell_dma)
#> # A tibble: 6 × 14
#>   Metabolite           Log2FC   p.adj   t.val HMDB  KEGG.ID KEGGCompound Pathway
#>   <chr>                 <dbl>   <dbl>   <dbl> <chr> <chr>   <chr>        <chr>  
#> 1 2-aminoadipic acid    0.153 2.38e-1 -7.56e5 HMDB… NA      NA           Not as…
#> 2 2-hydroxyglutarate    0.932 6.45e-5 -2.02e8 HMDB… C02630  2-Hydroxygl… Citrat…
#> 3 2-ketoglutarate       1.35  7.00e-6 -5.96e8 HMDB… C00026  2-Oxoglutar… Citrat…
#> 4 2/3-phosphoglycerate  0.699 9.46e-4 -1.93e7 HMDB… C00197  3-Phospho-D… Glycol…
#> 5 4-guanidinobutanoate -1.15  1.71e-4  2.37e6 HMDB… C01035  4-Guanidino… Argini…
#> 6 4-hydroxyphenyllact… -0.916 6.17e-8  1.91e6 HMDB… C03672  3-(4-Hydrox… Not as…
#> # ℹ 6 more variables: `786-M1A_1` <dbl>, `786-M1A_2` <dbl>, `786-M1A_3` <dbl>,
#> #   HK2_1 <dbl>, HK2_2 <dbl>, HK2_3 <dbl>
```
