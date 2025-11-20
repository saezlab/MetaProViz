# KEGG pathways

KEGG pathways

## Usage

``` r
metsigdb_kegg()
```

## Value

A data frame containing the KEGG pathways suitable for ORA.

## Examples

``` r
metsigdb_kegg()
#> # A tibble: 18,755 × 6
#>    Description MetaboliteID term               Metabolite pubchem compound_names
#>    <chr>       <chr>        <chr>              <chr>      <chr>   <list>        
#>  1 map00010    C00022       Glycolysis / Gluc… Pyruvate   3324    <chr [5]>     
#>  2 map00010    C00024       Glycolysis / Gluc… Acetyl-CoA 3326    <chr [2]>     
#>  3 map00010    C00031       Glycolysis / Gluc… D-Glucose  3333    <chr [5]>     
#>  4 map00010    C00033       Glycolysis / Gluc… Acetate    3335    <chr [3]>     
#>  5 map00010    C00036       Glycolysis / Gluc… Oxaloacet… 3338    <chr [6]>     
#>  6 map00010    C00068       Glycolysis / Gluc… Thiamin d… 3368    <chr [5]>     
#>  7 map00010    C00074       Glycolysis / Gluc… Phosphoen… 3374    <chr [3]>     
#>  8 map00010    C00084       Glycolysis / Gluc… Acetaldeh… 3384    <chr [2]>     
#>  9 map00010    C00085       Glycolysis / Gluc… D-Fructos… 3385    <chr [3]>     
#> 10 map00010    C00103       Glycolysis / Gluc… D-Glucose… 3403    <chr [4]>     
#> # ℹ 18,745 more rows
```
