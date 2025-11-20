# tissue_tvn_proteomics

The processed proteomics data was downloaded from the supplementary
table 3 of Mora & Schmidt et. al., which used the study from Clark et.
al. under Proteomics data Commons PDC000127. on their regulation
regulation in renal cancer, Genome Medicine 2024,
[doi:10.1186/s13073-024-01415-3](https://doi.org/10.1186/s13073-024-01415-3)
Clark et. al, Integrated proteogenomic characterization of clear cell
renal cell carcinoma, Cell 2019,
[doi:10.1016/j.cell.2019.10.007](https://doi.org/10.1016/j.cell.2019.10.007)

## Usage

``` r
tissue_tvn_proteomics
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 8769 rows and 10 columns.

## Examples

``` r
data(tissue_tvn_proteomics)
head(tissue_tvn_proteomics)
#> # A tibble: 6 × 10
#>   gene_name entrezgene_id Log2FC    p.val    p.adj t.val SiRCleCluster_RG2
#>   <chr>             <dbl>  <dbl>    <dbl>    <dbl> <dbl> <chr>            
#> 1 FGR                2268  0.594 6.45e-19 1.77e-18  9.98 MDE              
#> 2 PLXND1            23129  0.596 6.59e-38 4.17e-37 16.5  MDE              
#> 3 MPO                4353  0.504 1.22e- 9 2.23e- 9  6.42 MDE              
#> 4 IL32               9235  0.878 1.25e-12 2.63e-12  7.64 MDE              
#> 5 TRAF3IP3          80342  0.518 3.30e-15 7.78e-15  8.63 MDE              
#> 6 STAB1             23166  0.547 5.59e-16 1.36e-15  8.92 MDE              
#> # ℹ 3 more variables: `SiRCleThreshold_DNA-Methylation` <chr>,
#> #   SiRCleThreshold_RNA <chr>, SiRCleThreshold_Protein <chr>
```
