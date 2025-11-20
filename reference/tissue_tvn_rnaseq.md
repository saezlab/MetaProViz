# tissue_tvn_rnaseq

The processed transcriptomics data was downloaded from the supplementary
table 3 of Mora & Schmidt et. al., which used the study from Clark et.
al. under Proteomics data Commons PDC000127. on their regulation
regulation in renal cancer, Genome Medicine 2024,
[doi:10.1186/s13073-024-01415-3](https://doi.org/10.1186/s13073-024-01415-3)
Clark et. al, Integrated proteogenomic characterization of clear cell
renal cell carcinoma, Cell 2019,
[doi:10.1016/j.cell.2019.10.007](https://doi.org/10.1016/j.cell.2019.10.007)

## Usage

``` r
tissue_tvn_rnaseq
```

## Format

An object of class `spec_tbl_df` (inherits from `tbl_df`, `tbl`,
`data.frame`) with 29283 rows and 10 columns.

## Examples

``` r
data(tissue_tvn_rnaseq)
head(tissue_tvn_rnaseq)
#> # A tibble: 6 × 10
#>   gene_name entrezgene_id Log2FC    p.val    p.adj t.val SiRCleCluster_RG2
#>   <chr>             <dbl>  <dbl>    <dbl>    <dbl> <dbl> <chr>            
#> 1 FGR                2268   2.02 4.23e-90 8.96e-89 20.1  MDE              
#> 2 PLXND1            23129   1.54 1.27e-99 3.39e-98 21.2  MDE              
#> 3 MPO                4353   1.26 1.55e- 8 2.93e- 8  5.66 MDE              
#> 4 IL32               9235   1.13 1.77e-23 6.10e-23  9.99 MDE              
#> 5 TRAF3IP3          80342   1.46 1.02e-49 7.58e-49 14.8  MDE              
#> 6 STAB1             23166   1.61 2.68e-54 2.26e-53 15.5  MDE              
#> # ℹ 3 more variables: `SiRCleThreshold_DNA-Methylation` <chr>,
#> #   SiRCleThreshold_RNA <chr>, SiRCleThreshold_Protein <chr>
```
