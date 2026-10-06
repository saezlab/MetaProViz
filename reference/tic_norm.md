# Total ion count normalization

Total ion count normalization

## Usage

``` r
tic_norm(data, metadata_sample, metadata_info, tic = TRUE)
```

## Arguments

- data:

  DF which contains unique sample identifiers as row names and
  metabolite numerical values in columns with metabolite identifiers as
  column names. Use NA for metabolites that were not detected and
  consider converting any zeros to NA unless they are true zeros.

- metadata_sample:

  DF which contains information about the samples, which will be
  combined with the input data based on the unique sample identifiers
  used as rownames.

- metadata_info:

  Named vector containing the information about the names of the
  experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile").

- tic:

  *Optional:* If TRUE, total Ion Count normalization is performed. If
  FALSE, only RLA QC plots are returned. **Default = TRUE**

## Value

List with two elements: DF (including output table) and Plot (including
all plots generated)

This function can be used independently on numeric data with matching
sample metadata. Missing values and zero handling should already be in a
reasonable state for the intended normalization.

## Examples

``` r
data(intracell_raw)
Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
Res <- tic_norm(
    data =
        Intra[-c(49:58), -c(1:3)] %>%
        dplyr::mutate_all(
            ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
        ),
    metadata_sample = Intra[-c(49:58), c(1:3)],
    metadata_info = c(Conditions = "Conditions")
)
#> total Ion Count (tic) normalization: total Ion Count (tic) normalization is used to reduce the variation from non-biological sources, while maintaining the biological variation. REF: Wulff et. al., (2018), Advances in Bioscience and Biotechnology, 9, 339-351, doi:https://doi.org/10.4236/abb.2018.98022
#> Warning: Removed 39 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
#> Warning: Removed 39 rows containing non-finite outside the scale range
#> (`stat_boxplot()`).
```
