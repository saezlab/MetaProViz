# Normalize consumption release data

Normalize consumption release data

## Usage

``` r
core_norm(data, metadata_sample, metadata_info)
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
  experimental parameters. c(Conditions="ColumnName_Plot_SettingsFile",
  core_norm_factor = "Columnname_Input_SettingsFile", core_media =
  "Columnname_Input_SettingsFile"). Column core_norm_factor is used for
  normalization and core_media is used to specify the name of the media
  controls in the Conditions.

## Value

List with two elements: DF (including output table) and Plot (including
all plots generated)

This function is intended for standalone CoRe-style experiments with
control-media metadata available via `metadata_info`.

## Examples

``` r
data(medium_raw)
Media <-
    medium_raw %>%
    tibble::column_to_rownames("Code") %>%
    subset(!Conditions == "Pool") %>%
    dplyr::mutate_all(~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .))
Res <- core_norm(
    data = Media[, -c(1:3)],
    metadata_sample = Media[, c(1:3)],
    metadata_info = c(
        Conditions = "Conditions",
        core_norm_factor = "GrowthFactor",
        core_media = "blank"
    )
)
#> For Consumption Release experiment we are using the method from Jain M.  REF: Jain et. al, (2012), Science 336(6084):1040-4, doi: 10.1126/science.1218595.
#> 9 of variables have high variability (CV > 30) in the core_media control samples. Consider checking the pooled samples to decide whether to remove these metabolites or not.
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
#> Bin width defaults to 1/30 of the range of the data. Pick better value with
#> `binwidth`.
#> Warning: The core_media samples  MS51-04  were found to be different from the rest. They will not be included in the sum of the core_media samples.
#> core data are normalised by substracting mean (blank) from each sample and multiplying with the core_norm_factor
```
