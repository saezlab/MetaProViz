# outlier_detection

outlier_detection

## Usage

``` r
outlier_detection(
  data,
  metadata_sample,
  metadata_info,
  core = FALSE,
  hotellins_confidence = 0.99
)
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
  Biological_Replicates="ColumnName_Plot_SettingsFile", core_media =
  "Columnname_Input_SettingsFile"). Column "Conditions" with information
  about the sample conditions, Column "BiologicalReplicates" including
  numerical values and Column "Columnname_Input_SettingsFile" is used to
  specify the name of the media controls in the Conditions.

- core:

  *Optional:* If TRUE, a consumption-release experiment has been
  performed. If not normalised yet, provide information about control
  media sample names via metadata_info "core_media" samples. **Default =
  FALSE**

- hotellins_confidence:

  *Optional:* Confidence level for Hotellin's T2 test. **Default =
  0.99**

## Value

List with two elements: : DF (including output tables) and Plot
(including all plots generated)

This function can be run independently, but PCA-based outlier detection
is recommended on already preprocessed data such as filtered, imputed,
and where appropriate normalized input.

## Examples

``` r
data(intracell_raw)
Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
Res <- outlier_detection(
    data =
        Intra[-c(49:58), -c(1:3)] %>%
        dplyr::mutate_all(
            ~ ifelse(grepl("^0*(\\.0*)?$", as.character(.)), NA, .)
        ),
    metadata_sample = Intra[-c(49:58), c(1:3)],
    metadata_info = c(
        Conditions = "Conditions",
        Biological_Replicates = "Biological_Replicates"
    )
)
#> Outlier detection: Identification of outlier samples is performed using Hotellin's T2 test to define sample outliers in a mathematical way (Confidence = 0.99 ~ p.val < 0.01) (REF: Hotelling, H. (1931), Annals of Mathematical Statistics. 2 (3), 360-378, doi:https://doi.org/10.1214/aoms/1177732979). hotellins_confidence value selected: 0.99
#> No sample outliers were found
#> NA values are included in data that were set to 0 prior to performing PCA.
#> NA values are included in data that were set to 0 prior to performing PCA.
```
