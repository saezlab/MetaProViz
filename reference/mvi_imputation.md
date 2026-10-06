# Missing Value Imputation using half minimum value

Missing Value Imputation using half minimum value

## Usage

``` r
mvi_imputation(
  data,
  metadata_sample,
  metadata_info,
  core = FALSE,
  mvi_percentage = 50
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
  performed. Should not be normalised to media blank. Provide
  information about control media sample names via metadata_info
  "core_media" samples.**Default = FALSE**

- mvi_percentage:

  *Optional:* percentage 0-100 of imputed value based on the minimum
  value. **Default = 50**

## Value

DF with imputed values

This function can be used on raw or already filtered data, provided
sample metadata and condition information are available.

## Examples

``` r
data(intracell_raw)
Intra <- intracell_raw %>% tibble::column_to_rownames("Code")
Res <- mvi_imputation(
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
#> Missing Value Imputation: Missing value imputation is performed, as a complementary approach to address the missing value problem, where the missing values are imputing using the `half minimum value`. REF: Wei et. al., (2018), Reports, 8, 663, doi:https://doi.org/10.1038/s41598-017-19120-0
#> Some features had all NA values in all samples of certain conditions - hence no imputation was done and NA remains for: SAICAR
```
