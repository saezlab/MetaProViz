# Sample Metadata Analysis

## ![](Hexagon_MetaProViz.png)

  
Tissue metabolomics experiment is a standard metabolomics experiment
using tissue samples (e.g. from animals or patients).  
  
In this tutorial we showcase how to use **MetaProViz**:  

- to perform differential metabolite analysis (dma) to generate Log2FC
  and statistics and perform pathway analysis using Over Representation
  Analysis (ORA) on the results.  
- to do metabolite clustering analysis (MCA) to find clusters of
  metabolites with similar behaviors based on patients demographics like
  age, gender and tumour stage.  
- Find the main metabolite drivers that separate patients based on their
  demographics like age, gender and tumour stage.  
    
  First if you have not done yet, install the required dependencies and
  load the libraries:

``` r
# 1. Install Rtools if you haven’t done this yet, using the appropriate version (e.g.windows or macOS).
# 2. Install the latest development version from GitHub using devtools
# devtools::install_github("https://github.com/saezlab/MetaProViz")

library(MetaProViz)

# dependencies that need to be loaded:
library(magrittr)
library(dplyr)
library(rlang)
library(tidyr)
library(tibble)
library(stringr)

# Please install the Biocmanager Dependencies:
# BiocManager::install("EnhancedVolcano")
```

  
  

## 1. Loading the example data

  
Here we choose an example datasets, which is publicly available in the
[paper](https://www.cell.com/cancer-cell/comments/S1535-6108(15)00468-7#supplementaryMaterial)
“An Integrated Metabolic Atlas of Clear Cell Renal Cell Carcinoma”,
which includes metabolomic profiling on 138 matched clear cell renal
cell carcinoma (ccRCC)/normal tissue pairs (Hakimi et al. 2016).
Metabolomics was done using The company Metabolon, so this is untargeted
metabolomics. Here we use the median normalised data from the
supplementary table 2 of the paper. We have combined the metainformation
about the patients with the metabolite measurements and removed
unidentified metabolites. Lastly, we have added a column “Stage” where
Stage1 and Stage2 patients are summarised to “EARLY-STAGE” and Stage3
and Stage4 patients to “LATE-STAGE”. Moreover, we have added a column
“Age”, where patients with “AGE AT SURGERY” \<42 are defined as “Young”
and patients with AGE AT SURGERY \>58 as “Old” and the remaining
patients as “Middle”.  

\#As part of the **MetaProViz** package you can load the example data
into your global environment using the function `toy_data()`:  
`1.` Tissue experiment **(Intra)**  
We can access the built-in dataset `tissue_norm`, which includes columns
with Sample information and columns with the median normalised measured
metabolite integrated peaks.  

``` r
# Load the example data:
data(tissue_norm)

Tissue_Norm <- tissue_norm%>%
column_to_rownames("Code")
```

|            | TISSUE_TYPE | GENDER | AGE_AT_SURGERY | TYPE-STAGE     | STAGE       | AGE | 1,2-propanediol | 1,3-dihydroxyacetone |
|:-----------|:------------|:-------|---------------:|:---------------|:------------|:----|----------------:|---------------------:|
| DIAG-16076 | TUMOR       | Male   |        74.7778 | TUMOR-STAGE I  | EARLY-STAGE | Old |        0.710920 |             0.809182 |
| DIAG-16077 | NORMAL      | Male   |        74.7778 | NORMAL-STAGE I | EARLY-STAGE | Old |        0.339390 |             0.718725 |
| DIAG-16078 | TUMOR       | Male   |        77.1778 | TUMOR-STAGE I  | EARLY-STAGE | Old |        0.413386 |             0.276412 |
| DIAG-16079 | NORMAL      | Male   |        77.1778 | NORMAL-STAGE I | EARLY-STAGE | Old |        1.595697 |             4.332451 |
| DIAG-16080 | TUMOR       | Female |        59.0889 | TUMOR-STAGE I  | EARLY-STAGE | Old |        0.573787 |             0.646791 |

Preview of the DF `Tissue_Norm` including columns with sample
information and metabolite ids with their measured values.

`2.` Additional information mapping the trivial metabolite names to KEGG
IDs, HMDB IDs, etc. and selected pathways **(MappingInfo)**  

``` r
data(tissue_meta)

Tissue_MetaData <- tissue_meta%>%
column_to_rownames("Metabolite")
```

|                              | SUPER_PATHWAY | SUB_PATHWAY                                      | CAS         | PUBCHEM | KEGG   | Group_HMDB |
|:-----------------------------|:--------------|:-------------------------------------------------|:------------|--------:|:-------|:-----------|
| 1,2-propanediol              | Lipid         | Ketone bodies                                    | 57-55-6;    |      NA | C00583 | HMDB01881  |
| 1,3-dihydroxyacetone         | Carbohydrate  | Glycolysis, gluconeogenesis, pyruvate metabolism | 62147-49-3; |     670 | C00184 | HMDB01882  |
| 1,5-anhydroglucitol (1,5-AG) | Carbohydrate  | Glycolysis, gluconeogenesis, pyruvate metabolism | 154-58-5;   |      NA | C07326 | HMDB02712  |
| 10-heptadecenoate (17:1n7)   | Lipid         | Long chain fatty acid                            | 29743-97-3; | 5312435 | NA     | NA         |
| 10-nonadecenoate (19:1n9)    | Lipid         | Long chain fatty acid                            | 73033-09-7; | 5312513 | NA     | NA         |

Preview of the DF `Tissue_MetaData` including the trivial metabolite
identifiers used in the experiment as well as IDs and pathway
information.

  

## 2. Run MetaProViz Analysis

### Pre-processing

This has been done by the authors of the paper and we will use the
median normalized data. If you want to know how you can use the
**MetaProViz** pre-processing module, please check out the vignette:  
- [Standard metabolomics
data](https://saezlab.github.io/MetaProViz/articles/Standard%20Metabolomics.html)  
- [Consumption-Release (core) metabolomics data from cell culture
media](https://saezlab.github.io/MetaProViz/articles/core%20Metabolomics.html)  

### Metadata analysis

We can use the patient’s metadata to find the main metabolite drivers
that separate patients based on their demographics like age, gender,
etc.  
  
Here the metadata analysis is based on principal component analysis
(PCA), which is a dimensionality reduction method that reduces all the
measured features (=metabolites) of one sample into a few features in
the different principal components, whereby each principal component can
explain a certain percentage of the variance between the different
samples. Hence, this enables interpretation of sample clustering based
on the measured features (=metabolites).  
The [`metadata_analysis()`](../../reference/metadata_analysis.md)
function will perform PCA to extract the different PCs followed by
annova to find the main metabolite drivers that separate patients based
on their demographics.  
  

``` r
MetaRes <- metadata_analysis(data=Tissue_Norm[,-c(1:13)],
metadata_sample= Tissue_Norm[,c(2,4:5,12:13)],
scaling = TRUE,
percentage = 0.1,
cutoff_stat= 0.05,
cutoff_variance = 1)
#> The column names of the 'metadata_sample' contain special character that where removed.
```

![](sample-metadata_files/figure-html/code-3-1.png)  
Ultimately, this is leading to clusters of metabolites that are driving
the separation of the different demographics.  
  
We generated the general anova output DF:  

|      | PC    | tukeyHSD_Contrast      | term        |  anova_sumsq | anova_meansq | anova_statistic | anova_p.value | tukeyHSD_p.adjusted | Explained_Variance |
|:-----|:------|:-----------------------|:------------|-------------:|-------------:|----------------:|--------------:|--------------------:|-------------------:|
| 1    | PC1   | TUMOR-NORMAL           | TISSUE_TYPE | 7375.2212254 | 7375.2212254 |      90.1352800 |     0.0000000 |           0.0000000 |         19.0079573 |
| 2    | PC1   | Black-Asian            | RACE        |  156.7549451 |   52.2516484 |       0.4795311 |     0.6967838 |           0.9992824 |         19.0079573 |
| 1777 | PC232 | White-Asian            | RACE        |    0.2726952 |    0.0908984 |       0.9544580 |     0.4147472 |           0.5495634 |          0.0166997 |
| 1778 | PC232 | LATE-STAGE-EARLY-STAGE | STAGE       |    0.0360319 |    0.0360319 |       0.3776760 |     0.5393596 |           0.5393596 |          0.0166997 |
| 3191 | PC9   | White-Other            | RACE        |   20.2626222 |    6.7542074 |       0.7090398 |     0.5473277 |           0.8897908 |          1.6658973 |
| 3192 | PC9   | Young-Middle           | AGE         |   49.2030973 |   24.6015486 |       2.6213834 |     0.0745317 |           0.9979475 |          1.6658973 |

Preview of the DF MetaRes\[\[`res_aov`\]\] including the main metabolite
drivers that separate patients based on their demographics.

  
We generated the summarised results output DF, where each feature
(=metabolite) was assigned a main demographics parameter this feature is
separating:  

| feature                                     | term                                  | Sum(Explained_Variance)                                                                     | MainDriver                       | MainDriver_Term | MainDriver_Sum(VarianceExplained) |
|:--------------------------------------------|:--------------------------------------|:--------------------------------------------------------------------------------------------|:---------------------------------|:----------------|----------------------------------:|
| N2-methylguanosine                          | AGE, GENDER, RACE, STAGE, TISSUE_TYPE | 3.7598809871366, 2.75344828467363, 1.43932034747081, 25.2803273612931, 33.9484968841434     | FALSE, FALSE, FALSE, FALSE, TRUE | TISSUE_TYPE     |                           33.9485 |
| 5-methyltetrahydrofolate (5MeTHF)           | AGE, GENDER, RACE, STAGE, TISSUE_TYPE | 0.351143408196559, 0.294688079880209, 0.252515004077172, 19.4058160726611, 32.5442969233501 | FALSE, FALSE, FALSE, FALSE, TRUE | TISSUE_TYPE     |                           32.5443 |
| N-acetylalanine                             | AGE, GENDER, RACE, STAGE, TISSUE_TYPE | 0.381811888038144, 1.82507161869685, 2.97460134435356, 19.0079573465895, 32.5442969233501   | FALSE, FALSE, FALSE, FALSE, TRUE | TISSUE_TYPE     |                           32.5443 |
| N-acetyl-aspartyl-glutamate (NAAG)          | AGE, GENDER, RACE, STAGE, TISSUE_TYPE | 0.212976478603115, 2.75344828467363, 0.235952044704099, 20.3572056704699, 32.2825995362096  | FALSE, FALSE, FALSE, FALSE, TRUE | TISSUE_TYPE     |                           32.2826 |
| 1-heptadecanoylglycerophosphoethanolamine\* | AGE, GENDER, RACE, STAGE, TISSUE_TYPE | 4.33031140898824, 0.267851018697989, 0.437853800763463, 22.8424358182214, 30.8783995754163  | FALSE, FALSE, FALSE, FALSE, TRUE | TISSUE_TYPE     |                           30.8784 |
| 1-linoleoylglycerophosphoethanolamine\*     | AGE, STAGE, TISSUE_TYPE               | 4.20115004304429, 22.7555733683955, 30.8783995754163                                        | FALSE, FALSE, TRUE               | TISSUE_TYPE     |                           30.8784 |

Preview of the DF MetaRes\[\[`res_summary`\]\] including the metabolite
drivers in rows and list the patients demographics they can separate.

  

``` r
##1. Tissue_Type
TissueTypeList <- MetaRes[["res_summary"]]%>%
filter(MainDriver_Term == "TISSUE_TYPE")%>%
filter(`MainDriver_Sum(VarianceExplained)`>30)%>%
select(feature)%>%
pull()

# select columns tissue_norm that are in TissueTypeList if they exist
Input_Heatmap <- Tissue_Norm[ , names(Tissue_Norm) %in% TissueTypeList]#c("N1-methylguanosine", "N-acetylalanine", "lysylmethionine")

# Heatmap: Metabolites that separate the demographics, like here TISSUE_TYPE
viz_heatmap(data = Input_Heatmap,
metadata_sample = Tissue_Norm[,c(1:13)],
metadata_info = c(color_Sample = list("TISSUE_TYPE")),
scale ="column",
plot_name = "MainDrivers")
```

![](sample-metadata_files/figure-html/viz-heatmap-1.png)

### DMA

Here we use Differential Metabolite Analysis (`dma`) to compare two
conditions (e.g. Tumour versus Healthy) by calculating the Log2FC,
p-value, adjusted p-value and t-value.  
For more information please see the vignette:  
- [Standard metabolomics
data](https://saezlab.github.io/MetaProViz/articles/Standard%20Metabolomics.html)  
- [Consumption-Release (core) metabolomics data from cell culture
media](https://saezlab.github.io/MetaProViz/articles/core%20Metabolomics.html)  
  
We will perform multiple comparisons based on the different patient
demographics available: 1. Tumour versus Normal: All patients 2. Tumour
versus Normal: Subset of `Early Stage` patients 3. Tumour versus Normal:
Subset of `Late Stage` patients 4. Tumour versus Normal: Subset of
`Young` patients 5. Tumour versus Normal: Subset of `Old` patients  

``` r
# Prepare the different selections
EarlyStage <- Tissue_Norm %>%
filter(STAGE== "EARLY-STAGE")
LateStage <- Tissue_Norm %>%
filter(STAGE=="LATE-STAGE")
Old <- Tissue_Norm %>%
filter(AGE=="Old")
Young <- Tissue_Norm %>%
filter(AGE=="Young")

DFs <- list(
    "TissueType" = Tissue_Norm,
    "EarlyStage" = EarlyStage,
    "LateStage" = LateStage,
    "Old" = Old,
    "Young" = Young
)

# Run dma
ResList <- list()
for(item in names(DFs)){
    #Get the right DF:
    data <- DFs[[item]]

    message(paste("Running dma for", item))
    #Create folder for saving each comparison
    dir.create(paste(getwd(),"/MetaProViz_Results/dma/", sep=""), showWarnings = FALSE)
    dir.create(paste(getwd(),"/MetaProViz_Results/dma/", item, sep=""), showWarnings = FALSE)

    #Perform dma
    TvN <- dma(data =  data[,-c(1:13)],
    metadata_sample =  data[,c(1:13)],
    metadata_info = c(Conditions="TISSUE_TYPE", Numerator="TUMOR" , Denominator = "NORMAL"),
    shapiro=FALSE, #The data have been normalized by the company that provided the results and include metabolites with zero variance as they were all imputed with the same missing value.
    path = paste(getwd(),"/MetaProViz_Results/dma/", item, sep=""))

    #Add Results to list
    ResList[[item]] <- TvN
}
#> Running dma for TissueType
#> There are no NA/0 values
```

![](sample-metadata_files/figure-html/Run_not_Display-1.png)

    #> Running dma for EarlyStage
    #> There are no NA/0 values

![](sample-metadata_files/figure-html/Run_not_Display-2.png)

    #> Running dma for LateStage
    #> There are no NA/0 values

![](sample-metadata_files/figure-html/Run_not_Display-3.png)

    #> Running dma for Old
    #> There are no NA/0 values

![](sample-metadata_files/figure-html/Run_not_Display-4.png)

    #> Running dma for Young
    #> There are no NA/0 values

![](sample-metadata_files/figure-html/Run_not_Display-5.png)

  
  
  

We can see from the different Volcano plots have smaller p.adjusted
values and differences in Log2FC range.  
Here we can also use the `MetaproViz::viz_volcano()` function to plot
comparisons together on the same plot, such as Tumour versus Normal of
young and old patients:  

``` r
# Early versus Late Stage
viz_volcano(plot_types="Compare",
data=ResList[["EarlyStage"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
data2= ResList[["LateStage"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
name_comparison= c(data="EarlyStage", data2= "LateStage"),
plot_name= "EarlyStage-TUMOR_vs_NORMAL compared to LateStage-TUMOR_vs_NORMAL",
subtitle= "Results of dma" )
```

![](sample-metadata_files/figure-html/viz-volcano-1.png)

``` r

# Young versus Old
viz_volcano(plot_types="Compare",
data=ResList[["Young"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
data2= ResList[["Old"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
name_comparison= c(data="Young", data2= "Old"),
plot_name= "Young-TUMOR_vs_NORMAL compared to Old-TUMOR_vs_NORMAL",
subtitle= "Results of dma" )
```

![](sample-metadata_files/figure-html/viz-volcano-2.png)  
Here we can observe that Tumour versus Normal has lower significance
values for the Young patients compared to the Old patients. This can be
due to higher variance in the metabolite measurements from Young
patients compared to the Old patients.  
We can also check if the top changed metabolites comparing Tumour versus
Normal correlate with the main metabolite drivers that separate patients
based on their `TISSUE_TYPE`, which are Tumour or Normal.  

``` r
# Get the top changed metabolites
top_entries <- ResList[["TissueType"]][["dma"]][["TUMOR_vs_NORMAL"]] %>%
arrange(desc(t.val)) %>%
slice(1:25)%>%
select(Metabolite)%>%
pull()
bottom_entries <- ResList[["TissueType"]][["dma"]][["TUMOR_vs_NORMAL"]] %>%
arrange(desc(t.val)) %>%
slice((n()-24):n())%>%
select(Metabolite)  %>%
pull()

# Check if those overlap with the top demographics drivers
ggVennDiagram::ggVennDiagram(list(top = top_entries,
Bottom = bottom_entries,
TissueTypeList = TissueTypeList))+
ggplot2::scale_fill_gradient(low = "blue", high = "red")
```

![](sample-metadata_files/figure-html/viz-volcano-2-1.png)

``` r


MetaData_Metab <- merge(x=tissue_meta,
y= MetaRes[["res_summary"]][, c(1,5:6) ]%>%tibble::column_to_rownames("feature"),
by=0,
all.y=TRUE)%>%
column_to_rownames("Row.names")

# Make a Volcano plot:
viz_volcano(plot_types="Standard",
data=ResList[["TissueType"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
metadata_feature =  MetaData_Metab,
metadata_info = c(color = "MainDriver_Term"),
plot_name= "TISSUE_TYPE-TUMOR_vs_NORMAL",
subtitle= "Results of dma" )
```

![](sample-metadata_files/figure-html/viz-volcano-2-2.png)

### Metabolite ID QC

Before we perform enrichment analysis,it is important to check the
availability of metabolite IDs that we aim to use. We can first create
an overview plot:

    #> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    #> ℹ Please use tidy evaluation idioms with `aes()`.
    #> ℹ See also `vignette("ggplot2-in-packages")` for more information.
    #> ℹ The deprecated feature was likely used in the MetaProViz package.
    #>   Please report the issue at <https://github.com/saezlab/MetaProViz/issues>.
    #> This warning is displayed once every 8 hours.
    #> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    #> generated.
    #> Warning: themes$intersections_matrix is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: themes$intersections_matrix is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: selected_theme is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: selected_theme is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: themes$overall_sizes is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: themes$overall_sizes is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Note: Could not render upset plot to screen due to theme compatibility issues. Plot was saved to file successfully.

![](sample-metadata.Rmd_compare-pk.svg)  
Here we notice that 76 features have no metabolite ID assigned, yet have
a trivial name and a metabolite class assigned. For 135 metabolites we
only have a pubchem ID, yet no HMDB or KEGG ID. Next, we will try to
understand the missing IDs focusing on HMDb as an example:

``` r
Plot1_HMDB <- count_id(MetaboliteIDs, "HMDB")
```

![](sample-metadata_files/figure-html/code-6-1.png)  
Over 200 metabolites have no HMDB ID and if there is a HMDB ID assigned,
we only have one HMDB ID. Here we assume the experimental setup does not
account for stereoisomers and that other IDs are with different degrees
of ambiguity are possible (e.g. L-Alanine HMDB ID is present and we will
add D-Alanine HMDB ID, see PK vignette for more details). This enables
us to assign additional potential HMDB IDs per feature.

``` r
Input_HMDB <- MetaboliteIDs %>%
dplyr::filter(!is.na(HMDB)) %>% # ID in the measured data we want to use, hence we remove NA's
dplyr::select("Metabolite", "HMDB", "SUPER_PATHWAY") # only keep relevant columns

# Add equivalent IDs:
tissue_meta_AddIDs_HMDB <- equivalent_id(data= Input_HMDB,
metadata_info = c(InputID="HMDB"),# ID in the measured data, here we use the HMDB ID
from = "hmdb")
#> Warning in equivalent_id(data = Input_HMDB, metadata_info = c(InputID =
#> "HMDB"), : The following IDs are duplicated and removed: HMDB01859
#> chebi is used to find additional potential IDs for hmdb.
# Let's see how this has changed the number of entries:
tissue_meta_AddIDs_HMDB <- merge(x=MetaboliteIDs,
y=tissue_meta_AddIDs_HMDB,
by="HMDB",
all.x=TRUE)

# Plot 2 after equivalent IDs where added
Plot2_HMDB <- count_id(tissue_meta_AddIDs_HMDB, "AllIDs")
```

![](sample-metadata_files/figure-html/data-prep-1.png)  
Additionally, for the metabolites without HMDB ID we can check if a HMDB
ID is available using the other ID types. First we check if we have
cases in which we do have no HMDB ID, but other available IDs:

| Total_Cases | Not_NA_PUBCHEM | Not_NA_KEGG | Both_Not_NA |
|------------:|---------------:|------------:|------------:|
|         160 |            159 |          25 |          24 |

Overview of other ID types for metabolites without HMDB ID.

``` r
no_hmdb <- dplyr::filter(MetaboliteIDs, is.na(HMDB) & (!is.na(PUBCHEM) | !is.na(KEGG)))

tissue_meta_translated_HMDB <-
translate_id(no_hmdb, metadata_info = list(InputID = "PUBCHEM", grouping_variable = "SUPER_PATHWAY"), from = "pubchem", to = "hmdb") |>
extract2("TranslatedDF") |>
rename(hmdb_from_pubchem = hmdb) |>
translate_id(metadata_info = list(InputID = "KEGG", grouping_variable = "SUPER_PATHWAY"), from = "kegg", to = "hmdb") |>
extract2("TranslatedDF") |>
rename(hmdb_from_kegg = hmdb)
#> Warning in translate_id(no_hmdb, metadata_info = list(InputID = "PUBCHEM", :
#> The following IDs are duplicated within one group: 5.24601e+13, 6.99312e+13

# Here we combine the tables above, created by equivalent_id and translate_id:
Tissue_MetaData_HMDB <-
left_join(
    tissue_meta_AddIDs_HMDB |>
    select(Metabolite = Metabolite.x, SUPER_PATHWAY = SUPER_PATHWAY.x, SUB_PATHWAY ,COMP_ID, PLATFORM, RI, MASS, CAS, PUBCHEM, KEGG,  HMDB, hmdb_from_equivalentid = AllIDs),
    tissue_meta_translated_HMDB |>
    select(COMP_ID, hmdb_from_pubchem, hmdb_from_kegg) |> mutate(across(starts_with("hmdb_from"), ~na_if(., ""))),
    by = 'COMP_ID'
    ) |>
    rowwise() |>
    mutate(hmdb_combined = list(unique(na.omit(unlist(stringr::str_split(across(starts_with("hmdb_from")), ",")))))) |>
    mutate(hmdb_combined = paste0(hmdb_combined, collapse = ',')) |>  # we concatenate by "," only for the sake of printing in notebook
    ungroup()

# Plot 3:
Plot3_HMDB <- count_id(Tissue_MetaData_HMDB, "hmdb_combined")
```

![](sample-metadata_files/figure-html/data-prep-3-1.png)  
If we inspect the tables we can check that for example,
16-hydroxypalmitate certainly has an HMDB ID, yet was missing from the
original metadata table. We can check the structures associated to it on
[PubChem:10466](https://pubchem.ncbi.nlm.nih.gov/compound/10466) and
[KEGG:C18218](https://www.genome.jp/entry/C18218) are indeed identical,
and as we expect, the same structure is in HMDB too, under the ID
[HMDB0006294](https://hmdb.ca/metabolites/HMDB0006294).  
  
We can of course do the same for other IDs, such as the KEGG IDs:

``` r
Plot1_KEGG <- count_id(MetaboliteIDs, "KEGG")
```

![](sample-metadata_files/figure-html/data-prep-4-1.png)

``` r

##################################################################################################################
Input_KEGG <- Tissue_MetaData_HMDB %>%
dplyr::filter(!is.na(KEGG)) %>% # ID in the measured data we want to use, hence we remove NA's
dplyr::select("Metabolite", "KEGG", "SUPER_PATHWAY") # only keep relevant columns

# Add equivalent IDs:
tissue_meta_AddIDs_KEGG <- equivalent_id(data= Input_KEGG,
                                                metadata_info = c(InputID="KEGG"),# ID in the measured data, here we use the HMDB ID
                                                from = "kegg")
#> Warning in equivalent_id(data = Input_KEGG, metadata_info = c(InputID =
#> "KEGG"), : The following IDs are duplicated and removed: C01585, C00474,
#> C06804, C03594, C06429
#> chebi is used to find additional potential IDs for kegg.

# Let's see how this has changed the number of entries:
tissue_meta_AddIDs_KEGG <- merge(x=Tissue_MetaData_HMDB,
y=tissue_meta_AddIDs_KEGG,
by="KEGG",
all.x=TRUE)

Plot2_KEGG <-  count_id(tissue_meta_AddIDs_KEGG, "AllIDs")
```

![](sample-metadata_files/figure-html/data-prep-4-2.png)

``` r

###################################################################################################################
no_KEGG <- dplyr::filter(Tissue_MetaData_HMDB, is.na(KEGG) & (!is.na(PUBCHEM) | !is.na(HMDB)))

tissue_meta_translated_KEGG <-
translate_id(no_KEGG, metadata_info = list(InputID = "PUBCHEM", grouping_variable = "SUPER_PATHWAY"), from = "pubchem", to = "kegg") |>
extract2("TranslatedDF") |>
rename(kegg_from_pubchem = kegg) |>
translate_id(metadata_info = list(InputID = "HMDB", grouping_variable = "SUPER_PATHWAY"), from = "hmdb", to = "kegg") |>
extract2("TranslatedDF") |>
rename(kegg_from_hmdb = kegg)
#> Warning in translate_id(no_KEGG, metadata_info = list(InputID = "PUBCHEM", :
#> The following IDs are duplicated within one group: 5.24601e+13, 6.99312e+13

# here we combine the tables above, created by equivalent_id and translate_id:
Tissue_MetaData_KEGG <-
left_join(
    tissue_meta_AddIDs_KEGG |>
    select(Metabolite = Metabolite.x, SUPER_PATHWAY = SUPER_PATHWAY.x,SUB_PATHWAY, COMP_ID, PLATFORM, RI, MASS, CAS, PUBCHEM, HMDB=hmdb_combined , HMDB_Original = HMDB, hmdb_from_equivalentid, hmdb_from_kegg, hmdb_from_pubchem, , KEGG_Original = KEGG, kegg_from_equivalentid = AllIDs),
    tissue_meta_translated_KEGG |>
    select(COMP_ID, kegg_from_pubchem, kegg_from_hmdb) |> mutate(across(starts_with("kegg_from"), ~na_if(., ""))),
    by = 'COMP_ID'
    ) |>
    rowwise() |>
    mutate(KEGG= list(unique(na.omit(unlist(stringr::str_split(across(starts_with("kegg_from")), ",")))))) |>
    mutate(KEGG = paste0(KEGG, collapse = ',')) |>  # we concatenate by "," only for the sake of printing in notebook
    ungroup()

# Lets see the count now:
Plot3_KEGG <- count_id(Tissue_MetaData_KEGG, "KEGG")
```

![](sample-metadata_files/figure-html/data-prep-4-3.png)  
Lastly, we will have the new metadata table, which we can use for
enrichment analysis:

    #> Warning: themes$intersections_matrix is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: themes$intersections_matrix is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: selected_theme is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: selected_theme is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: themes$overall_sizes is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Warning: `legend.margin` must be specified using `margin()`
    #> ℹ For the old behavior use `legend.spacing`
    #> Warning: themes$overall_sizes is not a valid theme.
    #> Please use `theme()` to construct themes.
    #> Note: Could not render upset plot to screen due to theme compatibility issues. Plot was saved to file successfully.

![](sample-metadata.Rmd_compare-pk-2.svg)  
If we compare the results to the upset plot from the original metadata
we can see that initially 251 metabolites had a pubchem, HMDB and KEGG
ID, whilst after adding the equivalent IDs and translating the IDs, we
now have 285 metabolites with all three IDs. 135 metabolite had only a
Pubchem ID, whilst now we have found HMDB and/or KEGG IDs for half of
these metabolites with only 65 metabolites remaining with only a Pubchem
ID.

### ORA

We can perform Over Representation Analysis (`ORA`) using KEGG pathways
for each comparison and plot significant pathways. Noteworthy, since not
all metabolites have KEGG IDs, we will lose information.  

  
Given that in some cases we have multiple KEGG IDs for a measured
feature, we will check if this causes mapping to multiple, different
entries in the KEGG pathway-metabolite sets:

``` r
#Load Kegg pathways:
KEGG_Pathways <- metsigdb_kegg()

#check mapping with metadata

ccRCC_to_KEGGPathways <- checkmatch_pk_to_data(data = Tissue_MetaData_Extended,
input_pk = KEGG_Pathways,
metadata_info = c(InputID = "KEGG", PriorID = "MetaboliteID", grouping_variable = "term"))
#> Warning in checkmatch_pk_to_data(data = Tissue_MetaData_Extended, input_pk =
#> KEGG_Pathways, : 251 NA values were removed from column KEGG
#> Warning in checkmatch_pk_to_data(data = Tissue_MetaData_Extended, input_pk =
#> KEGG_Pathways, : 8 duplicated IDs were removed from column KEGG
#> data has multiple IDs per measurement = TRUE. input_pk has multiple IDs per entry = FALSE.
#> data has 318 unique entries with 328 unique KEGG IDs. Of those IDs, 249 match, which is 75.9146341463415%.
#> input_pk has 6541 unique entries with 6541 unique MetaboliteID IDs. Of those IDs, 249 are detected in the data, which is 3.80675737654793%.
#> Warning in checkmatch_pk_to_data(data = Tissue_MetaData_Extended, input_pk =
#> KEGG_Pathways, : There are cases where multiple detected IDs match to multiple
#> prior knowledge IDs of the same category

problems_terms <- ccRCC_to_KEGGPathways[["GroupingVariable_summary"]]%>%
filter(!Group_Conflict_Notes== "None")
```

| KEGG           | original_count | matches_count | MetaboliteID | term                                  |
|:---------------|---------------:|--------------:|:-------------|:--------------------------------------|
| C00221, C00031 |              2 |             2 | C00221       | Pentose phosphate pathway             |
| C00221, C00031 |              2 |             2 | C00221       | Biosynthesis of secondary metabolites |
| C00221, C00031 |              2 |             2 | C00031       | Glycolysis / Gluconeogenesis          |
| C00221, C00031 |              2 |             2 | C00221       | Glycolysis / Gluconeogenesis          |
| C00221, C00031 |              2 |             2 | C00221       | Metabolic pathways                    |
| C00221, C00031 |              2 |             2 | C00031       | Biosynthesis of secondary metabolites |
| C00221, C00031 |              2 |             2 | C00031       | Pentose phosphate pathway             |
| C00221, C00031 |              2 |             2 | C00031       | Metabolic pathways                    |
| C03460, C03722 |              2 |             2 | C03460       | Metabolic pathways                    |
| C03460, C03722 |              2 |             2 | C03722       | Metabolic pathways                    |
| C17737, C00695 |              2 |             2 | C17737       | Secondary bile acid biosynthesis      |
| C17737, C00695 |              2 |             2 | C00695       | Secondary bile acid biosynthesis      |

Terms in KEGG pathways where the same measured feature mapps with more
than one ID, which will inflate the enrichment analysis.

  
Dependent on the biological question and the organism and prior
knowledge, one can either maintain the metabolite ID of the more likely
metabolite (e.g. in human its moire likely that we have L-aminoacid than
D-aminoacid) or the metabolite ID that is represented in more/less
pathways (specificity). If we are looking into the cases where we do
have multiple IDs, for cases with no match to the prior knowledge we can
just maintain one ID, whilst for cases with exactly one match we should
maintain the ID that is found in the prior knowledge.

``` r
# Select the cases where a feature has multiple IDs
multipleIDs <- ccRCC_to_KEGGPathways[["data_summary"]]%>%
filter(original_count>1)
```

| KEGG                   | InputID_select | original_count | matches_count | matches        | Group_Conflict_Notes                                                                                                                     | ActionRequired | Action_Specific |
|:-----------------------|:---------------|---------------:|--------------:|:---------------|:-----------------------------------------------------------------------------------------------------------------------------------------|:---------------|:----------------|
| C00065, C00716         | NA             |              2 |             2 | C00065, C00716 | None                                                                                                                                     | Check          | KeepEachID      |
| C00155, C05330         | C00155         |              2 |             1 | C00155         | None                                                                                                                                     | None           | None            |
| C00221, C00031         | NA             |              2 |             2 | C00221, C00031 | None \|\| Pentose phosphate pathway \|\| Biosynthesis of secondary metabolites \|\| Glycolysis / Gluconeogenesis \|\| Metabolic pathways | Check          | KeepOneID       |
| C00671, C06008         | C00671         |              2 |             1 | C00671         | None                                                                                                                                     | None           | None            |
| C01835, C00721, C00420 | NA             |              3 |             2 | C01835, C00721 | None                                                                                                                                     | Check          | KeepEachID      |
| C01991, C00989         | C00989         |              2 |             1 | C00989         | None                                                                                                                                     | None           | None            |
| C02052, C00464         | C00464         |              2 |             1 | C00464         | None                                                                                                                                     | None           | None            |
| C03460, C03722         | NA             |              2 |             2 | C03460, C03722 | None \|\| Metabolic pathways                                                                                                             | Check          | KeepOneID       |
| C17737, C00695         | NA             |              2 |             2 | C17737, C00695 | Secondary bile acid biosynthesis \|\| None                                                                                               | Check          | KeepOneID       |

Terms in KEGG pathways where the same measured feature mapps with more
than one ID, which may inflate the enrichment analysis.

  
In case of `ActionRequired=="Check"`, we can look into the column
`Action_Specific` which contains additional information. In case of the
entry `KeepEachID`, multiple matches to the prior knowledge were found,
yet the features are in different pathways (=GroupingVariable). Yet, in
case of `KeepOneID`, the different IDs map to the same pathway in the
prior knowledge for at least one case and therefore keeping both would
inflate the enrichment analysis.  

``` r
SelectedIDs <- ccRCC_to_KEGGPathways[["data_summary"]]%>%
#Expand rows where Action == KeepEachID by splitting `matches`
dplyr::mutate(matches_split = if_else(Action_Specific == "KeepEachID", matches, NA_character_)) %>%
separate_rows(matches_split, sep = ",\\s*") %>%
mutate(InputID_select = if_else(Action_Specific  == "KeepEachID", matches_split, InputID_select)) %>%
select(-matches_split) %>%
#Select one ID for AcionSpecific==KeepOneID
dplyr::mutate(InputID_select = case_when(
    Action_Specific == "KeepOneID" & matches == "C03460, C03722" ~ "C03722", # 2-Methylprop-2-enoyl-CoA versus Quinolinate. No evidence, hence we keep the one present in more pathways ( C03722=7 pathways, C03460=2 pathway)
    Action_Specific == "KeepOneID" & matches ==  "C00221, C00031" ~ "C00031", # These are D- and L-Glucose. We have human samples, so in this conflict we will maintain L-Glucose
    Action_Specific == "KeepOneID" & matches ==  "C17737, C00695" ~ "C00695", # Allocholic acid versus Cholic acid. No evidence, hence we keep the one present in more pathways (C00695 = 4 pathways, C17737 = 1 pathway)
    Action_Specific == "KeepOneID" ~ InputID_select,  # Keep NA where not matched manually
    TRUE ~ InputID_select
    ))
```

  
Lastly, we need to add the column including our selected IDs to the
metadata table:

``` r
Tissue_MetaData_Extended <- merge(x= SelectedIDs %>%
                                        dplyr::select(KEGG, InputID_select),
                                        y= Tissue_MetaData_Extended,
                                        by= "KEGG",
                                        all.y=TRUE)
```

  
For results with p.adjusted value \< 0.1 and a minimum of 10% of the
pathway detected will be visualized as Volcano plots:  

``` r
# Since we have performed multiple comparisons (all_vs_HK2), we will run ORA for each of this comparison
DM_ORA_res<- list()

for(comparison in names(ResList)){#Res list includes the different comparisons we performed above <-
#Ensure that the Metabolite names match with KEGG IDs or KEGG trivial names.
dma_res <- merge(x= Tissue_MetaData_Extended,
y= ResList[[comparison]][["dma"]][["TUMOR_vs_NORMAL"]],
by="Metabolite",
all=TRUE)

#Ensure unique IDs and full background --> we include measured features that do not have a KEGG ID.
dma_res <- dma_res[,c(3,21:25)]%>%
    dplyr::mutate(InputID_select = if_else(
        is.na(InputID_select),
        paste0("NA_", cumsum(is.na(InputID_select))),
        InputID_select
        ))%>% #remove duplications and keep the higher Log2FC measurement
    group_by(InputID_select) %>%
    slice_max(order_by = Log2FC, n = 1, with_ties = FALSE) %>%
    ungroup()%>%
    remove_rownames()%>%
    tibble::column_to_rownames("InputID_select")

    #Perform ORA
    Res <- standard_ora(data= dma_res, #Input data requirements: column `t.val` and column `Metabolite`
    metadata_info=c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "MetaboliteID"),
    input_pathway=KEGG_Pathways,#Pathway file requirements: column `term`, `Metabolite` and `Description`. Above we loaded the Kegg_Pathways using Load_KEGG()
    pathway_name=paste0("KEGG_", comparison, sep=""),
    min_gssize=3,
    max_gssize=1000,
    cutoff_stat=0.01,
    cutoff_percentage=10)

    DM_ORA_res[[comparison]] <- Res

    #Select to plot:
    Res_Select <- Res[["ClusterGosummary"]]%>%
    filter(p.adjust<0.1)%>%
    #filter(pvalue<0.05)%>%
    filter(percentage_of_Pathway_detected>10)

    if(is.null(Res_Select)==FALSE){
        viz_volcano(plot_types="PEA",
        data= dma_res, #Must be the data you have used as an input for the pathway analysis
        data2=as.data.frame(Res_Select )%>%dplyr::rename("term"="ID"),
        metadata_info= c(PEA_Pathway="term",# Needs to be the same in both, metadata_feature and data2.
        PEA_stat="p.adjust",#Column data2
        PEA_score="GeneRatio",#Column data2
        PEA_Feature="MetaboliteID"),# Column metadata_feature (needs to be the same as row names in data)
        metadata_feature= KEGG_Pathways,#Must be the pathways used for pathway analysis
        plot_name= paste("KEGG_", comparison, sep=""),
        subtitle= "PEA" )
        }
}
```

  
  
  

### Biological regulated clustering

to understand which metabolites are changing independent of the patients
age, hence only due to tumour versus normal, and which metabolites
change independent of tumour versus normal, hence due to the different
age, we can use the [`mca_2cond()`](../../reference/mca_2cond.md)
function.  
Metabolite Clustering Analysis (`MCA`) enables clustering of metabolites
into groups based on logical regulatory rules. Here we set two different
thresholds, one for the differential metabolite abundance (Log2FC) and
one for the significance (e.g. p.adj). This will define if a feature (=
metabolite) is assigned into:  
1. “UP”, which means a metabolite is significantly up-regulated in the
underlying comparison.  
2. “DOWN”, which means a metabolite is significantly down-regulated in
the underlying comparison.  
3. “No Change”, which means a metabolite does not change significantly
in the underlying comparison and/or is not defined as
up-regulated/down-regulated based on the Log2FC threshold chosen.  
  
Thereby “No Change” is further subdivided into four states:  
1. “Not Detected”, which means a metabolite is not detected in the
underlying comparison.  
2. “Not Significant”, which means a metabolite is not significant in the
underlying comparison.  
3. “Significant positive”, which means a metabolite is significant in
the underlying comparison and the differential metabolite abundance is
positive, yet does not meet the threshold set for “UP” (e.g. Log2FC \>1
= “UP” and we have a significant Log2FC=0.8).  
4. “Significant negative”, which means a metabolite is significant in
the underlying comparison and the differential metabolite abundance is
negative, yet does not meet the threshold set for “DOWN”.  
  
For more information you can also check out the other vignettes.

``` r
MCAres <-  mca_2cond(data_c1=ResList[["Young"]][["dma"]][["TUMOR_vs_NORMAL"]],
data_c2=ResList[["Old"]][["dma"]][["TUMOR_vs_NORMAL"]],
metadata_info_c1=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
metadata_info_c2=c(ValueCol="Log2FC",StatCol="p.adj", cutoff_stat= 0.05, ValueCutoff=1),
feature = "Metabolite",
save_table = "csv",
method_background="C1&C2"#Most stringend background setting, only includes metabolites detected in both comparisons
)
```

  
Now we can use this information to colour code our volcano plot. We will
plot individual vocano plots for each metabolite pathway as defined by
the feature metadata provided as part of the data in (Hakimi et al.
2016).

``` r
# Add metabolite information such as KEGG ID or pathway to results
MetaData_Metab <- merge(x=Tissue_MetaData,
y= MCAres[["MCA_2Cond_Results"]][, c(1, 14:15)]%>%tibble::column_to_rownames("Metabolite"),
by=0,
all.y=TRUE)%>%
tibble::column_to_rownames("Row.names")%>%
dplyr::filter(!is.na(SUPER_PATHWAY))%>%
dplyr::filter(!is.na(SUB_PATHWAY))

viz_volcano(plot_types="Compare",
data=ResList[["Young"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
data2= ResList[["Old"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
name_comparison= c(data="Young", data2= "Old"),
metadata_feature =  MetaData_Metab,
plot_name= "Young-TUMOR_vs_NORMAL compared to Old-TUMOR_vs_NORMAL",
subtitle= "Results of dma",
metadata_info = c(individual = "SUPER_PATHWAY",
                                        color = "RG2_Significant"))

viz_volcano(plot_types="Compare",
data=ResList[["Young"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
data2= ResList[["Old"]][["dma"]][["TUMOR_vs_NORMAL"]]%>%tibble::column_to_rownames("Metabolite"),
name_comparison= c(data="Young", data2= "Old"),
metadata_feature =  MetaData_Metab,
plot_name= "Young-TUMOR_vs_NORMAL compared to Old-TUMOR_vs_NORMAL_Sub",
subtitle= "Results of dma",
metadata_info = c(individual = "SUB_PATHWAY",
                                        color = "RG2_Significant"))
```

  

![](sample-metadata_files/figure-html/viz-volcano-5-1.png)![](sample-metadata_files/figure-html/viz-volcano-5-2.png)![](sample-metadata_files/figure-html/viz-volcano-5-3.png)![](sample-metadata_files/figure-html/viz-volcano-5-4.png)![](sample-metadata_files/figure-html/viz-volcano-5-5.png)![](sample-metadata_files/figure-html/viz-volcano-5-6.png)![](sample-metadata_files/figure-html/viz-volcano-5-7.png)![](sample-metadata_files/figure-html/viz-volcano-5-8.png)

    #> Skipping viz_volcano compare plot for 'NA' because no rows are available after filtering.

![](sample-metadata_files/figure-html/viz-volcano-5-9.png)![](sample-metadata_files/figure-html/viz-volcano-5-10.png)![](sample-metadata_files/figure-html/viz-volcano-5-11.png)![](sample-metadata_files/figure-html/viz-volcano-5-12.png)![](sample-metadata_files/figure-html/viz-volcano-5-13.png)![](sample-metadata_files/figure-html/viz-volcano-5-14.png)![](sample-metadata_files/figure-html/viz-volcano-5-15.png)![](sample-metadata_files/figure-html/viz-volcano-5-16.png)![](sample-metadata_files/figure-html/viz-volcano-5-17.png)![](sample-metadata_files/figure-html/viz-volcano-5-18.png)![](sample-metadata_files/figure-html/viz-volcano-5-19.png)![](sample-metadata_files/figure-html/viz-volcano-5-20.png)![](sample-metadata_files/figure-html/viz-volcano-5-21.png)![](sample-metadata_files/figure-html/viz-volcano-5-22.png)![](sample-metadata_files/figure-html/viz-volcano-5-23.png)![](sample-metadata_files/figure-html/viz-volcano-5-24.png)![](sample-metadata_files/figure-html/viz-volcano-5-25.png)![](sample-metadata_files/figure-html/viz-volcano-5-26.png)![](sample-metadata_files/figure-html/viz-volcano-5-27.png)![](sample-metadata_files/figure-html/viz-volcano-5-28.png)![](sample-metadata_files/figure-html/viz-volcano-5-29.png)![](sample-metadata_files/figure-html/viz-volcano-5-30.png)![](sample-metadata_files/figure-html/viz-volcano-5-31.png)![](sample-metadata_files/figure-html/viz-volcano-5-32.png)![](sample-metadata_files/figure-html/viz-volcano-5-33.png)![](sample-metadata_files/figure-html/viz-volcano-5-34.png)![](sample-metadata_files/figure-html/viz-volcano-5-35.png)![](sample-metadata_files/figure-html/viz-volcano-5-36.png)![](sample-metadata_files/figure-html/viz-volcano-5-37.png)![](sample-metadata_files/figure-html/viz-volcano-5-38.png)![](sample-metadata_files/figure-html/viz-volcano-5-39.png)![](sample-metadata_files/figure-html/viz-volcano-5-40.png)![](sample-metadata_files/figure-html/viz-volcano-5-41.png)![](sample-metadata_files/figure-html/viz-volcano-5-42.png)![](sample-metadata_files/figure-html/viz-volcano-5-43.png)![](sample-metadata_files/figure-html/viz-volcano-5-44.png)![](sample-metadata_files/figure-html/viz-volcano-5-45.png)![](sample-metadata_files/figure-html/viz-volcano-5-46.png)![](sample-metadata_files/figure-html/viz-volcano-5-47.png)![](sample-metadata_files/figure-html/viz-volcano-5-48.png)![](sample-metadata_files/figure-html/viz-volcano-5-49.png)![](sample-metadata_files/figure-html/viz-volcano-5-50.png)![](sample-metadata_files/figure-html/viz-volcano-5-51.png)![](sample-metadata_files/figure-html/viz-volcano-5-52.png)![](sample-metadata_files/figure-html/viz-volcano-5-53.png)![](sample-metadata_files/figure-html/viz-volcano-5-54.png)![](sample-metadata_files/figure-html/viz-volcano-5-55.png)![](sample-metadata_files/figure-html/viz-volcano-5-56.png)![](sample-metadata_files/figure-html/viz-volcano-5-57.png)![](sample-metadata_files/figure-html/viz-volcano-5-58.png)![](sample-metadata_files/figure-html/viz-volcano-5-59.png)![](sample-metadata_files/figure-html/viz-volcano-5-60.png)![](sample-metadata_files/figure-html/viz-volcano-5-61.png)![](sample-metadata_files/figure-html/viz-volcano-5-62.png)![](sample-metadata_files/figure-html/viz-volcano-5-63.png)![](sample-metadata_files/figure-html/viz-volcano-5-64.png)![](sample-metadata_files/figure-html/viz-volcano-5-65.png)

    #> Skipping viz_volcano compare plot for 'NA' because no rows are available after filtering.

![](sample-metadata_files/figure-html/viz-volcano-5-66.png)![](sample-metadata_files/figure-html/viz-volcano-5-67.png)![](sample-metadata_files/figure-html/viz-volcano-5-68.png)![](sample-metadata_files/figure-html/viz-volcano-5-69.png)![](sample-metadata_files/figure-html/viz-volcano-5-70.png)![](sample-metadata_files/figure-html/viz-volcano-5-71.png)![](sample-metadata_files/figure-html/viz-volcano-5-72.png)![](sample-metadata_files/figure-html/viz-volcano-5-73.png)![](sample-metadata_files/figure-html/viz-volcano-5-74.png)![](sample-metadata_files/figure-html/viz-volcano-5-75.png)![](sample-metadata_files/figure-html/viz-volcano-5-76.png)![](sample-metadata_files/figure-html/viz-volcano-5-77.png)

  
  
  

### Pathway enrichment

Next, we perform Over Representation Analysis (`ORA`) using KEGG
pathways for each comparison.  

``` r
# Since we have performed multiple comparisons (all_vs_HK2), we will run ORA for each of this comparison
DM_ORA_res<- list()

KEGG_Pathways <- metsigdb_kegg()

for(comparison in names(ResList)){
    dma_res <- merge(x= Tissue_MetaData_Extended,
    y=ResList[[comparison]][["dma"]][["TUMOR_vs_NORMAL"]],
    by="Metabolite",
    all=TRUE)

    #Ensure unique IDs and full background --> we include measured features that do not have a KEGG ID.
    dma_res <- dma_res[,c(3,21:25)]%>%
    dplyr::mutate(InputID_select = if_else(
        is.na(InputID_select),
        paste0("NA_", cumsum(is.na(InputID_select))),
        InputID_select
        ))%>% #remove duplications and keep the higher Log2FC measurement
    group_by(InputID_select) %>%
    slice_max(order_by = Log2FC, n = 1, with_ties = FALSE) %>%
    ungroup()%>%
    remove_rownames()%>%
    tibble::column_to_rownames("InputID_select")

    #Perform ORA
    Res <- standard_ora(data= dma_res, #Input data requirements: column `t.val` and column `Metabolite`
    metadata_info=c(pvalColumn="p.adj", percentageColumn="t.val", PathwayTerm= "term", PathwayFeature= "MetaboliteID"),
    input_pathway=KEGG_Pathways,#Pathway file requirements: column `term`, `Metabolite` and `Description`. Above we loaded the Kegg_Pathways using Load_KEGG()
    pathway_name=paste0("KEGG_", comparison, sep=""),
    min_gssize=3,
    max_gssize=1000,
    cutoff_stat=0.01,
    cutoff_percentage=10)

    DM_ORA_res[[comparison]] <- Res

    #Select to plot:
    Res_Select <- Res[["ClusterGosummary"]]%>%
    filter(p.adjust<0.1)%>%
    #filter(pvalue<0.05)%>%
    filter(percentage_of_Pathway_detected>10)

    if(is.null(Res_Select)==FALSE){
        viz_volcano(plot_types="PEA",
        data= dma_res, #Must be the data you have used as an input for the pathway analysis
        data2=as.data.frame(Res_Select )%>%dplyr::rename("term"="ID"),
        metadata_info= c(PEA_Pathway="term",# Needs to be the same in both, metadata_feature and data2.
        PEA_stat="p.adjust",#Column data2
        PEA_score="GeneRatio",#Column data2
        PEA_Feature="MetaboliteID"),# Column metadata_feature (needs to be the same as row names in data)
        metadata_feature= KEGG_Pathways,#Must be the pathways used for pathway analysis
        plot_name= paste("KEGG_", comparison, sep=""),
        subtitle= "PEA" )
        }
}
```

  

![](sample-metadata_files/figure-html/viz-volcano-7-1.png)![](sample-metadata_files/figure-html/viz-volcano-7-2.png)![](sample-metadata_files/figure-html/viz-volcano-7-3.png)![](sample-metadata_files/figure-html/viz-volcano-7-4.png)![](sample-metadata_files/figure-html/viz-volcano-7-5.png)![](sample-metadata_files/figure-html/viz-volcano-7-6.png)![](sample-metadata_files/figure-html/viz-volcano-7-7.png)

  
  
  

## Session information

    #> R version 4.5.2 (2025-10-31)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.3 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    #>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    #>  [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    #> 
    #> time zone: Etc/UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] stringr_1.6.0      tibble_3.3.0       tidyr_1.3.1        rlang_1.1.6        dplyr_1.1.4        magrittr_2.0.4    
    #> [7] MetaProViz_3.99.26 BiocStyle_2.38.0  
    #> 
    #> loaded via a namespace (and not attached):
    #>   [1] RColorBrewer_1.1-3          rstudioapi_0.17.1           jsonlite_2.0.0              ggbeeswarm_0.7.2           
    #>   [5] farver_2.1.2                rmarkdown_2.30              fs_1.6.6                    ragg_1.5.0                 
    #>   [9] vctrs_0.6.5                 memoise_2.0.1               rstatix_0.7.3               htmltools_0.5.8.1          
    #>  [13] S4Arrays_1.10.0             progress_1.2.3              curl_7.0.0                  ComplexUpset_1.3.3         
    #>  [17] decoupleR_2.16.0            broom_1.0.10                cellranger_1.1.0            SparseArray_1.10.2         
    #>  [21] Formula_1.2-5               sass_0.4.10                 parallelly_1.45.1           bslib_0.9.0                
    #>  [25] htmlwidgets_1.6.4           desc_1.4.3                  plyr_1.8.9                  httr2_1.2.1                
    #>  [29] lubridate_1.9.4             cachem_1.1.0                igraph_2.2.1                lifecycle_1.0.4            
    #>  [33] pkgconfig_2.0.3             Matrix_1.7-4                R6_2.6.1                    fastmap_1.2.0              
    #>  [37] MatrixGenerics_1.22.0       selectr_0.5-0               digest_0.6.39               colorspace_2.1-2           
    #>  [41] patchwork_1.3.2             S4Vectors_0.48.0            textshaping_1.0.4           GenomicRanges_1.62.0       
    #>  [45] RSQLite_2.4.4               ggpubr_0.6.2                labeling_0.4.3              timechange_0.3.0           
    #>  [49] httr_1.4.7                  abind_1.4-8                 compiler_4.5.2              bit64_4.6.0-1              
    #>  [53] withr_3.0.2                 S7_0.2.1                    backports_1.5.0             BiocParallel_1.44.0        
    #>  [57] carData_3.0-5               DBI_1.2.3                   logger_0.4.1                OmnipathR_3.19.2           
    #>  [61] R.utils_2.13.0              ggsignif_0.6.4              cosmosR_1.18.0              MASS_7.3-65                
    #>  [65] rappdirs_0.3.3              DelayedArray_0.36.0         sessioninfo_1.2.3           scatterplot3d_0.3-44       
    #>  [69] gtools_3.9.5                tools_4.5.2                 vipor_0.4.7                 beeswarm_0.4.0             
    #>  [73] qcc_2.7                     zip_2.3.3                   R.oo_1.27.1                 glue_1.8.0                 
    #>  [77] grid_4.5.2                  checkmate_2.3.3             reshape2_1.4.5              generics_0.1.4             
    #>  [81] gtable_0.3.6                tzdb_0.5.0                  R.methodsS3_1.8.2           ggVennDiagram_1.5.4        
    #>  [85] hms_1.1.4                   xml2_1.5.0                  car_3.1-3                   XVector_0.50.0             
    #>  [89] BiocGenerics_0.56.0         ggrepel_0.9.6               pillar_1.11.1               vroom_1.6.6                
    #>  [93] limma_3.66.0                later_1.4.4                 splines_4.5.2               lattice_0.22-7             
    #>  [97] bit_4.6.0                   tidyselect_1.2.1            knitr_1.50                  gridExtra_2.3              
    #> [101] bookdown_0.45               IRanges_2.44.0              Seqinfo_1.0.0               SummarizedExperiment_1.40.0
    #> [105] svglite_2.2.2               stats4_4.5.2                xfun_0.54                   Biobase_2.70.0             
    #> [109] statmod_1.5.1               factoextra_1.0.7            matrixStats_1.5.0           pheatmap_1.0.13            
    #> [113] stringi_1.8.7               yaml_2.3.10                 kableExtra_1.4.0            evaluate_1.0.5             
    #> [117] codetools_0.2-20            tcltk_4.5.2                 qvalue_2.42.0               hash_2.2.6.3               
    #> [121] BiocManager_1.30.27         Polychrome_1.5.4            cli_3.6.5                   systemfonts_1.3.1          
    #> [125] jquerylib_0.1.4             EnhancedVolcano_1.13.2      Rcpp_1.1.0                  readxl_1.4.5               
    #> [129] XML_3.99-0.20               parallel_4.5.2              ggfortify_0.4.19            pkgdown_2.2.0              
    #> [133] ggplot2_4.0.1               readr_2.1.6                 blob_1.2.4                  prettyunits_1.2.0          
    #> [137] viridisLite_0.4.2           scales_1.4.0                writexl_1.5.4               inflection_1.3.7           
    #> [141] purrr_1.2.0                 crayon_1.5.3                rvest_1.0.5

## Bibliography

Hakimi, A Ari, Ed Reznik, Chung-Han Lee, Chad J Creighton, A Rose
Brannon, Augustin Luna, B Arman Aksoy, et al. 2016. “An Integrated
Metabolic Atlas of Clear Cell Renal Cell Carcinoma.” *Cancer Cell* 29
(1): 104–16. <https://doi.org/10.1016/j.ccell.2015.12.004>.
