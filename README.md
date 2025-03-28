
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaProViz <img src="vignettes/Hexagon_MetaProViz.png" align="right" width="200" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![GitHub
issues](https://img.shields.io/github/issues/saezlab/MetaProViz)](https://github.com/saezlab/MetaProViz/issues)
<!-- badges: end -->

## **Short Introduction**

**MetaProViz** (**Meta**bolomics **Pr**ocessing, functi**o**nal analysis
and **Vi**suali**z**ation), a free open-source R-package that provides
mechanistic hypotheses from metabolomics data by integrating prior
knowledge from literature with metabolomics. MetaProViz offers an
interactive framework consisting of five modules: Processing,
differential analysis, prior knwoledge access and refactoring,
functional analysis and visualization of both intracellular and
exometabolomics (=consumption-release (CoRe) data). Those modules and
their functions can be used independently from each other or in
combination (**Fig.1**).

  
<center>

<figure>
<img
src="https://github.com/saezlab/MetaProViz/blob/development/vignettes/Fig.1.png?raw=true"
style="width:90.0%" alt="Fig. 1: Overview of MetaProViz functions." />
<figcaption aria-hidden="true"><strong>Fig. 1:</strong> Overview of
MetaProViz functions.</figcaption>
</figure>

</center>

  
  
The first module, **MetaProViz**, `Processing`, allows the customized
processing of raw peak metabolomics data from different experimental
setups, including options to perform feature filtering due to
missingness, Total Ion Count (TIC) normalisation, Missing Value
Imputation (MVI) based on half-minimum and outlier detection based on
Hotellin’s T2. All of these pre-processing parameters can be customized
and combined as needed.  
The second module of **MetaProViz**,
`Differential Metabolite Analysis (DMA)`, allows the user to perform
differential analysis between two conditions (e.g. Tumour versus
Healthy) calculating the Log2FC, p-value, adjusted p-value and t-value,
whereby the user can choose all the test statistics. The input can
either be the output of the `Preprocessing` module or any DF including
metabolite values and information about the conditions that should be
compared.  
The third module of **MetaProViz**, `Functional Analysis`, includes
different methods to create clusters of metabolites based on their
distribution across the data using logical regulatory rules, prior
knowledge for enrichment analysis and functions to perform over
representation analysis (ORA). Here, the user can either input the
output of the `Processing` or `Differential Metabolite Analysis (DMA)`
module, or any other DF including Log2FC and statistics or metabolite
values.  
The fourth module of **MetaProViz**, `Visualization`, can easily create
customized visualizations of the output results of each of the other
**MetaProViz** modules or custom files. Here we not only enable overview
plots such as PCA, heatmap, Volcano plot, but also individual graphs of
each metabolite as bar graphs, box plots or violin plots. Moreover, the
user can provide additional information such as pathways the metabolites
correspond to, the clusters the metabolites where assigned to or any
other meta-information to customize the plots for color, shape or
selections, thus enabling biological interpretation of the results
otherwise missed in the data.

## Tutorials

We have generated several tutorials showcasing the different
functionalities MetaProViz offers using publicly available datasets,
which are included as example data within **MetaProViz**. You can find
those tutorial on the top under the “Tutorials” button, where you can
follow specific user case examples for different analysis. Otherwise,
you can also follow the links below:  
- [Standard metabolomics
data](https://saezlab.github.io/MetaProViz/articles/Standard%20Metabolomics.html)  
- [Consumption-Release (CoRe) metabolomics data from cell culture
media](https://saezlab.github.io/MetaProViz/articles/CoRe%20Metabolomics.html)  
- [Prior Knowledge Access &
Integration](https://saezlab.github.io/MetaProViz/articles/Prior%20Knowledge%20-%20Access%20&%20Integration.html)  
- [Sample Metadata
Analysis](https://saezlab.github.io/MetaProViz/articles/Sample%20Metadata%20Analysis.html)  
  
Here you will find a brief overview and information about the
installation of the package and its dependencies.

## Example Data

clear cell Renal Cell Carcinoma (ccRCC) patients data from Hakimi et. al
including 138 matched tumour and normal tissue pairs (Hakimi et al.
2016). Cell-lines data from intra- and extracellular metabolomics data
from cell culture media from [metabolomics workbench project
PR001418](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001418).
<img src="vignettesunnamed-chunk-2-1.png" width="100%" style="display: block; margin: auto auto auto 0;" />  
Additionally we also added transcriptomics and proteomics data of ccRCC
patients processed with SiRCle (Mora et al. 2024), originally from
[PDC000127](https://proteomic.datacommons.cancer.gov/pdc/study/PDC000127)
(Clark et al. 2019).

## Installation

**MetaProViz** is an R package and to install the package, start R and
enter:

``` r
devtools::install_github("https://github.com/saezlab/MetaProViz")
```

Now **MetaProViz** can be imported as:

``` r
library(MetaProViz)
```

### Windows specifications

Note if you are running Windows you might have an issue with long paths,
which you can resolve in the registry on Windows 10: Computer
Configuration \> Administrative Templates \> System \> Filesystem \>
Enable Win32 long paths (If you have a different version of Windows,
just google “Long paths fix” and your Windows version)

## Liscence

GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007

## Citation

    @Manual{,
      title = {MetaProViz: METabolomics pre-PRocessing, functiOnal analysis and VIZualisation},
      author = {Christina Schmidt, Dimitrios Prymidis, Macabe Daley, Denes Turei, Julio Saez-Rodriguez and Christian Frezza},
      year = {2023},
      note = {R package version 2.1.5},
    }

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Clark2019" class="csl-entry">

Clark, David J, Saravana M Dhanasekaran, Francesca Petralia, Jianbo Pan,
Xiaoyu Song, Yingwei Hu, Felipe da Veiga Leprevost, et al. 2019.
“Integrated Proteogenomic Characterization of Clear Cell Renal Cell
Carcinoma.” *Cell* 179 (4): 964–983.e31.
<https://doi.org/10.1016/j.cell.2019.10.007>.

</div>

<div id="ref-Hakimi2016" class="csl-entry">

Hakimi, A Ari, Ed Reznik, Chung-Han Lee, Chad J Creighton, A Rose
Brannon, Augustin Luna, B Arman Aksoy, et al. 2016. “An Integrated
Metabolic Atlas of Clear Cell Renal Cell Carcinoma.” *Cancer Cell* 29
(1): 104–16. <https://doi.org/10.1016/j.ccell.2015.12.004>.

</div>

<div id="ref-Mora_Schmidt2024" class="csl-entry">

Mora, Ariane, Christina Schmidt, Brad Balderson, Christian Frezza, and
Mikael Bodén. 2024. “SiRCle (Signature Regulatory Clustering) Model
Integration Reveals Mechanisms of Phenotype Regulation in Renal Cancer.”
*Genome Medicine* 16 (1): 144.
<https://doi.org/10.1186/s13073-024-01415-3>.

</div>

</div>
