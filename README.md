
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetaProViz <img src="vignettes/Hexagon_MetaProViz.png" align="right" width="200" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![GitHub
issues](https://img.shields.io/github/issues/saezlab/MetaProViz)](https://github.com/saezlab/MetaProViz/issues)
<!-- badges: end -->

## **Short Introduction**

**MetaProViz** enables the user to pre-process metabolomics data
including consumption-release (CoRe) data, to perform differential
analysis (DMA), do clustering based on regulatory rules (MCA), pathway
analysis (ORA) and contains different visualization methods to extract
biological interpretable graphs.

## Tutorials

We have generated several tutorials showcasing the different
functionalities MetaProViz offers using publicly available datasets,
which are included as example data within **MetaProViz**. You can find
those tutorial on the top under the “Tutorial” button, where you can
follow specific user case examples for different analysis. Otherwise,
you can also follow the links below:  
- [Standard metabolomics
data](https://github.com/saezlab/MetaProViz/docs/articles/Standard%20Metabolomics.html)  
- [Consumption-Release (CoRe) metabolomics data from cell culture
media](https://github.com/saezlab/MetaProViz/docs/articles/CoRe%20Metabolomics.html)  
  
Here you will find a brief overview and information about the
installation of the package and its dependencies.

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

### Dependencies

If you are using **MetaProViz** the following packages are require:  

    "tidyverse"
    "ggplot2"
    "factoextra"
    "qcc"
    "hash"
    "reshape"
    "gridExtra"
    "inflection"
    "patchwork"
    "clusterProfiler"
    "ggupset"
    "gtools"
    "EnhancedVolcano"
    "writexl"
    "pheatmap"
    "ggfortify"

  
While we have done our best to ensure all the dependencies are
documented, if they aren’t please let us know and we will try to resolve
them.

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
      author = {Christina Schmidt and Dimitrios Prymidis and Julio Saez-Rodriguez and Christian Frezza},
      year = {2023},
      note = {R package version 0.0.0.9000},
    }
