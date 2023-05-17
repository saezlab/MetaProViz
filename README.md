# MetaProViz - Metabolomics Processing and Vizualisation
## Overview
**`MetaProViz`** is an R package for the pre-processing, downstream analysis and visualisation of metabolomics data.



## Install
**`MetaProViz`** is an R package.
1. Install Rtools if you haven't done this yet, using the appropriate version (e.g.[windows](https://cran.r-project.org/bin/windows/Rtools/) or [macOS](https://cran.r-project.org/bin/macosx/tools/)).
2. Install the latest development version from GitHub with: **`SiRCleR package`** direcly in R:
    ```
    #install.packages("devtools")
    devtools::install_github()
    library(MetaProViz)
    ```
### Dependencies 
If you are using the visualisations you will need to install the following tools and cite them.\
1. CRAN packages
```

```
2. Biocmanager packages
```

```
While we have done our best to ensure all the dependencies are documented, if they aren't please let us know and we will try to resolve them.

### Windows specifications
Note if you are running **Windows** you might have an issue with long paths, which you can resolve in the registry on Windows 10:\
`Computer Configuration > Administrative Templates > System > Filesystem > Enable Win32 long paths`\
(If you have a different version of Windows, just google "Long paths fix" and your Windows version)

## Tutorial
Provide link to main tutorial

Add minitutorial here:

List the output files that can be generated:

Output preprocessing excel table:
Sheet_1 = Experimental_design (as provided by user)
Sheet_2 = Raw_data (as provided by user)
Sheet 3= processed_data (= all samples are there with the additional column of Outlier_detection and Experimental design, metabolites that were removed due to filtering are of course missing)
Sheet 4= processed_data_matrix_OutliersRemoved (input for downstream analysis)



## Citation
