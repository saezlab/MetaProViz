# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Important Note on Code Quality

This package was initially developed by a beginner programmer and
consequently has several areas requiring improvement:

- **Coding style** - Inconsistent formatting, naming conventions, and R
  idioms
- **Code quality** - Variable naming issues, redundant code, suboptimal
  patterns
- **Code organization** - Some functions are overly complex or poorly
  structured
- **API design** - Inconsistent parameter patterns, unclear function
  interfaces

Many of these issues can be addressed through standard refactoring and
application of R best practices. When working on this codebase,
prioritize: 1. **Correctness** - Fix bugs and ensure functionality works
as expected 2. **Consistency** - Apply R community conventions
(tidyverse style guide, etc.) 3. **Clarity** - Improve readability and
maintainability 4. **Simplicity** - Reduce complexity where possible

Be prepared for ongoing refactoring tasks to improve the overall code
quality.

## Project Overview

MetaProViz is a Bioconductor R package for **METabolomics
pre-PRocessing, functiOnal analysis and VIZualisation**. It processes
both intracellular metabolomics and exometabolomics
(consumption-release/CoRe) data through five integrated modules:

1.  **Processing** - Feature filtering, missing value imputation,
    normalization, outlier detection
2.  **Differential Metabolite Analysis (DMA)** - Statistical comparison
    between conditions
3.  **Functional Analysis** - Metabolite clustering and
    over-representation analysis (ORA)
4.  **Prior Knowledge** - Access and refactoring of KEGG, RaMP, and
    other metabolic databases
5.  **Visualization** - Publication-ready plots (PCA, heatmap, volcano,
    superplots, upset plots)

## Essential Commands

### Building and Checking

``` bash
# Build the package tarball
R CMD build .

# Check the package (standard R package checks)
R CMD check *.tar.gz --no-manual

# Run BiocCheck (Bioconductor-specific checks)
Rscript -e 'BiocManager::install("BiocCheck")'
Rscript -e 'BiocCheck::BiocCheck(list.files(pattern = "*.tar.gz"))'
```

### Testing

``` bash
# Run all tests
Rscript -e 'devtools::test()'

# Run specific test file
Rscript -e 'testthat::test_file("tests/testthat/test-checkmatch_pk_to_data.R")'

# Run tests with coverage
Rscript -e 'covr::package_coverage()'
```

### Installation

``` bash
# Install from local source
Rscript -e 'devtools::install()'

# Install with dependencies
Rscript -e 'devtools::install_deps(dependencies = TRUE)'

# Install OmnipathR from Bioconductor devel (required dependency)
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
Rscript -e 'BiocManager::install(version = "devel")'
Rscript -e 'BiocManager::install("OmnipathR")'
```

### Documentation

``` bash
# Build documentation
Rscript -e 'devtools::document()'

# Build vignettes
Rscript -e 'devtools::build_vignettes()'

# Build pkgdown site
Rscript -e 'pkgdown::build_site()'
```

### Loading and Debugging

``` bash
# Load package for development
Rscript -e 'devtools::load_all()'

# Check for code style issues
Rscript -e 'lintr::lint_package()'
```

## Code Architecture

### Module Organization

The codebase is organized into distinct functional modules in `/R`:

**Core Analysis Modules:** - `Processing.R` - Data preprocessing
pipeline (filtering, imputation, normalization, outlier detection) -
`DifferentialMetaboliteAnalysis.R` - [`dma()`](reference/dma.md)
function for statistical comparison - `MetaboliteClusteringAnalysis.R` -
[`mca_2cond()`](reference/mca_2cond.md), `mca_3cond()` for clustering by
regulatory rules - `MetaDataAnalysis.R` -
[`metadata_analysis()`](reference/metadata_analysis.md) for PCA and
ANOVA on sample metadata - `OverRepresentationAnalysis.R` -
[`cluster_ora()`](reference/cluster_ora.md) for pathway enrichment
analysis

**Prior Knowledge:** - `GetPriorKnowledge.R` - Functions to load KEGG,
RaMP pathways (`metsigdb_*` functions) - `RefactorPriorKnoweldge.R` -
[`translate_id()`](reference/translate_id.md) for ID translation between
databases

**Visualization (Viz\* files):** - `VizVolcano.R`, `VizHeatmap.R`,
`VizPCA.R` - Standard omics plots - `VizSuperplots.R` - Box/violin/bar
plots with statistical overlays - `VizStackedbar.R`, `VizUpset.R` -
Advanced visualization types - `HelperPlots.R` - Shared plotting
utilities and themes

**Infrastructure (Helper\* files):** - `HelperChecks.R` - Input
validation for all functions - `HelperLog.R` - Logging system (wraps
OmnipathR/logger) - `HelperSave.R` - File I/O for tables and plots -
`HelperOptions.R` - Configuration management - `HelperMisc.R` - Utility
functions

**Special files:** - `aaa_MetaProViz_Env.R` - Global environment for
state management (loaded first) - `zzz.R` - Package initialization hooks
(loaded last) - `enricher_internal.R` - Internal copy of DOSE package
enrichment logic

### Data Flow Pattern

All functions follow a standardized input/output pattern:

**Inputs:** - `data` - Data frame with samples as rows, metabolites as
columns - `metadata_sample` - Sample metadata with `Conditions` column
(required) - `metadata_feature` - Feature/metabolite metadata
(optional) - `metadata_info` - Named vector for flexible parameter
specification

**Outputs:**

``` r
list(
  DF = list(dataframe_1 = ..., dataframe_2 = ...),
  Plot = list(plot_1 = ggplot_object, plot_2 = ggplot_object)
)
```

### Function Execution Pattern

All public functions follow this structure:

``` r
function_name <- function(data, metadata_sample, metadata_info, ...) {
  # 1. Initialize logging
  metaproviz_init()

  # 2. General validation
  check_param(...)

  # 3. Function-specific validation
  check_param_<function_name>(...)

  # 4. Create output folders
  folder <- save_path(folder_name = "ModuleName", path = path)

  # 5. Main computation
  # ...

  # 6. Save results
  save_res(inputlist_df, inputlist_plot, ...)

  # 7. Return structured list
  return(list(DF = DFList, Plot = PlotList))
}
```

### Logging System

- **Initialization:** Call `metaproviz_init()` at the start of every
  public function
- **State Management:** Uses `metaproviz.env` isolated environment (not
  global variables)
- **Integration:** Wraps OmnipathR’s logging infrastructure (`logger`
  package)
- **Access:** Users can view logs via
  [`metaproviz_logfile()`](reference/metaproviz_logfile.md) and
  [`metaproviz_log()`](reference/metaproviz_log.md)
- **Levels:** Trace, Debug, Info, Success, Warn, Error

### Validation Strategy

- **Early validation:** All inputs checked before computation via
  `check_param()` functions
- **Reusable checks:** Validation logic is decoupled in `HelperChecks.R`
- **Informative errors:** Uses `logger` package for detailed error
  context
- **User-friendly messages:** Errors include suggested fixes

### Prior Knowledge Integration

- **Sources:** KEGG, RaMP (chemical classes), custom pathways
- **Format:** All PK files follow TERM2GENE/TERM2NAME pattern
- **Caching:** Data fetched from OmnipathR on demand, cached locally via
  [`rappdirs::user_cache_dir()`](https://rappdirs.r-lib.org/reference/user_cache_dir.html)
- **ID Translation:** [`translate_id()`](reference/translate_id.md)
  enables cross-database mapping (KEGG ↔︎ PubChem ↔︎ CHEBI ↔︎ HMDB)

## Important Development Patterns

### Naming Conventions

- **Public functions:** lowercase with underscores (e.g., `viz_volcano`,
  `dma`, `processing`)
- **Helper functions:** lowercase with underscores, prefixed (e.g.,
  `check_param_*`, `save_*`)
- **Internal functions:** Marked with `@noRd` in roxygen2 documentation
- **Files:** `Viz*.R` for visualization, `Helper*.R` for utilities

### State Management

- **DO NOT use global variables** - Use `metaproviz.env` environment
  instead
- **Package-wide state:** Accessed via `metaproviz.env$variable_name`
- **R Options:** Can be set via `metaproviz_*_config()` functions

### Output Organization

- **Default path:** `./MetaProViz_Results/` if not specified
- **Timestamped files:** All outputs include YYYYMMDD suffix for
  reproducibility
- **Modular saving:** Each function can independently save tables/plots
- **Formats:** Tables (CSV/XLSX/TXT), Plots (SVG/PDF/PNG)

### Adding New Functions

**For analysis functions:** 1. Start with `metaproviz_init()` +
validation 2. Use `logger::log_*()` for all major steps 3. Return
`list(DF = ..., Plot = ...)` structure 4. Use `save_res()` for output
management 5. Add roxygen docs with `@export`

**For visualization functions:** 1. Use utilities from `HelperPlots.R`
2. Return both raw ggplot2 object AND pre-sized version 3. Support
`metadata_info` for flexible coloring/shaping 4. Include automatic label
management for large datasets 5. Implement smart defaults but allow full
customization

**For prior knowledge sources:** 1. Add function to
`GetPriorKnowledge.R` 2. Follow TERM2GENE/TERM2NAME format 3. Cache
downloaded data in user cache directory 4. Integrate with
[`cluster_ora()`](reference/cluster_ora.md) and
[`translate_id()`](reference/translate_id.md)

## Dependency Notes

### Critical Dependencies

- **OmnipathR** (\>= 3.17.4) - Infrastructure for logging, config, and
  biological database access
- **limma** - Linear modeling for differential analysis
- **ggplot2** (\>= 3.3.5) - Visualization engine
- **ComplexUpset** (\>= 1.3.3) - Advanced upset plots

### Bioconductor Version

This package targets Bioconductor 3.22 and requires R \>= 4.4. OmnipathR
must be installed from Bioconductor devel.

### Installation Order

Some dependencies have specific installation requirements: 1. Install
BiocManager and set Bioc version 2. Install OmnipathR from devel 3.
Install EnhancedVolcano from GitHub 4. Install other dependencies 5.
Install ggplot2 from source (Linux/Ubuntu) 6. Install ComplexUpset from
source (after ggplot2)

## Testing Guidelines

- **Framework:** Uses testthat (\>= 3.1.4)
- **Location:** All tests in `tests/testthat/test-*.R`
- **Coverage:** Focus on validation functions and core analysis logic
- **Test data:** Use toy datasets from `ToyData.R`

## Continuous Integration

The package uses GitHub Actions for CI/CD: - **Build & Check:**
`.github/workflows/build-check.yaml` runs R CMD check on macOS, Ubuntu,
Windows - **pkgdown:** `.github/workflows/pkgdown-github-pages.yaml`
builds documentation site - **Triggers:** Push to main/devel branches,
daily scheduled builds

## Key Design Principles

1.  **Modular & Composable** - Functions can be chained; each returns
    standardized lists
2.  **Validation-First** - All public functions validate early using
    helper functions
3.  **Logging Everywhere** - Use `logger` package to track operations
    for debugging and reproducibility
4.  **Flexible Configuration** - Use `metadata_info` named vectors for
    parameter specification
5.  **Standardized I/O** - Data always flows: samples as rows,
    metabolites as columns
6.  **Publication-Ready Visualization** - Plots automatically sized with
    smart defaults
7.  **Integrated Infrastructure** - Built on OmnipathR for logging,
    config, and data access
8.  **Optional Pipeline Steps** - Users can skip preprocessing steps as
    needed
9.  **State Isolation** - Use `metaproviz.env` for package-wide state
    (not global variables)

## Common Patterns

### Statistical Testing

Functions like [`dma()`](reference/dma.md) and
[`viz_superplot()`](reference/viz_superplot.md) accept string
abbreviations for tests: - **One-vs-one:** t.test, wilcox.test,
cor.test, chisq.test, lmFit - **One-vs-all/all-vs-all:** aov, welch
(Welch ANOVA), kruskal.test, lmFit - **P-value adjustment:** BH, fdr,
bonferroni, holm, etc.

### Processing Pipeline

The [`processing()`](reference/processing.md) function is modular with
optional steps: - Feature filtering (Standard/Modified/NULL) - Missing
value imputation (TRUE/FALSE) - TIC normalization (TRUE/FALSE) - CoRe
normalization for exometabolomics (TRUE/FALSE) - Outlier detection
(always runs)

Each step can be called independently for debugging.

### Visualization Workflow

Every `viz_*()` function: 1. Returns interactive-ready ggplot2 objects
2. Supports flexible metadata annotations 3. Includes automatic sizing
logic for large datasets 4. Has consistent export options (svg/pdf/png)
5. Supports theme customization 6. Saves with timestamped filenames

## File References

When investigating specific functionality, refer to these key files:

- Input validation: `R/HelperChecks.R`
- Logging setup: `R/HelperLog.R`
- Plot utilities: `R/HelperPlots.R`
- Save/export logic: `R/HelperSave.R`
- Package state: `R/aaa_MetaProViz_Env.R`
- Example data: `R/ToyData.R`
