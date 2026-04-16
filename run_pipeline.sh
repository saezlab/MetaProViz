#!/usr/bin/env bash

# Exit immediately on error
set -e

echo "📦 Checking and initializing renv..."

# Step 0: Do you have R and sepcify which R to use?
# Define R installation path
R_PATH="/c/Program Files/R/R-4.x.x/bin/Rscript.exe"

# Minimal fallback for Git Bash when Rscript is not on PATH
if ! command -v Rscript >/dev/null 2>&1; then
    R_PATH=$(ls -1d /c/Program\ Files/R/R-*/bin/Rscript.exe 2>/dev/null | sort -V | tail -n 1)
    [ -x "$R_PATH" ] || { echo "❌ Rscript not found. Set R_PATH in run_pipeline.sh"; exit 1; }
    Rscript() { "$R_PATH" "$@"; }
fi


# Step 1: Install or load renv in R
# Initialize renv and install packages
Rscript - <<EOF
# Check if renv is installed, install if needed
if (!requireNamespace("renv", quietly = TRUE)) {
    message("📦 Installing renv package...")
    install.packages("renv")
}

# Ensure renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
    install.packages("renv")
}

# Define all required packages
packages <- c(
    "magrittr", "dplyr", "tibble", "rlang",
    "stringr", "tidyverse", "tidyr", "purrr", "rmarkdown",
    "remotes", "R.utils", "R.oo", "R.methodsS3", "lazyeval"
)

# Check if lockfile exists
if (!file.exists("renv.lock")) {
    message("🔧 No lockfile found. Initializing new environment...")
    renv::init(force = TRUE)

    # Install all CRAN packages at once
    message("📥 Installing required packages...")
    install.packages(packages, dependencies = NA)

    # Install GitHub packages
    message("📥 Installing MetaProViz from GitHub...")
    remotes::install_github("saezlab/MetaProViz")

    # Create initial snapshot
    message("📸 Creating initial snapshot...")
    renv::snapshot(force = TRUE)
} else {
    message("♻️ Loading existing environment...")
    renv::load()

    if (!renv::status()$synchronized) {
        message("⚠️ Restoring from lockfile...")
        renv::restore(prompt = FALSE)
    }
}

# Ensure runtime packages needed by MetaProViz are available
runtime_packages <- c("R.utils", "R.oo", "R.methodsS3", "lazyeval")
missing_runtime <- runtime_packages[!vapply(runtime_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_runtime) > 0) {
    install.packages(missing_runtime, dependencies = NA)
}
EOF

# Step 2: Render Rmd files in order
echo "📝 Rendering RMarkdown files..."
Rscript -e "renv::load(); rmarkdown::render('Cells_CoRe/core-metabolomics_revison.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_Intracellular/intra-metabolomics_revision.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('MetSigDB_Benchmark/metsigdb_benchmark.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Patients_Metadata/Fig2_FeatureMetadata.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Patients_Metadata/Fig5_EnrichmentAnalysis.Rmd')"

# Step 3: Run Source data tables
echo "📊 Creating source data tables..."
Rscript -e "renv::load(); rmarkdown::render('Supplementary_Tables.Rmd')"

# Step 4: Run supplementary tables last
echo "📊 Creating supplementary tables..."
Rscript -e "renv::load(); rmarkdown::render('SourceDataTables.Rmd')"

# Step 4.5: Remove stale package directories missing DESCRIPTION
echo "🧹 Cleaning stale package folders..."
Rscript -e "renv::load(); lib <- renv::paths\$library(); pkgs <- list.dirs(lib, full.names = TRUE, recursive = FALSE); stale <- pkgs[!file.exists(file.path(pkgs, 'DESCRIPTION'))]; if (length(stale) > 0) { unlink(stale, recursive = TRUE, force = TRUE); message('Removed stale folders: ', paste(basename(stale), collapse = ', ')) } else { message('No stale folders found.') }"

# Step 4.6: Align Bioconductor release for strict snapshot validation
echo "🧬 Aligning Bioconductor packages..."
Rscript -e "renv::load(); renv::settings\$bioconductor.version('3.20'); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); options(repos = BiocManager::repositories(version = '3.20')); cur <- tryCatch(as.character(utils::packageVersion('S4Vectors')), error = function(e) ''); if (cur != '0.44.0') BiocManager::install('S4Vectors', ask = FALSE, update = FALSE, force = TRUE)"


# Step 5: Snapshot the package environment
echo "📸 Saving package versions to renv.lock..."
Rscript -e "renv::load(); renv::settings\$bioconductor.version('3.20'); bioc_repos <- BiocManager::repositories(version = '3.20'); options(repos = bioc_repos, renv.config.repos.override = bioc_repos); cache_path <- file.path(renv::paths\$root(), 'cache'); if (dir.exists(cache_path)) unlink(list.files(cache_path, pattern = 'PACKAGES', full.names = TRUE, recursive = TRUE)); renv::snapshot(confirm = FALSE)"

echo "✅ Pipeline complete!"
