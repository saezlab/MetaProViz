#!/usr/bin/env bash

# Exit immediately on error
set -e

echo "📦 Checking and initializing renv..."

# Step 0: Do you have R and sepcify which R to use?
# Define R installation path
R_PATH="/c/Program Files/R/R-4.x.x/bin/Rscript.exe"


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
    "magrittr", "dplyr", "tibble", "rlang", "ggfortify",
    "stringr", "tidyverse", "tidyr", "purrr", "rmarkdown",
    "devtools", "easyalluvial"
)

# Check if lockfile exists
if (!file.exists("renv.lock")) {
    message("🔧 No lockfile found. Initializing new environment...")
    renv::init(force = TRUE)

    # Install all CRAN packages at once
    message("📥 Installing required packages...")
    install.packages(packages)

    # Install GitHub packages
    message("📥 Installing MetaProViz from GitHub...")
    devtools::install_github("saezlab/MetaProViz")

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
EOF

# Step 2: Render Rmd files in order
echo "📝 Rendering RMarkdown files..."

Rscript -e "renv::load(); rmarkdown::render('ExampleDataOverview/README_ExampleDataPlots.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_Intracellular/standard-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_CoRe/core-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Patients_Metadata/sample-metadata.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('PriorKnowledge/prior-knowledge.Rmd')"

# Step 3. Replace old dates in SupplementaryTables.Rmd with today's date
echo "🛠️ Replacing dated filenames in SupplementaryTables.Rmd..."
TODAY=$(date +%Y-%m-%d)

# Replace any date pattern like 20XX-XX-XX with today's date in file paths
sed -i.bak -E "s/([A-Za-z0-9_/.-]*)[0-9]{4}-[0-9]{2}-[0-9]{2}([A-Za-z0-9_/.-]*)/\1$TODAY\2/g" SupplementaryTables.Rmd

# Step 4: Run supplementary tables last
echo "📊 Creating supplementary tables..."
Rscript -e "renv::load(); rmarkdown::render('SupplementaryTables.Rmd')"

# Step 5: Snapshot the package environment
echo "📸 Saving package versions to renv.lock..."
Rscript -e "renv::snapshot(confirm = FALSE)"

echo "✅ Pipeline complete!"
