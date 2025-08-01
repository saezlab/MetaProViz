#!/usr/bin/env bash

# Exit immediately on error
set -e

echo "📦 Checking and initializing renv..."

# Step 0: Do you have R and sepcify which R to use?
# Define R installation path
R_PATH="/c/Program Files/R/R-4.x.x/bin/Rscript.exe"


# Step 1: Install or load renv in R
Rscript - <<EOF
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!file.exists("renv.lock")) {
  message("🔧 Initializing renv for the first time...")
  renv::init(bare = TRUE)
} else {
  renv::load()
}
EOF

# Step 2: Install MetaProViz R
echo "📝 MetaProViz install/loading..."

Rscript - <<EOF
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("MetaProViz", quietly = TRUE)) devtools::install_github("saezlab/MetaProViz")
EOF

# Step 2: Install Dependencies needed to run the analysis.rmd files
echo "📝 Script dependency install/loading..."

Rscript - <<EOF
# Install required CRAN packages if not present
packages <- c("magrittr", "dplyr", "tibble", "rlang", "ggfortify",
              "stringr", "tidyverse", "tidyr", "purrr", "easyalluvial")

for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

# Load all packages to verify installation
lapply(packages, library, character.only = TRUE)
EOF

# Step 4: Render Rmd files in order
echo "📝 Rendering RMarkdown files..."

Rscript - <<EOF
if (!requireNamespace("rmarkdown", quietly = TRUE)) install.packages("rmarkdown")
EOF

Rscript -e "renv::load(); rmarkdown::render('ExampleDataOverview/README_ExampleDataPlots.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_Intracellular/standard-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_CoRe/core-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Patients_Metadata/sample-metadata.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('PriorKnowledge/prior-knowledge.Rmd')"

# 5. Replace old dates in SupplementaryTables.Rmd with today's date
echo "🛠️ Replacing dated filenames in SupplementaryTables.Rmd..."
TODAY=$(date +%Y-%m-%d)

# Replace any date pattern like 20XX-XX-XX with today's date in file paths
sed -i.bak -E "s/([A-Za-z0-9_/.-]*)[0-9]{4}-[0-9]{2}-[0-9]{2}([A-Za-z0-9_/.-]*)/\1$TODAY\2/g" SupplementaryTables.Rmd

# Step 6: Run supplementary tables last
echo "📊 Creating supplementary tables..."
Rscript -e "renv::load(); rmarkdown::render('SupplementaryTables.Rmd')"

# Step 7: Snapshot the package environment
echo "📸 Saving package versions to renv.lock..."
Rscript -e "renv::snapshot(confirm = FALSE)"

echo "✅ Pipeline complete!"
