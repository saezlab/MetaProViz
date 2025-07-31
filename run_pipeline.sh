#!/usr/bin/env bash

# Exit immediately on error
set -e

echo "📦 Checking and initializing renv..."

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

# Step 2: Render Rmd files in order
echo "📝 Rendering RMarkdown files..."

Rscript -e "renv::load(); rmarkdown::render('ExampleDataOverview/README_ExampleDataPlots.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_Intracellular/standard-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Cells_CoRe/core-metabolomics.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('Patients_Metadata/sample-metadata.Rmd')"
Rscript -e "renv::load(); rmarkdown::render('PriorKnowledge/prior-knowledge.Rmd')"

# Step 3: Run supplementary tables last
echo "📊 Creating supplementary tables..."
Rscript -e "renv::load(); rmarkdown::render('SupplementaryTables.Rmd')"

# Step 4: Snapshot the package environment
echo "📸 Saving package versions to renv.lock..."
Rscript -e "renv::snapshot(confirm = FALSE)"

echo "✅ Pipeline complete!"
