library(MetaProViz)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "MetaProViz")

# Outputs
expected_dir <- system.file("testdata", "outputs", "ora", package = "MetaProViz")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_ora function -----------------------------------------------
test_that("test Preprocessing", {

})




