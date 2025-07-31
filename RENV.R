# RENV.R

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

if (!file.exists("renv.lock")) {
  message("Initializing renv for the first time...")
  renv::init(bare = TRUE)  # bare = don't snapshot yet
} else {
  message("Loading existing renv environment...")
  renv::load()
}
