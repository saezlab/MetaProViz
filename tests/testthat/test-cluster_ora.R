test_that("cluster_ora returns ORA tables for a small clustered dataset", {
  result <- cluster_ora(
    data = mpv_cluster_input_fixture(),
    metadata_info = c(
      ClusterColumn = "cluster",
      BackgroundColumn = "BG_method",
      PathwayTerm = "term",
      PathwayFeature = "Metabolite"
    ),
    remove_background = FALSE,
    input_pathway = mpv_simple_pathway_fixture(),
    min_gssize = 1,
    max_gssize = 10,
    save_table = NULL
  )

  expect_type(result, "list")
  expect_named(result, c("DF", "ClusterGo"))
  expect_true(all(vapply(result$DF, is.data.frame, logical(1))))
})

test_that("standard_ora returns the selected inputs and enrichment summary", {
  result <- standard_ora(
    data = mpv_standard_ora_input_fixture(),
    metadata_info = c(
      pvalColumn = "p.adj",
      percentageColumn = "t.val",
      PathwayTerm = "term",
      PathwayFeature = "Metabolite"
    ),
    cutoff_stat = 0.05,
    cutoff_percentage = 40,
    input_pathway = mpv_simple_pathway_fixture(),
    min_gssize = 1,
    max_gssize = 10,
    save_table = NULL
  )

  expect_type(result, "list")
  expect_true(all(c("InputSelection", "ClusterGosummary", "ClusterGo") %in% names(result)))
  expect_s3_class(result$ClusterGosummary, "data.frame")
})
