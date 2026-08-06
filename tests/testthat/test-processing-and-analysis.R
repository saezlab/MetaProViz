test_that("processing returns DF and Plot outputs on bundled data", {
  fx <- mpv_intracell_fixture()

  result <- processing(
    data = fx$data,
    metadata_sample = fx$metadata_sample,
    metadata_info = fx$metadata_info[c("Conditions", "Biological_Replicates")],
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("DF", "Plot") %in% names(result)))
  expect_true(is.list(result$DF))
})

test_that("replicate_sum aggregates analytical replicates", {
  fx <- mpv_intracell_fixture(8)
  input_data <- fx$data[c(1, 1, 2, 2), 1:4, drop = FALSE]
  rownames(input_data) <- paste0("S", seq_len(nrow(input_data)))
  metadata <- fx$metadata_sample[c(1, 1, 2, 2), , drop = FALSE]
  rownames(metadata) <- rownames(input_data)
  metadata$Analytical_Replicates <- c(1, 2, 1, 2)

  result <- replicate_sum(
    data = input_data,
    metadata_sample = metadata,
    metadata_info = c(
      Conditions = "Conditions",
      Biological_Replicates = "Biological_Replicates",
      Analytical_Replicates = "Analytical_Replicates"
    ),
    save_table = NULL
  )

  expect_s3_class(result, "data.frame")
  expect_true("n_AnalyticalReplicates_Summed" %in% colnames(result))
})

test_that("pool_estimation returns variability outputs for pool samples", {
  se <- mpv_intracell_se_fixture()

  result <- pool_estimation(
    data = se,
    metadata_info = c(PoolSamples = "Pool", Conditions = "Conditions"),
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("DF", "Plot") %in% names(result)))
})

test_that("dma returns differential analysis tables on a small subset", {
  fx <- mpv_intracell_fixture(10)

  result <- dma(
    data = fx$data,
    metadata_sample = fx$metadata_sample,
    metadata_info = c(Conditions = "Conditions"),
    shapiro = FALSE,
    bartlett = FALSE,
    vst = FALSE,
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(any(vapply(result, is.list, logical(1))))
})

test_that("metadata_analysis summarizes PCA and ANOVA outputs", {
  fx <- mpv_intracell_fixture(10)

  result <- metadata_analysis(
    data = fx$data,
    metadata_sample = fx$metadata_sample,
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("res_prcomp", "res_loadings", "res_topBottomFeatures", "res_aov", "res_summary") %in% names(result)))
})

test_that("meta_pk returns a metadata-derived prior knowledge table", {
  fx <- mpv_intracell_fixture(8)

  result <- meta_pk(
    data = fx$data,
    metadata_sample = fx$metadata_sample,
    metadata_info = c("Conditions"),
    save_table = NULL
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("SampleID", "term", "mor") %in% colnames(result)))
})

test_that("mca_2cond and mca_core return summary/result tables", {
  dma_df <- mpv_intracell_dma_fixture(12)
  cond1 <- dma_df[1:8, c("Metabolite", "Log2FC", "p.adj"), drop = FALSE]
  cond2 <- dma_df[5:12, c("Metabolite", "Log2FC", "p.adj"), drop = FALSE]

  result_2cond <- mca_2cond(
    data_c1 = cond1,
    data_c2 = cond2,
    save_table = NULL
  )

  expect_type(result_2cond, "list")
  expect_true(all(c("MCA_2Cond_summary", "MCA_2Cond_Results") %in% names(result_2cond)))

  result_core <- mca_core(
    data_intra = dma_df,
    data_core = mpv_core_dma_fixture(12),
    save_table = NULL
  )

  expect_type(result_core, "list")
  expect_true(any(grepl("summary", names(result_core), ignore.case = TRUE)))
})
