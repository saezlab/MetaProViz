test_that("viz_pca and viz_heatmap return plot containers", {
  fx <- mpv_intracell_fixture(10)

  pca_res <- viz_pca(
    data = fx$data,
    metadata_sample = fx$metadata_sample,
    metadata_info = c(color = "Conditions"),
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_type(pca_res, "list")
  expect_true(all(c("Plot", "Plot_Sized") %in% names(pca_res)))

  heatmap_res <- viz_heatmap(
    data = fx$data[, 1:6, drop = FALSE],
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_type(heatmap_res, "list")
  expect_true(all(c("Plot", "Plot_Sized") %in% names(heatmap_res)))
})

test_that("viz_superplot returns boxplot output for one metabolite", {
  fx <- mpv_intracell_fixture(10)

  result <- viz_superplot(
    data = fx$data[, 1, drop = FALSE],
    metadata_sample = fx$metadata_sample,
    metadata_info = c(Conditions = "Conditions"),
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("Plot", "Plot_Sized") %in% names(result)))
})

test_that("viz_graph builds a graph plot from a similarity matrix", {
  sim <- matrix(
    c(1, 0.8, 0.1, 0.8, 1, 0.2, 0.1, 0.2, 1),
    nrow = 3,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  result <- viz_graph(
    similarity_matrix = sim,
    clusters = c(A = "cluster1", B = "cluster1", C = "cluster2"),
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_true(inherits(result, "ggplot") || is.null(result))
})

test_that("viz_volcano returns plot objects from bundled DMA results", {
  dma_df <- mpv_intracell_dma_fixture(12)
  rownames(dma_df) <- dma_df$Metabolite

  result <- viz_volcano(
    data = dma_df,
    metadata_info = c(pval = "p.adj", log2fc = "Log2FC", feature_name = "Metabolite"),
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("Plot", "Plot_Sized") %in% names(result)))
})
