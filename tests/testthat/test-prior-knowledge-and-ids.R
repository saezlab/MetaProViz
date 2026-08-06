test_that("translate_id reshapes translated columns using mocked ID translation", {
  testthat::local_mocked_bindings(
    translate_ids = mpv_mock_translate_ids,
    id_types = function() data.frame(in_a = "hmdb", in_b = "kegg", in_c = "chebi", in_d = "pubchem"),
    .package = "MetaProViz"
  )

  result <- translate_id(
    data = mpv_translation_fixture(),
    from = "kegg",
    to = c("hmdb", "chebi"),
    summary = FALSE,
    save_table = NULL
  )

  expect_type(result, "list")
  expect_true("TranslatedDF" %in% names(result))
  expect_true(all(c("hmdb", "chebi") %in% colnames(result$TranslatedDF)))
})

test_that("equivalent_id adds alternative IDs using mocked translation", {
  testthat::local_mocked_bindings(
    translate_ids = mpv_mock_translate_ids,
    id_types = function() data.frame(in_a = "hmdb", in_b = "kegg", in_c = "chebi", in_d = "pubchem"),
    count_id = function(...) list(result = data.frame(), plot = ggplot2::ggplot()),
    .package = "MetaProViz"
  )

  result <- equivalent_id(
    data = data.frame(MetaboliteID = c("HMDB0000001", "HMDB0000002"), stringsAsFactors = FALSE),
    metadata_info = c(InputID = "MetaboliteID"),
    from = "hmdb",
    save_table = NULL
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("PotentialAdditionalIDs", "AllIDs") %in% colnames(result)))
})

test_that("mapping_ambiguity validates inputs and seed compatibility returns full output", {
  expect_error(
    mapping_ambiguity(
      data = data.frame(x = "A"),
      from = "missing_from",
      to = "x",
      save_table = NULL
    ),
    "column was not found in data"
  )

  result <- seed_id_compatibility_check(
    data = mpv_seed_input_data(),
    edge_table = mpv_seed_edge_table(),
    delimiter = ";",
    handle_partially_compatible = TRUE,
    handle_completely_incompatible = TRUE
  )

  expect_true(all(c(
    "ID_pair_compatibility",
    "data_with_compatibility",
    "feature_compatibility_summary",
    "data_after_handling",
    "ID_pair_compatibility_after_handling",
    "handling_summary_text",
    "handling_summary_metrics"
  ) %in% names(result)))
})

test_that("traverse_ids expands rows using a mocked edge table builder", {
  testthat::local_mocked_bindings(
    build_id_edges_bidirectional = function(selected_types, verbose = FALSE) {
      mpv_seed_edge_table()
    },
    .package = "MetaProViz"
  )

  input <- data.frame(
    HMDB = c("HMDB0000001", NA),
    KEGG = c(NA, "C00001"),
    CHEBI = c(NA, "CHEBI:1"),
    PUBCHEM = c(NA, NA),
    stringsAsFactors = FALSE
  )

  result <- traverse_ids(
    data = input,
    delimiter = ";",
    save_table = NULL
  )

  expect_type(result, "list")
  expect_true(all(c("ExpandedDF", "ID_pair_compatibility", "ID_Edges_prior_knowledge") %in% names(result)))
})

test_that("cluster_pk, compare_pk and count_id return structured outputs", {
  cluster_input <- data.frame(
    term = c("T1", "T1", "T2", "T2", "T3"),
    MetaboliteID = c("HMDB0000001", "HMDB0000002", "HMDB0000002", "HMDB0000003", "HMDB0000004"),
    stringsAsFactors = FALSE
  )

  cluster_result <- cluster_pk(
    data = cluster_input,
    metadata_info = c(metabolite_column = "MetaboliteID", pathway_column = "term"),
    threshold = 0.1,
    save_plot = NULL,
    print_plot = FALSE
  )

  expect_type(cluster_result, "list")
  expect_true(all(c("data", "cluster_summary", "clusters", "similarity_matrix", "distance_matrix", "graph_plot") %in% names(cluster_result)))

  compare_result <- compare_pk(
    data = mpv_compare_pk_fixture(),
    metadata_info = list(KEGG = "MetaboliteID", Reactome = "MetaboliteID"),
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(compare_result, "list")
  expect_true(all(c("summary_table", "upset_plot") %in% names(compare_result)))

  count_result <- count_id(
    data = data.frame(HMDB = c(NA, "HMDB1", "HMDB1, HMDB2"), stringsAsFactors = FALSE),
    column = "HMDB",
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE
  )

  expect_type(count_result, "list")
  expect_true(all(c("Table", "Plot", "Plot_Sized") %in% names(count_result)))
})

test_that("get_exclusion_metabolites and make_gene_metab_set return expected structures", {
  testthat::local_mocked_bindings(
    translate_ids = mpv_mock_translate_ids,
    .package = "MetaProViz"
  )

  exclusion <- get_exclusion_metabolites()
  expect_s3_class(exclusion, "data.frame")
  expect_true(all(c("metabolite_id", "class", "id_type") %in% colnames(exclusion)))

  hallmarks <- mpv_load_dataset("hallmarks")
  result <- make_gene_metab_set(
    input_pk = hallmarks[seq_len(min(10, nrow(hallmarks))), , drop = FALSE],
    save_table = NULL
  )

  expect_type(result, "list")
  expect_true(all(c("GeneMetabSet", "MetabSet") %in% names(result)))
})
