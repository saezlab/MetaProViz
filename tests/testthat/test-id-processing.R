test_that("id_processing default run returns staged QC object with automatic compatibility handling", {
  testthat::local_mocked_bindings(
    compare_pk = function(data, ...) {
      list(
        summary_table = data.frame(stage_metric = nrow(data[[1]]), stringsAsFactors = FALSE),
        upset_plot = ggplot2::ggplot()
      )
    },
    count_id = function(data, column, delimiter = ";", ...) {
      values <- data[[column]]
      split_pattern <- if (delimiter == ";") ";\\s*" else ",\\s*"
      entry_count <- vapply(values, function(cell) {
        if (is.na(cell) || cell == "") {
          0L
        } else {
          length(strsplit(as.character(cell), split_pattern)[[1]])
        }
      }, integer(1))
      id_label <- ifelse(entry_count == 0L, "No ID", ifelse(entry_count == 1L, "Single ID", "Multiple IDs"))
      list(
        Table = data.frame(entry_count = entry_count, id_label = id_label, stringsAsFactors = FALSE),
        Plot = ggplot2::ggplot(),
        Plot_Sized = ggplot2::ggplot()
      )
    },
    .package = "MetaProViz"
  )

  result <- id_processing(
    data = mpv_seed_input_data(),
    delimiter = ";",
    edge_table = mpv_seed_edge_table(),
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_true(all(c("Data", "Plot", "Workflow") %in% names(result)))
  expect_true(all(c("input", "after_compatibility", "id_count_summary", "compatibility") %in% names(result$Data)))
  expect_true("compatibility_check" %in% result$Workflow$steps_run)
  expect_true(is.list(result$Data$compatibility))
  expect_null(result$Workflow$compatibility_result)
  expect_identical(names(result$Plot), c("compare_pk_by_stage", "count_id_by_stage"))
  expect_identical(result$Data$after_compatibility$HMDB[[2]], "HMDB0000002")
  expect_identical(result$Data$after_compatibility$CHEBI[[4]], "CHEBI:3")
  expect_true(all(c("stage", "namespace", "delta_total_ids_vs_input") %in% colnames(result$Data$id_count_summary)))
})

test_that("id_processing rejects simultaneous translation and traversal", {
  expect_error(
    id_processing(
      data = mpv_seed_input_data(),
      run_translation = TRUE,
      translation_from = "KEGG",
      translation_to = "HMDB",
      run_traversal = TRUE,
      save_plot = NULL,
      save_table = NULL,
      print_plot = FALSE,
      verbose = FALSE
    ),
    "mutually exclusive"
  )
})

test_that("id_processing wraps translate_id and merges translated IDs back into wide metadata", {
  testthat::local_mocked_bindings(
    compare_pk = function(data, ...) {
      list(
        summary_table = data.frame(stage_metric = nrow(data[[1]]), stringsAsFactors = FALSE),
        upset_plot = ggplot2::ggplot()
      )
    },
    count_id = function(data, column, delimiter = ";", ...) {
      values <- data[[column]]
      split_pattern <- if (delimiter == ";") ";\\s*" else ",\\s*"
      entry_count <- vapply(values, function(cell) {
        if (is.na(cell) || cell == "") {
          0L
        } else {
          length(strsplit(as.character(cell), split_pattern)[[1]])
        }
      }, integer(1))
      id_label <- ifelse(entry_count == 0L, "No ID", ifelse(entry_count == 1L, "Single ID", "Multiple IDs"))
      list(
        Table = data.frame(entry_count = entry_count, id_label = id_label, stringsAsFactors = FALSE),
        Plot = ggplot2::ggplot(),
        Plot_Sized = ggplot2::ggplot()
      )
    },
    translate_id = function(data, metadata_info, from, to, summary = FALSE, save_table = NULL, path = NULL) {
      out <- data
      out$hmdb <- ifelse(out$workflow_input_id == "C00001", "HMDB0000001", "HMDB0000002")
      list(TranslatedDF = out)
    },
    .package = "MetaProViz"
  )

  input <- data.frame(
    Metabolite = c("m1", "m2"),
    KEGG = c("C00001", "C00002"),
    HMDB = c(NA, NA),
    stringsAsFactors = FALSE
  )

  result <- id_processing(
    data = input,
    id_types = c("HMDB", "KEGG"),
    delimiter = ";",
    run_compatibility_check = FALSE,
    run_translation = TRUE,
    translation_from = "KEGG",
    translation_to = "HMDB",
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE,
    verbose = FALSE
  )

  expect_true("translation" %in% result$Workflow$steps_run)
  expect_true(is.list(result$Data$translation))
  expect_null(result$Workflow$translation_result)
  expect_identical(result$Data$after_translation$HMDB[[1]], "HMDB0000001")
  expect_identical(result$Data$after_translation$HMDB[[2]], "HMDB0000002")
})

test_that("id_processing wraps traverse_ids and keeps traversal artifacts", {
  testthat::local_mocked_bindings(
    compare_pk = function(data, ...) {
      list(
        summary_table = data.frame(stage_metric = nrow(data[[1]]), stringsAsFactors = FALSE),
        upset_plot = ggplot2::ggplot()
      )
    },
    count_id = function(data, column, delimiter = ";", ...) {
      values <- data[[column]]
      split_pattern <- if (delimiter == ";") ";\\s*" else ",\\s*"
      entry_count <- vapply(values, function(cell) {
        if (is.na(cell) || cell == "") {
          0L
        } else {
          length(strsplit(as.character(cell), split_pattern)[[1]])
        }
      }, integer(1))
      id_label <- ifelse(entry_count == 0L, "No ID", ifelse(entry_count == 1L, "Single ID", "Multiple IDs"))
      list(
        Table = data.frame(entry_count = entry_count, id_label = id_label, stringsAsFactors = FALSE),
        Plot = ggplot2::ggplot(),
        Plot_Sized = ggplot2::ggplot()
      )
    },
    traverse_ids = function(data, id_types, delimiter, save_table = "csv", path = NULL, verbose = FALSE) {
      expanded <- data
      expanded$HMDB <- NA_character_
      expanded$HMDB_translated <- c("HMDB0000001", "HMDB0000002")
      list(
        ExpandedDF = expanded,
        ID_pair_compatibility = data.frame(row_id = 1:2, stringsAsFactors = FALSE),
        ID_Edges_prior_knowledge = data.frame(id1 = character(), type1 = character(), id2 = character(), type2 = character(), stringsAsFactors = FALSE)
      )
    },
    .package = "MetaProViz"
  )

  input <- data.frame(
    Metabolite = c("m1", "m2"),
    KEGG = c("C00001", "C00002"),
    stringsAsFactors = FALSE
  )

  result <- id_processing(
    data = input,
    id_types = c("HMDB", "KEGG"),
    delimiter = ";",
    run_compatibility_check = FALSE,
    run_traversal = TRUE,
    save_plot = NULL,
    save_table = NULL,
    print_plot = FALSE,
    verbose = FALSE
  )

  expect_true("traversal" %in% result$Workflow$steps_run)
  expect_true(is.list(result$Data$traversal))
  expect_null(result$Workflow$traversal_result)
  expect_identical(result$Data$after_traversal$HMDB[[2]], "HMDB0000002")
  expect_true(all(c("pair_compatibility", "prior_knowledge_edges") %in% names(result$Data$traversal)))
})
