library(testthat)
library(MetaProViz)
library(tibble)

make_edge_table <- function() {
  tibble(
    id1 = c(
      "HMDB0000001", "C00001",
      "C00001", "1",
      "HMDB0000002", "C00002"
    ),
    type1 = c(
      "HMDB", "KEGG",
      "KEGG", "CHEBI",
      "HMDB", "KEGG"
    ),
    id2 = c(
      "C00001", "HMDB0000001",
      "1", "C00001",
      "C00002", "HMDB0000002"
    ),
    type2 = c(
      "KEGG", "HMDB",
      "CHEBI", "KEGG",
      "KEGG", "HMDB"
    )
  )
}

make_input_data <- function() {
  data.frame(
    metabolite = c("full", "partial", "complete_hmdb", "complete_chebi"),
    HMDB = c("HMDB1", "HMDB2; HMDB3", "HMDB4; HMDB5", NA),
    KEGG = c("C00001", "C00002", NA, "C00003"),
    CHEBI = c("CHEBI:1", NA, NA, "CHEBI:3; CHEBI:4"),
    PUBCHEM = NA_character_,
    stringsAsFactors = FALSE
  )
}

test_that("seed_id_compatibility_check preserves raw compatibility behavior by default", {
  res <- seed_id_compatibility_check(
    data = make_input_data(),
    edge_table = make_edge_table(),
    delimiter = ";"
  )

  expect_true(all(c(
    "ID_pair_compatibility",
    "data_with_compatibility",
    "feature_compatibility_summary",
    "data_after_handling",
    "ID_pair_compatibility_after_handling"
  ) %in% names(res)))
  expect_false("handling_summary_text" %in% names(res))
  expect_false("handling_summary_metrics" %in% names(res))

  classes <- res$feature_compatibility_summary$compatibility_class
  expect_identical(
    classes,
    c("fully_compatible", "partially_compatible", "completely_incompatible", "completely_incompatible")
  )

  expect_identical(res$data_after_handling$HMDB[[2]], "HMDB2; HMDB3")
  expect_identical(res$data_after_handling$CHEBI[[4]], "CHEBI:3; CHEBI:4")
})

test_that("seed_id_compatibility_check handles partially compatible features", {
  res <- seed_id_compatibility_check(
    data = make_input_data(),
    edge_table = make_edge_table(),
    delimiter = ";",
    handle_partially_compatible = TRUE
  )

  expect_true("handling_summary_text" %in% names(res))
  expect_true("handling_summary_metrics" %in% names(res))

  expect_identical(res$data_after_handling$HMDB[[2]], "HMDB0000002")
  expect_identical(res$data_after_handling$KEGG[[2]], "C00002")
  expect_true(is.na(res$data_after_handling$CHEBI[[2]]))

  partial_pairs <- subset(res$ID_pair_compatibility_after_handling, row_id == 2)
  expect_equal(sum(partial_pairs$pair_compatible %in% TRUE, na.rm = TRUE), 1)
  expect_equal(sum(partial_pairs$pair_compatible %in% FALSE, na.rm = TRUE), 0)
})

test_that("seed_id_compatibility_check handles completely incompatible features with default priority", {
  res <- seed_id_compatibility_check(
    data = make_input_data(),
    edge_table = make_edge_table(),
    delimiter = ";",
    handle_completely_incompatible = TRUE
  )

  expect_identical(res$data_after_handling$HMDB[[3]], "HMDB0000004")
  expect_true(is.na(res$data_after_handling$KEGG[[3]]))
  expect_identical(res$data_after_handling$CHEBI[[4]], "CHEBI:3")
  expect_true(is.na(res$data_after_handling$KEGG[[4]]))

  summary_metrics <- res$handling_summary_metrics
  expect_true(summary_metrics$user_choices$handle_completely_incompatible)
  expect_true(grepl("priority=HMDB > CHEBI > PUBCHEM > KEGG", res$handling_summary_text, fixed = TRUE))
})

test_that("seed_id_compatibility_check respects custom complete handling priority", {
  res <- seed_id_compatibility_check(
    data = make_input_data(),
    edge_table = make_edge_table(),
    delimiter = ";",
    handle_completely_incompatible = TRUE,
    completely_incompatible_priority = c("KEGG", "CHEBI", "PUBCHEM", "HMDB")
  )

  expect_identical(res$data_after_handling$KEGG[[4]], "C00003")
  expect_true(is.na(res$data_after_handling$CHEBI[[4]]))
})
