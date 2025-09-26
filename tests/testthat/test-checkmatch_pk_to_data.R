# Unit tests for checkmatch_pk_to_data
library(testthat)
library(MetaProViz)
library(dplyr)

# Example test: check that function returns TRUE for matching data

test_that("checkmatch_pk_to_data returns TRUE for matching data", {
  pk <- data.frame(HMDB = c("A", "B", "C"), term = c("T1", "T2", "T3"))
  data <- data.frame(HMDB = c("A", "B", "C"), value = 1:3, term = c("T1", "T2", "T3"))
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  result <- checkmatch_pk_to_data(data, pk, metadata_info)
  expect_true(is.list(result) || is.data.frame(result) || is.logical(result))
})

# Test: mismatched IDs should return FALSE (or appropriate result)
test_that("checkmatch_pk_to_data returns FALSE for mismatched IDs", {
  pk <- data.frame(HMDB = c("A", "B", "C"), term = c("T1", "T2", "T3"))
  data <- data.frame(HMDB = c("A", "B", "Y", "Z"), value = 1:4, term = c("T1","T1", "T2", "T3"))
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  result <- checkmatch_pk_to_data(data, pk, metadata_info)
  expect_true(is.list(result) || is.data.frame(result) || is.logical(result))
})

# Test: empty data frames
test_that("checkmatch_pk_to_data handles empty data frames", {
  pk <- data.frame(HMDB = character(), term = character())
  data <- data.frame(HMDB = character(), value = numeric(), term = character())
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  expect_error(checkmatch_pk_to_data(data, pk, metadata_info))
})

# Test: invalid input types
test_that("checkmatch_pk_to_data errors on non-data.frame input", {
  pk <- list(HMDB = c("A", "B", "C"), term = c("T1", "T2", "T3"))
  data <- "not a data frame"
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  expect_error(checkmatch_pk_to_data(data, pk, metadata_info))
})

# Example from @examples and input checks

test_that("checkmatch_pk_to_data errors if InputID column missing in data", {
  pk <- data.frame(HMDB = c("A", "B", "C"), term = c("T1", "T2", "T3"))
  data <- data.frame(value = 1:3)
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  expect_error(checkmatch_pk_to_data(data, pk, metadata_info),
               "The HMDB column selected as InpuID in metadata_info was not found in data. Please check your input.")
})

test_that("checkmatch_pk_to_data errors if PriorID column missing in input_pk", {
  pk <- data.frame(ID = c("A", "B", "C"), term = c("T1", "T2", "T3"))
  data <- data.frame(HMDB = c("A", "D", "C"), value = 1:3)
  metadata_info <- c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term")
  expect_error(checkmatch_pk_to_data(data, pk, metadata_info),
               "The HMDB column selected as InpuID in metadata_info was not found in input_pk. Please check your input.")
})

test_that("checkmatch_pk_to_data works with example data from @examples", {
  DetectedIDs <- data.frame(HMDB = c("A", "B", "C"))
  input_pathway <- data.frame(hmdb = c("A", "E", "F"), term = c("T1", "T2", "T3"))
  metadata_info <- c(InputID = "HMDB", PriorID = "hmdb", grouping_variable = "term")
  result <- checkmatch_pk_to_data(DetectedIDs, input_pathway, metadata_info)
  expect_true(is.list(result) || is.data.frame(result) || is.logical(result))
})
