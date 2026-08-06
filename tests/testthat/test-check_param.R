test_that("check_param validates numeric data and aligned sample metadata", {
  fx <- mpv_intracell_fixture()
  non_numeric <- fx$data
  non_numeric[[1]] <- as.character(non_numeric[[1]])

  expect_error(
    MetaProViz:::check_param(non_numeric),
    "needs to be of class numeric"
  )

  bad_metadata <- tibble::remove_rownames(fx$metadata_sample)
  expect_error(
    MetaProViz:::check_param(
      fx$data,
      data_num = FALSE,
      metadata_sample = bad_metadata
    ),
    "row.names data need to match row.names metadata_sample"
  )
})

test_that("internal delimiter and ID normalization helpers behave deterministically", {
  expect_identical(MetaProViz:::normalize_delimiter(";"), ";")
  expect_identical(MetaProViz:::normalize_delimiter("comma"), ",")
  expect_identical(MetaProViz:::delimiter_to_pattern(";"), ";\\s*")
  expect_identical(MetaProViz:::normalize_hmdb(c("1544", "HMDB0000190")), c("HMDB0001544", "HMDB0000190"))
  expect_identical(MetaProViz:::split_ids("A; B;C", type = "HMDB"), c("A", "B", "C"))
})
