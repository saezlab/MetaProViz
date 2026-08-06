test_that("checkmatch_pk_to_data returns detailed match summaries", {
  pk <- data.frame(
    HMDB = c("A", "B", "C"),
    term = c("T1", "T2", "T3"),
    stringsAsFactors = FALSE
  )
  detected <- data.frame(
    HMDB = c("A", "B", "Y", NA),
    value = 1:4,
    stringsAsFactors = FALSE
  )

  result <- checkmatch_pk_to_data(
    data = detected,
    input_pk = pk,
    metadata_info = c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term"),
    save_table = NULL
  )

  expect_type(result, "list")
  expect_true(any(grepl("summary", names(result), ignore.case = TRUE)))
  expect_true(any(vapply(result, is.data.frame, logical(1))))
})

test_that("checkmatch_pk_to_data validates required columns", {
  pk <- data.frame(HMDB = c("A", "B"), term = c("T1", "T2"))

  expect_error(
    checkmatch_pk_to_data(
      data = data.frame(value = 1:2),
      input_pk = pk,
      metadata_info = c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term"),
      save_table = NULL
    ),
    "selected as InpuID in metadata_info was not found in data"
  )

  expect_error(
    checkmatch_pk_to_data(
      data = data.frame(HMDB = c("A", "B")),
      input_pk = data.frame(ID = c("A", "B"), term = c("T1", "T2")),
      metadata_info = c(InputID = "HMDB", PriorID = "HMDB", grouping_variable = "term"),
      save_table = NULL
    ),
    "selected as InpuID in metadata_info was not found in input_pk"
  )
})
