test_that("InputChecks", {
  d <- intracell_raw%>%tibble::column_to_rownames("Code")
  expect_error(MetaProViz:::check_param(d), 'needs to be of class numeric')

  d_meta <-d[,c(1:3)]%>%
    tibble::remove_rownames()
  expect_error(MetaProViz:::check_param(d, data_num=FALSE, metadata_sample=d_meta), 'row.names InputData need to match row.names SettingsFile_Sample.')
})
