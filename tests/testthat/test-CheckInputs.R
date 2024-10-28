test_that("InputChecks", {
  d <- ToyData('IntraCells_Raw')
  expect_error(CheckInput(d), 'needs to be of class numeric')

  d_meta <-d[,c(1:3)]%>%
    tibble::remove_rownames()
  expect_error(CheckInput(d, InputData_Num=FALSE, SettingsFile_Sample=d_meta), 'row.names InputData need to match row.names SettingsFile_Sample.')
})
