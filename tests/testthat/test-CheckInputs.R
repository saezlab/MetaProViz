test_that("multiplication works", {
  d <- ToyData('IntraCells_Raw')
  expect_error(CheckInput(d), 'needs to be of class numeric')
})
