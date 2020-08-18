context("2-dimensional cubical complexes")
library("ripserr")

test_that("basic function works", {
  # reproducibility
  set.seed(42)
  
  # create data
  test_data <- matrix(rnorm(100, mean = 1000, sd = 150),
                      nrow = 10, ncol = 10)
  
  # create cubical complex
  cub_comp <- ripserr::cubical(test_data)
  
  # test cubical complex
  expect_equal(ncol(cub_comp), 3)
  expect_equal(nrow(cub_comp), 28)
  
  counts <- base::table(cub_comp[, 1])
  names(counts) <- NULL
  counts <- as.numeric(counts)
  
  expect_equal(counts[1], 1)
  expect_equal(counts[2], 20)
  expect_equal(counts[3], 7)
})
