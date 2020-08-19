context("cubical")
library("ripserr")

test_that("basic 2-dim cubical works", {
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
  
  expect_equal(0, sum(cub_comp[, 2] > cub_comp[, 3]))
})

test_that("basic 4-dim cubical works", {
  # reproducibility
  set.seed(42)
  
  # create data
  test_data <- rnorm(5 ^ 4)
  dim(test_data) <- rep(5, 4)
  
  # create cubical complex
  cub_comp <- ripserr::cubical(test_data)
  
  # test cubical complex
  expect_equal(ncol(cub_comp), 3)
  expect_equal(nrow(cub_comp), 1592)
  
  counts <- base::table(cub_comp[, 1])
  names(counts) <- NULL
  counts <- as.numeric(counts)
  
  expect_equal(counts[1], 870)
  expect_equal(counts[2], 74)
  expect_equal(counts[3], 72)
  expect_equal(counts[4], 573)
  expect_equal(counts[5], 3)
  
  # expect_equal(0, sum(cub_comp[, 2] > cub_comp[, 3]))
})