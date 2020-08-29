context("cubical error-checks")

# make sure error checking is in place
# appropriately increases code coverage
test_that("errors are detected appropriately by cubical", {
  test_data <- rnorm(10 ^ 2)
  dim(test_data) <- rep(10, 2)
  
  # non-numeric threshold
  expect_error(cubical(test_data, threshold = "9999"))
  
  # dataset of type list
  expect_error(cubical(list(test_data)))
  
  # non-numeric matrix dataset
  expect_error(cubical(as.character(test_data)))
  
  # invalid numeric method
  expect_error(cubical(test_data, method = 2))
  
  # invalid method class
  expect_error(cubical(test_data, method = "0"))
  
  # invalid standardize class
  expect_error(cubical(test_data, standardize = "TRUE"))
  
  skip_on_cran()
  
  # too large dataset (2-dim)
  test_data_large <- rnorm(1500 ^ 2)
  dim(test_data_large) <- rep(1500, 2)
  expect_error(cubical(test_data_large))
  
  # too small dataset (2-dim)
  test_data_small <- numeric()
  dim(test_data_small) <- c(0, 0)
  expect_error(cubical(test_data_small))
  
  # too large dataset (3-dim)
  test_data_large <- rnorm(515 * 10 * 10)
  dim(test_data_large) <- c(515, 10, 10)
  expect_error(cubical(test_data_large))
  
  # too small dataset (3-dim)
  test_data_small <- numeric()
  dim(test_data_small) <- c(0, 0)
  expect_error(cubical(test_data_small))
  
  # too large dataset (4-dim)
  test_data_large <- rnorm(75 * 10 * 10 * 10)
  dim(test_data_large) <- c(75, 10, 10, 10)
  expect_error(cubical(test_data_large))
  
  # too small dataset (4-dim)
  test_data_small <- numeric()
  dim(test_data_small) <- c(0, 0)
  expect_error(cubical(test_data_small))
})