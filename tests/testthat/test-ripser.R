context("vietoris-rips")

test_that("basic vietoris-rips works", {
  # setup vals (including reproducibility)
  set.seed(42)
  INPUT_SIZE <- 10
  
  # check dimensions 2 through 5
  for (curr_dim in 2:5) {
    curr_data <- runif(50 * curr_dim)
    dim(curr_data) <- c(50, curr_dim)
    
    curr_vr <- vietoris_rips(curr_data, dim = curr_dim - 1)
    
    # check df dimensions
    expect_true(nrow(curr_vr) > 0)
    expect_equal(3, ncol(curr_vr))
    
    # make sure there's at least 1 feature from each dimension
    expected <- 0:(curr_dim - 1)
    actual <- sort(unique(curr_vr$dimension))
    expect_equal(expected, actual)
  }
  
})