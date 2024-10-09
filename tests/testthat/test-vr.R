context("vietoris_rips")

test_that("basic vietoris-rips works", {
  # setup vals (including reproducibility)
  set.seed(42)
  INPUT_SIZE <- 10
  
  MIN_DIM <- 2
  MAX_DIM <- 4
  
  # check dimensions 2 through 4
  for (curr_dim in MIN_DIM:MAX_DIM) {
    curr_data <- runif(50 * curr_dim)
    dim(curr_data) <- c(50, curr_dim)
    
    curr_vr <- vietoris_rips(curr_data, max_dim = curr_dim - 1)
    
    # check df dimensions
    expect_true(nrow(curr_vr) > 0)
    expect_equal(3, ncol(curr_vr))
    
    # make sure there's at least 1 feature from each dimension
    expected <- 0:(curr_dim - 1)
    actual <- sort(unique(curr_vr$dimension))
    expect_equal(expected, actual)
  }
  
})

# generate dataset (set seed for reproducibility, altho new one would be fine)
set.seed(42)
angles <- runif(25, min = 0, max = 2 * pi)
circle_mat <- cbind(cos(angles), sin(angles))
circle_df <- as.data.frame(circle_mat)
circle_dist <- dist(circle_mat)

test_that("consistency across generic methods for tidy data", {
  # calculate persistent homology for each class
  mat_phom <- vietoris_rips(circle_mat)
  df_phom <- vietoris_rips(circle_df)
  dist_phom <- vietoris_rips(circle_dist)
  
  # compare persistent homology across classes
  expect_equal(mat_phom, df_phom)
  expect_equal(mat_phom, dist_phom)
})

test_that("`dim` deprecation warns and replaces", {
  # use data above, print warnings
  expect_warning(vietoris_rips(circle_mat, dim = 1L), "max_dim")
  expect_warning(vietoris_rips(circle_dist, dim = 1L), "max_dim")
  # use data above, prefer `max_dim`
  circle_vr <- vietoris_rips(circle_dist, dim = 1L, max_dim = 0L)
  expect_lt(max(circle_vr$dimension), 1L)
})

test_that("consistency across generic methods for time series", {
  # generate dataset (set seed for reproducibility, altho new one would be fine)
  set.seed(42)
  val_num <- runif(50, min = 0, max = 1)
  val_ts <- ts(val_num)
  
  # calculate persistent homology
  num_phom <- vietoris_rips(val_num)
  ts_phom <- vietoris_rips(val_ts)
  
  # compare persistent homology across classes
  expect_equal(num_phom, ts_phom)
})
