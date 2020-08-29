context("cubical 2-dim")
library("ripserr")

# setup vars
INPUT_SIZE <- 10
DIM <- 2
set.seed(42)
test_data <- rnorm(INPUT_SIZE ^ DIM)
dim(test_data) <- rep(INPUT_SIZE, DIM)

test_that("basic 2-dim cubical works", {
  # create cubical complex
  cub_comp <- cubical(test_data)
  
  # test cubical complex frequency/counts
  expect_equal(ncol(cub_comp), 3)
  expect_true(nrow(cub_comp) > 0)
  
  counts <- table(cub_comp$dimension)
  names(counts) <- NULL
  counts <- as.numeric(counts)
  
  # at least 1 feature from each dimension
  expect_true(counts[1] > 0)
  expect_true(counts[2] > 0)
  expect_true(counts[3] > 0)
  
  # make sure no births after deaths
  expect_equal(0, sum(cub_comp$birth > cub_comp$death))
  
  # can return a matrix or a data frame (both equivalent)
  expect_equal(cub_comp, 
               as.data.frame(cubical(test_data, return_format = "mat")))
})

# these tests use example data + original code from Github: CubicalRipser/Cubical_2dim
#   to validate accuracy
test_that("2-dim cubical returns same values as validated tests", {
  
  # read validated input and output data
  input_data <- readRDS("input_2dim.rds")
  output_data <- readRDS("output_2dim.rds")
  
  # re-calculate output w/ ripserr
  THRESH <- 9999
  test_output <- cubical(input_data, threshold = THRESH)
  
  # filter out threshold value features to avoid spurious differences in equality
  output_data <- subset(output_data, death < THRESH)
  test_output <- subset(test_output, death < THRESH)
  
  # ensure no NAs
  expect_equal(0, sum(is.na(output_data)))
  expect_equal(0, sum(is.na(test_output)))
  
  # make sure # of features is close enough
  expect_equal(nrow(test_output), nrow(output_data), tolerance = 5)
  
  # check means of births and deaths to ensure close enough
  expect_equal(mean(test_output$birth), mean(output_data$birth), tolerance = 0.025)
  expect_equal(mean(test_output$death), mean(output_data$death), tolerance = 0.025)
})