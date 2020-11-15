context("PHom S3 class")

test_that("S3 class functionality", {
  # create data frame with sample persistence data and convert to PHom
  df <- data.frame(dimension = c(0, 0, 1, 1, 1, 2),
                   birth = 0, death = 1)
  df_phom <- as.PHom(df)
  
  # confirm that is.PHom works
  expect_true(is.PHom(df_phom))
  expect_false(is.PHom(df))
  
  # confirm that internal data was preserved
  expect_equal(df, as.data.frame(df_phom))
})
