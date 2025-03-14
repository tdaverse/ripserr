context("printing")

# arbitrary 'PHom' object
x <- cbind(
  x = runif(n = 24L),
  y = runif(n = 24L)
)
x_vr <- vietoris_rips(x, dimension = 2L)

test_that("`print()` method returns 'PHom' object invisibly", {
  # detect printed output
  expect_output(out <- print(x_vr), "PHom object")
  # confirm common class
  expect_equal(class(x_vr), class(out))
  # test equality
  expect_equal(x_vr, out)
})
