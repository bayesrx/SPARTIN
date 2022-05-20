library(spatstat)

test_that("reject malformed", {
  expect_error(IsHole(c(1,2,3), c(1,2)), "length of x must equal length of y")
  expect_error(IsHole(c("a", "b", "c"), c(1,2,3)), "x and y must be numeric vectors")
  expect_error(IsHole(c(1,2,3), c("a", "b", "c")), "x and y must be numeric vectors")
})

test_that("warn malformed", {
  expect_warning(IsHole(c(1,2), c(1,2)), "length of point list is less than three")
})

test_that("non-holes", {
  expect_equal(IsHole(c(0, 1, 1, 0), c(0, 0, 1, 1)), FALSE)
  expect_equal(IsHole(c(0, 1, 0.5), c(0, 0, 1)), FALSE)
})

test_that("holes", {
  expect_equal(IsHole(c(0, 0, 1, 1), c(0, 1, 1, 0)), TRUE)
  expect_equal(IsHole(c(0, 0.5, 1), c(0, 1, 1)), TRUE)
})
