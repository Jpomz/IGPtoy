test_that("B_extinct output dimensions", {
  N <- matrix(rpois(15, lambda = 100),
              nrow = 3,
              ncol = 15)
  N[1, 1:2] <- 0
  out <- B_extinct(N)
  expect_equal(dim(out), dim(N))
})

test_that("B_extinct output colSums <= input", {
  N <- matrix(rpois(15, lambda = 100),
              nrow = 3,
              ncol = 15)
  N[1, 1:2] <- 0
  out <- B_extinct(N)
  expect_equal(all(colSums(out) <= colSums(N)), TRUE)
})
