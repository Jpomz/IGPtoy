test_that("get_fcl returns length(vector) = ncol(N)", {
  N <- matrix(rpois(3 * 10, lambda = 1), nrow = 3, ncol = 10)
  N[,N[1,] == 0] <- 0
  expect_length(get_fcl(N = N), ncol(N))
})

test_that("get_fcl returns value from 0-3", {
  N <- matrix(rpois(3 * 10, lambda = 1), nrow = 3, ncol = 10)
  N[,N[1,] == 0] <- 0
  out <- get_fcl(N = N)
  expect_gte(min(out), 0)
  expect_lte(max(out), 3)
})
