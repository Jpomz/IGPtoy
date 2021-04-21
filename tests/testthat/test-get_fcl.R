test_that("get_fcl returns vector of 0, 1, 2, 3", {
  N = sim_N(fcl_5_states = TRUE)
  expect_equal(get_fcl(N), c(0, 1, 2, 2, 3))
})

test_that("get_fcl returns 0 when 1st row N is 0", {
  N = sim_N(n_patch = 3)
  N[1,] <- 0
  expect_equal(get_fcl(N), c(0, 0, 0))
})

test_that("get_fcl returns vector of length(n_patch)", {
  n_patch = sample(10:100, 1)
  N = sim_N(n_patch = n_patch)
  expect_equal(length(get_fcl(N)), n_patch)
})

test_that("get_fcl returns value between (0,3)", {
  N = matrix(rpois(3*100, lambda = 0.5), nrow = 3, ncol = 100)
  output <- get_fcl(N)
  expect_lte(max(output), 3)
  expect_gte(min(output), 0)
})
