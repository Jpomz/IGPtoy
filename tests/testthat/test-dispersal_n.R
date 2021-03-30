test_that("dispersal_n output dim", {
  n_patch = sample(3:10, 1)
  n_sp = 3

  x_coord = runif(n_patch, 0, 5)
  y_coord = runif(n_patch, 0, 5)
  dist_mat = data.matrix(dist(cbind(x_coord,
                                    y_coord)))
  N <- matrix(rpois(n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
  N[,N[1,] == 0] <- 0

  out <- dispersal_n(N = N, dist_mat = dist_mat, theta = 1,
                     v_p_dispersal = c(0.1, 0.1, 0.1))

  expect_length(out, length(N))
})

test_that("dispersal_n output not NA", {
  n_patch = sample(3:10, 1)
  n_sp = 3

  x_coord = runif(n_patch, 0, 5)
  y_coord = runif(n_patch, 0, 5)
  dist_mat = data.matrix(dist(cbind(x_coord,
                                    y_coord)))
  N <- matrix(rpois(n_sp*n_patch, lambda = 100), nrow = n_sp, ncol = n_patch)
  N[,N[1,] == 0] <- 0

  out <- dispersal_n(N = N, dist_mat = dist_mat, theta = 1,
                     v_p_dispersal = c(0.1, 0.1, 0.1))

  out_na <- all(!is.na(out))
  out_gte_0 <- all(out >=0)
  expect_equal(out_na, TRUE)
  expect_equal(out_gte_0, TRUE)
})


