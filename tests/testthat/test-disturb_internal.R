test_that("disturb_type = NULL, no disturbance", {
  N <- matrix(rpois(3*5, lambda = 100),
              nrow = 3, ncol = 5)
  out <- disturb_internal(N = N, disturb_type = NULL)
  expect_equal(out$N, N)
  expect_equal(length(out$patch_extinction), ncol(N))
})


test_that("extreme disturbance does not result in -N", {
n_patch = 5
N <- matrix(c(50, 50, 50, #patch 1, at bottom
              100, 100, 100, #next three patches are source
              100, 50, 10,
              50, 50, 50,
              100, 100, 100), # patch 5, confluence
            ncol = n_patch)

dist_mat <- matrix(
  c(0, 1, 2, 2, 1,
    1, 0, 3, 3, 2,
    2, 3, 0, 2, 1,
    2, 3, 2, 0, 1,
    1, 2, 1, 1, 0),
  ncol = n_patch,
  dimnames = list(
    c("patch1", "patch2", "patch3", "patch4", "patch5"),
    c("patch1", "patch2", "patch3", "patch4", "patch5")))

out <- disturb_internal(N = N, disturb_type = "point-source",
                 adjacency_matrix = NULL,
                 dist_mat = dist_mat,
                 disturb_mag = 1,
                 disturb_p = 1, disturb_rho = 1)
expect_equal(out$N, matrix(0, nrow = 3, ncol = n_patch))
expect_equal(length(out$patch_extinction), ncol(N))}
)





test_that("regional disturbance reduces N correctly when env_val known",{
  n_patch = 5
  N <- matrix(100, # patch 5, confluence
              ncol = n_patch,
              nrow = 3)
  environment_value = c(-2, -1, 0, 1, 2)
  # this turns into 0.11, 0.27, 0.5, 0.73, 0.88 on inv.logit

out <- disturb_internal(N = N, environment_value = environment_value,
                 disturb_type = "regional",
                 disturb_p = 1, disturb_mag = 1)
expect_equal(round(colSums(out$N) / colSums(N), 3),
             c(0.881, 0.731, 0.500, 0.269, 0.119))
expect_equal(length(out$patch_extinction), ncol(N))
}
)

