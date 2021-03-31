test_that("env_to_k returns positive", {
  out <- env_to_k(c(-100, 0, 100), k_base = 0.001)
  expect_identical(all(out >= 0), TRUE)
})


test_that("env_to_k returns k_base when length(env) == 1", {
  k_base = 150
  out <- env_to_k(0, k_base = k_base)
  expect_identical(out, k_base)
})
