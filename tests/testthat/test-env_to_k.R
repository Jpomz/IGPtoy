test_that("env_to_k returns positive", {
  out <- env_to_k(c(-100, 0, 100), k_base = 0.001)
  expect_identical(all(out >= 0), TRUE)
})
