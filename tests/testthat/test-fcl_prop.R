test_that("fcl_prop output names", {
  names.out <- names(
    fcl_prop(
      rep(c(0, 1, 2, 2.5, 3),
          each = 3)))
  expect_identical(names.out,
                   c("p0", "p1", "p2", "p2.5", "p3"))
})

test_that("fcl_prop output values", {
  out <- fcl_prop(rep(c(0, 1, 2, 2.5, 3), each = 3))
  expect_equal(out, c(0.2, 0.2, 0.2, 0.2, 0.2),
               ignore_attr = TRUE)
})
