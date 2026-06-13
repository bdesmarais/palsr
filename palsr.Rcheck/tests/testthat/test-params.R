test_that("pals_params builds a four-parameter set", {
  p <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
  expect_s3_class(p, "pals_params")
  expect_equal(p$model, "four")
  expect_equal(p$alpha, 0.9)
  expect_equal(p$beta, 0.2)
})

test_that("pals_params builds a one-parameter set with pi fixed at 0", {
  p <- pals_params(alpha = 0.5, model = "one")
  expect_equal(p$model, "one")
  expect_equal(p$alpha, 0.5)
  expect_output(print(p), "pi fixed at 0")
})

test_that("pals_params validates inputs", {
  expect_error(pals_params(alpha = -1), "non-negative")
  expect_error(pals_params(alpha = c(1, 2)), "single")
  expect_error(pals_params(alpha = 1, beta = -1), "non-negative")
  expect_error(pals_params(alpha = 1, gamma = NA), "number")
})

test_that("as_pals_params extracts params from a fit", {
  p <- pals_params(alpha = 0.9, model = "one")
  expect_identical(palsr:::as_pals_params(p), p)
  expect_error(palsr:::as_pals_params(list()), "Expected")
})
