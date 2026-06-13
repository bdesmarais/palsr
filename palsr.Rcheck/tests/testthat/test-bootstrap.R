test_that("bootstrap_pals returns a pals_boot with replicate estimates", {
  ev <- make_sim_events()
  bt <- bootstrap_pals(ev, R = 6, model = "one", seed = 1)
  expect_s3_class(bt, "pals_boot")
  expect_equal(bt$R, 6L)
  expect_true(nrow(bt$estimates) >= 1)
  expect_true("alpha" %in% names(bt$estimates))
  expect_s3_class(bt$estimate, "pals_fit")
})

test_that("bootstrap_pals is reproducible with a seed", {
  ev <- make_sim_events()
  a <- bootstrap_pals(ev, R = 5, model = "one", seed = 42)
  b <- bootstrap_pals(ev, R = 5, model = "one", seed = 42)
  expect_equal(a$estimates$alpha, b$estimates$alpha)
})

test_that("bootstrap_pals produces projections when predict_time is given", {
  ev <- make_sim_events()
  bt <- bootstrap_pals(ev, R = 4, model = "one", seed = 1,
                       predict_time = as.Date("2012-01-01"))
  expect_false(is.null(bt$projections))
  expect_true(all(c("actor", "lon", "lat", "replicate") %in% names(bt$projections)))
})

test_that("bootstrap summary / print / coef methods work", {
  ev <- make_sim_events()
  bt <- bootstrap_pals(ev, R = 6, model = "one", seed = 1)
  expect_output(print(bt), "pals_boot")
  tab <- summary(bt)
  expect_true(all(c("estimate", "boot_se", "lower", "upper") %in% names(tab)))
  expect_type(coef(bt), "double")
})

test_that("bootstrap_pals validates its inputs", {
  ev <- make_sim_events()
  expect_error(bootstrap_pals(ev, R = 0), "positive")
  expect_error(bootstrap_pals(list()), "pal_events")
})

test_that("pool_rubin pools estimates with Rubin's Rules", {
  q <- c(1.10, 0.95, 1.20, 1.05, 0.98)
  u <- c(0.04, 0.05, 0.045, 0.038, 0.052)
  out <- pool_rubin(q, u)
  expect_equal(out$qbar, mean(q))
  expect_equal(out$ubar, mean(u))
  expect_equal(out$b, var(q))
  expect_equal(out$t, out$ubar + (1 + 1/5) * out$b)
  expect_equal(out$se, sqrt(out$t))
  expect_true(out$fmi > 0 && out$fmi < 1)
})

test_that("pool_rubin returns df and p-value when requested", {
  q <- c(1.10, 0.95, 1.20, 1.05, 0.98)
  u <- c(0.04, 0.05, 0.045, 0.038, 0.052)
  out <- pool_rubin(q, u, df = TRUE, dfcom = 100)
  expect_true(all(c("df", "p.value") %in% names(out)))
  expect_true(out$df > 0)
  expect_true(out$p.value >= 0 && out$p.value <= 1)
})

test_that("pool_rubin validates its inputs", {
  expect_error(pool_rubin(1:3, 1:2), "same length")
  expect_error(pool_rubin(1, 1), "at least 2")
  expect_error(pool_rubin(c(1, 2), c(-1, 1)), "non-negative")
})
