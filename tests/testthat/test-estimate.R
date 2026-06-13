test_that("estimate_pals fits the one-parameter model", {
  ev  <- make_sim_events()
  fit <- estimate_pals(ev, model = "one")
  expect_s3_class(fit, "pals_fit")
  expect_equal(fit$model, "one")
  cf <- coef(fit)
  expect_named(cf, "alpha")
  expect_true(cf[["alpha"]] > 0)
  expect_true(is.finite(fit$objective))
  expect_true(fit$n_used > 0)
})

test_that("estimate_pals fits the four-parameter model", {
  ev  <- make_sim_events(n_events = 200)
  fit <- estimate_pals(ev, model = "four",
                       control = list(maxit = 80))
  expect_equal(fit$model, "four")
  cf <- coef(fit)
  expect_named(cf, c("alpha", "beta", "gamma", "eta"))
  expect_true(all(is.finite(cf)))
})

test_that("estimate_pals errors when fit_events lacks required columns", {
  ev <- make_sim_events()
  bad <- data.frame(actor1 = "G01", actor2 = "G02", time = Sys.Date())
  expect_error(estimate_pals(ev, fit_events = bad), "must contain")
})

test_that("fit print / summary / coef methods work", {
  ev  <- make_sim_events()
  fit <- estimate_pals(ev, model = "one")
  expect_output(print(fit), "pals_fit")
  expect_output(summary(fit), "alter weighting")
  expect_type(coef(fit), "double")
})

test_that("predict.pals_fit dispatches on type", {
  ev  <- make_sim_events()
  fit <- estimate_pals(ev, model = "one")

  pal <- predict(fit, predict_time = as.Date("2012-01-01"), type = "pal")
  expect_true(all(c("actor", "lon", "lat") %in% names(pal)))

  evt <- predict(fit, type = "event")
  expect_true(all(c("pred_lon", "pred_lat") %in% names(evt)))

  expect_error(predict(fit, type = "pal"), "predict_time")
})

test_that("the objective improves on a naive alpha", {
  ev  <- make_sim_events()
  fit <- estimate_pals(ev, model = "one")
  units <- pals:::.build_units(ev, ev)
  naive <- pals:::.objective(log(1), units, "one", FALSE, 0.01, 6371.0088, "mean")
  expect_lte(fit$objective, naive + 1e-6)
})
