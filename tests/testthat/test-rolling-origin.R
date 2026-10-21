test_that("rolling_origin_pals evaluates genuinely held-out events", {
  ev <- make_sim_events()
  origins <- as.Date(c("2008-06-01", "2010-06-01"))
  ro <- rolling_origin_pals(ev, origins, model = "one")
  expect_s3_class(ro, "pals_rolling")
  expect_true(all(c("n_train", "n_holdout", "mean_error_km") %in% names(ro$by_origin)))
  expect_true(ro$pooled$n > 0)
  expect_output(print(ro), "pals_rolling")
})

test_that("rolling_origin_pals origins use only strictly earlier training data", {
  ev <- make_sim_events()
  origin <- as.Date("2009-01-01")
  ro <- rolling_origin_pals(ev, origin, model = "one")
  expect_true(all(ev$time[ev$time < origin] <= origin))
  expect_equal(ro$by_origin$n_train, sum(ev$time < origin))
})
