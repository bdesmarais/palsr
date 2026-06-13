test_that("haversine matches known great-circle distances", {
  # One degree of latitude is ~111 km.
  expect_equal(haversine(0, 0, 0, 1), 111.195, tolerance = 1e-2)
  # Distance from a point to itself is zero.
  expect_equal(haversine(7.4, 9.1, 7.4, 9.1), 0)
  # Symmetry.
  expect_equal(haversine(0, 0, 10, 20), haversine(10, 20, 0, 0))
})

test_that("haversine recycles arguments and propagates NA", {
  d <- haversine(0, 0, c(0, 0, 0), c(1, 2, 3))
  expect_length(d, 3)
  expect_true(all(diff(d) > 0))
  expect_true(is.na(haversine(NA, 0, 0, 1)))
})

test_that("haversine respects the radius argument", {
  km <- haversine(0, 0, 0, 1)
  mi <- haversine(0, 0, 0, 1, radius = 3958.7613)
  expect_equal(mi / km, 3958.7613 / 6371.0088, tolerance = 1e-6)
})

test_that("predict_event_locations returns predictions and error scoring", {
  ev <- make_sim_events()
  fit <- estimate_pals(ev, model = "one")
  tg <- ev[ev$time > as.Date("2010-01-01"), ]
  out <- predict_event_locations(ev, tg, fit)
  expect_true(all(c("pred_lon", "pred_lat", "error_km") %in% names(out)))
  expect_equal(nrow(out), nrow(tg))
  expect_true(all(out$error_km[!is.na(out$error_km)] >= 0))
})

test_that("pal_distance returns non-negative distances and a log transform", {
  ev <- make_sim_events()
  p  <- pals_params(alpha = 0.9, model = "one")
  dy <- data.frame(actor1 = "G01", actor2 = "G02",
                   time = as.Date("2012-12-01"))
  d <- pal_distance(ev, dy, p)
  expect_true("pal_distance" %in% names(d))
  expect_true(d$pal_distance >= 0)

  dlog <- pal_distance(ev, dy, p, transform = "log")
  expect_true("pal_log_distance" %in% names(dlog))
  expect_equal(dlog$pal_log_distance, log(dlog$pal_distance + 0.01))
})
