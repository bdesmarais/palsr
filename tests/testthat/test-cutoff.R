test_that("cutoff granularity controls the history window", {
  ev <- pal_events(data.frame(
    actor1 = c("A", "A", "A"), actor2 = c("B", "B", "B"),
    time = as.Date(c("2009-06-01", "2010-03-01", "2010-10-01")),
    lon = c(1, 2, 3), lat = c(1, 2, 3), stringsAsFactors = FALSE),
    drop_self = FALSE)
  p <- pals_params(alpha = 0.5, model = "one")

  # day: all events strictly before 2010-12-01 -> 3
  expect_equal(project_pal(ev, "A", as.Date("2010-12-01"), p, cutoff = "day")$n_focal, 3)
  # month: events before 2010-10-01 -> the 2009 and 2010-03 events -> 2
  expect_equal(project_pal(ev, "A", as.Date("2010-10-15"), p, cutoff = "month")$n_focal, 2)
  # year: events before 2010-01-01 -> only the 2009 event -> 1
  expect_equal(project_pal(ev, "A", as.Date("2010-12-01"), p, cutoff = "year")$n_focal, 1)
})

test_that("default cutoff is day (strict date), and bad values error", {
  ev <- simulate_conflict_events(n_actors = 6, n_events = 80, seed = 3)
  p <- pals_params(alpha = 0.6, model = "one")
  pt <- as.Date("2010-12-01")
  expect_identical(
    project_pals(ev, predict_time = pt, params = p),
    project_pals(ev, predict_time = pt, params = p, cutoff = "day"))
  expect_error(project_pals(ev, predict_time = pt, params = p, cutoff = "decade"))
})
