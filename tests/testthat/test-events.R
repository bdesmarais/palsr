test_that("pal_events builds a validated, sorted pal_events object", {
  ev <- make_toy_events()
  expect_s3_class(ev, "pal_events")
  expect_s3_class(ev, "data.frame")
  expect_named(ev, c("actor1", "actor2", "time", "lon", "lat"))
  expect_true(all(diff(ev$time) >= 0))          # sorted by time
  expect_type(ev$actor1, "character")
  expect_s3_class(ev$time, "Date")
})

test_that("pal_events errors on missing columns", {
  df <- data.frame(actor1 = "A", actor2 = "B", time = Sys.Date())
  expect_error(pal_events(df), "not found")
})

test_that("pal_events drops self-dyads and missing rows", {
  df <- data.frame(
    actor1 = c("A", "B", "A"),
    actor2 = c("A", "C", "B"),                 # first row is a self-dyad
    time   = as.Date(c("2001-01-01", "2001-02-01", "2001-03-01")),
    lon    = c(1, 2, NA),                      # third row has a missing coord
    lat    = c(1, 2, 3)
  )
  expect_message(expect_message(
    ev <- pal_events(df), "self-dyad"), "missing")
  expect_equal(nrow(ev), 1L)
  expect_equal(ev$actor1, "B")
})

test_that("pal_events rejects out-of-range coordinates", {
  df <- data.frame(actor1 = "A", actor2 = "B", time = Sys.Date(),
                   lon = 999, lat = 0)
  expect_error(pal_events(df), "lon")
})

test_that("print and summary methods run without error", {
  ev <- make_toy_events()
  expect_output(print(ev), "pal_events")
  expect_output(summary(ev), "actors")
  expect_invisible(print(ev))
})
