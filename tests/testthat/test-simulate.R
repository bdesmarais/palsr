test_that("simulate_conflict_events is reproducible with a seed", {
  a <- simulate_conflict_events(n_actors = 6, n_events = 50, seed = 99)
  b <- simulate_conflict_events(n_actors = 6, n_events = 50, seed = 99)
  expect_equal(a, b)
})

test_that("simulate_conflict_events respects requested dimensions", {
  ev <- simulate_conflict_events(n_actors = 10, n_events = 120, seed = 7)
  expect_s3_class(ev, "pal_events")
  expect_lte(nrow(ev), 120L)                    # self-dyads may be dropped
  expect_lte(length(unique(c(ev$actor1, ev$actor2))), 10L)
})

test_that("simulated coordinates fall inside the bounding box", {
  bbox <- c(2.7, 14.7, 4.0, 13.9)
  ev <- simulate_conflict_events(n_actors = 8, n_events = 100,
                                 bbox = bbox, seed = 3)
  expect_true(all(ev$lon >= bbox[1] & ev$lon <= bbox[2]))
  expect_true(all(ev$lat >= bbox[3] & ev$lat <= bbox[4]))
})

test_that("simulation does not leave the global RNG state altered", {
  set.seed(123)
  before <- .Random.seed
  invisible(simulate_conflict_events(n_actors = 5, n_events = 30, seed = 1))
  expect_identical(.Random.seed, before)
})

test_that("bundled nigeria_sim dataset is available and well-formed", {
  data(nigeria_sim, package = "pals")
  expect_s3_class(nigeria_sim, "pal_events")
  expect_equal(nrow(nigeria_sim), 1500L)
  expect_equal(length(unique(c(nigeria_sim$actor1, nigeria_sim$actor2))), 25L)
})
