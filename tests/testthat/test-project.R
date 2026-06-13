test_that("project_pal returns one row per prediction time", {
  ev <- make_toy_events()
  p  <- pals_params(alpha = 0.9, model = "one")
  times <- as.Date(c("2001-07-01", "2002-02-01"))
  out <- project_pal(ev, actor = "A", predict_time = times, params = p)
  expect_equal(nrow(out), 2L)
  expect_named(out, c("actor", "time", "lon", "lat",
                      "n_focal", "n_alter", "has_history"))
  expect_true(all(out$has_history))
})

test_that("project_pal returns NA when the actor has no prior history", {
  ev <- make_toy_events()
  p  <- pals_params(alpha = 0.9, model = "one")
  # Before any events, A has no history.
  out <- project_pal(ev, actor = "A", predict_time = as.Date("2000-01-01"),
                     params = p)
  expect_true(is.na(out$lon))
  expect_false(out$has_history)
})

test_that("project_pals covers the actor-by-time grid", {
  ev <- make_sim_events()
  p  <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
  times <- as.Date(c("2008-01-01", "2012-01-01"))
  out <- project_pals(ev, predict_time = times, params = p)
  n_actors <- length(unique(c(ev$actor1, ev$actor2)))
  expect_equal(nrow(out), n_actors * length(times))
})

test_that("the compiled kernel agrees with the pure-R reference", {
  ev <- make_sim_events()
  p  <- pals_params(alpha = 0.8, beta = 0.3, gamma = -5, eta = -3)
  tt <- as.Date("2011-06-01")
  for (actor in c("G01", "G02", "G03")) {
    h <- pals:::.history(ev, actor, tt)
    r_out <- pals:::project_one_r(h, p$alpha, p$beta, p$gamma, p$eta,
                                  pi_zero = FALSE)
    c_out <- pals:::.project_actor_time(ev, actor, tt, p)
    expect_equal(unname(r_out[["lon"]]), c_out[[1]], tolerance = 1e-8)
    expect_equal(unname(r_out[["lat"]]), c_out[[2]], tolerance = 1e-8)
  }
})

test_that("one-parameter projection ignores alters (pi = 0)", {
  ev <- make_sim_events()
  tt <- as.Date("2011-06-01")
  p_one  <- pals_params(alpha = 0.8, model = "one")
  # With pi forced to 0, the four-parameter projection at gamma -> -Inf matches.
  h <- pals:::.history(ev, "G01", tt)
  r_one <- pals:::project_one_r(h, 0.8, 0, 0, 0, pi_zero = TRUE)
  c_one <- pals:::.project_actor_time(ev, "G01", tt, p_one)
  expect_equal(unname(r_one[["lon"]]), c_one[[1]], tolerance = 1e-8)
})
