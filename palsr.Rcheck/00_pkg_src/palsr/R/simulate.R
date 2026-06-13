# Simulated example data ---------------------------------------------------------

#' Simulate dyadic conflict events between moving actors
#'
#' Generates a synthetic dataset of dyadic interaction events with the qualitative
#' structure PALS is designed for: a set of armed-group-like actors, each following a
#' slowly drifting spatial trajectory, that interact preferentially with nearby actors.
#' Useful for examples, tests, and the package vignette. The geographic frame
#' approximates Nigeria, echoing the application in Kim, Liu and Desmarais (2023).
#'
#' @param n_actors Number of actors (default 30).
#' @param n_events Number of dyadic events (default 2000).
#' @param years Integer vector of years to span (default `2000:2016`).
#' @param drift Standard deviation (decimal degrees per year) of each actor's directional
#'   drift; larger values make actors more mobile (default `0.18`).
#' @param jitter Standard deviation (decimal degrees) of event-location noise around the
#'   midpoint of the two actors' current locations (default `0.25`).
#' @param decay Spatial interaction scale (degrees): partners are chosen with probability
#'   proportional to `exp(-distance / decay)`, so nearby actors interact more (default `2`).
#' @param bbox Bounding box `c(lon_min, lon_max, lat_min, lat_max)` for actor home
#'   locations (default approximately Nigeria).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A [pal_events] object with columns `actor1`, `actor2`, `time`, `lon`, `lat`.
#'
#' @examples
#' ev <- simulate_conflict_events(n_actors = 12, n_events = 400, seed = 42)
#' ev
#' summary(ev)
#'
#' @export
simulate_conflict_events <- function(n_actors = 30, n_events = 2000,
                                     years = 2000:2016, drift = 0.18,
                                     jitter = 0.25, decay = 2,
                                     bbox = c(2.7, 14.7, 4.0, 13.9),
                                     seed = NULL) {
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  stopifnot(n_actors >= 2, n_events >= 1, length(bbox) == 4)

  ids <- sprintf("G%02d", seq_len(n_actors))
  start <- as.Date(sprintf("%d-01-01", min(years)))
  end   <- as.Date(sprintf("%d-12-31", max(years)))
  span_days <- as.numeric(end - start)
  span_years <- span_days / 365.25

  # Actor trajectories: home location + linear drift velocity (deg / year).
  home_lon <- stats::runif(n_actors, bbox[1], bbox[2])
  home_lat <- stats::runif(n_actors, bbox[3], bbox[4])
  vel_lon  <- stats::rnorm(n_actors, 0, drift)
  vel_lat  <- stats::rnorm(n_actors, 0, drift)
  activity <- stats::rgamma(n_actors, shape = 2, rate = 2)  # heterogeneous activity

  loc_at <- function(t_years) {
    cbind(
      .clamp(home_lon + vel_lon * t_years, bbox[1], bbox[2]),
      .clamp(home_lat + vel_lat * t_years, bbox[3], bbox[4])
    )
  }

  # Event times, sorted.
  day_off <- sort(sample.int(span_days + 1L, n_events, replace = TRUE) - 1L)
  t_years <- day_off / 365.25

  a1 <- character(n_events); a2 <- character(n_events)
  lon <- numeric(n_events);  lat <- numeric(n_events)

  for (e in seq_len(n_events)) {
    L <- loc_at(t_years[e])
    i <- sample.int(n_actors, 1, prob = activity)
    d <- sqrt((L[, 1] - L[i, 1])^2 + (L[, 2] - L[i, 2])^2)
    w <- exp(-d / decay) * activity
    w[i] <- 0
    j <- sample.int(n_actors, 1, prob = w)
    a1[e] <- ids[i]; a2[e] <- ids[j]
    lon[e] <- .clamp((L[i, 1] + L[j, 1]) / 2 + stats::rnorm(1, 0, jitter), bbox[1], bbox[2])
    lat[e] <- .clamp((L[i, 2] + L[j, 2]) / 2 + stats::rnorm(1, 0, jitter), bbox[3], bbox[4])
  }

  pal_events(
    data.frame(actor1 = a1, actor2 = a2, time = start + day_off,
               lon = lon, lat = lat, stringsAsFactors = FALSE)
  )
}

#' @noRd
.clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

#' @noRd
.Random.seed_safe <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL
}

#' @noRd
.restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
}
