# Projection ---------------------------------------------------------------------

#' @noRd
# Extract the focal and alter event histories for one actor strictly before a time.
.history <- function(events, actor, predict_time) {
  before <- events$time < predict_time
  is_focal <- before & (events$actor1 == actor | events$actor2 == actor)
  foc <- events[is_focal, , drop = FALSE]

  # alters = distinct partners that shared an event with the focal actor (before t)
  partners <- ifelse(foc$actor1 == actor, foc$actor2, foc$actor1)
  alters <- unique(partners)

  if (length(alters)) {
    is_alter <- before & (events$actor1 %in% alters | events$actor2 %in% alters)
    alt <- events[is_alter, , drop = FALSE]
  } else {
    alt <- events[0, , drop = FALSE]
  }

  list(
    focal_age = as.numeric(predict_time - foc$time),
    focal_lon = foc$lon, focal_lat = foc$lat,
    alter_age = as.numeric(predict_time - alt$time),
    alter_lon = alt$lon, alter_lat = alt$lat
  )
}

#' @noRd
# Pure-R reference implementation of the projection kernel (mirrors project_one_cpp);
# used in tests to validate the compiled kernel.
project_one_r <- function(h, alpha, beta, gamma, eta,
                          pi_zero = FALSE, alter_legacy = FALSE, eps = 0.01) {
  n_i <- length(h$focal_age); n_k <- length(h$alter_age)
  if (n_i == 0) return(c(lon = NA_real_, lat = NA_real_, n_focal = 0, n_alter = n_k))
  wi <- 1 / (h$focal_age^alpha + eps); wi <- wi / sum(wi)
  foc_lon <- sum(wi * h$focal_lon); foc_lat <- sum(wi * h$focal_lat)
  pi <- 0
  if (!pi_zero && n_k > 0) {
    v <- (n_i / n_k)^(1 / sqrt(n_k))
    pi <- stats::plogis(gamma + eta * v)
  }
  alt_lon <- 0; alt_lat <- 0
  if (pi > 0 && n_k > 0) {
    if (alter_legacy) {
      scal <- 1 / sum(1 / (h$alter_age^beta + eps))
      alt_lon <- scal * sum(h$alter_lon); alt_lat <- scal * sum(h$alter_lat)
    } else {
      wk <- 1 / (h$alter_age^beta + eps); wk <- wk / sum(wk)
      alt_lon <- sum(wk * h$alter_lon); alt_lat <- sum(wk * h$alter_lat)
    }
  }
  c(lon = (1 - pi) * foc_lon + pi * alt_lon,
    lat = (1 - pi) * foc_lat + pi * alt_lat,
    n_focal = n_i, n_alter = n_k)
}

#' @noRd
.project_actor_time <- function(events, actor, predict_time, p,
                                alter_legacy = FALSE, eps = 0.01) {
  h <- .history(events, actor, predict_time)
  out <- project_one_cpp(h$focal_age, h$focal_lon, h$focal_lat,
                         h$alter_age, h$alter_lon, h$alter_lat,
                         alpha = p$alpha, beta = p$beta,
                         gamma = p$gamma, eta = p$eta,
                         pi_zero = (p$model == "one"),
                         alter_legacy = alter_legacy, eps = eps)
  out
}

#' Project the location of a single actor
#'
#' Computes the Projected Actor Location (PAL) for one actor at one or more prediction
#' times, given a parameter set.
#'
#' @param events A [pal_events] object.
#' @param actor The focal actor id (length-1 character).
#' @param predict_time A `Date` (or vector of `Date`s) at which to project. Only events
#'   strictly earlier than each prediction time are used.
#' @param params A [pals_params] object or a fitted [estimate_pals] object.
#' @param alter_weight `"normalized"` (default) uses the paper's normalized alter
#'   weights; `"legacy"` reproduces the original replication code (an un-normalized sum
#'   of alter coordinates). See `DECISIONS.md`.
#' @param eps Numerical offset inside each age weight (default `0.01`).
#'
#' @return A `data.frame` with one row per prediction time and columns `actor`, `time`,
#'   `lon`, `lat`, `n_focal`, `n_alter`, `has_history`. `lon`/`lat` are `NA` when the
#'   actor has no prior events.
#'
#' @examples
#' ev <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
#' p  <- pals_params(alpha = 0.9, model = "one")
#' project_pal(ev, actor = "G01", predict_time = as.Date("2010-12-01"), params = p)
#'
#' @export
project_pal <- function(events, actor, predict_time, params,
                        alter_weight = c("normalized", "legacy"), eps = 0.01) {
  stopifnot(inherits(events, "pal_events"))
  alter_weight <- match.arg(alter_weight)
  p <- as_pals_params(params)
  predict_time <- .as_date(predict_time)
  legacy <- identical(alter_weight, "legacy")

  rows <- lapply(predict_time, function(tt) {
    o <- .project_actor_time(events, actor, tt, p, alter_legacy = legacy, eps = eps)
    data.frame(actor = actor, time = tt,
               lon = o[[1]], lat = o[[2]],
               n_focal = o[[3]], n_alter = o[[4]],
               has_history = o[[3]] > 0,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Project locations for multiple actors
#'
#' Computes Projected Actor Locations for a set of actors at one or more prediction
#' times (an actor-by-time grid).
#'
#' @inheritParams project_pal
#' @param actors Character vector of actor ids. Defaults to all actors in `events`.
#' @param predict_time A `Date` or vector of `Date`s.
#'
#' @return A `data.frame` with columns `actor`, `time`, `lon`, `lat`, `n_focal`,
#'   `n_alter`, `has_history` (one row per actor-time combination).
#'
#' @examples
#' ev <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
#' p  <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
#' pal <- project_pals(ev, predict_time = as.Date("2010-12-01"), params = p)
#' head(pal)
#'
#' @export
project_pals <- function(events, actors = NULL, predict_time, params,
                         alter_weight = c("normalized", "legacy"), eps = 0.01) {
  stopifnot(inherits(events, "pal_events"))
  alter_weight <- match.arg(alter_weight)
  if (is.null(actors)) actors <- .actors(events)
  parts <- lapply(actors, function(a)
    project_pal(events, a, predict_time, params,
                alter_weight = alter_weight, eps = eps))
  out <- do.call(rbind, parts)
  out <- out[order(out$time, out$actor), , drop = FALSE]
  rownames(out) <- NULL
  out
}
