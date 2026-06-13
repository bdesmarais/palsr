# Distances and event-location prediction ----------------------------------------

#' Great-circle (Haversine) distance
#'
#' Vectorized great-circle distance between longitude/latitude points. Arguments are
#' recycled to a common length, so any may be length 1.
#'
#' @param lon1,lat1,lon2,lat2 Numeric vectors of coordinates in decimal degrees.
#' @param radius Sphere radius in the desired output units. Defaults to the mean Earth
#'   radius, `6371.0088` km, so distances are returned in kilometres.
#'
#' @return A numeric vector of distances. `NA` in any coordinate gives `NA`.
#'
#' @examples
#' haversine(0, 0, 0, 1)          # ~111 km per degree of latitude
#' haversine(7.4, 9.1, 8.5, 12.0) # Abuja-ish to Kano-ish
#'
#' @export
haversine <- function(lon1, lat1, lon2, lat2, radius = 6371.0088) {
  haversine_cpp(as.numeric(lon1), as.numeric(lat1),
                as.numeric(lon2), as.numeric(lat2), radius = radius)
}

#' @noRd
# Compute a PAL for each unique (actor, time) pair; return a named-lookup data.frame.
.pal_lookup <- function(events, actor, time, p, legacy = FALSE, eps = 0.01) {
  key <- paste(actor, as.character(time), sep = "\r")
  uk <- !duplicated(key)
  ua <- actor[uk]; ut <- time[uk]; ukey <- key[uk]
  lon <- numeric(length(ukey)); lat <- numeric(length(ukey))
  for (i in seq_along(ukey)) {
    o <- .project_actor_time(events, ua[i], ut[i], p, alter_legacy = legacy, eps = eps)
    lon[i] <- o[[1]]; lat[i] <- o[[2]]
  }
  idx <- match(key, ukey)
  list(lon = lon[idx], lat = lat[idx])
}

#' Predict dyadic event locations
#'
#' For each dyadic target (a pair of actors at a time), predicts the interaction
#' location as the mean of the two actors' Projected Actor Locations. Optionally scores
#' the prediction against an observed location.
#'
#' @param events A [pal_events] object supplying the histories.
#' @param targets A `data.frame` with columns `actor1`, `actor2`, `time`, and optionally
#'   `lon`/`lat` giving the observed event location (for error scoring).
#' @param params A [pals_params] or fitted [estimate_pals] object.
#' @param alter_weight,eps See [project_pal].
#'
#' @return `targets` augmented with `pred_lon`, `pred_lat`, and (if observed `lon`/`lat`
#'   were supplied) `error_km`, the Haversine distance between predicted and observed
#'   locations. Predictions are `NA` when both actors lack usable history.
#'
#' @examples
#' ev  <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
#' fit <- estimate_pals(ev, model = "one")
#' tg  <- ev[ev$time > as.Date("2012-01-01"), ]
#' head(predict_event_locations(ev, tg, fit))
#'
#' @export
predict_event_locations <- function(events, targets, params,
                                     alter_weight = c("normalized", "legacy"),
                                     eps = 0.01) {
  stopifnot(inherits(events, "pal_events"))
  alter_weight <- match.arg(alter_weight)
  p <- as_pals_params(params)
  legacy <- identical(alter_weight, "legacy")
  tm <- .as_date(targets$time)

  a <- .pal_lookup(events, targets$actor1, tm, p, legacy, eps)
  b <- .pal_lookup(events, targets$actor2, tm, p, legacy, eps)

  # mean of the two PALs, na.rm = TRUE (fall back to the available one)
  pred_lon <- .rowmean2(a$lon, b$lon)
  pred_lat <- .rowmean2(a$lat, b$lat)

  out <- targets
  out$pred_lon <- pred_lon
  out$pred_lat <- pred_lat
  if (all(c("lon", "lat") %in% names(targets))) {
    out$error_km <- haversine(out$lon, out$lat, pred_lon, pred_lat)
  }
  out
}

#' @noRd
.rowmean2 <- function(x, y) {
  m <- cbind(x, y)
  rowMeans(m, na.rm = TRUE) * ifelse(rowSums(!is.na(m)) == 0, NA, 1)
}

#' Dyadic distance between Projected Actor Locations
#'
#' Builds the dyadic distance covariate used to model interaction likelihood: the
#' Haversine distance between the two actors' Projected Actor Locations.
#'
#' @param events A [pal_events] object.
#' @param dyads A `data.frame` with columns `actor1`, `actor2`, `time`.
#' @param params A [pals_params] or fitted [estimate_pals] object.
#' @param transform `"none"` (default) returns distance in km; `"log"` returns
#'   `log(distance + offset)`, as used for interstate-conflict-style specifications.
#' @param offset Offset added before logging (default `0.01`).
#' @param alter_weight,eps See [project_pal].
#'
#' @return `dyads` augmented with `pal_distance` (and, for `transform = "log"`,
#'   `pal_log_distance`).
#'
#' @examples
#' ev  <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
#' fit <- estimate_pals(ev, model = "one")
#' dy  <- data.frame(actor1 = "G01", actor2 = "G02",
#'                   time = as.Date("2012-12-01"))
#' pal_distance(ev, dy, fit)
#'
#' @export
pal_distance <- function(events, dyads, params,
                         transform = c("none", "log"), offset = 0.01,
                         alter_weight = c("normalized", "legacy"), eps = 0.01) {
  stopifnot(inherits(events, "pal_events"))
  transform <- match.arg(transform)
  alter_weight <- match.arg(alter_weight)
  p <- as_pals_params(params)
  legacy <- identical(alter_weight, "legacy")
  tm <- .as_date(dyads$time)

  a <- .pal_lookup(events, dyads$actor1, tm, p, legacy, eps)
  b <- .pal_lookup(events, dyads$actor2, tm, p, legacy, eps)
  d <- haversine(a$lon, a$lat, b$lon, b$lat)

  out <- dyads
  out$pal_distance <- d
  if (transform == "log") out$pal_log_distance <- log(d + offset)
  out
}
