# Parameter estimation -----------------------------------------------------------

#' @noRd
# Pre-extract and flatten the focal/alter histories for both actors of every target,
# so the optimizer objective is a single compiled call per evaluation.
.build_units <- function(events, targets) {
  n <- nrow(targets)
  actor <- c(targets$actor1, targets$actor2)        # 2n: first n = side A
  time  <- rep(.as_date(targets$time), 2L)

  key <- paste(actor, as.character(time), sep = "\r")
  uk  <- unique(key)
  parts <- strsplit(uk, "\r", fixed = TRUE)
  hists <- vector("list", length(uk))
  for (i in seq_along(uk)) {
    hists[[i]] <- .history(events, parts[[i]][1], as.Date(parts[[i]][2]))
  }
  idx <- match(key, uk)

  flatten <- function(which_age, which_lon, which_lat) {
    lens <- vapply(idx, function(j) length(hists[[j]][[which_age]]), integer(1))
    off  <- c(0L, cumsum(lens)[-length(lens)])
    age  <- unlist(lapply(idx, function(j) hists[[j]][[which_age]]), use.names = FALSE)
    lon  <- unlist(lapply(idx, function(j) hists[[j]][[which_lon]]), use.names = FALSE)
    lat  <- unlist(lapply(idx, function(j) hists[[j]][[which_lat]]), use.names = FALSE)
    if (is.null(age)) age <- numeric(0)
    if (is.null(lon)) lon <- numeric(0)
    if (is.null(lat)) lat <- numeric(0)
    list(age = age, lon = lon, lat = lat, off = as.integer(off), len = as.integer(lens))
  }

  list(
    n = n,
    focal = flatten("focal_age", "focal_lon", "focal_lat"),
    alter = flatten("alter_age", "alter_lon", "alter_lat"),
    obs_lon = targets$lon, obs_lat = targets$lat
  )
}

#' @noRd
.theta_to_params <- function(theta, model) {
  if (model == "one") {
    pals_params(alpha = exp(theta[1]), model = "one")
  } else {
    pals_params(alpha = exp(theta[3]), beta = exp(theta[4]),
                gamma = theta[1], eta = theta[2], model = "four")
  }
}

#' @noRd
.objective <- function(theta, units, model, legacy, eps, radius, aggregate) {
  p <- .theta_to_params(theta, model)
  M <- project_batch_cpp(units$focal$age, units$focal$lon, units$focal$lat,
                         units$focal$off, units$focal$len,
                         units$alter$age, units$alter$lon, units$alter$lat,
                         units$alter$off, units$alter$len,
                         alpha = p$alpha, beta = p$beta, gamma = p$gamma, eta = p$eta,
                         pi_zero = (model == "one"), alter_legacy = legacy, eps = eps)
  n <- units$n
  lon <- .rowmean2(M[seq_len(n), 1], M[n + seq_len(n), 1])
  lat <- .rowmean2(M[seq_len(n), 2], M[n + seq_len(n), 2])
  d <- haversine_cpp(units$obs_lon, units$obs_lat, lon, lat, radius)
  if (aggregate == "sum") sum(d, na.rm = TRUE) else mean(d, na.rm = TRUE)
}

#' Estimate PALS parameters
#'
#' Estimates the PALS smoothing parameters by minimizing the mean great-circle
#' (Haversine) distance between observed event locations and the locations predicted
#' from each event's preceding history ("marching forward": every prediction uses only
#' events strictly earlier than the event being predicted).
#'
#' @param events A [pal_events] object providing the actor histories.
#' @param fit_events Optional `data.frame` of target events to fit against (needs
#'   `actor1`, `actor2`, `time`, `lon`, `lat`). Defaults to all of `events`; events with
#'   no usable history contribute nothing and are ignored.
#' @param model `"four"` (full: focal + alter histories, estimates `alpha`, `beta`,
#'   `gamma`, `eta`) or `"one"` (focal-only; estimates `alpha`, with `pi = 0`).
#' @param start Optional numeric starting vector on the optimizer's scale
#'   (`c(gamma, eta, log_alpha, log_beta)` for the four-parameter model, `log_alpha` for
#'   the one-parameter model). Sensible defaults are used if `NULL`.
#' @param method Optimizer method passed to [stats::optim] (`"Nelder-Mead"` for the
#'   four-parameter model, `"Brent"` for the one-parameter model by default).
#' @param aggregate `"mean"` (default) or `"sum"` of per-event distances.
#' @param alter_weight,eps See [project_pal].
#' @param radius Sphere radius for the Haversine objective (km).
#' @param control A list of control parameters for [stats::optim].
#'
#' @return An object of class `pals_fit` with components `params` (estimated
#'   [pals_params]), `model`, `objective` (minimized mean/sum distance), `n_used`
#'   (events contributing), `convergence`, `optim` (raw optimizer output), `events`,
#'   `settings`, and `call`. Methods: [print()], [summary()], [coef()], [predict()],
#'   [plot()].
#'
#' @seealso [project_pals()], [predict_event_locations()], [bootstrap_pals()].
#'
#' @examples
#' ev  <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
#' fit <- estimate_pals(ev, model = "one")
#' fit
#' coef(fit)
#'
#' @export
estimate_pals <- function(events, fit_events = NULL,
                          model = c("four", "one"),
                          start = NULL, method = NULL,
                          aggregate = c("mean", "sum"),
                          alter_weight = c("normalized", "legacy"),
                          eps = 0.01, radius = 6371.0088,
                          control = list()) {
  stopifnot(inherits(events, "pal_events"))
  model <- match.arg(model)
  aggregate <- match.arg(aggregate)
  alter_weight <- match.arg(alter_weight)
  legacy <- identical(alter_weight, "legacy")
  if (is.null(fit_events)) fit_events <- events
  need <- c("actor1", "actor2", "time", "lon", "lat")
  if (!all(need %in% names(fit_events)))
    stop("`fit_events` must contain columns: ", paste(need, collapse = ", "),
         call. = FALSE)

  units <- .build_units(events, fit_events)

  obj <- function(theta)
    .objective(theta, units, model, legacy, eps, radius, aggregate)

  if (model == "one") {
    if (is.null(start)) start <- 0          # log_alpha = 0  => alpha = 1
    if (is.null(method)) method <- "Brent"
    op <- stats::optim(par = start, fn = obj, method = method,
                       lower = -8, upper = 5, control = control)
  } else {
    if (is.null(start)) start <- c(-10, -10, -0.5, -5)  # gamma, eta, log a, log b
    if (is.null(method)) method <- "Nelder-Mead"
    op <- stats::optim(par = start, fn = obj, method = method, control = control)
  }

  params <- .theta_to_params(op$par, model)
  # how many events actually contributed (non-NA predicted distance at the optimum)?
  M <- project_batch_cpp(units$focal$age, units$focal$lon, units$focal$lat,
                         units$focal$off, units$focal$len,
                         units$alter$age, units$alter$lon, units$alter$lat,
                         units$alter$off, units$alter$len,
                         alpha = params$alpha, beta = params$beta,
                         gamma = params$gamma, eta = params$eta,
                         pi_zero = (model == "one"), alter_legacy = legacy, eps = eps)
  n <- units$n
  predlon <- .rowmean2(M[seq_len(n), 1], M[n + seq_len(n), 1])
  n_used <- sum(!is.na(haversine_cpp(units$obs_lon, units$obs_lat,
                                     predlon, predlon, radius)))

  structure(
    list(params = params, model = model, objective = op$value,
         n_used = n_used, convergence = op$convergence, optim = op,
         events = events,
         settings = list(alter_weight = alter_weight, eps = eps,
                         radius = radius, aggregate = aggregate),
         call = match.call()),
    class = "pals_fit"
  )
}

#' @export
print.pals_fit <- function(x, ...) {
  cat(sprintf("<pals_fit> %s-parameter PALS model\n", x$model))
  print(x$params)
  cat(sprintf("  objective (%s Haversine km): %.4f over %d events\n",
              x$settings$aggregate, x$objective, x$n_used))
  if (x$convergence != 0)
    cat(sprintf("  NOTE: optimizer convergence code %d\n", x$convergence))
  invisible(x)
}

#' @export
summary.pals_fit <- function(object, ...) {
  print(object)
  cat(sprintf("  alter weighting: %s; eps: %g; radius: %g km\n",
              object$settings$alter_weight, object$settings$eps,
              object$settings$radius))
  invisible(object)
}

#' @export
coef.pals_fit <- function(object, ...) {
  p <- object$params
  if (object$model == "one") c(alpha = p$alpha)
  else c(alpha = p$alpha, beta = p$beta, gamma = p$gamma, eta = p$eta)
}

#' Project locations from a fitted PALS model
#'
#' @param object A `pals_fit` from [estimate_pals].
#' @param newdata Optional [pal_events] (for `type = "pal"`) or target `data.frame` with
#'   `actor1`, `actor2`, `time` (for `type = "event"`). Defaults to the fitted events.
#' @param predict_time For `type = "pal"`, a `Date` (or vector) at which to project.
#' @param type `"pal"` projects per-actor locations; `"event"` predicts dyadic event
#'   locations (mean of the two actors' PALs).
#' @param actors For `type = "pal"`, which actors to project (default: all).
#' @param ... Unused.
#'
#' @return A `data.frame` of projections (see [project_pals] / [predict_event_locations]).
#'
#' @examples
#' ev  <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
#' fit <- estimate_pals(ev, model = "one")
#' predict(fit, predict_time = as.Date("2013-12-01"), type = "pal")[1:5, ]
#'
#' @export
predict.pals_fit <- function(object, newdata = NULL, predict_time = NULL,
                             type = c("pal", "event"), actors = NULL, ...) {
  type <- match.arg(type)
  events <- if (is.null(newdata) || identical(type, "event")) object$events else newdata
  aw <- object$settings$alter_weight; eps <- object$settings$eps
  if (type == "pal") {
    if (is.null(predict_time))
      stop("Provide `predict_time` for type = 'pal'.", call. = FALSE)
    ev <- if (inherits(newdata, "pal_events")) newdata else object$events
    project_pals(ev, actors = actors, predict_time = predict_time,
                 params = object$params, alter_weight = aw, eps = eps)
  } else {
    targets <- if (is.null(newdata)) object$events else newdata
    predict_event_locations(object$events, targets, object$params,
                            alter_weight = aw, eps = eps)
  }
}
