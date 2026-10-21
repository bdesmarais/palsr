# Genuine out-of-sample (rolling-origin) evaluation ------------------------------

#' Rolling-origin evaluation of a PALS model
#'
#' Estimates PALS at each of a set of origin dates using only events strictly
#' before the origin, then scores the model's predictions for events in the
#' following horizon window, which did not contribute to the fitted parameters.
#' This is a genuine out-of-sample evaluation, unlike the "marching forward"
#' one-step-ahead (prequential) criterion [estimate_pals] minimizes, in which
#' every event's own location contributes to the single fitted parameter vector
#' even though its prediction uses only earlier history.
#'
#' @param events A [pal_events] object.
#' @param origins A `Date` vector of origin dates; parameters are re-estimated at
#'   each one from events strictly before it.
#' @param horizon Width of the evaluation window after each origin, passed to
#'   [seq.Date] (default `"1 year"`, i.e. the origin date through 364 days later).
#' @param model `"four"` or `"one"` (passed to [estimate_pals]).
#' @param ... Further arguments passed to [estimate_pals].
#'
#' @return An object of class `pals_rolling` with components `by_origin` (a
#'   `data.frame`, one row per origin, with `n_train`, `n_holdout`, `mean_error_km`,
#'   `median_error_km`, and the fit's `convergence` code) and `pooled` (the mean and
#'   median error pooling every origin's held-out events). Methods: [print()].
#'
#' @examples
#' ev <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
#' origins <- as.Date(c("2008-01-01", "2010-01-01"))
#' ro <- rolling_origin_pals(ev, origins, model = "one")
#' print(ro)
#'
#' @export
rolling_origin_pals <- function(events, origins, horizon = "1 year",
                                model = c("four", "one"), ...) {
  stopifnot(inherits(events, "pal_events"))
  model <- match.arg(model)
  origins <- .as_date(origins)

  horizon_end <- function(o) seq(o, by = horizon, length.out = 2)[2] - 1

  rows <- vector("list", length(origins))
  pooled_err <- numeric(0)
  for (i in seq_along(origins)) {
    o <- origins[i]
    train <- events[events$time < o, , drop = FALSE]
    class(train) <- c("pal_events", "data.frame")
    holdout <- events[events$time >= o & events$time <= horizon_end(o), , drop = FALSE]
    if (nrow(holdout) == 0 || nrow(train) == 0) next

    fit <- tryCatch(estimate_pals(train, model = model, ...), error = function(e) NULL)
    if (is.null(fit)) next
    pred <- predict_event_locations(events, holdout, fit$params)

    rows[[i]] <- data.frame(
      origin = o, n_train = nrow(train), n_holdout = nrow(holdout),
      mean_error_km = mean(pred$error_km, na.rm = TRUE),
      median_error_km = median(pred$error_km, na.rm = TRUE),
      convergence = fit$convergence
    )
    pooled_err <- c(pooled_err, pred$error_km)
  }
  by_origin <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  rownames(by_origin) <- NULL

  structure(
    list(by_origin = by_origin,
         pooled = data.frame(n = sum(!is.na(pooled_err)),
                             mean_error_km = mean(pooled_err, na.rm = TRUE),
                             median_error_km = median(pooled_err, na.rm = TRUE)),
         model = model, origins = origins, horizon = horizon, call = match.call()),
    class = "pals_rolling"
  )
}

#' @export
print.pals_rolling <- function(x, ...) {
  cat(sprintf("<pals_rolling> %s-parameter PALS, %d origins, horizon = %s\n",
              x$model, length(x$origins), x$horizon))
  print(x$by_origin)
  cat(sprintf("Pooled: n=%d  mean=%.1f km  median=%.1f km\n",
              x$pooled$n, x$pooled$mean_error_km, x$pooled$median_error_km))
  invisible(x)
}
