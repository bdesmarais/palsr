# Bootstrap uncertainty + multiple-imputation pooling --------------------------

#' Nonparametric bootstrap for PALS estimates and projections
#'
#' Quantifies uncertainty in PALS parameter estimates and in projected actor
#' locations by resampling events with replacement and re-estimating the model on
#' each bootstrap replicate, following Kim, Liu and Desmarais (2023). Each replicate
#' yields a parameter vector and (optionally) a set of Projected Actor Locations; the
#' collection of replicate PAL sets can be treated as multiple imputations and pooled
#' with Rubin's Rules (see [pool_rubin()]).
#'
#' @param events A [pal_events] object.
#' @param R Number of bootstrap replicates (default `50`; the paper uses `10`).
#' @param model `"four"` or `"one"` (passed to [estimate_pals]).
#' @param predict_time Optional `Date` (or vector of `Date`s). When supplied, every
#'   replicate also projects all actors at these times, so the spread of projected
#'   coordinates across replicates is available for confidence regions and pooling.
#' @param actors For projection, which actors to project (default: all in `events`).
#' @param seed Optional integer seed; replicate `r` uses `seed + r` so the whole run
#'   is reproducible.
#' @param ... Further arguments passed to [estimate_pals] (e.g. `aggregate`,
#'   `alter_weight`, `method`, `control`).
#'
#' @return An object of class `pals_boot` with components:
#'   \describe{
#'     \item{`estimates`}{An `R`-row `data.frame` of replicate parameter estimates.}
#'     \item{`estimate`}{The point estimate on the full sample (an [estimate_pals] fit).}
#'     \item{`projections`}{If `predict_time` was given, a `data.frame` of projected
#'       `lon`/`lat` for every actor-time-replicate combination; otherwise `NULL`.}
#'     \item{`R`, `model`, `call`}{Bookkeeping.}
#'   }
#'   Methods: [print()], [summary()] (bootstrap SEs / percentile intervals), and
#'   [coef()] (the full-sample point estimate).
#'
#' @details Resampling is over rows of `events` (the nonparametric event bootstrap).
#'   Duplicated events are kept as ordinary repeated events. Replicates whose
#'   optimizer fails to converge are retained but flagged via the `convergence`
#'   column of `estimates`.
#'
#' @seealso [estimate_pals()], [pool_rubin()].
#'
#' @examples
#' ev <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
#' bt <- bootstrap_pals(ev, R = 10, model = "one", seed = 1)
#' summary(bt)
#'
#' @export
bootstrap_pals <- function(events, R = 50, model = c("four", "one"),
                           predict_time = NULL, actors = NULL,
                           seed = NULL, ...) {
  stopifnot(inherits(events, "pal_events"))
  model <- match.arg(model)
  if (!is.numeric(R) || length(R) != 1 || R < 1)
    stop("`R` must be a positive integer.", call. = FALSE)
  R <- as.integer(R)
  if (!is.null(predict_time)) predict_time <- .as_date(predict_time)

  # Full-sample point estimate.
  point <- estimate_pals(events, model = model, ...)

  n <- nrow(events)
  est_list <- vector("list", R)
  proj_list <- if (is.null(predict_time)) NULL else vector("list", R)

  for (r in seq_len(R)) {
    if (!is.null(seed)) {
      old <- .Random.seed_safe(); on.exit(.restore_seed(old), add = TRUE)
      set.seed(seed + r)
    }
    idx <- sample.int(n, n, replace = TRUE)
    bev <- events[idx, , drop = FALSE]
    rownames(bev) <- NULL
    class(bev) <- c("pal_events", "data.frame")

    fit <- tryCatch(estimate_pals(bev, model = model, ...), error = function(e) NULL)
    if (is.null(fit)) {
      est_list[[r]] <- NULL
      next
    }
    cf <- coef(fit)
    est_list[[r]] <- data.frame(
      replicate = r, as.list(cf),
      objective = fit$objective, convergence = fit$convergence,
      check.names = FALSE, stringsAsFactors = FALSE
    )

    if (!is.null(predict_time)) {
      # Project from the ORIGINAL histories using the replicate's parameters, so the
      # spread reflects estimation uncertainty in the smoothing parameters.
      pj <- project_pals(events, actors = actors, predict_time = predict_time,
                         params = fit$params,
                         alter_weight = point$settings$alter_weight,
                         eps = point$settings$eps)
      pj$replicate <- r
      proj_list[[r]] <- pj
    }
  }

  estimates <- do.call(rbind, est_list[!vapply(est_list, is.null, logical(1))])
  rownames(estimates) <- NULL
  projections <- if (is.null(predict_time)) NULL else
    do.call(rbind, proj_list[!vapply(proj_list, is.null, logical(1))])

  structure(
    list(estimates = estimates, estimate = point, projections = projections,
         R = R, model = model, call = match.call()),
    class = "pals_boot"
  )
}

#' @export
print.pals_boot <- function(x, ...) {
  cat(sprintf("<pals_boot> %d replicates of the %s-parameter PALS model\n",
              x$R, x$model))
  cat("  point estimate:\n  ")
  print(x$estimate$params)
  ok <- sum(x$estimates$convergence == 0)
  cat(sprintf("  replicates that converged: %d / %d\n", ok, nrow(x$estimates)))
  invisible(x)
}

#' @export
coef.pals_boot <- function(object, ...) coef(object$estimate)

#' @export
summary.pals_boot <- function(object, conf = 0.95, ...) {
  pars <- names(coef(object$estimate))
  est <- object$estimates
  a <- (1 - conf) / 2
  tab <- do.call(rbind, lapply(pars, function(p) {
    v <- est[[p]]
    data.frame(
      term = p,
      estimate = coef(object$estimate)[[p]],
      boot_mean = mean(v),
      boot_se = stats::sd(v),
      lower = stats::quantile(v, a, names = FALSE),
      upper = stats::quantile(v, 1 - a, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(tab) <- NULL
  cat(sprintf("<pals_boot> summary (%d replicates, %.0f%% percentile intervals)\n",
              nrow(est), 100 * conf))
  print(tab, digits = 4)
  invisible(tab)
}

# Rubin's Rules ----------------------------------------------------------------

#' Pool estimates across imputations with Rubin's Rules
#'
#' Combines per-imputation point estimates and variances of a scalar quantity into a
#' single pooled estimate with a variance that accounts for both within- and
#' between-imputation uncertainty (Rubin, 1987). Use it to pool estimands computed on
#' each bootstrap/imputation replicate of [bootstrap_pals()] — for example a
#' regression coefficient from a dyadic model fit on each replicate's PAL distances.
#'
#' @param estimates Numeric vector of per-imputation point estimates \eqn{Q_j}.
#' @param variances Numeric vector of per-imputation variances \eqn{U_j}
#'   (the squared standard errors), the same length as `estimates`.
#' @param df Logical; if `TRUE`, also return the Barnard-Rubin adjusted degrees of
#'   freedom and a corresponding two-sided p-value for `H0: Q = 0`. Default `FALSE`
#'   reproduces the source code's normal-based pooling.
#' @param dfcom Complete-data degrees of freedom, used only when `df = TRUE`
#'   (default `Inf`, the large-sample limit).
#'
#' @return A one-row `data.frame` with the pooled estimate `qbar`, within-imputation
#'   variance `ubar`, between-imputation variance `b`, total variance `t`, standard
#'   error `se`, fraction of missing information `fmi`, and (if `df = TRUE`) `df` and
#'   `p.value`.
#'
#' @details With \eqn{m} imputations,
#'   \deqn{\bar Q = \tfrac1m \sum_j Q_j,\quad \bar U = \tfrac1m \sum_j U_j,\quad
#'         B = \tfrac{1}{m-1}\sum_j (Q_j-\bar Q)^2,}
#'   and total variance \eqn{T = \bar U + (1 + 1/m) B}. The fraction of missing
#'   information is \eqn{(1 + 1/m)B / T}. When `df = TRUE`, the Barnard-Rubin (1999)
#'   small-sample degrees of freedom are used.
#'
#' @references
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
#'
#' Barnard, J. and Rubin, D. B. (1999). Small-sample degrees of freedom with multiple
#' imputation. *Biometrika*, 86(4), 948-955.
#'
#' @seealso [bootstrap_pals()].
#'
#' @examples
#' # Five imputations of a coefficient and its variance.
#' q <- c(1.10, 0.95, 1.20, 1.05, 0.98)
#' u <- c(0.04, 0.05, 0.045, 0.038, 0.052)
#' pool_rubin(q, u)
#' pool_rubin(q, u, df = TRUE, dfcom = 100)
#'
#' @export
pool_rubin <- function(estimates, variances, df = FALSE, dfcom = Inf) {
  estimates <- as.numeric(estimates)
  variances <- as.numeric(variances)
  if (length(estimates) != length(variances))
    stop("`estimates` and `variances` must have the same length.", call. = FALSE)
  m <- length(estimates)
  if (m < 2) stop("Need at least 2 imputations to pool.", call. = FALSE)
  if (any(variances < 0, na.rm = TRUE))
    stop("`variances` must be non-negative.", call. = FALSE)

  qbar <- mean(estimates)
  ubar <- mean(variances)
  b    <- stats::var(estimates)              # denominator m - 1
  t    <- ubar + (1 + 1 / m) * b
  se   <- sqrt(t)
  fmi  <- ((1 + 1 / m) * b) / t

  out <- data.frame(qbar = qbar, ubar = ubar, b = b, t = t, se = se, fmi = fmi)

  if (df) {
    # Barnard-Rubin (1999) adjusted degrees of freedom.
    lambda <- ((1 + 1 / m) * b) / t
    lambda <- max(lambda, 1e-04)
    dfold <- (m - 1) / lambda^2
    dfobs <- if (is.finite(dfcom))
      (dfcom + 1) / (dfcom + 3) * dfcom * (1 - lambda) else Inf
    nu <- if (is.finite(dfobs)) dfold * dfobs / (dfold + dfobs) else dfold
    out$df <- nu
    out$p.value <- 2 * stats::pt(-abs(qbar / se), df = nu)
  }
  out
}
