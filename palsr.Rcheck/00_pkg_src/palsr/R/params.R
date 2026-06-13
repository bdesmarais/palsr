# PALS parameter set ------------------------------------------------------------

#' Create a PALS parameter set
#'
#' A lightweight container for the four PALS smoothing parameters. Use it to project
#' actor locations with known parameters (e.g. values reported in a paper), without
#' estimating them from data.
#'
#' @param alpha Time-decay of the focal actor's own event history (`>= 0`). Larger
#'   values down-weight older events more steeply.
#' @param beta Time-decay of the alters' event histories (`>= 0`). Ignored in the
#'   one-parameter model.
#' @param gamma Intercept of the logistic mixing weight \eqn{\pi}. Higher values place
#'   more weight on the alters' average location. Ignored in the one-parameter model.
#' @param eta Slope of the logistic mixing weight on the event-count ratio. Ignored in
#'   the one-parameter model.
#' @param model Either `"four"` (full model: focal + alter histories) or `"one"`
#'   (focal-only; \eqn{\pi = 0}, so only `alpha` is used).
#'
#' @return An object of class `pals_params`.
#'
#' @details The mixing weight is \eqn{\pi = \mathrm{plogis}(\gamma + \eta v)}, where
#'   \eqn{v = (n_i / n_k)^{1/\sqrt{n_k}}} compares the number of focal events
#'   (\eqn{n_i}) with the number of alter events (\eqn{n_k}). The projected location is
#'   \eqn{(1-\pi)} times the recency-weighted mean of the focal actor's own past event
#'   locations plus \eqn{\pi} times the recency-weighted mean of its alters' locations.
#'   See [project_pals()] and the package vignette.
#'
#' @examples
#' p <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
#' p
#' pals_params(alpha = 0.9, model = "one")
#'
#' @export
pals_params <- function(alpha, beta = 0, gamma = 0, eta = 0,
                        model = c("four", "one")) {
  model <- match.arg(model)
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0)
    stop("`alpha` must be a single non-negative number.", call. = FALSE)
  if (model == "four") {
    if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta < 0)
      stop("`beta` must be a single non-negative number.", call. = FALSE)
    for (nm in c("gamma", "eta")) {
      v <- get(nm)
      if (!is.numeric(v) || length(v) != 1 || is.na(v))
        stop(sprintf("`%s` must be a single number.", nm), call. = FALSE)
    }
  }
  structure(
    list(alpha = alpha, beta = beta, gamma = gamma, eta = eta, model = model),
    class = "pals_params"
  )
}

#' @export
print.pals_params <- function(x, ...) {
  cat(sprintf("<pals_params> (%s-parameter model)\n", x$model))
  if (x$model == "one") {
    cat(sprintf("  alpha = %.4g  (pi fixed at 0)\n", x$alpha))
  } else {
    cat(sprintf("  alpha = %.4g, beta = %.4g, gamma = %.4g, eta = %.4g\n",
                x$alpha, x$beta, x$gamma, x$eta))
  }
  invisible(x)
}

#' @noRd
as_pals_params <- function(x) {
  if (inherits(x, "pals_params")) return(x)
  if (inherits(x, "pals_fit")) return(x$params)
  stop("Expected a `pals_params` or `pals_fit` object.", call. = FALSE)
}
