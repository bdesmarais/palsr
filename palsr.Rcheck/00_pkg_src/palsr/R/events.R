# Event container ---------------------------------------------------------------

#' Construct a validated dyadic-event table
#'
#' `pal_events()` builds the core data object used throughout \pkg{palsr}: a table of
#' dyadic interaction events, each involving two actors at a known time and location.
#' Projected Actor Locations are computed from these histories.
#'
#' @param data A `data.frame` with one row per dyadic event.
#' @param actor1,actor2 Column names (length-1 character) identifying the two actors
#'   involved in each event. The pair is treated as unordered.
#' @param time Column name of the event time. Must be a `Date`, or coercible to one via
#'   [as.Date()].
#' @param lon,lat Column names of the event longitude and latitude, in decimal degrees.
#' @param drop_self Logical; drop events whose two actors are identical (default `TRUE`).
#'
#' @return An object of class `pal_events` (a `data.frame` subclass) with canonical
#'   columns `actor1`, `actor2`, `time`, `lon`, `lat`, sorted by `time`.
#'
#' @details Longitudes must lie in `[-180, 180]` and latitudes in `[-90, 90]`. Rows with
#'   missing actor, time, or coordinate values are dropped with a message.
#'
#' @examples
#' df <- data.frame(
#'   from = c("A", "A", "B"),
#'   to   = c("B", "C", "C"),
#'   when = as.Date(c("2001-01-01", "2001-06-01", "2002-01-01")),
#'   x    = c(7.1, 8.0, 7.5),
#'   y    = c(9.0, 9.4, 10.1)
#' )
#' ev <- pal_events(df, actor1 = "from", actor2 = "to",
#'                  time = "when", lon = "x", lat = "y")
#' ev
#'
#' @export
pal_events <- function(data, actor1 = "actor1", actor2 = "actor2",
                       time = "time", lon = "lon", lat = "lat",
                       drop_self = TRUE) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  cols <- c(actor1, actor2, time, lon, lat)
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols)) {
    stop("Column(s) not found in `data`: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  out <- data.frame(
    actor1 = as.character(data[[actor1]]),
    actor2 = as.character(data[[actor2]]),
    time   = .as_date(data[[time]]),
    lon    = as.numeric(data[[lon]]),
    lat    = as.numeric(data[[lat]]),
    stringsAsFactors = FALSE
  )

  complete <- stats::complete.cases(out)
  if (any(!complete)) {
    message(sprintf("Dropping %d event(s) with missing values.", sum(!complete)))
    out <- out[complete, , drop = FALSE]
  }

  if (any(out$lon < -180 | out$lon > 180, na.rm = TRUE))
    stop("`lon` must be within [-180, 180].", call. = FALSE)
  if (any(out$lat < -90 | out$lat > 90, na.rm = TRUE))
    stop("`lat` must be within [-90, 90].", call. = FALSE)

  if (drop_self) {
    self <- out$actor1 == out$actor2
    if (any(self)) {
      message(sprintf("Dropping %d self-dyad event(s).", sum(self)))
      out <- out[!self, , drop = FALSE]
    }
  }

  if (nrow(out) == 0) stop("No valid events remain after validation.", call. = FALSE)

  out <- out[order(out$time), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("pal_events", "data.frame")
  out
}

#' @noRd
.as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(x)
}

#' @noRd
.actors <- function(events) sort(unique(c(events$actor1, events$actor2)))

#' @export
print.pal_events <- function(x, ...) {
  cat("<pal_events>\n")
  cat(sprintf("  events: %d\n", nrow(x)))
  cat(sprintf("  actors: %d\n", length(.actors(x))))
  if (nrow(x)) {
    cat(sprintf("  time:   %s to %s\n",
                format(min(x$time)), format(max(x$time))))
    cat(sprintf("  bbox:   lon [%.2f, %.2f], lat [%.2f, %.2f]\n",
                min(x$lon), max(x$lon), min(x$lat), max(x$lat)))
  }
  invisible(x)
}

#' @export
summary.pal_events <- function(object, ...) {
  acts <- .actors(object)
  tab <- table(c(object$actor1, object$actor2))
  cat("<pal_events> summary\n")
  cat(sprintf("  %d events, %d actors, %s to %s\n",
              nrow(object), length(acts),
              format(min(object$time)), format(max(object$time))))
  cat(sprintf("  events per actor: min %d, median %g, max %d\n",
              min(tab), stats::median(tab), max(tab)))
  invisible(object)
}
