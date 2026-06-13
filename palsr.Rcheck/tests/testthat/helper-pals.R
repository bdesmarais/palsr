# Small, fast fixtures shared across test files.

# A tiny hand-built event table with a known structure.
make_toy_events <- function() {
  df <- data.frame(
    from = c("A", "A", "B", "B", "C", "A"),
    to   = c("B", "C", "C", "A", "B", "B"),
    when = as.Date(c("2001-01-01", "2001-03-01", "2001-06-01",
                     "2001-09-01", "2002-01-01", "2002-03-01")),
    x    = c(7.0, 7.5, 8.0, 7.2, 8.4, 7.1),
    y    = c(9.0, 9.4, 10.1, 9.2, 10.5, 9.1)
  )
  pal_events(df, actor1 = "from", actor2 = "to",
             time = "when", lon = "x", lat = "y")
}

# A small simulated set that is quick to estimate on.
make_sim_events <- function(n_actors = 8, n_events = 150, seed = 1) {
  simulate_conflict_events(n_actors = n_actors, n_events = n_events, seed = seed)
}
