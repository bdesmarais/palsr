# Generate the figure used in the JOSS paper (paper/figure-trajectories.png).
#
# Illustrates the core idea of PALS: actors are *mobile*, so their projected
# locations trace a trajectory through space over time. We fit the one-parameter
# model to the bundled `nigeria_acled` ACLED data and project the most active
# actors at yearly intervals, drawing each actor's projected path over the cloud
# of observed conflict events.
#
# Run from the package root with the package loaded:
#   Rscript paper/make_figure.R

suppressMessages({
  if (requireNamespace("palsr", quietly = TRUE)) {
    library(palsr)
  } else {
    devtools::load_all(quiet = TRUE)   # development fallback
  }
  library(ggplot2)
})

data(nigeria_acled)
fit <- estimate_pals(nigeria_acled, model = "one")

# Four of the most active actors; shorten the long ACLED labels for the legend.
tab    <- sort(table(c(nigeria_acled$actor1, nigeria_acled$actor2)), decreasing = TRUE)
actors <- names(tab)[1:4]
short  <- function(a) sub(" \\(Nigeria\\)$", "", a)

dates <- as.Date(sprintf("%d-12-01", seq(2009, 2016, by = 1)))
traj  <- project_pals(nigeria_acled, actors = actors, predict_time = dates, params = fit)
traj  <- traj[!is.na(traj$lon), ]
traj$actor <- factor(short(traj$actor), levels = short(actors))

p <- ggplot() +
  geom_point(data = nigeria_acled, aes(lon, lat),
             colour = "grey80", size = 0.5, alpha = 0.5) +
  geom_path(data = traj, aes(lon, lat, colour = actor),
            linewidth = 0.8, lineend = "round",
            arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed")) +
  geom_point(data = traj, aes(lon, lat, colour = actor), size = 1.6) +
  scale_colour_brewer(palette = "Dark2", name = "Actor") +
  labs(x = "Longitude", y = "Latitude") +
  coord_quickmap() +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("paper/figure-trajectories.png", p, width = 7.2, height = 4.2,
       dpi = 150, bg = "white")
cat("wrote paper/figure-trajectories.png; actors:",
    paste(short(actors), collapse = ", "), "\n")
