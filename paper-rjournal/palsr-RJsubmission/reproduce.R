# Reproduces the numerical results and figure in the article
#   "palsr: Projecting the Locations of Moving Actors for Spatial Models of
#    Dyadic Interaction"
#
# Requires the palsr package (on CRAN) and ggplot2.
#   install.packages("palsr")

library(palsr)
library(ggplot2)

data(nigeria_acled)

## One-parameter model -------------------------------------------------------
fit1 <- estimate_pals(nigeria_acled, model = "one")
coef(fit1)                 # alpha ~ 1.241
fit1$objective             # mean Haversine error ~ 147 km

## Full four-parameter model -------------------------------------------------
fit4 <- estimate_pals(nigeria_acled, model = "four")
coef(fit4)                 # gamma is large and negative => pi ~ 0

## Projected trajectories (Figure 1) -----------------------------------------
actors <- c("Boko Haram", "Fulani Ethnic Militia (Nigeria)",
            "Ijaw Ethnic Militia (Nigeria)",
            "MEND: Movement for the Emancipation of the Niger Delta")
dates  <- as.Date(sprintf("%d-01-01", 2005:2016))
traj   <- project_pals(nigeria_acled, actors = actors,
                       predict_time = dates, params = fit1)
traj   <- traj[!is.na(traj$lon), ]
ends   <- do.call(rbind, lapply(split(traj, traj$actor),
                                function(d) d[which.max(d$time), ]))

ggplot() +
  geom_point(data = nigeria_acled, aes(lon, lat),
             colour = "grey80", size = 0.5, alpha = 0.5) +
  geom_path(data = traj, aes(lon, lat, colour = actor), linewidth = 0.8,
            arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed")) +
  geom_point(data = traj, aes(lon, lat, colour = actor), size = 1.5) +
  coord_quickmap() +
  theme_minimal()

## Dyadic distance covariate -------------------------------------------------
dyads <- data.frame(actor1 = "Boko Haram",
                    actor2 = "Fulani Ethnic Militia (Nigeria)",
                    time   = as.Date("2014-06-01"))
pal_distance(nigeria_acled, dyads, fit1, transform = "log")  # ~445.6 km

## Bootstrap uncertainty -----------------------------------------------------
bt <- bootstrap_pals(nigeria_acled, R = 20, model = "one", seed = 1)
summary(bt)                # 95% interval for alpha ~ [1.08, 1.45]
