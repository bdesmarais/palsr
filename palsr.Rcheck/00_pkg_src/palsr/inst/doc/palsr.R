## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5,
  dpi = 96
)
set.seed(1)

## ----load---------------------------------------------------------------------
library(palsr)

## ----data---------------------------------------------------------------------
data(nigeria_sim)
nigeria_sim
summary(nigeria_sim)

## ----build--------------------------------------------------------------------
raw <- data.frame(
  from = c("A", "A", "B"),
  to   = c("B", "C", "C"),
  when = as.Date(c("2001-01-01", "2001-06-01", "2002-01-01")),
  x    = c(7.1, 8.0, 7.5),
  y    = c(9.0, 9.4, 10.1)
)
pal_events(raw, actor1 = "from", actor2 = "to",
           time = "when", lon = "x", lat = "y")

## ----fit-one------------------------------------------------------------------
fit1 <- estimate_pals(nigeria_sim, model = "one")
fit1
coef(fit1)

## ----fit-four-----------------------------------------------------------------
fit4 <- estimate_pals(nigeria_sim, model = "four",
                      control = list(maxit = 60))
coef(fit4)

## ----project------------------------------------------------------------------
pal_2015 <- project_pals(nigeria_sim,
                         predict_time = as.Date("2015-01-01"),
                         params = fit1)
head(pal_2015)

## ----map, fig.alt = "Projected actor locations on 2015-01-01"-----------------
library(ggplot2)
ggplot(pal_2015, aes(lon, lat)) +
  geom_point(colour = "#2b6cb0", size = 2) +
  geom_text(aes(label = actor), vjust = -0.8, size = 3) +
  labs(title = "Projected actor locations, 2015-01-01",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

## ----predict-events-----------------------------------------------------------
targets <- nigeria_sim[nigeria_sim$time > as.Date("2014-01-01"), ]
scored  <- predict_event_locations(nigeria_sim, targets, fit1)
summary(scored$error_km)

## ----distance-----------------------------------------------------------------
dyads <- data.frame(actor1 = "G01", actor2 = "G02",
                    time = as.Date("2014-06-01"))
pal_distance(nigeria_sim, dyads, fit1, transform = "log")

## ----bootstrap----------------------------------------------------------------
bt <- bootstrap_pals(nigeria_sim, R = 10, model = "one", seed = 1)
summary(bt)

## ----rubin--------------------------------------------------------------------
q <- c(1.10, 0.95, 1.20, 1.05, 0.98)   # per-replicate estimates
u <- c(0.04, 0.05, 0.045, 0.038, 0.052) # per-replicate variances
pool_rubin(q, u, df = TRUE, dfcom = 100)

