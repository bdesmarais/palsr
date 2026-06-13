# pals

<!-- badges: start -->
[![R-CMD-check](https://github.com/bdesmarais/pals/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bdesmarais/pals/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**pals** implements the **Projected Actor Location (PALS)** method for spatial
modeling of dyadic interactions between geographically mobile actors, as
described in Kim, Liu and Desmarais (2023),
[*Spatial modeling of dyadic geopolitical interactions between moving actors*](https://doi.org/10.1017/psrm.2022.6),
*Political Science Research and Methods* **11**(3).

Actors such as armed groups, firms, or diplomats have no fixed location: they
move through space over time, and *where* two such actors interact is itself an
outcome worth modeling. PALS projects where a mobile actor "is" at any moment
from the spatiotemporal history of its past interactions, using
exponential-smoothing weights that favour recent events, and combines a focal
actor's own history with that of its interaction partners ("alters").

## Installation

```r
# install.packages("remotes")
remotes::install_github("bdesmarais/pals")
```

The package compiles a small amount of C++ via Rcpp, so a working toolchain
(Rtools on Windows, Xcode command-line tools on macOS) is required.

## The method

For a focal actor *i* at time *t*, the projected location is a convex
combination of a recency-weighted mean of *i*'s own past event locations and a
recency-weighted mean of its alters' event locations:

```
g_i(t) = (1 - π) Σ W_i(e) g(e)  +  π Σ W_k(e) g(e),   π = logistic(γ + η·v)
```

with four smoothing parameters:

| Parameter | Role                                                |
|-----------|-----------------------------------------------------|
| `alpha`   | time-decay of the focal actor's own history         |
| `beta`    | time-decay of the alters' histories                 |
| `gamma`   | intercept of the focal-vs-alter mixing weight `π`   |
| `eta`     | dependence of `π` on relative event activity `v`    |

A reduced **one-parameter** model fixes `π = 0` (focal history only) and
estimates `alpha` alone.

## Quick start

```r
library(pals)

# Bundled simulated dataset: 1,500 dyadic events, 25 mobile actors, 2000-2016.
data(nigeria_sim)

# Estimate the one-parameter model (marching forward in time, minimizing
# great-circle prediction error).
fit <- estimate_pals(nigeria_sim, model = "one")
fit
#> <pals_fit> one-parameter PALS model
#> <pals_params> (one-parameter model)
#>   alpha = 0.874  (pi fixed at 0)
#>   objective (mean Haversine km): 103.37 over 1500 events

# Project where each actor is on a given date.
project_pals(nigeria_sim, predict_time = as.Date("2015-01-01"), params = fit)

# Predict interaction locations and score them in kilometres.
predict_event_locations(nigeria_sim, nigeria_sim, fit)

# Build the dyadic distance covariate (optionally log-transformed).
dyads <- data.frame(actor1 = "G01", actor2 = "G02", time = as.Date("2014-06-01"))
pal_distance(nigeria_sim, dyads, fit, transform = "log")
```

The full four-parameter model adds the alter component:

```r
fit4 <- estimate_pals(nigeria_sim, model = "four")
coef(fit4)
```

## Uncertainty

`bootstrap_pals()` resamples events with replacement and re-estimates the model
on each replicate, giving bootstrap standard errors and percentile intervals.
Downstream estimands computed per replicate can be pooled with Rubin's Rules
via `pool_rubin()`:

```r
bt <- bootstrap_pals(nigeria_sim, R = 10, model = "one", seed = 1)
summary(bt)

# Pool five per-replicate estimates (q) and their variances (u).
pool_rubin(q = c(1.10, 0.95, 1.20, 1.05, 0.98),
           u = c(0.04, 0.05, 0.045, 0.038, 0.052))
```

## Learn more

See the vignette for a full worked example:

```r
vignette("pals")
```

## Key functions

| Function                     | Purpose                                            |
|------------------------------|----------------------------------------------------|
| `pal_events()`               | Build a validated dyadic-event table               |
| `simulate_conflict_events()` | Simulate mobile-actor interaction data             |
| `estimate_pals()`            | Estimate parameters by Haversine error minimization|
| `project_pal()` / `project_pals()` | Project actor locations                      |
| `predict_event_locations()`  | Predict dyadic interaction locations               |
| `pal_distance()`             | Dyadic distance covariate between projected locations |
| `bootstrap_pals()`           | Nonparametric bootstrap uncertainty                |
| `pool_rubin()`               | Pool replicate estimates with Rubin's Rules        |

## Citation

If you use **pals**, please cite the method paper:

> Kim, S., Liu, H., & Desmarais, B. A. (2023). Spatial modeling of dyadic
> geopolitical interactions between moving actors. *Political Science Research
> and Methods*, 11(3), 617–635. https://doi.org/10.1017/psrm.2022.6

## License

MIT © the package authors.
