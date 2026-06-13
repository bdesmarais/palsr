---
title: 'palsr: Projected Actor Locations for spatial modeling of dyadic interactions between moving actors'
tags:
  - R
  - political science
  - spatial statistics
  - conflict
  - dyadic data
  - exponential smoothing
authors:
  - name: Bruce A. Desmarais
    corresponding: true
    affiliation: 1
  - name: Sangyeon Kim
    affiliation: 2
  - name: Howard Liu
    affiliation: 3
affiliations:
  - name: Department of Political Science, Pennsylvania State University, USA
    index: 1
  - name: Pennsylvania State University, USA
    index: 2
  - name: University of Southern Denmark, Denmark
    index: 3
date: 13 June 2026
bibliography: paper.bib
---

# Summary

Many actors studied in the social sciences have no fixed geographic location.
Armed groups roam across a theater of conflict, firms relocate, diplomats and
international organizations shift their operational focus, and migrating
populations move in response to events. For such *moving actors*, the location
at which two of them interact is itself a meaningful, modelable outcome, but
standard spatial-interaction tools assume actors sit at fixed points and so
cannot be applied directly. The Projected Actor Location (PALS) method, introduced
by @kim2023pals, fills this gap by estimating where a mobile actor effectively
"is" at any point in time from the spatiotemporal history of its past
interactions, and using these projected locations to model the distance between
--- and the probability of interaction among --- pairs of actors.

`palsr` is an R package that implements the complete PALS workflow: constructing
and validating dyadic event data, estimating the smoothing parameters by
minimizing great-circle prediction error, projecting actor locations at
arbitrary times, building the dyadic distance covariates used in downstream
interaction models, and quantifying uncertainty through a nonparametric
bootstrap with multiple-imputation (Rubin's Rules) pooling. Performance-critical
kernels --- the Haversine distance and the projection smoother --- are
implemented in C++ through `Rcpp` [@rcpp], so that parameter estimation, which
repeatedly re-projects every actor across a grid of times, runs quickly even on
data sets with thousands of events.

# Statement of need

The empirical study of geopolitical interaction --- alliances, conflict,
cooperation, trade --- is overwhelmingly *dyadic*: the unit of analysis is a
pair of actors, and a central predictor of whether and how two actors interact
is the distance between them [@gleditsch2001; @ward2007]. When actors are
states with fixed capitals or centroids, this distance is trivial to compute.
But a large and growing share of substantively important actors --- rebel
groups, militias, terrorist organizations, peacekeeping deployments,
non-governmental organizations --- are mobile, and treating them as
fixed-location entities discards exactly the spatial dynamics that drive their
interactions. Prior to PALS there was no general, estimable method for assigning
time-varying locations to such actors in a way that is optimized for predicting
their interactions.

PALS addresses this need, and the original article [@kim2023pals] demonstrates
that PALS-based distances substantially improve the prediction of subnational
conflict relative to naive location measures. However, that contribution was
accompanied only by replication scripts specific to a single application.
`palsr` turns the method into reusable, documented, and tested research software.
It lets applied researchers in political science, conflict studies, economics,
and geography apply PALS to their own dyadic event data without re-deriving the
estimator, provides a principled bootstrap procedure for propagating estimation
uncertainty into downstream models, and supplies a simulated example data set so
that the full workflow can be learned and reproduced without access to
restricted conflict data. To our knowledge it is the first general-purpose
software implementation of projected-location modeling for moving actors.

# Method

For a focal actor $i$ at prediction time $t$, PALS forms a recency-weighted mean
of the locations of $i$'s own past events (the *focal* component) and a
recency-weighted mean of the locations of events involving $i$'s past
interaction partners, or *alters* (the *alter* component). The projected
location is a convex combination of the two,

$$
g_i(t) = (1-\pi)\sum_e W_i(e)\, g(e) \;+\; \pi \sum_e W_k(e)\, g(e),
\qquad \pi = \mathrm{logistic}(\gamma + \eta\, v),
$$

where event weights decay with age --- governed by $\alpha$ for the focal actor
and $\beta$ for the alters --- and the mixing weight $\pi$ depends through
$\gamma$ and $\eta$ on the focal actor's activity relative to its alters,
summarized by the event-count statistic $v$. The four parameters are estimated
by "marching forward" through time: every event is predicted using only events
strictly preceding it, and parameters are chosen to minimize the mean Haversine
distance between predicted and observed interaction locations. The package also
exposes a parsimonious one-parameter variant that fixes $\pi = 0$ and estimates
only $\alpha$, which is fast and frequently competitive.

# Key features

- **Validated data container** (`pal_events()`) for dyadic, time-stamped,
  georeferenced events.
- **Parameter estimation** (`estimate_pals()`) for both the full four-parameter
  and the reduced one-parameter models, via Haversine-error minimization with a
  compiled objective.
- **Projection** (`project_pal()`, `project_pals()`) of actor locations at any
  set of times.
- **Prediction and covariate construction** (`predict_event_locations()`,
  `pal_distance()`) for interaction locations and dyadic distances, including a
  log transform for conflict-style specifications.
- **Uncertainty quantification** (`bootstrap_pals()`, `pool_rubin()`) via a
  nonparametric event bootstrap and Rubin's Rules [@rubin1987].
- **Reproducible example data** (`nigeria_sim`, `simulate_conflict_events()`),
  a vignette, and a full `testthat` suite that also checks the C++ kernels
  against pure-R reference implementations.

# Example

```r
library(palsr)
data(nigeria_sim)                       # 1,500 dyadic events, 25 mobile actors

fit <- estimate_pals(nigeria_sim, model = "one")   # estimate alpha
coef(fit)

# Project where each actor is on a given date.
project_pals(nigeria_sim, predict_time = as.Date("2015-01-01"), params = fit)

# Build the dyadic distance covariate for an interaction model.
dyads <- data.frame(actor1 = "G01", actor2 = "G02", time = as.Date("2014-06-01"))
pal_distance(nigeria_sim, dyads, fit, transform = "log")

# Bootstrap uncertainty in the smoothing parameter.
summary(bootstrap_pals(nigeria_sim, R = 10, model = "one", seed = 1))
```

# Acknowledgements

We thank the editors and reviewers of *Political Science Research and Methods*
for feedback on the method, and the developers of `Rcpp` for the tooling that
underpins the package's performance.

# References
