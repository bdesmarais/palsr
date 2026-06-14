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
  - name: Sangyeon Kim
    affiliation: 1
  - name: Howard Liu
    affiliation: 2
  - name: Bruce A. Desmarais
    corresponding: true
    affiliation: 3
affiliations:
  - name: Division of Communication and Media, Ewha Womans University, Republic of Korea
    index: 1
  - name: Department of Political Science, University of South Carolina, USA
    index: 2
  - name: Department of Political Science, Pennsylvania State University, USA
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

`palsr` is an R package [@rcoreteam] that implements the complete PALS workflow:
constructing and validating dyadic event data, estimating the smoothing parameters by
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
non-governmental organizations --- are mobile. The fine-grained, geocoded event
data that now document them, such as the Armed Conflict Location and Event Data
project [@acled], record *where each interaction happened* but assign actors no
persistent coordinates. Treating these actors as fixed-location entities
discards exactly the spatial dynamics that drive their interactions.

Existing spatial software for R does not fill this gap. General geospatial
toolkits such as `sf` [@sf] and spatial point-process libraries such as
`spatstat` [@spatstat] are powerful, but they operate on objects with given
coordinates; they provide no way to *infer* a moving actor's effective location,
let alone one tuned to predict future interactions. Conversely, dyadic
interaction models in political science and economics typically take inter-actor
distance as an exogenous input. Prior to PALS there was thus no general,
estimable method for assigning time-varying locations to mobile actors in a way
that is optimized for predicting their interactions.

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

where event weights decay with age in the manner of exponential smoothing
[@gardner2006] --- governed by $\alpha$ for the focal actor and $\beta$ for the
alters --- and the mixing weight $\pi$ depends through $\gamma$ and $\eta$ on the
focal actor's activity relative to its alters, summarized by the event-count
statistic $v$. The four parameters are estimated by "marching forward" through
time: every event is predicted using only events strictly preceding it, and
parameters are chosen to minimize the mean Haversine distance between predicted
and observed interaction locations. The package also exposes a parsimonious
one-parameter variant that fixes $\pi = 0$ and estimates only $\alpha$, which is
fast and frequently competitive.

Because the projected location is recomputed as time advances, each actor traces
a *trajectory* through space rather than occupying a fixed point.
\autoref{fig:traj} illustrates this for four actors in the bundled simulated
data: PALS projects each actor's location at yearly intervals, and the resulting
paths drift through the cloud of observed events.

![Projected locations of four actors in the simulated `nigeria_sim` data, fit with the one-parameter model and projected at yearly intervals from 2005 to 2016. Grey points are the observed dyadic events; coloured paths trace each actor's PALS-projected trajectory (arrows point forward in time). \label{fig:traj}](figure-trajectories.png){ width=85% }

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

# Availability and quality control

`palsr` is developed openly on GitHub. The package ships with documentation for
every exported function, an introductory vignette that walks through the full
workflow on `nigeria_sim`, and a `testthat` suite covering data validation,
estimation, projection, distance construction, and uncertainty quantification.
The tests additionally verify the compiled C++ kernels against independent
pure-R reference implementations, and continuous integration runs `R CMD check`
across Linux, macOS, and Windows on every change.

# Acknowledgements

We thank the editors and reviewers of *Political Science Research and Methods*
for feedback on the method, and the developers of `Rcpp` [@rcpp] for the tooling
that underpins the package's performance.

# References
