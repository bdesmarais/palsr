# palsr 0.1.1

* New `rolling_origin_pals()`: a genuine out-of-sample (rolling-origin)
  evaluation, re-estimating the model at each of a set of origin dates from
  strictly earlier events and scoring predictions on the following holdout
  window. Complements the one-step-ahead ("marching forward") criterion that
  `estimate_pals()` minimizes, in which every event's own location
  contributes to the single fitted parameter vector.
* `bootstrap_pals()` gains a `block` argument for a block bootstrap
  (`block = "year"`, or a user-supplied grouping vector) that preserves
  within-block temporal/actor dependence, as an alternative to the default
  i.i.d. event bootstrap.
* README: fixed the `pool_rubin()` quick-start example, which errored with
  "unused arguments" when called with `q =`/`u =` (the function's arguments
  are named `estimates`/`variances`); updated install instructions to lead
  with CRAN and to add `build_vignettes = TRUE` to the GitHub install
  example, so `vignette("palsr")` is available after following the README.

# palsr 0.1.0

* Initial release implementing the Projected Actor Location (PALS) method of
  Kim, Liu and Desmarais (2023) <doi:10.1017/psrm.2022.6>.
* Core workflow: `pal_events()`, `estimate_pals()` (one- and four-parameter
  models), `project_pal()` / `project_pals()`, `predict_event_locations()`, and
  `pal_distance()`.
* Uncertainty quantification via `bootstrap_pals()` and Rubin's Rules pooling
  with `pool_rubin()`.
* Fast great-circle distance and projection kernels in C++ via Rcpp.
* Bundled simulated example dataset `nigeria_sim` and
  `simulate_conflict_events()`.
* Introductory vignette and full function documentation.
