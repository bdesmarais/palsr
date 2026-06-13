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
