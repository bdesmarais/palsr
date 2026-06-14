# Design Decisions & Rationale

Numbered decisions, each with the question, the choice, and why. Referenced from
`ALGORITHM.md` and `DESIGN.md`.

1. **Alter-weight normalization.** The Dataverse code overwrites the alter weight vector
   with a scalar (`weights_k <- 1/sum(weights_k)`), computing an un-normalized sum of
   alter coordinates rather than the paper's weighted average `W_k(e)=w_k/Σw_k`.
   → **Implement the normalized weighted average as the default** (`alter_weight =
   "normalized"`); provide `alter_weight = "legacy"` for byte-exact reproduction.
   *Why:* the normalized form is what the paper's equations specify and is statistically
   correct; estimated `pi≈0` in the application makes the two numerically equivalent
   there, so we lose nothing and gain correctness + an explicit replication switch.

2. **Objective: mean vs sum of distances.** Paper says "sum"; code minimizes the mean km.
   → **Default `aggregate = "mean"` (km)**, `"sum"` available. *Why:* matches released
   code and is scale-stable when the count of non-NA events varies across parameters.

3. **Optimizer.** Source used `optimParallel` L-BFGS-B across 16 cores ("hill-climbing").
   → **Use base `stats::optim`** (Nelder–Mead for 4-param, Brent for 1-param) with a
   `method` override and optional restarts. *Why:* zero extra dependencies, reproducible,
   adequate for the low-dimensional smooth objective; users can pass `method="L-BFGS-B"`.

4. **Positivity of alpha/beta.** → Optimize on **log scale** (`alpha=exp(.)`), exactly as
   the source code, so no box constraints are needed and estimates match.

5. **Prediction date + history cutoff.** Source fixes Dec-1 of the target year and
   truncates history at the calendar-year boundary (`filter(YEAR <= year-1)`). → **Expose
   `predict_time` as a free argument** plus a **`cutoff = c("day","month","year")`** option
   on the history filter: `"day"` (default) = events strictly before `predict_time`;
   `"year"` = events before Jan 1 of the prediction year, reproducing the source's
   convention (verified to machine precision against the Dataverse `actorYearDF` outputs,
   571/571 actor-years). *Why:* generality without losing replicability; `cutoff` subsumes
   the originally-planned `estimate_pals_yearly()` helper.

6. **Package name collision.** CRAN has an unrelated `pals` (colour palettes). → **Keep
   `pals`** per the project brief; JOSS does not require CRAN. Documented prominently.
   If CRAN submission is later desired, rename (e.g. `palsr`) is a one-shot change.
   *Why:* honor the brief; flag the constraint rather than silently renaming.

7. **Age units / timescale.** → **Days**, as in source. Internally `as.numeric(t - time)`
   on Dates. `alpha/beta` are therefore on a per-day decay scale (document the
   interpretation, incl. the paper's "2-yr event ≈ 0.53× a 1-yr event" example).

8. **Zero-alter fallback.** `v=(n_i/0)^…` is undefined. → When `n_k==0`, set `pi=0`
   (focal-only) and warn at most once. *Why:* avoids `NaN` PALs; the only sensible blend.

9. **No-history → NA, not (0,0).** → Track a `has_history` flag; return `NA` coordinates
   when `n_i==0`. *Why:* the source maps predicted-0 to NA, which would wrongly null a
   genuine (0,0); an explicit flag is robust.

10. **Haversine in C++.** → Implement great-circle distance in Rcpp (mean radius
    6371.0088 km) instead of depending on `geosphere`. *Why:* it is the hot inner loop of
    estimation/covariate building; a tiny compiled kernel keeps deps light and is the
    package's substantive compiled contribution. A pure-R reference is kept for tests.
    (Source scored with `geosphere::distHaversine`, default radius 6378.137 km; our mean
    radius shifts reported distances ~0.11%, but as a constant scale on the objective it
    leaves the estimated parameters unchanged. Pass `radius=6378.137` to match exactly.)

11. **Scope: AMEN + Nigeria covariate plumbing.** The downstream AMEN network model and
    the bespoke missing-distance regression imputation are **not re-implemented** in the
    core package. → The package outputs PAL coordinates and dyadic PAL **distances**
    (linear and `log(d+0.01)`); a vignette shows feeding them to a simple dyadic model on
    simulated data. *Why:* keep the package focused, dependency-light, and broadly useful;
    re-bundling AMEN (a separate package, non-CRAN git ref) is out of scope for JOSS.

12. **Example data: real ACLED + a simulator.** → Ship the real ACLED Nigeria data used in
    the paper (`nigeria_acled`, from the public Dataverse replication, doi:10.7910/DVN/NLWWPE)
    as the featured example, plus `simulate_conflict_events()` + `nigeria_sim` for
    deterministic, license-free tests/examples. *Why:* the real data is already public in the
    replication archive and makes the paper concrete; the simulator keeps tests
    self-contained and seeded. (ACLED's terms of use still apply to the bundled data; this
    reverses the earlier simulated-only decision at the author's request.)

13. **Rubin pooling completeness.** → Provide the source's normal-based
    `T=W+(1+1/m)B` plus optional Barnard–Rubin df. *Why:* exact replication by default,
    with a more complete option for general use; cross-checked against `mice::pool` in a
    (suggested-pkg-guarded) test.

14. **Alter-event multiplicity.** The Dataverse code builds the alter history with a
    per-partner loop, so an event between two of the focal actor's alters is appended
    **twice** (once per partner). → The package counts each qualifying event **once** (a
    set filter over events involving any alter). *Why:* an event shouldn't be weighted up
    merely because two of your partners both attended it; the single-count is the cleaner,
    intended behavior. Numerically moot in the application (`pi≈0`); a deliberate
    divergence, now described accurately in `ALGORITHM.md`.

15. **Dyadic-prediction NA handling.** Source computes the event prediction as
    `rowMeans(cbind(PAL_A, PAL_B))` with `na.rm=FALSE`, so an event is **dropped** when
    *either* actor lacks history. → The package uses `na.rm=TRUE`, falling back to the one
    available PAL (so it predicts one-sided-history events the source drops). *Why:* a
    one-sided projection is still informative and increases coverage. Documented so the
    scored-event set / reported distances are understood to differ from the source; a
    `dyad_na="drop"` switch for byte-exact replication is a one-line add if wanted.
