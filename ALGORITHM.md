# PALS Algorithm Specification (implementation contract)

This is the precise specification the R package implements, distilled from
`docs-research/01-pals-method.md` (which was reconstructed from the published article
**and** the authors' Dataverse replication code). Where the released code and the
package diverge, the divergence is intentional and documented in `DECISIONS.md`.

Reference: Kim, Liu & Desmarais (2023) "Spatial modeling of dyadic geopolitical
interactions between moving actors," *PSRM* 11(3):633–644, doi:10.1017/psrm.2022.6.

---

## 1. Inputs

A set of **dyadic events**. Each event `e` records:

- `actor1`, `actor2` — the two actors involved (character/factor ids).
- `time` — a `Date` (or numeric day index).
- `lon`, `lat` — event location in decimal degrees.

Derived quantities (computed internally, never asked of the user):

- **Focal history** of actor `i` at time `t`: all events with `time < t` in which
  `i == actor1` or `i == actor2`. Call this `E_i(t)`.
- **Alters** of `i` (through time `t`): the set of distinct actors that share at least
  one event with `i` (i.e. the other actor in any of `i`'s events with `time < t`).
- **Alter history** `E_k(t)`: the **distinct** events with `time < t` that involve *any*
  alter of `i` — each such event counted **once**, however many of `i`'s alters it
  involves. (The Dataverse code loops per partner and so counts an event between two
  alters twice; the package counts it once by design — see `DECISIONS.md` #14.)

`a(e) = (t - time(e))` measured in **days** (the smoothing timescale). The prediction
time `t` defaults, in the yearly workflow, to **December 1 of the target year**, but is
a free argument in the package (`predict_time`).

---

## 2. Projection for one actor at one time

Given parameters `(alpha, beta, gamma, eta)` with `alpha, beta > 0`:

**Focal weights** over `e in E_i(t)`:
```
w_i(e) = 1 / ( a(e)^alpha + 0.01 )
W_i(e) = w_i(e) / sum_r w_i(r)          # normalized (sums to 1)
```

**Alter weights** over `e in E_k(t)`:
```
w_k(e) = 1 / ( a(e)^beta + 0.01 )
W_k(e) = w_k(e) / sum_r w_k(r)          # normalized — DEFAULT (see note)
```

**Relevance / event-count ratio** (`n_i = |E_i(t)|`, `n_k = |E_k(t)|`):
```
v = ( n_i / n_k ) ^ ( 1 / sqrt(n_k) )
```

**Mixing weight** (logistic; `pi` = weight on alters):
```
pi = plogis(gamma + eta * v) = 1 / (1 + exp(-(gamma + eta*v)))
```

**Projected location** (applied independently to lon and lat):
```
g_i(t) = (1 - pi) * sum_e W_i(e) g(e)  +  pi * sum_e W_k(e) g(e)
```

### Fixed numerical constants
- `EPS_AGE = 0.01` — offset inside every age weight so age-0 events give a finite
  weight (`1/0.01 = 100`) instead of `Inf`. **Do not remove.**

### Model variants
- **Four-parameter ("full")**: estimate `(alpha, beta, gamma, eta)`; blends focal+alter.
- **One-parameter ("focal" / alpha-only)**: force `pi = 0`; only `alpha` is active and
  estimated. Equivalent to `gamma -> -Inf`. The projection is a pure recency-weighted
  average of the focal actor's own past event locations.

### Edge cases (must be guarded)
- **No focal history** (`n_i == 0`): PAL is `NA` (NOT (0,0)). Carry a `has_history`
  flag; never coerce missing to the origin.
- **No alters** (`n_k == 0`): `v`/`pi` undefined → fall back to `pi = 0` (focal-only).
- **Duplicate events** (from bootstrap): treated as ordinary repeated events; no dedup.
- **Shared dyadic events**: an event involving `i` also appears in `E_i`; events
  involving `i`'s alters (possibly including events `i` also attended) appear in `E_k`.
  No de-duplication *across* the two sets. *Within* `E_k`, an event involving two alters
  is counted once (the source counts it per-alter; deliberate change — `DECISIONS.md` #14).

### Alter-weight normalization (the one substantive correction)
The released code line `weights_k <- 1/(sum(weights_k))` overwrites the alter weight
vector with a **scalar**, so the code computes an un-normalized *sum* of alter
coordinates, not a weighted average. The package implements the **paper's intended
normalized** `W_k(e)` as the default, and exposes `alter_weight = "legacy"` to
reproduce the Dataverse behavior exactly. Because estimated `pi` is ~0 in the
application, the two agree numerically there. See `DECISIONS.md` #1.

---

## 3. Dyadic event-location prediction

For an event between actors A and B at time `t`, the predicted event location is the
**arithmetic mean** of the two actors' PALs:
```
pred_lon = mean(PAL_lon(A,t), PAL_lon(B,t))   # na.rm = TRUE: fall back to the present one
pred_lat = mean(PAL_lat(A,t), PAL_lat(B,t))
```
If both PALs are `NA`, the event has no prediction and is dropped from the objective.

---

## 4. Parameter estimation

**Objective**: minimize the **mean Haversine (great-circle) distance in kilometres**
between predicted and observed event locations, over a set of fit events, each predicted
from history strictly preceding its own time ("marching forward"):
```
J(theta) = mean_e  haversine_km( pred(e; theta), observed(e) )     # na.rm = TRUE
```
- Optimize over the unconstrained vector `theta`:
  - four-param: `(gamma, eta, log_alpha, log_beta)`, with `alpha = exp(log_alpha)`,
    `beta = exp(log_beta)` (guarantees positivity, no box constraints needed).
  - one-param: `log_alpha` only.
- Optimizer: base R `stats::optim`. Default method `"Nelder-Mead"` for the 4-param
  model and `"Brent"` (1-D) for the 1-param model; allow `method=` override (the paper
  used L-BFGS-B / "hill-climbing"). Multiple random restarts optional.
- Earth radius for Haversine: `R = 6371.0088` km (mean radius).

Return a fitted object storing: estimated params (natural scale), model type,
objective value, optimizer diagnostics, the events used, and `predict_time`.

---

## 5. Bootstrap + multiple imputation (Rubin's Rules)

**Bootstrap** (nonparametric, over events): resample all events with replacement,
re-estimate parameters, re-project PALs. Repeat `R` times (paper: `R = 10`, fixed
seeds). Returns the distribution of parameter estimates and of projected locations —
the basis for confidence bands and for treating each PAL set as one imputation.

**Rubin pooling** across `m` imputations, for a scalar estimand with per-imputation
point estimates `Q_j` and variances `U_j`:
```
Qbar = mean_j Q_j                              # pooled estimate
Ubar = mean_j U_j                              # within-imputation variance
B    = var_j Q_j  (denominator m-1)            # between-imputation variance
T    = Ubar + (1 + 1/m) * B                    # total variance
SE   = sqrt(T)
```
Provide a `pool_rubin(estimates, variances)` helper returning `Qbar`, `Ubar`, `B`,
`T`, `SE`, and (optionally) the Barnard–Rubin degrees of freedom. The source code uses
the simple normal-based `T = W + (1+1/m)B`; the package reproduces that and additionally
offers the df correction.

---

## 6. What the package deliberately generalizes

- `predict_time` is a free argument (source code fixes Dec-1).
- Estimation works on any user-supplied set of fit events / time windows. A
  `cutoff = c("day","month","year")` argument sets the history boundary: `"day"` (default)
  = strict date; `"year"` reproduces the source's calendar-year convention (events before
  Jan 1 of the prediction year). Ages are always measured from `predict_time`.
- Alter-weight normalization corrected by default (`legacy` available).
- Objective is mean-km by default; `aggregate = "sum"` available.
- The downstream AMEN network model and the Nigeria-specific missing-distance
  regression imputation are **out of scope** for the core package (shown in a vignette
  using simulated data, not re-implemented). PAL distances are exposed so users can feed
  them to any dyadic model.
