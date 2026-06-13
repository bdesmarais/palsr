pkgname <- "palsr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('palsr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bootstrap_pals")
### * bootstrap_pals

flush(stderr()); flush(stdout())

### Name: bootstrap_pals
### Title: Nonparametric bootstrap for PALS estimates and projections
### Aliases: bootstrap_pals

### ** Examples

ev <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
bt <- bootstrap_pals(ev, R = 10, model = "one", seed = 1)
summary(bt)




cleanEx()
nameEx("estimate_pals")
### * estimate_pals

flush(stderr()); flush(stdout())

### Name: estimate_pals
### Title: Estimate PALS parameters
### Aliases: estimate_pals

### ** Examples

ev  <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
fit <- estimate_pals(ev, model = "one")
fit
coef(fit)




cleanEx()
nameEx("haversine")
### * haversine

flush(stderr()); flush(stdout())

### Name: haversine
### Title: Great-circle (Haversine) distance
### Aliases: haversine

### ** Examples

haversine(0, 0, 0, 1)          # ~111 km per degree of latitude
haversine(7.4, 9.1, 8.5, 12.0) # Abuja-ish to Kano-ish




cleanEx()
nameEx("nigeria_sim")
### * nigeria_sim

flush(stderr()); flush(stdout())

### Name: nigeria_sim
### Title: Simulated subnational conflict events (Nigeria-like)
### Aliases: nigeria_sim
### Keywords: datasets

### ** Examples

data(nigeria_sim)
nigeria_sim
fit <- estimate_pals(nigeria_sim, model = "one")
coef(fit)



cleanEx()
nameEx("pal_distance")
### * pal_distance

flush(stderr()); flush(stdout())

### Name: pal_distance
### Title: Dyadic distance between Projected Actor Locations
### Aliases: pal_distance

### ** Examples

ev  <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
fit <- estimate_pals(ev, model = "one")
dy  <- data.frame(actor1 = "G01", actor2 = "G02",
                  time = as.Date("2012-12-01"))
pal_distance(ev, dy, fit)




cleanEx()
nameEx("pal_events")
### * pal_events

flush(stderr()); flush(stdout())

### Name: pal_events
### Title: Construct a validated dyadic-event table
### Aliases: pal_events

### ** Examples

df <- data.frame(
  from = c("A", "A", "B"),
  to   = c("B", "C", "C"),
  when = as.Date(c("2001-01-01", "2001-06-01", "2002-01-01")),
  x    = c(7.1, 8.0, 7.5),
  y    = c(9.0, 9.4, 10.1)
)
ev <- pal_events(df, actor1 = "from", actor2 = "to",
                 time = "when", lon = "x", lat = "y")
ev




cleanEx()
nameEx("pals_params")
### * pals_params

flush(stderr()); flush(stdout())

### Name: pals_params
### Title: Create a PALS parameter set
### Aliases: pals_params

### ** Examples

p <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
p
pals_params(alpha = 0.9, model = "one")




cleanEx()
nameEx("pool_rubin")
### * pool_rubin

flush(stderr()); flush(stdout())

### Name: pool_rubin
### Title: Pool estimates across imputations with Rubin's Rules
### Aliases: pool_rubin

### ** Examples

# Five imputations of a coefficient and its variance.
q <- c(1.10, 0.95, 1.20, 1.05, 0.98)
u <- c(0.04, 0.05, 0.045, 0.038, 0.052)
pool_rubin(q, u)
pool_rubin(q, u, df = TRUE, dfcom = 100)




cleanEx()
nameEx("predict.pals_fit")
### * predict.pals_fit

flush(stderr()); flush(stdout())

### Name: predict.pals_fit
### Title: Project locations from a fitted PALS model
### Aliases: predict.pals_fit

### ** Examples

ev  <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
fit <- estimate_pals(ev, model = "one")
predict(fit, predict_time = as.Date("2013-12-01"), type = "pal")[1:5, ]




cleanEx()
nameEx("predict_event_locations")
### * predict_event_locations

flush(stderr()); flush(stdout())

### Name: predict_event_locations
### Title: Predict dyadic event locations
### Aliases: predict_event_locations

### ** Examples

ev  <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
fit <- estimate_pals(ev, model = "one")
tg  <- ev[ev$time > as.Date("2012-01-01"), ]
head(predict_event_locations(ev, tg, fit))




cleanEx()
nameEx("project_pal")
### * project_pal

flush(stderr()); flush(stdout())

### Name: project_pal
### Title: Project the location of a single actor
### Aliases: project_pal

### ** Examples

ev <- simulate_conflict_events(n_actors = 8, n_events = 200, seed = 1)
p  <- pals_params(alpha = 0.9, model = "one")
project_pal(ev, actor = "G01", predict_time = as.Date("2010-12-01"), params = p)




cleanEx()
nameEx("project_pals")
### * project_pals

flush(stderr()); flush(stdout())

### Name: project_pals
### Title: Project locations for multiple actors
### Aliases: project_pals

### ** Examples

ev <- simulate_conflict_events(n_actors = 10, n_events = 300, seed = 1)
p  <- pals_params(alpha = 0.9, beta = 0.2, gamma = -10, eta = -10)
pal <- project_pals(ev, predict_time = as.Date("2010-12-01"), params = p)
head(pal)




cleanEx()
nameEx("simulate_conflict_events")
### * simulate_conflict_events

flush(stderr()); flush(stdout())

### Name: simulate_conflict_events
### Title: Simulate dyadic conflict events between moving actors
### Aliases: simulate_conflict_events

### ** Examples

ev <- simulate_conflict_events(n_actors = 12, n_events = 400, seed = 42)
ev
summary(ev)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
