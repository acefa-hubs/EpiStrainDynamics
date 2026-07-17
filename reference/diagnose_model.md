# Diagnose model convergence and fit

Checks MCMC convergence diagnostics —
[R-hat](https://mc-stan.org/docs/reference-manual/analysis.html) and
[effective sample
size](https://mc-stan.org/docs/reference-manual/analysis.html) — against
user-specified thresholds, and reports any parameters that fail either
check.

## Usage

``` r
diagnose_model(fitted_model, rhat_threshold = 1.1, eff_sample_threshold = 100)
```

## Arguments

- fitted_model:

  A fitted model object of class `EpiStrainDynamics.fit`, as returned by
  [`fit_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/fit_model.md).

- rhat_threshold:

  R-hat threshold for convergence (default 1.1). Values above this
  threshold indicate chains have not mixed well.

- eff_sample_threshold:

  Effective sample size threshold (default 100). Values below this
  threshold indicate the posterior samples for a parameter are too
  autocorrelated to reliably estimate its distribution.

## Value

An object of class `list`, returned invisibly, containing:

- convergence:

  Logical; `TRUE` if no parameter fails either threshold

- rhat_issues:

  Character vector of parameter names with R-hat above `rhat_threshold`

- eff_sample_issues:

  Character vector of parameter names with effective sample size below
  `eff_sample_threshold`

- max_rhat:

  The largest R-hat value across all parameters

- min_neff:

  The smallest effective sample size across all parameters

- summary:

  The full `rstan::summary()` table the above are derived from

## Examples

``` r
if (FALSE) { # interactive()
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date))
  fit <- fit_model(mod)
  diagnose_model(fit)
}
```
