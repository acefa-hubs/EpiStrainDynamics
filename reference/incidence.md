# Generic Method for Incidence Analysis

Computes epidemiological incidence, defined as the number of new cases
occurring at a specific time point, derived by exponentiating the
log-incidence estimates from the fitted model: \$\$I_t =
\exp(\log\text{-incidence}\_t)\$\$

## Usage

``` r
incidence(fitted_model, dow, ...)

# S3 method for class 'ps'
incidence(fitted_model, dow, ...)

# S3 method for class 'rw'
incidence(fitted_model, dow, ...)

# S3 method for class 'ps_single'
incidence(fitted_model, dow, ...)

# S3 method for class 'rw_single'
incidence(fitted_model, dow, ...)
```

## Arguments

- fitted_model:

  Fitted model object with class `EpiStrainDynamics.fit`

- dow:

  Logical whether or not to include day-of-week in incidence calc

- ...:

  Additional arguments passed to metrics calculation

## Value

named list of class `EpiStrainDynamics.metric` containing a dataframe of
the calculated metric outcome (`$measure`), the fit object (`$fit`), and
the constructed model object (`$constructed_model`). The `measure` data
frame contains the median of the epidemiological quantity (`y`), the 50%
credible interval of the quantity (`lb_50` & `ub_50`), the 95% credible
interval (`lb_95` & `ub_95`), the proportion greater than a defined
threshold value (`prop`), the pathogen name (`pathogen`), and the time
label (`time`).

## Details

This metric quantifies the **absolute number of new cases** at each time
point, where it is:

- Always positive (since it's an exponentiated value)

- Represents the expected case count at time t

- Can be adjusted for day-of-week effects when modeled

- Provides uncertainty quantification through posterior credible
  intervals

**Day-of-week adjustment**: When day-of-week effects are included in the
model, the incidence is further adjusted as: \$\$I_t^{adj} = I_t \times
\text{week\\effect} \times \text{dow\\simplex}\[\text{DOW}(t)\]\$\$

This accounts for systematic variations in case reporting or
transmission patterns across different days of the week (e.g., lower
weekend reporting, higher weekday transmission).

This metric function can be run directly on the fitted model output.

## See also

Other metrics:
[`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
[`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
[`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date))

  fit <- fit_model(mod)
  inc <- incidence(fit, dow = TRUE)
} # }
```
