# Generic Method for Growth Rate Analysis

Computes the epidemiological growth rate, defined as the instantaneous
rate of change in log-incidence over time. Mathematically, it
represents: \$\$r_t = \log(I_t) - \log(I\_{t-1}) =
\log\left(\frac{I_t}{I\_{t-1}}\right)\$\$ where \\I_t\\ is incidence at
time t.

## Usage

``` r
growth_rate(fitted_model, ...)

# S3 method for class 'ps'
growth_rate(fitted_model, ...)

# S3 method for class 'rw'
growth_rate(fitted_model, ...)

# S3 method for class 'ps_single'
growth_rate(fitted_model, ...)

# S3 method for class 'rw_single'
growth_rate(fitted_model, ...)
```

## Arguments

- fitted_model:

  Fitted model object with class `EpiStrainDynamics.fit`

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

This metric quantifies the **proportional change** in disease incidence
from one time period to the next on a logarithmic scale, where:

- Positive values (\\r_t \> 0\\) indicate exponential growth

- Negative values (\\r_t \< 0\\) indicate exponential decline

- Values near zero (\\r_t \approx 0\\) indicate stable incidence

- The magnitude indicates the rate of exponential change

For example:

- A growth rate of 0.1 means incidence increased by approximately 10.5%
  (\\e^{0.1} - 1\\)

- A growth rate of -0.05 means incidence decreased by approximately 4.9%

- A growth rate of 0 means no change in incidence

This metric function can be run directly on the fitted model output.

## See also

Other metrics:
[`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
[`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md),
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
  gr <- growth_rate(fit)
} # }
```
