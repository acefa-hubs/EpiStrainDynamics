# Generic Method for plotting metrics calculation outputs

S3 generic for plotting

## Usage

``` r
plot(df, xlab = "Time", ...)

# S3 method for class 'incidence'
plot(df, xlab = "Time", ...)

# S3 method for class 'growth_rate'
plot(df, xlab = "Time", ...)

# S3 method for class 'Rt'
plot(df, xlab = "Time", ...)

# S3 method for class 'proportion'
plot(df, xlab = "Time", ...)
```

## Arguments

- df:

  Metrics calculation output of class `EpiStrainDynamics.metric` from
  either
  [`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md),
  [`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
  [`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
  or
  [`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md).

- xlab:

  Time label for x axis, defaults to "Time"

- ...:

  Additional arguments passed to plot

## Value

ggplot2 plot output

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- construct_model(
    method = random_walk(),
    pathogen_structure = single(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date))

  fit <- fit_model(mod)
  gr <- growth_rate(mod)
  plot(gr)
} # }
```
