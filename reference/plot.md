# Generic Method for plotting metrics calculation outputs

S3 generic for plotting

## Usage

``` r
plot(df, ...)

# S3 method for class 'incidence'
plot(df)

# S3 method for class 'growth_rate'
plot(df)

# S3 method for class 'Rt'
plot(df)

# S3 method for class 'proportion'
plot(df)
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
