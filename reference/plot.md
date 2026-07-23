# Plot metrics calculation outputs

S3 methods for plotting metrics calculated from the output sof either
[`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md),
[`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
[`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
or
[`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md).

## Usage

``` r
# S3 method for class 'incidence'
plot(x, xlab = "Time", ...)

# S3 method for class 'growth_rate'
plot(x, xlab = "Time", ...)

# S3 method for class 'Rt'
plot(x, xlab = "Time", ...)

# S3 method for class 'proportion'
plot(x, xlab = "Time", ...)
```

## Arguments

- x:

  Metrics calculation output of class `EpiStrainDynamics.metric`

- xlab:

  Time label for x axis, defaults to "Time"

- ...:

  Additional arguments passed to plot

## Value

ggplot2 plot output

## Examples

``` r
if (FALSE) { # interactive()
mod <- construct_model(
  method = random_walk(),
  pathogen_structure = single(
    case_timeseries = sarscov2$cases,
    time = sarscov2$date
  )
)

fit <- fit_model(mod)
gr <- growth_rate(mod)
plot(gr)
}
```
