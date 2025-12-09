# Construct model

Construct model

## Usage

``` r
construct_model(
  method,
  pathogen_structure,
  smoothing_params = smoothing_structure(),
  dispersion_params = dispersion_structure(),
  pathogen_noise = FALSE,
  dow_effect = FALSE
)
```

## Arguments

- method:

  either
  [`random_walk()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/random_walk.md)
  or
  [`p_spline()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/p_spline.md)

- pathogen_structure:

  either
  [`single()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/single.md),
  [`multiple()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/multiple.md),
  or
  [`subtyped()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/subtyped.md)

- smoothing_params:

  argument is optional and defines the structure of the smoothing terms
  including optionally setting the smoothing prior tau. Created with
  [`smoothing_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/smoothing_structure.md).
  NULL option defaults to 'shared' smoothing structure and default
  priors.

- dispersion_params:

  argument is optional and defines priors for the overdispersion
  parameter of the negative binomial likelihood for the case timeseries.
  Created using
  [`dispersion_structure()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/dispersion_structure.md).
  NULL option uses default priors for phi.

- pathogen_noise:

  logical whether individual pathogen counts have additional
  gamma-distributed noise. Default is FALSE. Models with `single`
  pathogen structure will be set to FALSE.

- dow_effect:

  logical whether to incorporate a day of week model.

## Value

a list containing the data, the model parameters, and pathogen names of
class `EpiStrainDynamics.model`

## Examples

``` r
mod <- construct_model(

  method = p_spline(),

  pathogen_structure = multiple(
    data = sarscov2,
    case_timeseries = 'cases',
    time = 'date',
    component_pathogen_timeseries = c('alpha', 'delta', 'omicron', 'other')
    ),

   smoothing_params = smoothing_structure(
      'independent', tau_mean = c(0, 0.1, 0.3, 0), tau_sd = rep(1, times = 4)),
   dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
   pathogen_noise = FALSE,
   dow_effect = TRUE
)
#> Registered S3 method overwritten by 'tsibble':
#>   method               from 
#>   as_tibble.grouped_df dplyr

```
