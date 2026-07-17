# Create Dispersion Structure Specification

This function creates a dispersion structure object that specifies
priors for the overdispersion parameter of the negative binomial
likelihood for the case timeseries. `phi` controls how much the observed
counts vary around the expected trend beyond Poisson noise: smaller
values allow more overdispersion (noisier counts relative to the mean),
larger values approach Poisson-like variance.

## Usage

``` r
dispersion_structure(phi_mean = NULL, phi_sd = NULL)
```

## Arguments

- phi_mean:

  Numeric scalar specifying the prior mean for the negative binomial
  dispersion parameter. Must be positive. Optional - if not provided, no
  priors will be set.

- phi_sd:

  Numeric scalar specifying the prior standard deviation for the
  dispersion parameter. Must be positive. Optional - if not provided, no
  priors will be set.

## Value

An object of class `EpiStrainDynamics.dispersion` containing:

- mean:

  The prior mean for phi

- sd:

  The prior standard deviation for phi

- priors_provided:

  Integer flag passed to the Stan model: `1` if no priors were supplied
  (Stan's built-in default prior is used), `2` if priors were supplied
  (the `phi_mean`/`phi_sd` values are used as the prior)

## Examples

``` r
# Create dispersion structure
disp_struct <- dispersion_structure(phi_mean = 2.0, phi_sd = 0.5)
```
