# Create Dispersion Structure Specification

This function creates a dispersion structure object that specifies
priors for the overdispersion parameter of the negative binomial
likelihood for the case timeseries.

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

An object of class `EpiStrainDynamics.dispersion` containing prior
specifications for phi parameter

## Examples

``` r
# Create dispersion structure
disp_struct <- dispersion_structure(phi_mean = 2.0, phi_sd = 0.5)
```
