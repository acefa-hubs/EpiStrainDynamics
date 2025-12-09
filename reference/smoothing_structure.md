# Create Smoothing Structure Specification with Priors

This function creates a standardized smoothing structure object that
specifies both the smoothing structure and associated priors for
EpiStrainDynamics models.

## Usage

``` r
smoothing_structure(smoothing_type = "shared", tau_mean = NULL, tau_sd = NULL)
```

## Arguments

- smoothing_type:

  Character string specifying the smoothing type:

  - `"shared"`: All pathogens have the same smoothing parameter
    (equivalent to `tau[1]`). By default a model with a single pathogen
    will have `shared` smoothing type.

  - `"independent"`: Independent smoothing per pathogen (equivalent to
    `tau[number of pathogens]`)

  - `"correlated"`: Correlated smoothing type (equivalent to
    `Sigma[number of pathogens, number of pathogens]`)

- tau_mean:

  Optional numeric vector specifying the prior mean(s) for tau
  parameter. Can be provided for `shared` (single value) and
  `independent` smoothing types (can provide a single value which will
  be repeated for each pathogen or can provide a unique prior for each
  pathogen). Prior for tau for `correlated` smoothing type is not
  currently supported.

- tau_sd:

  Numeric vector specifying the prior standard deviation(s) for tau
  parameter. Can be provided for `shared` (single value) and
  `independent` smoothing types (can provide a single value which will
  be repeated for each pathogen or can provide a unique prior for each
  pathogen). Prior for tau for `correlated` smoothing type is not
  currently supported.

## Value

An object of class `EpiStrainDynamics.smoothing` containing:

- smoothing_type:

  The specified smoothing structure type

- tau_priors:

  Prior specifications for tau

## Examples

``` r
# Shared smoothing with scalar priors
shared_smooth <- smoothing_structure("shared", tau_mean = 0, tau_sd = 1)

# Independent smoothing with vector priors
indep_smooth <- smoothing_structure("independent",
                                    tau_mean = c(0, 0, 0),
                                    tau_sd = c(1, 1, 1))
```
