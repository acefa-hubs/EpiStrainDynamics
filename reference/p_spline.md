# Specify p_spline method

Penalised splines is one of two optional Bayesian smoothing prior
methods that can be selected and used in the model definition with
`EpiStrainDynamics`. The main benefit of selecting the penalised spline
method over random walks is that it can capture dynamical effects (fine
enough temporal resolution) while not being too computationally
expensive.

## Usage

``` r
p_spline(spline_degree = 3, days_per_knot = 3)
```

## Arguments

- spline_degree:

  polynomial degree of the individual spline segments used to construct
  the overall curve (must be a positive whole number)

- days_per_knot:

  number of days for each knot (must be a positive whole number)

## Value

list with method and model parameters of class
`EpiStrainDynamics.method`

## See also

Other method:
[`random_walk()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/random_walk.md)

## Examples

``` r
# Valid usage
p_spline(spline_degree = 2L, days_per_knot = 5L)
#> $method
#> [1] "p-spline"
#> 
#> $model_params
#> $model_params$spline_degree
#> [1] 2
#> 
#> $model_params$days_per_knot
#> [1] 5
#> 
#> 
#> attr(,"class")
#> [1] "EpiStrainDynamics.method"
p_spline(spline_degree = 3, days_per_knot = 7)
#> $method
#> [1] "p-spline"
#> 
#> $model_params
#> $model_params$spline_degree
#> [1] 3
#> 
#> $model_params$days_per_knot
#> [1] 7
#> 
#> 
#> attr(,"class")
#> [1] "EpiStrainDynamics.method"

# These will produce validation errors (as intended):
# \donttest{
# Non-positive values
try(p_spline(spline_degree = 0, days_per_knot = 5))
#> Error in validate_positive_whole_number(spline_degree, "spline_degree") : 
#>   Argument spline_degree must be a positive number
try(p_spline(spline_degree = 3, days_per_knot = -1))
#> Error in validate_positive_whole_number(days_per_knot, "days_per_knot") : 
#>   Argument days_per_knot must be a positive number

# Non-whole numbers
try(p_spline(spline_degree = 2.5, days_per_knot = 5))
#> Error in validate_positive_whole_number(spline_degree, "spline_degree") : 
#>   Argument spline_degree must be a whole number
try(p_spline(spline_degree = 3, days_per_knot = 4.2))
#> Error in validate_positive_whole_number(days_per_knot, "days_per_knot") : 
#>   Argument days_per_knot must be a whole number

# Non-numeric values
try(p_spline(spline_degree = "invalid", days_per_knot = 5))
#> Error in validate_positive_whole_number(spline_degree, "spline_degree") : 
#>   Argument spline_degree must be numeric
# }
```
