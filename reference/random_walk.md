# Specify random walk method

Random walk is one of two optional Bayesian smoothing prior methods that
can be selected and used in the model definition with
`EpiStrainDynamics`. The random walk is a stochastic process that
describes a path of a series of random steps on a mathematical space,
and each next step's direction only depends on the current position, not
the previous path. No additional arguments must be supplied to define
the random walk method.

## Usage

``` r
random_walk()
```

## Value

list with method identified as random walk of class
`EpiStrainDynamics.method`

## See also

Other method:
[`p_spline()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/p_spline.md)

## Examples

``` r
random_walk()
#> $method
#> [1] "random-walk"
#> 
#> attr(,"class")
#> [1] "EpiStrainDynamics.method"
```
