# Print method for constructed EpiStrainDynamics models

Prints a concise summary instead of dumping the full data and standata
list to the console.

## Usage

``` r
# S3 method for class 'EpiStrainDynamics.model'
print(x, ...)
```

## Arguments

- x:

  an `EpiStrainDynamics.model` object, as returned by
  [`construct_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/construct_model.md)

- ...:

  further arguments passed to or from other methods (unused)

## Value

`x`, invisibly

## Examples

``` r
mod <- construct_model(
  method = random_walk(),
  pathogen_structure = single(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date"
  )
)
print(mod)
#> <EpiStrainDynamics model>
#> Method:             Random walk
#> Pathogen structure: single
#> Pathogens:          cases
#> Observations:       830
#> Day-of-week effect: FALSE
#> 
#> Data (head):
#> # A tsibble: 6 x 2 [1D]
#>   date       case_timeseries
#>   <date>               <dbl>
#> 1 2020-09-23            7096
#> 2 2020-09-24            7682
#> 3 2020-09-25            7318
#> 4 2020-09-26            6810
#> 5 2020-09-27            7023
#> 6 2020-09-28            9983
```
