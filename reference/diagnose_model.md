# Diagnose model convergence and fit

Diagnose model convergence and fit

## Usage

``` r
diagnose_model(fitted_model, rhat_threshold = 1.1, eff_sample_threshold = 100)
```

## Arguments

- fitted_model:

  fitted model object

- rhat_threshold:

  R-hat threshold for convergence (default 1.1)

- eff_sample_threshold:

  effective sample size threshold

## Value

list of diagnostic information
