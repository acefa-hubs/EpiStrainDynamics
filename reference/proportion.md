# Generic Method for proportion Analysis

Computes epidemiological proportion, defined as the relative fraction of
cases attributable to specific pathogen(s) or strain(s) at a given time
point, calculated as the ratio of incidence from selected pathogen(s) to
a reference group: \$\$P_t =
\frac{I\_{\text{numerator}}(t)}{I\_{\text{denominator}}(t)}\$\$

## Usage

``` r
proportion(
  fitted_model,
  numerator_combination = NULL,
  denominator_combination = NULL,
  ...
)

# S3 method for class 'ps'
proportion(
  fitted_model,
  numerator_combination = NULL,
  denominator_combination = NULL,
  ...
)

# S3 method for class 'rw'
proportion(
  fitted_model,
  numerator_combination = NULL,
  denominator_combination = NULL,
  ...
)
```

## Arguments

- fitted_model:

  Fitted model object with class `EpiStrainDynamics.fit` with `multiple`
  or `subtyped` pathogen structure.

- numerator_combination:

  Named pathogens or subtypes to be included in proportion numerator.

- denominator_combination:

  Named pathogens or subtypes to be included in proportion denominator,
  or 'all'.

- ...:

  Additional arguments passed to metrics calculation

## Value

named list of class `EpiStrainDynamics.metric` containing a dataframe of
the calculated metric outcome (`$measure`), the fit object (`$fit`), and
the constructed model object (`$constructed_model`). The `measure` data
frame contains the median of the epidemiological quantity (`y`), the 50%
credible interval of the quantity (`lb_50` & `ub_50`), the 95% credible
interval (`lb_95` & `ub_95`), the proportion greater than a defined
threshold value (`prop`), the pathogen name (`pathogen`), and the time
label (`time`).

## Details

Where incidences are derived from the exponential of log-incidence
estimates: \$\$P_t = \frac{\sum\_{i \in \text{numerator}}
\exp(\log\text{-incidence}\_{i,t})}{\sum\_{j \in \text{denominator}}
\exp(\log\text{-incidence}\_{j,t})}\$\$

This metric quantifies the **relative contribution** of specific
pathogen(s) or strain(s) to the total disease burden. Key
characteristics:

- Values between 0 and 1 (\\0 \leq P_t \leq 1\\) when denominator
  includes numerator components

- Values can exceed 1 (\\P_t \> 1\\) when denominator excludes numerator
  components

- Represents the fractional share of cases at each time point

- Time-varying to capture changing pathogen/strain dynamics

**Flexible combinations**:

- Individual proportions: Each pathogen relative to all pathogens
  (default)

- Custom numerator: Specific pathogen(s) of interest (e.g., variant of
  concern)

- Custom denominator: Either 'all' pathogens or a specified subset

- Subset comparisons: Compare specific groups (e.g., Alpha vs. Delta +
  Omicron)

This metric function can be run directly on the fitted model output.

## See also

Other metrics:
[`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md),
[`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md),
[`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  mod <- construct_model(
    method = p_spline(),
    pathogen_structure = multiple(
      case_timeseries = sarscov2$cases,
      time = sarscov2$date,
      component_pathogen_timeseries = list(
        alpha = sarscov2$alpha,
        delta = sarscov2$delta,
        omicron = sarscov2$omicron,
        other = sarscov2$other))
  )

  fit <- fit_model(mod)
  prop <- proportion(fit)

  # or a unique combination, compared to all pathogens
  prop2 <- proportion(fit,
    numerator_combination = c('alpha', 'delta', 'omicron'),
    denominator_combination = 'all'
  )

  # or a user-specified combination in both numerator and denominator
  prop3 <- proportion(fit,
    numerator_combination = 'alpha',
    denominator_combination = c('alpha', 'delta', 'omicron')
  )
} # }
```
