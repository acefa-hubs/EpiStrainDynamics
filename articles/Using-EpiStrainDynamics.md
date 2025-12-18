# Using EpiStrainDynamics

EpiStrainDynamics extends an existing statistical modelling framework
capable of inferring the trends of up to two pathogens. The modeling
framework has been extended here to handle any number of pathogens; fit
to time series data of counts (eg, daily number of cases); incorporate
influenza testing data in which the subtype for influenza A samples may
be undetermined; account for day-of-the-week effects in daily data;
include options for fitting penalized splines or random walks; support
additional (optional) correlation structures in the parameters
describing the smoothness of the penalized splines (or random walks);
and account for additional (optional) sources of noise in the
observation process.

## Step 1: Construct model

These modelling specifications are specified using the `construct_model`
function. The correct stan model is then applied based on the
specifications provided.
[`construct_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/construct_model.md)
requires specification of the `method` and `pathogen_structure`, and
allows the user to additionally specify four arguments related to model
design: `smoothing_params`, `dispersion_params`, `pathogen_noise`, and
`dow_effect`. We will break these down one by one.

### Method

EpiStrainDynamics has pre-compiled stan models that fit either with
bayesian penalised splines or random walks. These are specified using
the `method` argument of
[`construct_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/construct_model.md)
as functions, either with
[`random_walk()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/random_walk.md)
or
[`p_spline()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/p_spline.md).
The penalised spline model has two further options to specify:
`spline_degree` is the polynomial degree of the individual spline
segments used to construct the overall curve (must be a positive whole
number) and `days_per_knot`, which is the number of days for each knot
(must also be a positive whole number).

So we may specify: \`\`\` method = random_walk(),

## OR

method = p_spline(spline_degree = 3, \# example value for spline degree
days_per_knot = 2) \# example value for days per knot

    ## Pathogen structure

    There are three main types of pathogen structure available to model: `single()`, `multiple()`, and `subtyped()`.
    These functions require the name of the dataset itself and the column names for different data elements.

    The `single()` pathogen structure is the simplest and models a single pathogen timeseries. The name of the dataframe is passed to argument `data`, the name of the column with total case data is passed to `case_timeseries`, and the name of the column of time data is passed to `time`. It can be specified as follows, illustrated using data provided with the package `sarscov2`:

pathogen_structure = single( data = sarscov2, \# dataframe
case_timeseries = ‘cases’, \# timeseries of case data time = ‘date’ \#
date or time variable )

    The `multiple()` pathogen structure allows modelling of different component pathogens. In addition to specifying `data`, `case_timeseries`, and `time`, these additional pathogens are specified as a vector of column names with the argument `component_pathogen_timeseries`. Example pathogen structure specification for multiple pathogens model using the `sarscov2` dataset:

pathogen_structure = multiple( data = sarscov2, \# dataframe
case_timeseries = ‘cases’, \# timeseries of case data time = ‘date’, \#
date or time variable labels

component_pathogen_timeseries = c( \# vector of column names of ‘alpha’,
‘delta’, ‘omicron’, ‘other’ \# component pathogens ) )

    The `subtyped()` pathogen structure enables additional complexity specifically for an influenza modelling scenario by allowing the user to incorporate testing data for influenza A subtypes. The unsubtyped column is specified with `influenzaA_unsubtyped_timeseries`, and the subtyped data are specified with vector of column names of provided to `influenzaA_subtyped_timeseries`. Additional pathogens are provided as a vector of column names to `other_pathogen_timeseries`. Example pathogen structure specification for subtyped model using the `influenza` dataset included in the package:

pathogen_structure = subtyped( data = influenza, \# dataframe
case_timeseries = ili, \# timeseries of case data time = week, \# date
or time variable labels

influenzaA_unsubtyped_timeseries = ‘inf_A’, \# unsubtyped influenzaA
influenzaA_subtyped_timeseries = c( \# subtyped influenzaA ‘inf_H3N2’,
‘inf_H1N1’ ), other_pathogen_timeseries = c( \# other pathogens ‘inf_B’,
‘other’ ) )

    ## Smoothing parameters

    The argument `smoothing_params` allows users to modify the correlation structures in the parameters describing smoothness and set related priors. These are specified with the function `smoothing_structure()` that requires the user to specify a `smoothing_type` that is either `shared` (all pathogens have the same smoothness), `independent` (each pathogen has completely independent smoothing structure), or `correlated` (smoothing structure is correlated among pathogens). For `shared` or `independent` smoothing types parameters for the mean and standard deviation of the prior on tau can also be specified, as below:

smoothing_params = smoothing_structure( smoothing_type = ‘independent’,
tau_mean = c(0, 0.1, 0.3, 0), tau_sd = rep(1, times = 4) )

    ## Dispersion parameters

    The argument `dispersion_params` allows users to set a prior for the overdispersion parameter of the negative binomial likelihood for the case timeseries. It is specified using `dispersion_structure()` as below:

dispersion_params = dispersion_structure( phi_mean = 0, phi_sd = 1 )

    ## Pathogen noise

    Whether to include noise between individual pathogens as well as the observation noise is specified as a logical (`TRUE` or `FALSE`) to the argument `pathogen_noise`.

    ## Day of week effect

    Day of week effect is specified as a logical (`TRUE` or `FALSE`) to the `dow_effect` argument. In plotting the day of week effect can be selectively removed.

    ## Worked example

    A full worked example using a subtyped structure:


    ``` r
    mod <- construct_model(

      method = p_spline(),

      pathogen_structure = subtyped(
        data = influenza,
        case_timeseries = 'ili',
        time = 'week',
        influenzaA_unsubtyped_timeseries = 'inf_A',
        influenzaA_subtyped_timeseries = c('inf_H3N2', 'inf_H1N1'),
        other_pathogen_timeseries = c('inf_B', 'other')
      ),

      smoothing_params = smoothing_structure(
        'independent', tau_mean = c(0, 0.1, 0.3, 0), tau_sd = rep(1, times = 4)),
      dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
      pathogen_noise = FALSE,
      dow_effect = TRUE
    )

### Output

The constructed model object is a named list containing input data,
accessed with `$data`, parameter values (such as smoothing structure,
observation noise, and penalised spline parameters, if appropriate),
accessed with `model_params`, names of provided pathogens, accessed with
`pathogen_names`, and record of whether day of week effect has been
selected, accessed with `dow_effect`.

## Step 2: Fit model

The model estimates the expected value of the time series (eg, a
smoothed trend in the daily number of cases accounting for noise) for
each individual pathogen. Model parameterisation decisions specified
when constructing the model in step 1 mean the correct stan model will
be applied at this stage by simply calling
[`fit_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/fit_model.md)
onto the constructed model object.

``` r
fit <- fit_model(
  mod,
  n_iter = 2000,
  n_warmup = 1000,
  verbose = FALSE
)
```

The output of this function is a list with the stan fit object, accessed
with `$fit`, and the elements of the constructed model object from the
previous step, accessed with `$constructed_model`.

### Step 3: Check fit

A summary of model convergence diagnostics can be viewed using:

``` r
diagnose_model(fit)
#> Model Convergence Diagnostics
#> =============================
#> Overall convergence: GOOD 
#> Maximum R-hat: 1.006 
#> Minimum n_eff: 1018
```

Beyond this summary, the stan fit object can be interrogated with any
package that works with stan outputs. For example, one could visualise
posteriors with the package `bayesplot`:

``` r
bayesplot::mcmc_areas(as.matrix(fit$fit), pars = 'tau[1]', prob = 0.8)
```

![plot of chunk
unnamed-chunk-5](articles/figures/Using-EpiStrainDynamics-unnamed-chunk-5-1.png)

plot of chunk unnamed-chunk-5

The `bayesplot` package can also be used to evaluate convergence:

``` r
bayesplot::mcmc_trace(rstan::extract(fit$fit, permuted = FALSE), pars = c('a[1,1]'))
```

![plot of chunk
unnamed-chunk-6](articles/figures/Using-EpiStrainDynamics-unnamed-chunk-6-1.png)

plot of chunk unnamed-chunk-6

Or else the application `shinystan` is a great tool to visualise the fit
and diagnose issues:

    library(shinystan)
    launch_shinystan(fit$fit)

Model not converging? You can find more information about divergent
transitions
[here](https://mc-stan.org/docs/2_19/reference-manual/divergent-transitions).
Some quick tips:

- Increase the `adapt_delta` parameter in fit_model(), which governs how
  big the steps are in the MCMC.
- Try reducing the complexity of the model. For example, by using
  `pathogen_noise = FALSE` or a shared smoothing structure.
- Increase the number of warmup samples with the n_warmup parameter.

## Step 3: Calculate and explore epidemiological quantities

The package provides helper functions to calculate a number of useful
epidemiological quantities. The output of these methods functions are a
named list containing a data frame of the outcome quantity (`$measure`),
the fit object (`$fit`), and the constructed model object
(`$constructed_model`). `measure` is a data frame containing the median
of the epidemiological quantity (`y`), the 50% credible interval of the
quantity (`lb_50` & `ub_50`), the 95% credible interval (`lb_95` &
`ub_95`), the proportion greater than a defined threshold value
(`prop`), the pathogen name (`pathogen`), and the time label (`time`).

Calculate epidemic growth rate with
[`growth_rate()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/growth_rate.md):

``` r
gr <- growth_rate(fit)
head(gr$measure)
#>   pathogen pathogen_idx          y         lb_50      ub_50       lb_95      ub_95    prop       time
#> 1    Total           NA 0.04083607 -0.0069918329 0.08851913 -0.09467077 0.18206126 0.71150 2012-01-09
#> 2    Total           NA 0.04080373  0.0006121294 0.07873632 -0.06878750 0.15687638 0.75425 2012-01-16
#> 3    Total           NA 0.03748100  0.0082013129 0.06783071 -0.04555991 0.13255478 0.80475 2012-01-23
#> 4    Total           NA 0.03330234  0.0090557712 0.05835547 -0.03434353 0.10992005 0.82600 2012-01-30
#> 5    Total           NA 0.03060279  0.0065369216 0.05356837 -0.03843168 0.10196366 0.80675 2012-02-06
#> 6    Total           NA 0.02843207  0.0066621189 0.04939690 -0.03564196 0.09146533 0.81300 2012-02-13
plot(gr)
```

![plot of chunk
growth_rate](articles/figures/Using-EpiStrainDynamics-growth_rate-1.png)

plot of chunk growth_rate

Calculate effective reproduction number over time with
[`Rt()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/Rt.md)
(requiring specification of a generation interval distribution):

``` r
rt <- Rt(fit, gi_dist = function(x) 4*x*exp(-2*x))
plot(rt)
```

![plot of chunk Rt](articles/figures/Using-EpiStrainDynamics-Rt-1.png)

plot of chunk Rt

Calculate incidence with or without a day of week effect with
[`incidence()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/incidence.md):

``` r
inc_dow <- incidence(fit, dow = TRUE)
plot(inc_dow)
```

![plot of chunk
incidence](articles/figures/Using-EpiStrainDynamics-incidence-1.png)

plot of chunk incidence

Calculate proportions of different combinations of cases attributable to
different pathogens/subtypes using
[`proportion()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/proportion.md).
By default, the function will return a dataframe with proportions of
each pathogen or subtype out of all pathogens/subtypes. Alternatively,
one can specify a selection of pathogens/subtypes by their names in the
named lists provided to
[`construct_model()`](https://acefa-hubs.github.io/EpiStrainDynamics/reference/construct_model.md):

    prop <- proportion(
      fit,
      numerator_combination = c('inf_H3N2', 'inf_H1N1'),
      denominator_combination = c('inf_H3N2', 'inf_H1N1', 'inf_B', 'other')
    )

``` r
prop <- proportion(fit)
plot(prop)
```

![plot of chunk
proportion](articles/figures/Using-EpiStrainDynamics-proportion-1.png)

plot of chunk proportion
