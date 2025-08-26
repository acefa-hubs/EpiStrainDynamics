# EpiStrainDynamics
[![codecov](https://codecov.io/gh/acefa-hubs/EpiStrainDynamics/graph/badge.svg)](https://app.codecov.io/gh/acefa-hubs/EpiStrainDynamics)
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

This repository contains code corresponding to the R package `EpiStrainDynamics`. 
For code corresponding to the analyses in the paper: _Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data_, see branch [`paper_analysis`](https://github.com/acefa-hubs/EpiStrainDynamics/tree/paper_analysis). 

Estimating the temporal trends in infectious disease activity is crucial for monitoring disease spread and the impact of interventions. 
Surveillance indicators routinely collected to monitor these trends are often a composite of multiple pathogens. 
For example, ‘influenza-like illness’ — routinely monitored as a proxy for influenza infections — is a symptom definition that could be caused by a wide range of pathogens, including multiple subtypes of influenza, SARS-CoV-2, and RSV. 
Inferred trends from such composite time series may not reflect the trends of any one of the component pathogens, each of which can exhibit distinct dynamics. 
Although many surveillance systems routinely test a subset of individuals contributing to a surveillance indicator — providing information on the relative contribution of the component pathogens — trends may be obscured by time-varying testing rates or substantial noise in the observation process. 

EpiStrainDynamics extends an existing statistical modelling framework capable of inferring the trends of up to two pathogens. 
The modeling framework has been extended here to handle any number of pathogens; fit to time series data of counts (eg, daily number of cases); incorporate influenza testing data in which the subtype for influenza A samples may be undetermined; account for day-of-the-week effects in daily data; include options for fitting penalized splines or random walks; support additional (optional) correlation structures in the parameters describing the smoothness of the penalized splines (or random walks); and account for additional (optional) sources of noise in the observation process. 

## Step 1: Construct model

Modelling specifications are provided using the `construct_model` function. The correct stan model is then applied based on the specifications provided. `construct_model()` takes three arguments: `method`, `pathogen_structure`, and `dow_effect`. 

### Method

EpiStrainDynamics has pre-compiled stan models that fit either with bayesian penalised splines or random walks. These are specified using the `method` argument of `construct_model()` as functions, either with `random_walk()` or `p_spline()`. The penalised spline model has two further options to specify: `spline_degree` is the polynomial degree of the individual spline segments used to construct the overall curve (must be a positive whole number) and `days_per_knot`, which is the number of days for each knot (must also be a positive whole number). 

So we may specify: 
```
mod <- construct_model(
  method = random_walk(), 
  ...
)
# OR

mod <- construct_model(
  method = p_spline(spline_degree = 3, days_per_knot = 2),  # example options
  ...
)
```

### Pathogen structure

There are three main types of pathogen structure available to model: `single()`, `multiple()`, and `subtyped()`. 

The `single()` pathogen structure is the simplest and models a single pathogen timeseries. A vector is the case timeseries is provided to the argument `case_timeseries`, a vector of date or time labels is provided to `time`, and optionally a pathogen name can be assigned with `pathogen_name`. It can be specified as follows, illustrated using data provided with the package `sarscov2`: 

```
single(
  case_timeseries = sarscov2$cases,           # timeseries of case data
  time = sarscov2$date,                       # date or time variable labels
  pathogen_name = 'SARS-COV-2'                # optional name of pathogen 
)
```

The `multiple()` and `subtyped()` pathogen structures both also have the `case_timeseries` and `time` arguments, as in the `single()` structure. 

In addition to these two arguments, instead of a single argument for `pathogen_name`, `multiple()` and `subtyped()` have one or more arguments that allow the user to define the names and data for the different components or subtypes to be modelled. 
For `multiple()`, these are specified as a named list with the argument `component_pathogen_timeseries`. The names in this list will be using in subsequent plotting. 
For `subtyped()`, which allows the user to incorporate testing data for influenza A subtypes, a vector of total unsubtyped influenza A cases is specified with `influenzaA_unsubtyped_timeseries`, a named list of the subtyped influenza A timeseries is provided to `influenzaA_subtyped_timeseries`, and further pathogens are provided in a named list to `other_pathogen_timeseries`. 

In addition to specifying the data and pathogen names, `multiple()` and `subtyped()` both allow the user to modify the correlation structures in the parameters describing the smoothness (with argument `smoothing_structure`) and account for additional sources of noise in the observation process (`observation_noise`). `smoothing_structure` is either `'shared'` (all pathogens have the same smoothness), `'independent'` (each pathogen has completely independent smoothing structure), or `'correlated'` (smoothing structure is correlated among pathogens). 
`observation_noise` is either `'observation_noise_only'` (only includes observation noise - the same between pathogens) or `'pathogen_specific_noise'` (includes noise in individual pathogens as well). 

Example pathogen structure specification for multiple pathogens model:
```
multiple(
   case_timeseries = sarscov2$cases,        # timeseries of case data
   time = sarscov2$date,                    # date or time variable labels
   
   component_pathogen_timeseries = list(    # named list of component pathogens
     alpha = sarscov2$alpha,
     delta = sarscov2$delta,
     omicron = sarscov2$omicron,
     other = sarscov2$other
   ),
   
   smoothing_structure = 'independent',             # correlation structures
   observation_noise = 'observation_noise_only'     # observation noise
 )
```

Example pathogen structure specification for subtyped model:
```
subtyped(
   case_timeseries = influenza$ili,         # timeseries of case data
   time = influenza$week,                   # date or time variable labels
   
   influenzaA_unsubtyped_timeseries = influenza$inf_A,  # unsubtyped influenzaA
   influenzaA_subtyped_timeseries = list(       # named list of subtyped infA
     H3N2 = influenza$inf_H3N2,
     H1N1 = influenza$inf_H1N1
   ),
   other_pathogen_timeseries = list(            # named list of other pathogens
     influenzaB = influenza$inf_B,
     other = influenza$num_spec - influenza$inf_all
   ),
   
   smoothing_structure = 'correlated',            # correlation structures
   observation_noise = 'pathogen_specific_noise'  # observation noise
 )
```

### Day of week effect

Day of week effect is specified as a logical (`TRUE` or `FALSE`) to the `dow_effect` argument. In plotting the day of week effect can be selectively removed. 

```
mod <- construct_model(
  method = random_walk(), 
  pathogen_structure = ...,
  dow_effect = TRUE
)
```

## Step 2: fit model

The model estimates the expected value of the time series (eg, a smoothed trend in the daily number of cases accounting for noise) for each individual pathogen. 
Model parameterisation decisions specified when constructing the model in step 1 mean the correct stan model will be applied at this stage by simply calling `fit_model()` onto the constructed model object. 

## Step 3: Calculate and visualise epidemiological quantities

Calculate epidemic growth rate with `growth_rate(fit)`, effective reproduction number over time with `Rt(fit, gi_dist = X)` (requiring specification of a generation interval distribution X), incidence with or without a day of week effect with `incidence(fit, dow_effect = TRUE)`, and proportions of different combinations of cases attributable to different pathogens/subtypes using `proportion()`, where the pathogens/subtypes are specified by their names in the named lists provided to `construct_model()`:
```
prop <- proportion(
  fit, 
  numerator_combination = c('H3N2', 'H1N1', 'influenzaB'),
  denominator_combination = c('H3N2', 'H1N1', 'influenzaB', 'other')
)
```

These quantities can be visualised using `plot()`. For a more detailed discussion, check out the [vignette](https://acefa-hubs.github.io/EpiStrainDynamics/articles/Using-EpiStrainDyamics.html).  


## Citation
Oliver Eales, Saras M Windecker, James M McCaw, Freya M Shearer, Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data, American Journal of Epidemiology, 2025;, kwaf119, https://doi.org/10.1093/aje/kwaf119

## Contribution
This is a work in progress. 
If you see any mistakes in the package (branch `main`) or in the code associated with the analyses for the preprint (branch `preprint`), let us know by logging a bug on the [issues](https://github.com/acefa-hubs/EpiStrainDynamics/issues) page. 

## Code of Conduct
Please note that the EpiStrainDynamics project is released with a [Contributor Code of Conduct](https://acefa-hubs.github.io/EpiStrainDynamics/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
