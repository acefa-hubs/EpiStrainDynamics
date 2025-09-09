# EpiStrainDynamics 
[![codecov](https://codecov.io/gh/acefa-hubs/EpiStrainDynamics/graph/badge.svg)](https://app.codecov.io/gh/acefa-hubs/EpiStrainDynamics)
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/acefa-hubs/EpiStrainDynamics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/acefa-hubs/EpiStrainDynamics/actions/workflows/R-CMD-check.yaml)

`EpiStrainDynamics` is a statistical modelling framework capable of inferring trends of multiple pathogens. 
Estimating the temporal trends in infectious disease activity is crucial for monitoring disease spread and the impact of interventions. 
Surveillance indicators routinely collected to monitor these trends are often a composite of multiple pathogens. 
For example, ‘influenza-like illness’ — routinely monitored as a proxy for influenza infections — is a symptom definition that could be caused by a wide range of pathogens, including multiple subtypes of influenza, SARS-CoV-2, and RSV. 
Inferred trends from such composite time series may not reflect the trends of any one of the component pathogens, each of which can exhibit distinct dynamics. 
Although many surveillance systems routinely test a subset of individuals contributing to a surveillance indicator — providing information on the relative contribution of the component pathogens — trends may be obscured by time-varying testing rates or substantial noise in the observation process. 

`EpiStrainDynamics` builds on existing modelling frameworks built to handle two pathogens, and extends them be to able to:

- infer the trends of any number of pathogens
- fit to time series data of counts (eg, daily number of cases)
- incorporate influenza testing data in which the subtype for influenza A samples may be undetermined
- account for day-of-the-week effects in daily data
- include options for fitting penalized splines or random walks 
- support additional (optional) correlation structures in the parameters describing the smoothness of the penalized splines (or random walks), and 
- account for additional (optional) sources of noise in the observation process. 

## Installation

You can install `EpiStrainDynamics` with: 

```
devtools::install_github('acefa-hubs/EpiStrainDynamics')
```

## Using `EpiStrainDynamics`
Detailed instructions can be found on the vignette on the [website](https://acefa-hubs.github.io/EpiStrainDynamics/articles/Using-EpiStrainDyamics.html). 
Here we provide a short overview. 

### Step 1: Construct model

Modelling specifications are provided using the `construct_model` function. The correct stan model is then applied based on the specifications provided. 
`construct_model()` takes three arguments: `method`, `pathogen_structure`, and `dow_effect`. 

#### `method`

EpiStrainDynamics has pre-compiled stan models that fit either with bayesian penalised splines or random walks. These are specified using the `method` argument of `construct_model()` as functions, either with `random_walk()` or `p_spline()`. The penalised spline model has two further options to specify: `spline_degree` is the polynomial degree of the individual spline segments used to construct the overall curve (must be a positive whole number) and `days_per_knot`, which is the number of days for each knot (must also be a positive whole number). 

So we may specify: 
```
method = random_walk(), 

# OR
method = p_spline(spline_degree = 3, days_per_knot = 2),  # example options
```

#### `pathogen_structure`

There are three main types of pathogen structure available to model: `single()`, `multiple()`, and `subtyped()`. 

The `single()` pathogen structure is the simplest and models a single pathogen timeseries. 
A vector of the case timeseries is provided to the argument `case_timeseries`, a vector of date or time labels is provided to `time`, and optionally a pathogen name can be assigned with `pathogen_name`. 
The `multiple()` and `subtyped()` pathogen structures both also have the `case_timeseries` and `time` arguments, as in the `single()` structure. 
In addition to these two arguments, instead of a single argument for `pathogen_name`, `multiple()` and `subtyped()` have one or more arguments that allow the user to define the names and data for the different components or subtypes to be modelled. 
For `multiple()`, these are specified as a named list with the argument `component_pathogen_timeseries`. 
The names in this list will be using in subsequent plotting. 
For `subtyped()`, which allows the user to incorporate testing data for influenza A subtypes, a vector of total unsubtyped influenza A cases is specified with `influenzaA_unsubtyped_timeseries`, a named list of the subtyped influenza A timeseries is provided to `influenzaA_subtyped_timeseries`, and further pathogens are provided in a named list to `other_pathogen_timeseries`. 

In addition to specifying the data and pathogen names, `multiple()` and `subtyped()` both allow the user to modify the correlation structures in the parameters describing the smoothness (with argument `smoothing_structure`) and account for additional sources of noise in the observation process (`observation_noise`). `smoothing_structure` is either `'shared'` (all pathogens have the same smoothness), `'independent'` (each pathogen has completely independent smoothing structure), or `'correlated'` (smoothing structure is correlated among pathogens). 
`observation_noise` is either `'observation_noise_only'` (only includes observation noise - the same between pathogens) or `'pathogen_specific_noise'` (includes noise in individual pathogens as well). 

#### `dow_effect`

Day of week effect is specified as a logical (`TRUE` or `FALSE`) to the `dow_effect` argument. In plotting the day of week effect can be selectively removed. 

Altogether, an example constructed model for a random walk model with multiple pathogens might look as follows, illustrated using data provided with the package `sarscov2`: 

```
mod <- construct_model(
  method = random_walk(),                   # specify `random_walk` method
  
  pathogen_structure = multiple(            # specify `multiple` pathogen structure
   case_timeseries = sarscov2$cases,        # timeseries of case data
   time = sarscov2$date,                    # date or time variable labels
   
   component_pathogen_timeseries = list(    # named list of component pathogens
     alpha = sarscov2$alpha,
     delta = sarscov2$delta,
     omicron = sarscov2$omicron,
     other = sarscov2$other
   ),
   
   smoothing_structure = 'independent',             # correlation structure
   observation_noise = 'observation_noise_only'     # observation noise
 ),
  dow_effect = TRUE                         # day of week effect
)

```

### Step 2: fit model

The model estimates the expected value of the time series (eg, a smoothed trend in the daily number of cases accounting for noise) for each individual pathogen. 
Model parameterisation decisions specified when constructing the model in step 1 mean the correct stan model will be applied at this stage by simply calling:

```
fit <- fit_model(mod)
```

### Step 3: Calculate and visualise epidemiological quantities

Calculate epidemic growth rate with `growth_rate(fit)`, effective reproduction number over time with `Rt(fit, gi_dist = X)` (requiring specification of a generation interval distribution X), incidence with or without a day of week effect with `incidence(fit, dow_effect = TRUE)`, and proportions of different combinations of cases attributable to different pathogens/subtypes using `proportion()`.
```
prop <- proportion(fit)
plot(prop)
```


For a more detailed discussion, check out the [vignette](https://acefa-hubs.github.io/EpiStrainDynamics/articles/Using-EpiStrainDyamics.html).  



## Citation
When using this package, please cite both the package and the research article underlying the statistical model developments: 

Windecker S, Eales O (2025). EpiStrainDynamics: Infer temporal trends of multiple pathogens. R package version 0.0.0.9000, https://acefa-hubs.github.io/EpiStrainDynamics/.
Oliver Eales, Saras M Windecker, James M McCaw, Freya M Shearer, Inferring temporal trends of multiple pathogens, variants, subtypes or serotypes from routine surveillance data, American Journal of Epidemiology, 2025;, kwaf119, https://doi.org/10.1093/aje/kwaf119
For code corresponding to this paper, see branch [`paper_analysis`](https://github.com/acefa-hubs/EpiStrainDynamics/tree/paper_analysis). 

## Contribution
This is a work in progress. 
If you see any mistakes in the package (branch `main`) or in the code associated with the analyses for the paper (branch `paper_analysis`), let us know by logging a bug on the [issues](https://github.com/acefa-hubs/EpiStrainDynamics/issues) page. 

## Code of Conduct
Please note that the `EpiStrainDynamics` project is released with a [Contributor Code of Conduct](https://acefa-hubs.github.io/EpiStrainDynamics/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## Support
This project was supported by the Australia-Aotearoa Consortium for Epidemic Forecasting and Analytics.

<a href="https://acefa-hubs.github.io/"><img src="man/figures/ACEFA.png" align = "center" height="150" alt="EpiStrainDynamics website" /></a>
