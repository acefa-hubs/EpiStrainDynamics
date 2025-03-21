# example analysis with new fit_model function

# Loading required packages
library(ggplot2)
library(rstan)
library(EpiStrainDynamics)
# source('R/fit_model.R')
source('R/get_pathogen_info.R')
source('R/get_model_type.R')
source('R/get_knots.R')
source('R/growth_rate.R')
source('R/predict_B_true.R')
source('R/construct_model_data.R')

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# Example 1 "influenza" models
# create data objects
# Load aus data
influenza_data <- read.csv('tests/test-workflow/example_data/aus_influenza_data.csv')

# Set limits on dates to consider
min_date <- as.Date("2011-08-29")
max_date <- as.Date("2020-03-01")

influenza_data$week <- as.Date(influenza_data$week)
influenza_data <- influenza_data[influenza_data$week < max_date &
                                   influenza_data$week >= min_date, ]
influenza_data <- influenza_data[order(influenza_data$week), ]
influenza_data$time <- seq(1, nrow(influenza_data))

#option 1
influenza_data_list <- construct_model_data(
  total_cases = influenza_data$ili,
  time = influenza_data$time,
  component_pathogens = list(
    influenzaA = influenza_data$inf_A,
    influenzaB = influenza_data$inf_B,
    other = influenza_data$num_spec - influenza_data$inf_all
  ),
  influenzaA_subtypes = list(
    influenzaA.H3N2 = influenza_data$inf_H3N2,
    influenzaA.H1N1 = influenza_data$inf_H1N1
  )
)


source('R/test-fit_model.R')
option1_fit <- fit_model(
  influenza_data_list,
  options = model_options(
    method = 'random-walk',
    pathogen_structure = 'subtyped',
    spline_degree = NULL,
    days_per_knot = NULL,
    dow_effect = TRUE,
    smoothing_structure = 'single',
    observation_noise = 'observation_noise_only',
    iter = 1000,
    warmup = 100,
    chains = 2
  )
)

#option 2
source('R/test-fit_model-2.R')
my_model <- construct_model(
  method = random_walk(), # p_spline(spline_degree = 5, days_per_knot = 5, time_data = time)
  pathogen_structure = subtyped(
    total_case_data = influenza_data$ili,
    influenzaA_data = influenza_data$inf_A,
    other_component_pathogen_data = list(
      influenzaB = influenza_data$inf_B,
      other = influenza_data$num_spec - influenza_data$inf_all
    ),
    influenzaA_subtype_data = list(
      influenzaA.H3N2 = influenza_data$inf_H3N2,
      influenzaA.H1N1 = influenza_data$inf_H1N1
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  ),
  dow_data = influenza_data$time
)
option2_fit <- fit_model(
  constructed_model,
  iter = 2000,
  warmup = 1000,
  chains = 3
)






#### NOT YET CHANGED
# fit
rw_influenza <- fit_model(
  influenza_data_list,
  method = 'random_walk',
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_influenza <- fit_model(
  influenza_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_influenza_gr <- ps_growth_rate(
  ps_influenza,
  influenza_data_list$time,
  time_labels = influenza_data_list$time
)

# Example 2 "mp" models
# Fitting the model 'mid-season', only include first 140 days
# Load some simulated data
sim_data_raw <- read.csv('tests/test-workflow/example_data/simulated_data1.csv')
sim_data <- sim_data_raw[sim_data_raw$t < 140, ]

sim_data_list <- list(
  cases = sim_data$y,
  time = sim_data$t,
  component_pathogens = list(
    influenzaA.H3N2 = sim_data$H3N2,
    influenzaA.H1N1 = sim_data$H1N1,
    influenzaB = sim_data$B
  )
)

rw_mp <- fit_model(
  sim_data_list,
  method = 'random_walk',
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_mp <- fit_model(
  sim_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)

ps_mp_gr <- ps_growth_rate(
  ps_mp,
  sim_data_list$time,
  time_labels = sim_data_list$time
)


# Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated
# dataset to demonstrate utility during flu-season/pandemic
# Load some COVID data from the UK
covid_data <- read.csv('tests/test-workflow/example_data/ukhsa-covid-data.csv')
covid_data <- covid_data[covid_data$geography == "England", ]

covid_data <- covid_data[order(covid_data$date), ]
covid_data$time <- seq(1, nrow(covid_data))

covid_data_list <- list(
  cases = covid_data$metric_value,
  time = covid_data$time
)

# very very slow?
rw_single <- fit_model(
  covid_data_list,
  method = 'random_walk',
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_single <- fit_model(
  covid_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  smoothing_structure = 'independent',
  iter = 500,
  warmup = 300,
  chains = 3
)

ps_single_gr <- ps_single_growth_rate(
  ps_single,
  covid_data_list$time,
  time_labels = covid_data_list$time
)

