# example analysis with new fit_model function

# Loading required packages
library(ggplot2)
library(rstan)
library(EpiStrainDynamics)
source('R/fit_model.R')
source('R/get_pathogen_info.R')
source('R/get_model_info.R')
source('R/get_knots.R')

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)


# Example 1 "influenza" models
# create data objects
# Load aus data
influenza_data <- read.csv('example_data/aus_influenza_data.csv')

# Set limits on dates to consider
min_date <- as.Date("2011-08-29")
max_date <- as.Date("2020-03-01")

influenza_data$week <- as.Date(influenza_data$week)
influenza_data <- influenza_data[influenza_data$week < max_date &
                                   influenza_data$week >= min_date, ]
influenza_data <- influenza_data[order(influenza_data$week), ]
influenza_data$time <- seq(1, nrow(influenza_data))

main_pathogens <- list(
  influenzaA = influenza_data$inf_A,
  influenzaB = influenza_data$inf_B,
  other = influenza_data$num_spec - influenza_data$inf_all
)

influenzaA_subtypes <- list(
  influenzaA.H3N2 = influenza_data$inf_H3N2,
  influenzaA.H1N1 = influenza_data$inf_H1N1
)

influenza_data_list <- list(
  cases = influenza_data$ili,
  time = influenza_data$time,
  main_pathogens = main_pathogens,
  influenzaA_subtypes = influenzaA_subtypes
)

# fit
rw_influenza <- fit_model(
  influenza_data_list,
  method = 'random_walk',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_influenza <- fit_model(
  influenza_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  iter = 500,
  warmup = 300,
  chains = 3
)

# Example 2 "mp" models
# Fitting the model 'mid-season', only include first 140 days
# Load some simulated data
sim_data_raw <- read.csv('example_data/simulated_data1.csv')
sim_data <- sim_data_raw[sim_data_raw$t < 140, ]

main_pathogens <- list(
  influenzaA.H3N2 = sim_data$H3N2,
  influenzaA.H1N1 = sim_data$H1N1,
  influenzaB = sim_data$B
)

sim_data_list <- list(cases = sim_data$y,
                      time = sim_data$t,
                      main_pathogens = main_pathogens)

rw_mp <- fit_model(
  sim_data_list,
  method = 'random_walk',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_mp <- fit_model(
  sim_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  iter = 500,
  warmup = 300,
  chains = 3
)

# Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated
# dataset to demonstrate utility during flu-season/pandemic
# Load some COVID data from the UK
covid_data <- read.csv('example_data/ukhsa-covid-data.csv')
covid_data <- covid_data[covid_data$geography == "England", ]

covid_data <- covid_data[order(covid_data$date), ]
covid_data$time <- seq(1, nrow(covid_data))

covid_data_list <- list(cases = covid_data$metric_value,
                        time = covid_data$time)

rw_single <- fit_model(
  covid_data_list,
  method = 'random_walk',
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_single <- fit_model(
  covid_data_list,
  method = 'p-spline',
  spline_degree = 3,
  days_per_knot = 5,
  iter = 500,
  warmup = 300,
  chains = 3
)

