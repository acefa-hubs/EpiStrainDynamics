# example analysis with new fit_model function

# Loading required packages
library(ggplot2)
library(rstan)
library(EpiStrainDynamics)
source('R/fit_model.R')
source('R/construct_model.R')
source('R/get_model_type.R')
source('R/get_knots.R')
# source('R/growth_rate.R')
# source('R/predict_B_true.R')

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

rw_subtyped_mod <- construct_model(
  method = random_walk(),
  pathogen_structure = subtyped(
    case_timeseries = influenza_data$ili,
    influenzaA_unsubtyped_timeseries = influenza_data$inf_A,
    other_pathogen_timeseries = list(
      influenzaB = influenza_data$inf_B,
      other = influenza_data$num_spec - influenza_data$inf_all
    ),
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza_data$inf_H3N2,
      influenzaA.H1N1 = influenza_data$inf_H1N1
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)
rw_subtyped_fit <- fit_model(
  rw_subtyped_mod,
  iter = 2000,
  warmup = 1000,
  chains = 3
)

ps_subtyped_mod <- construct_model(
  method = p_spline(),
  pathogen_structure = subtyped(
    case_timeseries = influenza_data$ili,
    influenzaA_unsubtyped_timeseries = influenza_data$inf_A,
    influenzaA_subtyped_timeseries = list(
      influenzaA.H3N2 = influenza_data$inf_H3N2,
      influenzaA.H1N1 = influenza_data$inf_H1N1
    ),
    other_pathogen_timeseries = list(
      influenzaB = influenza_data$inf_B,
      other = influenza_data$num_spec - influenza_data$inf_all
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)
ps_subtyped_fit <- fit_model(
  ps_subtyped_mod,
  iter = 2000,
  warmup = 1000,
  chains = 3
)
# ps_subtyped_gr <- ps_growth_rate(
#   ps_influenza,
#   influenza_data_list$time,
#   time_labels = influenza_data_list$time
# )


# Example 2 "mp" models

# knots <- get_knots(df$t, days_per_knot = 5, spline_degree = 3)
#
# ps_mp_data <- list(
#   num_data = nrow(df),
#   num_path = 11,
#   num_knots = length(knots),
#   knots = knots,
#   spline_degree = 3,
#   Y = df$cases,
#   X = df$t,
#   P = t(df[3:13]),
#   # Might have to fix this line
#   week_effect = 7,
#   DOW = (df$t %% 7) + 1,
#   cov_structure = 0,
#   noise_structure = 0
# )
#
# ps_mp_fit <- sampling(
#   ps_mp_mod,
#   iter = 5000,
#   warmup = 1000,
#   chains = 4,
#   data = ps_mp_data
# )

ps_multiple_mod <- construct_model(
  method = p_spline(spline_degree = 3, days_per_knot = 5),
  pathogen_structure = multiple(
    case_timeseries = sarscov2$cases,
    component_pathogen_timeseries = list(
      B.1.177 = sarscov2$B.1.177,
      B.1.1.7 = sarscov2$B.1.1.7,
      B.1.617.2 = sarscov2$B.1.617.2,
      BA.1 = sarscov2$BA.1,
      BA.2 = sarscov2$BA.2,
      BA.2.75 = sarscov2$BA.2.75,
      BA.4 = sarscov2$BA.4,
      BA.5 = sarscov2$BA.5,
      BQ.1 = sarscov2$BQ.1,
      XBB = sarscov2$XBB,
      Other = sarscov2$Other
    ),
    smoothing_structure = 'single',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)

# Fitting the model 'mid-season', only include first 140 days
# Load some simulated data
sim_data_raw <- read.csv('tests/test-workflow/example_data/simulated_data1.csv')
sim_data <- sim_data_raw[sim_data_raw$t < 140, ]

rw_multiple_mod <- construct_model(
  method = random_walk(),
  pathogen_structure = multiple(
    case_timeseries = sim_data$y,
    component_pathogen_timeseries = list(
      influenzaA.H3N2 = sim_data$H3N2,
      influenzaA.H1N1 = sim_data$H1N1,
      influenzaB = sim_data$B
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)
rw_multiple_fit <- fit_model(
  rw_multiple_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)

ps_multiple_mod <- construct_model(
  method = p_spline(spline_degree = 3, days_per_knot = 5),
  pathogen_structure = multiple(
    case_timeseries = sim_data$y,
    component_pathogen_timeseries = list(
      influenzaA.H3N2 = sim_data$H3N2,
      influenzaA.H1N1 = sim_data$H1N1,
      influenzaB = sim_data$B
    ),
    smoothing_structure = 'independent',
    observation_noise = 'observation_noise_only'
  ),
  dow_effect = TRUE
)
ps_multiple_fit <- fit_model(
  ps_multiple_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)
# ps_mp_gr <- ps_growth_rate(
#   ps_mp,
#   sim_data_list$time,
#   time_labels = sim_data_list$time
# )

# Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated
# dataset to demonstrate utility during flu-season/pandemic
# Load some COVID data from the UK
covid_data <- read.csv('tests/test-workflow/example_data/ukhsa-covid-data.csv')
covid_data <- covid_data[covid_data$geography == "England", ]

covid_data <- covid_data[order(covid_data$date), ]
covid_data$time <- seq(1, nrow(covid_data))

rw_single_mod <- construct_model(
  method = random_walk(),
  pathogen_structure = single(
    case_timeseries = covid_data$metric_value
  ),
  dow_effect = TRUE
)
#THIS MODEL GETS STUCK
rw_single_fit <- fit_model(
  rw_single_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)

ps_single_mod <- construct_model(
  method = p_spline(),
  pathogen_structure = single(
    case_timeseries = covid_data$metric_value
  ),
  dow_effect = TRUE
)
ps_single_fit <- fit_model(
  ps_single_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)
# ps_single_gr <- ps_single_growth_rate(
#   ps_single,
#   covid_data_list$time,
#   time_labels = covid_data_list$time
# )
