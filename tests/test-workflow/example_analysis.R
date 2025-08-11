# example analysis with new fit_model function

# Loading required packages
library(ggplot2)
library(rstan)
library(EpiStrainDynamics)
source('../for_comparison/previous_growth_rate.R')
source('../for_comparison/Rt_old.R')
source('../for_comparison/incidence_old.R')

gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)
rw_gi_dist <- function(x) gammaDist(b=1, n=1, x)
ps_gi_dist <- function(x) gammaDist(b=2, n=2, x)

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

rw_subtyped_mod <- construct_model(
  method = random_walk(),
  pathogen_structure = subtyped(
    case_timeseries = influenza_data$ili,
    time = influenza_data$week,
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
rw_subtyped_gr <- growth_rate(rw_subtyped_fit)
rw_subtyped_rt <- Rt(rw_subtyped_fit,
                     gi_dist = rw_gi_dist)
rw_subtyped_inc <- incidence(rw_subtyped_fit, dow = FALSE)
plot(rw_subtyped_inc)
plot(rw_subtyped_gr)
plot(rw_subtyped_rt)

rw_subtyped_inc_dow <- incidence(rw_subtyped_fit, dow = TRUE)
rw_subtyped_prop <- proportion(rw_subtyped_fit,
                               numerator_combination = c('influenzaA.H3N2',
                                                         'influenzaA.H1N1',
                                                         'influenzaB'),
                               denominator_combination = c('influenzaA.H3N2',
                                                           'influenzaA.H1N1',
                                                           'influenzaB',
                                                           'other'))



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
ps_subtyped_gr_orig <- ps_growth_rate_orig(ps_subtyped_fit) |>
  dplyr::arrange(pathogen)
ps_subtyped_gr <- ps_growth_rate(ps_subtyped_fit) |>
  dplyr::arrange(pathogen) |>
  as.data.frame()
ps_subtyped_gr <- ps_subtyped_gr[,colnames(ps_subtyped_gr_orig)]

ps_subtyped_rt_orig <- ps_Rt_orig(ps_subtyped_fit,
                                  gi_dist = ps_gi_dist) |>
  dplyr::arrange(pathogen)
ps_subtyped_rt <- ps_Rt(ps_subtyped_fit,
                        gi_dist = ps_gi_dist) |>
  dplyr::arrange(pathogen) |>
  as.data.frame()
ps_subtyped_rt <- ps_subtyped_rt[,colnames(ps_subtyped_rt_orig)]

ps_subtyped_incidence <- ps_incidence(ps_subtyped_fit)
ps_subtyped_DOW_incidence <- ps_incidence_dow(ps_subtyped_fit)

# ps_sub_prop <- proportion(ps_subtyped_fit$fit,
#                         X = NULL,
#                         method = c('rw', 'ps'),
#                         num_days = length(X),
#                         time_labels,
#                         days_per_knot = 5,
#                         spline_degree = 3,
#                         num_path = 4,
#                         comb_num = list(c(1,2,3), c(1,2), c(3), c(1), c(2)),
#                         comb_den = list(c(1,2,3,4), c(1,2,3,4), c(1,2,3,4), c(1,2), c(1,2) ),
#                         comb_names = c("Influenza", "Influenza A", "Influenza B", "H3N2", "H1N1")


# Example 2 "mp" models
# Fitting the model 'mid-season', only include first 140 days
# Load some simulated data
# use for the sarscov2
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
rw_multiple_gr_orig <- rw_growth_rate_orig(rw_multiple_fit) |>
  dplyr::arrange(pathogen)
rw_multiple_gr <- growth_rate(rw_multiple_fit) |>
  dplyr::arrange(pathogen) |>
  as.data.frame()
rw_multiple_gr <- rw_multiple_gr[,colnames(rw_multiple_gr_orig)]

rw_multiple_incidence <- rw_incidence(rw_subtyped_fit)

# sarscov2 <-
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
ps_multiple_gr_orig <- ps_growth_rate_orig(ps_multiple_fit) |>
  dplyr::arrange(pathogen)
ps_multiple_gr <- ps_growth_rate(ps_multiple_fit) |>
  dplyr::arrange(pathogen) |>
  as.data.frame()
ps_multiple_gr <- ps_multiple_gr[,colnames(ps_multiple_gr_orig)]

ps_multiple_incidence <- ps_incidence(ps_multiple_fit)
ps_multiple_DOW_incidence <- ps_incidence_dow(ps_multiple_fit)

# Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated
# dataset to demonstrate utility during flu-season/pandemic
# Load some COVID data from the UK
covid_data <- read.csv('tests/test-workflow/example_data/ukhsa-covid-data.csv')
covid_data <- covid_data[covid_data$geography == "England", ]

covid_data <- covid_data[order(covid_data$date), ]
covid_data$time <- seq(1, nrow(covid_data))

# rw_single_mod <- construct_model(
#   method = random_walk(),
#   pathogen_structure = single(
#     case_timeseries = covid_data$metric_value
#   ),
#   dow_effect = TRUE
# )
#THIS MODEL GETS STUCK
# rw_single_fit <- fit_model(
#   rw_single_mod,
#   iter = 500,
#   warmup = 300,
#   chains = 3
# )

ps_single_mod <- construct_model(
  method = p_spline(),
  pathogen_structure = single(
    case_timeseries = covid_data$metric_value
  )#,
  # dow_effect = FALSE
)
ps_single_fit <- fit_model(
  ps_single_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_single_mod_dow <- construct_model(
  method = p_spline(),
  pathogen_structure = single(
    case_timeseries = covid_data$metric_value
  ),
  dow_effect = TRUE
)
ps_single_fit_dow <- fit_model(
  ps_single_mod,
  iter = 500,
  warmup = 300,
  chains = 3
)
ps_single_gr_orig <- ps_single_growth_rate_orig(ps_single_fit)
ps_single_gr <- ps_single_growth_rate(ps_single_fit) |>
  as.data.frame()
ps_single_gr <- ps_single_gr[,colnames(ps_single_gr_orig)]

ps_single_rt_orig <- ps_single_Rt_orig(ps_single_fit, gi_dist = ps_gi_dist)
ps_single_rt <- ps_single_Rt(ps_single_fit, gi_dist = ps_gi_dist) |>
  as.data.frame()
ps_single_rt <- ps_single_rt[,colnames(ps_single_rt_orig)]

ps_single_incidence_orig <- ps_single_incidence_orig(ps_single_fit)
ps_single_incidence_dow_orig <- ps_single_incidence_dow(ps_single_fit_dow)
ps_single_incidence <- ps_single_incidence(ps_single_fit) |>
  dplyr::select(-prop) |>
  as.data.frame()
ps_single_incidence_dow <- ps_single_incidence(ps_single_fit_dow) |>
  dplyr::select(-prop) |>
  as.data.frame()

