# example analysis with new fit_model function

# Loading required packages
library(ggplot2)
library(rstan)
library(EpiStrainDynamics)
source('R/fit_model.R')
source('R/get_pathogen_info.R')

# Load aus data
df1 <- read.csv('example_data/aus_influenza_data.csv')

# Set limits on dates to consider
min_date <- as.Date("2011-08-29")
max_date <- as.Date("2020-03-01")

df1$week <- as.Date(df1$week)

df1 <- df1[df1$week<max_date & df1$week>=min_date,]

df1 <- df1[order(df1$week),]

df1$time <- seq(1, nrow(df1))

# Load some simulated data
dfS1 <- read.csv('example_data/simulated_data1.csv')

# Load some COVID data from the UK
dfC1 <- read.csv('example_data/ukhsa-covid-data.csv')
dfC1 <- dfC1[dfC1$geography=="England",]

dfC1 <- dfC1[order(dfC1$date),]
dfC1$time <- seq(1, nrow(dfC1))

# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)


# Example 1
# create data objects
cases <- df1$ili

main_pathogens <- list(
  influenzaA = df1$inf_A,
  influenzaB = df1$inf_B,
  influenzaOther = df1$num_spec - df1$inf_all)

influenzaA_subtypes <- list(
  influenzaA.H3N2 = df1$inf_H3N2,
  influenzaA.H1N1 = df1$inf_H1N1)


data <- list(cases = cases,
             main_pathogens = main_pathogens,
             influenzaA_subtypes = influenzaA_subtypes)

# fit
out1 <- fit_model(data,
                 method = 'random_walk',
                 iter = 500,
                 warmup = 300,
                 chains = 3)
