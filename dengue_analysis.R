# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)
library(RColorBrewer)
library(DescTools)
library(lubridate)
library(tidyverse)

# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')
###############################################################################

# Read in Taiwan dengue data and format it into new data frame df
df1 <- read.csv('Data/Dengue_Daily.csv')

df1 <- df1[c(1, 3, 22)]
colnames(df1) <- c("Onset_date", "Notification_date", "Serotype")
df1$Onset_date <- as.Date(df1$Onset_date)
df1$Notification_date <- as.Date(df1$Notification_date)
table(df1$Serotype)
df1$Serotype <- factor(df1$Serotype)
levels(df1$Serotype) <- c("No data",
                          "Serotype 1",
                          "Serotype 3",
                          "Serotype 2",
                          "Serotype 4")

df <- data.frame()
onset_dates <- seq.Date(min(df1$Onset_date), max(df1$Onset_date), by = 1)

for (i in 1:length(onset_dates)) {
  df_tmp <- df1[df1$Onset_date == onset_dates[i], ]

  table(df_tmp$Serotype)

  row_df <- data.frame(
    date = onset_dates[i],
    cases = nrow(df_tmp),
    sero1 = nrow(df_tmp[df_tmp$Serotype == "Serotype 1", ]),
    sero2 = nrow(df_tmp[df_tmp$Serotype == "Serotype 2", ]),
    sero3 = nrow(df_tmp[df_tmp$Serotype == "Serotype 3", ]),
    sero4 = nrow(df_tmp[df_tmp$Serotype == "Serotype 4", ]),
    unsero = nrow(df_tmp[df_tmp$Serotype == "No data", ])
  )

  df <- rbind(df, row_df)

}

dfR <- df

###############################################################################
# Subset data into specific periods to fit model to
################################################################
## Single season analysis (2023 season)
# Set limits on dates to consider
min_date <- as.Date("2023-04-01")#as.Date("2006-01-01")
max_date <- as.Date("2024-04-01")#as.Date("2012-06-01")

dfS <- df[df$date < max_date & df$date >= min_date, ]
dfS$time <- seq(1, nrow(dfS))
################################################################
## Many season analysis 2006-2015 seasons
# Set limits on dates to consider
min_date <- as.Date("2006-04-01")#as.Date("2006-01-01")
max_date <- as.Date("2016-04-01")#as.Date("2016-04-01")

df <- df[df$date < max_date & df$date >= min_date, ]
df$time <- seq(1, nrow(df))

###############################################################################
# Set stan options
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Loading stan model used for dengue analysis
ps_mp_modV2 <- stan_model('stan/ps_mp_finalV2.stan')
ps_mp_modS <- stan_model('stan/ps_mp_final - CopyWithPrior.stan')

###############################################################################
# Fit model to 2006-2015 seasons

knots <- get_knots(df$time, days_per_knot = 5, spline_degree = 3)

ps_mp_data <- list(
  num_data = nrow(df),
  num_path = 4,
  num_knots = length(knots),
  knots = knots,
  spline_degree = 3,
  Y = df$cases,
  X = df$time,
  P = t(matrix(
    data = c(df$sero1, df$sero2, df$sero3, df$sero4),
    ncol = 4
  )),
  week_effect = 1,
  DOW = (df$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 0
)

ps_mp_fit <- sampling(
  ps_mp_modV2,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data
)
saveRDS(ps_mp_fit, "FitStanModels/dengue_fit4.rds")

###############################################################################
# Extract posterior distributions to include as priors in 2023 season analysis

# Re-reading model fit saved above
ps_mp_fit <- readRDS("FitStanModels/dengue_fit4.rds")
# Extracting posterior
post <- rstan::extract(ps_mp_fit)

# Get mean and standard deviation of tau and phi parameters
tau_mn <- apply(post$tau, 2, mean)
tau_sd <- apply(post$tau, 2, sd)
phi_mn <- c(mean(post$phi))
phi_sd <- c(sd(post$phi))

###############################################################################
# Fit model to 2023 season using informed priors

knotsS <- get_knots(dfS$time, days_per_knot = 5, spline_degree = 3)

ps_mp_dataS <- list(
  num_data = nrow(dfS),
  num_path = 4,
  num_knots = length(knotsS),
  knots = knotsS,
  spline_degree = 3,
  Y = dfS$cases,
  X = dfS$time,
  P = t(matrix(
    data = c(dfS$sero1, dfS$sero2, dfS$sero3, dfS$sero4),
    ncol = 4
  )),
  week_effect = 1,
  DOW = (dfS$time %% 1) + 1,
  cov_structure = 1,
  tau_mn = tau_mn,
  tau_sd = tau_sd,
  phi_mn = phi_mn,
  phi_sd = phi_sd
)

ps_mp_fitS <- sampling(
  ps_mp_modS,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_dataS
)
saveRDS(ps_mp_fitS, "FitStanModels/dengue_fitS.rds")

###############################################################################
# Making the figures
###############################################################################
# Figure 4
###############################################################################
# Read in model
ps_mp_fit <- readRDS("FitStanModels/dengue_fit4.rds")

# Get incidence curves
mod_inc <- ps_incidence(
  ps_mp_fit,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Get time-varying growth rates
mod_gr <- ps_growth_rate(
  ps_mp_fit,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Add mask to values for periods when modelled cases < 0.5
mod_gr$mask <- 0
mod_inc$mask <- 0
mod_gr$group <- 0

pathogens <- c("1", "2", "3", "4")
for (i in 1:length(pathogens)) {
  index <- mod_inc[mod_inc$pathogen == pathogens[i], ]$y > 0.5
  dates <- mod_inc[mod_inc$pathogen == pathogens[i], ][index, ]$time
  mod_gr[mod_gr$pathogen == pathogens[i] &
           mod_gr$time %in% dates, ]$mask <- 1
  mod_inc[index, ]$mask <- 1
  print(pathogens[i])
  #mod_gr<-mod_gr[ (!(mod_gr$pathogen==pathogens[i]))
  # | (mod_gr$time>min_date & mod_gr$time<max_date) ,]
  grp <- 1
  for (j in 1:(length(mod_gr[mod_gr$pathogen == pathogens[i] &
                             mod_gr$time %in% dates, ]$mask) - 1)) {
    mod_gr[mod_gr$pathogen == pathogens[i] &
             mod_gr$time %in% dates, ]$group[j] <- grp

    if (mod_gr[mod_gr$pathogen == pathogens[i] &
               mod_gr$time %in% dates, ]$time[j + 1] -
        mod_gr[mod_gr$pathogen == pathogens[i] &
               mod_gr$time %in% dates, ]$time[j] != 1) {
      grp <- grp + 1
    }
  }
  mod_gr[mod_gr$pathogen == pathogens[i] &
           mod_gr$time %in% dates, ]$group[j + 1] <- grp
}

# Dividing modelled values after 2014-07-01 by 10 - this is only for
# plotting purposes (second axes is added)
mod_inc$cat <- 0
mod_inc[mod_inc$time >= as.Date("2014-07-01"), ]$cat <- 1
mod_inc_tmp <- mod_inc
mod_inc_tmp[mod_inc_tmp$cat == 1, ]$y <-
  mod_inc_tmp[mod_inc_tmp$cat == 1, ]$y / 10
mod_inc_tmp[mod_inc_tmp$cat == 1, ]$lb_50 <-
  mod_inc_tmp[mod_inc_tmp$cat == 1, ]$lb_50 / 10
mod_inc_tmp[mod_inc_tmp$cat == 1, ]$ub_50 <-
  mod_inc_tmp[mod_inc_tmp$cat == 1, ]$ub_50 / 10
mod_inc_tmp[mod_inc_tmp$cat == 1, ]$lb_95 <-
  mod_inc_tmp[mod_inc_tmp$cat == 1, ]$lb_95 / 10
mod_inc_tmp[mod_inc_tmp$cat == 1, ]$ub_95 <-
  mod_inc_tmp[mod_inc_tmp$cat == 1, ]$ub_95 / 10

df$cases_tmp <- df$cases
df[df$date >= as.Date("2014-07-01"), ]$cases_tmp <-
  df[df$date >= as.Date("2014-07-01"), ]$cases / 10

# Creating subplots for figure
plt1 <- ggplot(mod_inc_tmp[mod_inc_tmp$pathogen == "Total", ]) +
  geom_line(aes(x = time, y = y)) +
  geom_point(data = df, aes(x = date, y = cases_tmp), size = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95
  ), alpha = 0.2) +
  theme_bw(base_size = 14) +
  ylab("Cases") +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  geom_vline(xintercept = as.Date("2014-07-01")) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 60,
    label = "600",
    hjust = 1.4
  ) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 40,
    label = "400",
    hjust = 1.4
  ) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 20,
    label = "200",
    hjust = 1.4
  ) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 60,
    y = 60
  )) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 40,
    y = 40
  )) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 20,
    y = 20
  )) +
  coord_cartesian(xlim = c(min_date + 80, max_date - 80))

plt2 <- ggplot(mod_inc_tmp[mod_inc_tmp$pathogen != "Total", ]) +
  geom_line(aes(x = time, y = y, color = pathogen)) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = pathogen
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = pathogen
  ), alpha = 0.2) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_fill_brewer("Serotype", palette = "Dark2") +
  theme_bw(base_size = 14) +
  ylab("Modelled cases") +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  geom_vline(xintercept = as.Date("2014-07-01")) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 60,
    label = "600",
    hjust = 1.4
  ) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 40,
    label = "400",
    hjust = 1.4
  ) +
  annotate(
    'text',
    x = as.Date("2014-07-01"),
    y = 20,
    label = "200",
    hjust = 1.4
  ) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 60,
    y = 60
  )) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 40,
    y = 40
  )) +
  geom_segment(aes(
    x = as.Date("2014-07-01") - 20,
    xend = as.Date("2014-07-01"),
    yend = 20,
    y = 20
  )) +
  coord_cartesian(xlim = c(min_date + 80, max_date - 80))

plt3 <- ggplot(mod_gr[mod_gr$mask == 1, ]) +
  geom_line(aes(
    x = time,
    y = y,
    color = pathogen,
    group = interaction(pathogen, group)
  )) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = pathogen,
    group = interaction(pathogen, group)
  ),
  alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = pathogen,
    group = interaction(pathogen, group)
  ),
  alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_fill_brewer("Serotype", palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab("Date") +
  ylab("Growth rate") +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  coord_cartesian(xlim = c(min_date + 80, max_date - 80))

plt1 <- plt1 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

plt2 <- plt2 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.background = element_rect(color = "black"),
    legend.position = c(0.6, 0.8),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))

# Plot figure and save
plt1 / plt2 / plt3
ggsave("Figures/Dengue1.png", width = 10, height = 10)


###############################################################################
# Figure 5
###############################################################################
ps_mp_fitS <- readRDS("FitStanModels/dengue_fitS.rds")

# Get incidence curves
mod_incS <- ps_incidence(
  ps_mp_fitS,
  dfS$time,
  num_days = nrow(dfS),
  time_labels = dfS$date,
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Get time-varying growth rates
mod_grS <- ps_growth_rate(
  ps_mp_fitS,
  dfS$time,
  num_days = nrow(dfS),
  time_labels = dfS$date,
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Assign mask variable for periods in which modelled cases is < 0.5
mod_grS$mask <- 0
mod_incS$mask <- 0

pathogens <- c("1", "2", "3", "4")
for (i in 1:length(pathogens)) {
  print(i)

  index <- mod_incS[mod_incS$pathogen == pathogens[i], ]$y > 0.5
  dates <- mod_incS[mod_incS$pathogen == pathogens[i], ][index, ]$time
  if (TRUE %in% index) {
    mod_grS[mod_grS$pathogen == pathogens[i] &
              mod_grS$time %in% dates, ]$mask <- 1
    mod_incS[index, ]$mask <- 1
  }

}

# Set colour scheme
cols <- c(brewer.pal(4, "Dark2"), "black")

# Plot all the subfigures
plt1 <- ggplot(mod_incS) +
  geom_line(aes(x = time, y = y, color = pathogen)) +
  geom_point(data = dfS, aes(x = date, y = cases), size = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = pathogen
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = pathogen
  ), alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_color_manual("Serotype", values = cols) +
  scale_fill_manual("Serotype", values = cols) +
  ylab("Cases") +
  xlab("Date") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(color = "black")
  ) +
  coord_cartesian(xlim = c(as.Date("2023-04-01") + 12, as.Date("2024-04-01") -
                             12))

plt2 <- ggplot(mod_grS[mod_grS$mask == 1, ]) +
  geom_line(aes(x = time, y = y, color = pathogen)) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = pathogen
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = pathogen
  ), alpha = 0.2) +
  theme_bw(base_size = 14) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_fill_brewer("Serotype", palette = "Dark2") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  xlab("Date") +
  ylab("Growth rate") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  coord_cartesian(xlim = c(as.Date("2023-04-01") + 12, as.Date("2024-04-01") -
                             12))

plt1 <- plt1 +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = c(0.8, 0.8)
  )

plt2 <- plt2 +
  theme(
    legend.background = element_rect(color = "black"),
    legend.position = "none",
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  guides(fill = guide_legend(nrow = 1), color = guide_legend(nrow = 1))

# Plot the main figure and save
plt1 / plt2
ggsave("Figures/DengueS.png", width = 10, height = 10)

###############################################################################
# S Fig 14
###############################################################################

source('R/SFig14_functions.R')

# Extract posterior of smoothed values for each serotype (a)
post <- rstan::extract(ps_mp_fit)
a <- ps_get_a(
  post,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2006-04-01"),
  max_time = as.Date("2007-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

postS <- rstan::extract(ps_mp_fitS)
aS <- ps_get_a(
  postS,
  dfS$time,
  num_days = nrow(dfS),
  time_labels = dfS$date,
  min_time = as.Date("2006-04-01"),
  max_time = as.Date("2007-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Calculating cumulative incidence for each serotype for each season
cum_inc1 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2006-04-01"),
  max_time = as.Date("2007-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc2 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2007-04-01"),
  max_time = as.Date("2008-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc3 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2008-04-01"),
  max_time = as.Date("2009-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc4 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2009-04-01"),
  max_time = as.Date("2010-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc5 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2010-04-01"),
  max_time = as.Date("2011-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc6 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2011-04-01"),
  max_time = as.Date("2012-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc7 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2012-04-01"),
  max_time = as.Date("2013-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc8 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2013-04-01"),
  max_time = as.Date("2014-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc9 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2014-04-01"),
  max_time = as.Date("2015-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc10 <- ps_cumulative_incidence(
  a,
  df$time,
  num_days = nrow(df),
  time_labels = df$date,
  min_time = as.Date("2015-04-01"),
  max_time = as.Date("2016-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)
cum_inc11 <- ps_cumulative_incidence(
  aS,
  dfS$time,
  num_days = nrow(dfS),
  time_labels = dfS$date,
  min_time = as.Date("2023-04-01"),
  max_time = as.Date("2024-03-31"),
  num_path = 4,
  days_per_knot = 5,
  pathogen_names = c("1", "2", "3", "4")
)

# Extract estimate of cumulative incidence by serotype for each season
# using raw data instead
raw_inc1 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2006-04-01"),
                                     max_time = as.Date("2007-03-31"))
raw_inc2 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2007-04-01"),
                                     max_time = as.Date("2008-03-31"))
raw_inc3 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2008-04-01"),
                                     max_time = as.Date("2009-03-31"))
raw_inc4 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2009-04-01"),
                                     max_time = as.Date("2010-03-31"))
raw_inc5 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2010-04-01"),
                                     max_time = as.Date("2011-03-31"))
raw_inc6 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2011-04-01"),
                                     max_time = as.Date("2012-03-31"))
raw_inc7 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2012-04-01"),
                                     max_time = as.Date("2013-03-31"))
raw_inc8 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2013-04-01"),
                                     max_time = as.Date("2014-03-31"))
raw_inc9 <- raw_cumulative_incidence(df,
                                     min_time = as.Date("2014-04-01"),
                                     max_time = as.Date("2015-03-31"))
raw_inc10 <- raw_cumulative_incidence(df,
                                      min_time = as.Date("2015-04-01"),
                                      max_time = as.Date("2016-03-31"))
raw_inc11 <- raw_cumulative_incidence(dfS,
                                      min_time = as.Date("2023-04-01"),
                                      max_time = as.Date("2024-03-31"))

# Labelling and merging data frames
cum_inc1$year <- "2006"
cum_inc2$year <- 2007
cum_inc3$year <- 2008
cum_inc4$year <- 2009
cum_inc5$year <- 2010
cum_inc6$year <- 2011
cum_inc7$year <- 2012
cum_inc8$year <- 2013
cum_inc9$year <- 2014
cum_inc10$year <- 2015
cum_inc11$year <- 2023

cum_inc <- rbind(
  cum_inc1,
  cum_inc2,
  cum_inc3,
  cum_inc4,
  cum_inc5,
  cum_inc6,
  cum_inc7,
  cum_inc8,
  cum_inc9,
  cum_inc10,
  cum_inc11
)

raw_inc1$year <- "2006"
raw_inc2$year <- 2007
raw_inc3$year <- 2008
raw_inc4$year <- 2009
raw_inc5$year <- 2010
raw_inc6$year <- 2011
raw_inc7$year <- 2012
raw_inc8$year <- 2013
raw_inc9$year <- 2014
raw_inc10$year <- 2015
raw_inc11$year <- 2023

raw_inc <- rbind(
  raw_inc1,
  raw_inc2,
  raw_inc3,
  raw_inc4,
  raw_inc5,
  raw_inc6,
  raw_inc7,
  raw_inc8,
  raw_inc9,
  raw_inc10,
  raw_inc11
)

ggplot(cum_inc[cum_inc$label == "proportion", ]) +
  geom_point(
    aes(
      x = year,
      y = y,
      color = pathogen),
    position = position_dodge(width = .5)) +
  geom_errorbar(
    aes(
      x = year,
      y = y,
      ymin = lb_95,
      ymax = ub_95,
      color = pathogen
    ),
    position = position_dodge(width = .5),
    width = 0
  ) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_y_continuous("Proportion of cases") +
  xlab("Dengue season") +
  theme_bw() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5,
                            6.5, 7.5, 8.5, 9.5, 10.5)) +
  coord_cartesian(ylim = c(0.045, 0.955)) +
  geom_hline(yintercept = 0) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

ggplot(raw_inc[raw_inc$label == "proportion", ]) +
  geom_point(
    aes(
      x = year,
      y = y,
      color = pathogen),
    position = position_dodge(width = .5)) +
  geom_errorbar(
    aes(
      x = year,
      y = y,
      ymin = lb_95,
      ymax = ub_95,
      color = pathogen
    ),
    position = position_dodge(width = .5),
    width = 0
  ) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_y_continuous("Proportion of cases") +
  xlab("Dengue season") +
  theme_bw() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5,
                            6.5, 7.5, 8.5, 9.5, 10.5)) +
  coord_cartesian(ylim = c(0.045, 0.955)) +
  geom_hline(yintercept = 0) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

cum_inc$method <- "Multi-pathogen Bayesian P-spline"
raw_inc$method <- "Naive method"
test <- rbind(cum_inc[cum_inc$label == "proportion", ],
              raw_inc[raw_inc$label == "proportion", ])

ggplot(test) +
  geom_point(aes(
    x = year,
    y = y,
    color = pathogen,
    shape = method
  ),
  position = position_dodge(width = .5)) +
  geom_errorbar(
    aes(
      x = year,
      y = y,
      ymin = lb_95,
      ymax = ub_95,
      color = pathogen,
      shape = method,
      linetype = method
    ),
    position = position_dodge(width = .5),
    width = 0
  ) +
  scale_color_brewer("Serotype", palette = "Dark2") +
  scale_shape_discrete("Method") +
  scale_linetype_discrete("Method") +
  scale_y_continuous("Proportion of cases") +
  xlab("Dengue season") +
  theme_bw() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5,
                            6.5, 7.5, 8.5, 9.5, 10.5)) +
  coord_cartesian(ylim = c(0.045, 0.955)) +
  geom_hline(yintercept = 0) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.background = element_rect(color = "black")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2, byrow = TRUE))

ggsave("Figures/DengueP.png", width = 12, height = 6)

###############################################################################
# Sup Figure 5
###############################################################################

df_new <- pivot_longer(dfR, cols = colnames(dfR)[3:6])

# Set colour scheme
cols <- brewer.pal(4, "Dark2")

# Plot subplots for daily raw data
raw1 <- ggplot(data = df_new[df_new$date >= as.Date("2006-04-01") &
                               df_new$date < as.Date("2016-04-01"), ]) +
  geom_col(aes(x = date, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(as.Date("2006-04-01") + 70,
                           as.Date("2016-04-01") - 70)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual("Serotype", values = cols) +
  scale_y_continuous("Serotypes detected per day") +
  theme(
    legend.position = "none",
    legend.background = element_rect(color = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

raw2 <- ggplot(data = df_new[df_new$date >= as.Date("2023-04-01") &
                               df_new$date < as.Date("2024-04-01"), ]) +
  geom_col(aes(x = date, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(as.Date("2023-04-01") + 10,
                           as.Date("2024-04-01") - 10)) +
  scale_x_date(date_breaks = "4 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual(values = cols) +
  scale_y_continuous("Serotypes detected per day") +
  theme(
    legend.position = "none",
    legend.background = element_rect(color = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
  )

# Aggregate data into weekly values for second pair of subplots
dfT <- rbind(df[, 1:8], dfS[, 1:8])
weekno <- as.numeric(dfT$date - dfT$date[1]) %/% 7
Week <- as.Date(dfT$date[1] + 7 * weekno)
abc1 <- aggregate(sero1 ~ Week , dfT, sum)
abc2 <- aggregate(sero2 ~ Week , dfT, sum)
abc3 <- aggregate(sero3 ~ Week , dfT, sum)
abc4 <- aggregate(sero4 ~ Week , dfT, sum)

abc1$name <- "1"
abc2$name <- "2"
abc3$name <- "3"
abc4$name <- "4"

colnames(abc1)[2] <- "value"
colnames(abc2)[2] <- "value"
colnames(abc3)[2] <- "value"
colnames(abc4)[2] <- "value"

df_week <- rbind(abc1, abc2, abc3, abc4)

# Plot second pair of subplots
raw3 <- ggplot(data = df_week[df_week$Week >= as.Date("2006-04-01") &
                                df_week$Week < as.Date("2016-04-01"), ]) +
  geom_col(aes(x = Week, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(as.Date("2006-04-01") + 70,
                           as.Date("2016-04-01") - 70)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual("Serotype", values = cols) +
  scale_y_continuous("Serotypes detected per week") +
  theme(
    legend.position = c(0.9, 0.7),
    legend.background = element_rect(color = "black")
  )

raw4 <- ggplot(data = df_week[df_week$Week >= as.Date("2023-04-01") &
                                df_week$Week < as.Date("2024-04-01"), ]) +
  geom_col(aes(x = Week, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(as.Date("2023-04-01") + 10,
                           as.Date("2024-04-01") - 10)) +
  scale_x_date(date_breaks = "4 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual(values = cols) +
  scale_y_continuous("Serotypes detected per week") +
  theme(
    legend.position = "none",
    legend.background = element_rect(color = "black"),
    axis.title.y = element_blank()
  )

# Plot overall figure and save
raw1 + raw2 + raw3 + raw4 + plot_layout(nrow = 2, widths = c(10, 1.5))

ggsave("Figures/Dengue_RAW.png",
       width = 12,
       height = 8)
