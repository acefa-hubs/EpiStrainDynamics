
# Loading required packages
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(patchwork)
library(lubridate)
library(tidyverse)

# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')

###############################################################################
## Reading in and cleaning data

## Loading the influenza data
df1 <- read.csv('Data/aus_influenza_data.csv')
df2 <- read.csv('Data/usa_influenza_data.csv')
df3 <- read.csv('Data/sin_influenza_data.csv')

## Selecting subset of data (2012-2023)
# Set limits on dates to consider
min_date <- as.Date("2012-01-01")
max_date <- as.Date("2023-12-31")

df1$week <- as.Date(df1$week)
df2$week <- as.Date(df2$week)
df3$week <- as.Date(df3$week)

df1 <- df1[df1$week < max_date & df1$week >= min_date, ]
df2 <- df2[df2$week < max_date & df2$week >= min_date, ]
df3 <- df3[df3$week < max_date & df3$week >= min_date, ]

df1 <- df1[order(df1$week), ]
df2 <- df2[order(df2$week), ]
df3 <- df3[order(df3$week), ]

df1$time <- seq(1, nrow(df1))
df2$time <- seq(1, nrow(df2))
df3$time <- seq(1, nrow(df3))

########################################################################################################################################################
## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Loading stan model used for influenza analysis
rw_influenza_mod <- stan_model('stan/rw_influenza_finalV2.stan')

#############################################################################################################################################
## Models fit to Australia (for main analysis)
rw_influenza_data1 <- list(
  num_data = nrow(df1),
  num_path = 4,
  Y = df1$ili,
  P1 = t(matrix(
    data = c(df1$inf_A, df1$inf_B, df1$num_spec - df1$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df1$inf_H3N2, df1$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df1$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 0
)

rw_influenza_fit1 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data1
)

saveRDS(rw_influenza_fit1, "FitStanModels/rw_influenza_fit1_N0.rds")

# Model fit to Australia with additional noise (for supplementary analysis)
rw_influenza_data1_N1 <- list(
  num_data = nrow(df1),
  num_path = 4,
  Y = df1$ili,
  P1 = t(matrix(
    data = c(df1$inf_A, df1$inf_B, df1$num_spec - df1$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df1$inf_H3N2, df1$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df1$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 1
)

rw_influenza_fit1_N1 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data1_N1
)

saveRDS(rw_influenza_fit1_N1,
        "FitStanModels/rw_influenza_fit1_N1.rds")

#############################################################################################################################################
## Model fit to United States of America (for main analysis)
rw_influenza_data2 <- list(
  num_data = nrow(df2),
  num_path = 4,
  Y = df2$ili,
  P1 = t(matrix(
    data = c(df2$inf_A, df2$inf_B, df2$num_spec - df2$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df2$inf_H3N2, df2$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df2$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 0
)

rw_influenza_fit2 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data2
)

saveRDS(rw_influenza_fit2, "FitStanModels/rw_influenza_fit2_N0.rds")

## Model fit to United States of America with additional noise (for supplementary analysis)
rw_influenza_data2_N1 <- list(
  num_data = nrow(df2),
  num_path = 4,
  Y = df2$ili,
  P1 = t(matrix(
    data = c(df2$inf_A, df2$inf_B, df2$num_spec - df2$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df2$inf_H3N2, df2$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df2$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 1
)

rw_influenza_fit2_N1 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data2_N1
)

saveRDS(rw_influenza_fit2_N1, "FitStanModels/rw_influenza_fit2_N1.rds")

#############################################################################################################################################
## Model fit to Singapore (for main analysis)
rw_influenza_data3 <- list(
  num_data = nrow(df3),
  num_path = 4,
  Y = df3$ili,
  P1 = t(matrix(
    data = c(df3$inf_A, df3$inf_B, df3$num_spec - df3$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df3$inf_H3N2, df3$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df3$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 0
)

rw_influenza_fit3 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data3
)

## Model fit to Singapore with noise term (for supplementary analysis)
saveRDS(rw_influenza_fit3, "FitStanModels/rw_influenza_fit3_N0.rds")

rw_influenza_data3_N1 <- list(
  num_data = nrow(df3),
  num_path = 4,
  Y = df3$ili,
  P1 = t(matrix(
    data = c(df3$inf_A, df3$inf_B, df3$num_spec - df3$inf_all),
    ncol = 3
  )),
  P2 = t(matrix(
    data = c(df3$inf_H3N2, df3$inf_H1N1), ncol = 2
  )),
  week_effect = 1,
  DOW = (df3$time %% 1) + 1,
  cov_structure = 1,
  noise_structure = 1
)

rw_influenza_fit3_N1 <- sampling(
  rw_influenza_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = rw_influenza_data3_N1
)

saveRDS(rw_influenza_fit3_N1,
        "FitStanModels/rw_influenza_fit3_N1.rds")

######################################################################################################################################################
## Main figures
#########################################################################################################################################################
# Read in the stan models that were saved
rw_influenza_fit1 <- readRDS("FitStanModels/rw_influenza_fit1_N0.rds")
rw_influenza_fit2 <- readRDS("FitStanModels/rw_influenza_fit2_N0.rds")
rw_influenza_fit3 <- readRDS("FitStanModels/rw_influenza_fit3_N0.rds")

# Calculate the incidence curves
rw_mod_inc1 <- rw_incidence(rw_influenza_fit1,
                            num_days = nrow(df1),
                            time_labels = df1$week)
rw_mod_inc2 <- rw_incidence(rw_influenza_fit2,
                            num_days = nrow(df2),
                            time_labels = df2$week)
rw_mod_inc3 <- rw_incidence(rw_influenza_fit3,
                            num_days = nrow(df3),
                            time_labels = df3$week)

# Re-level pathogen names
rw_mod_inc1$pathogen <- factor(
  rw_mod_inc1$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)
rw_mod_inc2$pathogen <- factor(
  rw_mod_inc2$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)
rw_mod_inc3$pathogen <- factor(
  rw_mod_inc3$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)

# Color scheme for figures
cols <- c("navy", "red4", "gold3","green4", "black")

###############################################################################
# Figure 1
##############################################################################
## Making all the subplots
plt1a <- ggplot(rw_mod_inc1) +
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
  geom_point(data = df1, aes(x = week, y = ili)) +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Cases of influenza-like illness") +
  theme_bw() +
  theme(legend.position = "none")

plt1b <- ggplot(rw_mod_inc2) +
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
  geom_point(data = df2, aes(x = week, y = ili)) +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Cases of influenza-like illness") +
  theme_bw() +
  theme(legend.position = "none")

plt1c <- ggplot(rw_mod_inc3) +
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
  geom_point(data = df3, aes(x = week, y = ili)) +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Cases of influenza-like illness") +
  theme_bw() +
  theme(legend.position = "none")

plt1a <- plt1a +
  labs(tag = "C") +
  theme(plot.tag.position = c(0.01, 0.97)) +
  geom_label(aes(
    label = "Australia",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)
plt1b <- plt1b +
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.01, 0.97),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  geom_label(aes(
    label = "United States of America",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)
plt1c <- plt1c +
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.01, 0.97),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.93, 0.7),
    legend.background = element_rect(color = "black"),
    legend.title = element_blank()
  ) +
  geom_label(aes(
    label = "Singapore",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)



plt1b / plt1c / plt1a

ggsave("Figures/Influenza1.png",
       width = 14,
       height = 10)

###############################################################################
# Supplementary Figure 6
##############################################################################


plt2a <- ggplot(rw_mod_inc1[rw_mod_inc1$pathogen %in%
                              c("Influenza A H3N2",
                                "Influenza A H1N1",
                                "Influenza B"), ]) +
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
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180, as.Date("2024-01-01") -
                             180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Modelled ILI+") +
  theme_bw() +
  theme(legend.position = "none")

plt2b <- ggplot(rw_mod_inc2[rw_mod_inc2$pathogen %in%
                              c("Influenza A H3N2",
                                "Influenza A H1N1",
                                "Influenza B"), ]) +
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
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Modelled ILI+") +
  theme_bw() +
  theme(legend.position = "none")

plt2c <- ggplot(rw_mod_inc3[rw_mod_inc3$pathogen %in%
                              c("Influenza A H3N2",
                                "Influenza A H1N1",
                                "Influenza B"), ]) +
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
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_manual("", values = cols) +
  scale_fill_manual("", values = cols) +
  xlab("Date") +
  ylab("Modelled ILI+") +
  theme_bw() +
  theme(legend.position = "none")


plt2a <- plt2a +
  labs(tag = "C") +
  theme(plot.tag.position = c(0.01, 0.97)) +
  geom_label(aes(
    label = "Australia",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)

plt2b <- plt2b +
  labs(tag = "A") +
  theme(
    plot.tag.position = c(0.01, 0.97),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

plt2c <- plt2c +
  labs(tag = "B") +
  theme(
    plot.tag.position = c(0.01, 0.97),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = c(0.93, 0.7),
    legend.background = element_rect(color = "black"),
    legend.title = element_blank()
  ) +
  geom_label(aes(
    label = "Singapore",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)



plt2a <- plt2a +
  annotate(
    "rect",
    xmin = as.Date("2012-05-01"),
    xmax = as.Date("2012-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2013-05-01"),
    xmax = as.Date("2013-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2014-05-01"),
    xmax = as.Date("2014-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2015-05-01"),
    xmax = as.Date("2015-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2016-05-01"),
    xmax = as.Date("2016-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2017-05-01"),
    xmax = as.Date("2017-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2018-05-01"),
    xmax = as.Date("2018-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2019-05-01"),
    xmax = as.Date("2019-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2020-05-01"),
    xmax = as.Date("2020-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2021-05-01"),
    xmax = as.Date("2021-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2022-05-01"),
    xmax = as.Date("2022-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2023-05-01"),
    xmax = as.Date("2023-10-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  geom_label(aes(
    label = "Australia",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)


plt2b <- plt2b +
  annotate(
    "rect",
    xmin = as.Date("2011-11-01"),
    xmax = as.Date("2012-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2012-11-01"),
    xmax = as.Date("2013-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2013-11-01"),
    xmax = as.Date("2014-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2014-11-01"),
    xmax = as.Date("2015-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2015-11-01"),
    xmax = as.Date("2016-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2016-11-01"),
    xmax = as.Date("2017-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2017-11-01"),
    xmax = as.Date("2018-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2018-11-01"),
    xmax = as.Date("2019-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2019-11-01"),
    xmax = as.Date("2020-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2020-11-01"),
    xmax = as.Date("2021-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2021-11-01"),
    xmax = as.Date("2022-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2022-11-01"),
    xmax = as.Date("2023-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  annotate(
    "rect",
    xmin = as.Date("2023-11-01"),
    xmax = as.Date("2024-04-30"),
    ymin = -Inf,
    ymax = Inf,
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  ) +
  geom_label(aes(
    label = "United States of America",
    x = as.Date("2018-01-01"),
    y = Inf
  ), vjust = 1.2)


plt2b / plt2c / plt2a


ggsave("Figures/Influenza2.png",
       width = 14,
       height = 10)

###############################################################################
# Supplementary Figure 7
##############################################################################

# Calculate the time-varying growth rates for all countries
rw_mod_gr1 <- rw_growth_rate(rw_influenza_fit1,
                             num_days = nrow(df1),
                             time_labels = df1$week)
rw_mod_gr2 <- rw_growth_rate(rw_influenza_fit2,
                             num_days = nrow(df2),
                             time_labels = df1$week)
rw_mod_gr3 <- rw_growth_rate(rw_influenza_fit3,
                             num_days = nrow(df3),
                             time_labels = df1$week)

# Label countries and merge dataframes
rw_mod_gr1$country <- "Australia"
rw_mod_gr2$country <- "United States of America"
rw_mod_gr3$country <- "Singapore"

rw_mod_gr <- rbind(rw_mod_gr1, rw_mod_gr2, rw_mod_gr3)


# Plot
ggplot(rw_mod_gr[rw_mod_gr$pathogen %in%
                   c("Influenza A H3N2",
                     "Influenza A H1N1",
                     "Influenza B"), ]) +
  facet_wrap(. ~ pathogen, nrow = 3) +
  geom_line(aes(x = time, y = y, color = country)) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = country
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = country
  ), alpha = 0.2) +
  scale_fill_brewer("Country", palette = "Dark2") +
  scale_color_brewer("Country", palette = "Dark2") +
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Growth rate") +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  theme(legend.position = "bottom",
        strip.background = element_rect(color = "black", fill = NA))

ggsave("Figures/Influenza3.png",
       width = 14,
       height = 10)

###############################################################################
# Supplementary Figures 16 & 17
##############################################################################
# Read in the models with additional noise
rw_influenza_fit1_N1 <- readRDS("FitStanModels/rw_influenza_fit1_N1.rds")
rw_influenza_fit2_N1 <- readRDS("FitStanModels/rw_influenza_fit2_N1.rds")
rw_influenza_fit3_N1 <- readRDS("FitStanModels/rw_influenza_fit3_N1.rds")

# Get incidence curves
rw_mod_inc1_N1 <- rw_incidence(rw_influenza_fit1_N1,
                               num_days = nrow(df1),
                               time_labels = df1$week)
rw_mod_inc2_N1 <- rw_incidence(rw_influenza_fit2_N1,
                               num_days = nrow(df2),
                               time_labels = df2$week)
rw_mod_inc3_N1 <- rw_incidence(rw_influenza_fit3_N1,
                               num_days = nrow(df3),
                               time_labels = df3$week)

# Reorder pathogen labels
rw_mod_inc1_N1$pathogen <- factor(
  rw_mod_inc1_N1$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)
rw_mod_inc2_N1$pathogen <- factor(
  rw_mod_inc2_N1$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)
rw_mod_inc3_N1$pathogen <- factor(
  rw_mod_inc3_N1$pathogen,
  levels = c(
    "Influenza A H1N1",
    "Influenza A H3N2",
    "Influenza B",
    "Unknown",
    "Total"
  )
)

# Label models (extra noise and main analyses) and merge
rw_mod_inc1_N1$label <- "Including pathogen-specific noise"
rw_mod_inc2_N1$label <- "Including pathogen-specific noise"
rw_mod_inc3_N1$label <- "Including pathogen-specific noise"

rw_mod_inc1$label <- "Observation noise only"
rw_mod_inc2$label <- "Observation noise only"
rw_mod_inc3$label <- "Observation noise only"

rw_mod1 <- rbind(rw_mod_inc1, rw_mod_inc1_N1)
rw_mod2 <- rbind(rw_mod_inc2, rw_mod_inc2_N1)
rw_mod3 <- rbind(rw_mod_inc3, rw_mod_inc3_N1)

# Create first plot and save
ggplot(rw_mod1) +
  geom_line(aes(x = time, y = y, color = label), size = 0.2) +
  facet_wrap(. ~ pathogen, nrow = 5, scales = "free_y") +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = label
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = label
  ), alpha = 0.2) +
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 150)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_brewer("Model", palette = "Dark2") +
  scale_fill_brewer("Model", palette = "Dark2") +
  xlab("Date") +
  ylab("Modelled cases") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_rect(color = "black", fill = NA))

ggsave("Figures/InfluenzaNoise1.png",
       width = 10,
       height = 12)

# Create second plot and save
ggplot(rw_mod2) +
  geom_line(aes(x = time, y = y, color = label), size = 0.2) +
  facet_wrap(. ~ pathogen, nrow = 5, scales = "free_y") +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = label
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = label
  ), alpha = 0.2) +
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180,
                           as.Date("2024-01-01") - 150)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_y_continuous(
    labels = function(x)
      format(x, scientific = FALSE)
  ) +
  scale_color_brewer("Model", palette = "Dark2") +
  scale_fill_brewer("Model", palette = "Dark2") +
  xlab("Date") +
  ylab("Modelled cases") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_rect(color = "black", fill = NA))

ggsave("Figures/InfluenzaNoise2.png",
       width = 10,
       height = 12)

# Can also plot for Singapore, but model does not converge for the two noise parameters (see sup materials and below)
# so this was not included in the paper
ggplot(rw_mod3) +
  geom_line(aes(x = time, y = y, color = label), size = 0.2) +
  facet_wrap(. ~ pathogen, nrow = 5, scales = "free_y") +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = label
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = label
  ), alpha = 0.2) +
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180, as.Date("2024-01-01") -
                             150)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_brewer("Model", palette = "Dark2") +
  scale_fill_brewer("Model", palette = "Dark2") +
  xlab("Date") +
  ylab("Modelled cases") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(color = "black", fill = NA))

ggsave("Figures/InfluenzaNoise3.png",
       width = 6,
       height = 8)

# Compare the noise terms for models fit to each data

post1 <- rstan::extract(rw_influenza_fit1)
post2 <- rstan::extract(rw_influenza_fit2)
post3 <- rstan::extract(rw_influenza_fit3)

post1_N1 <- rstan::extract(rw_influenza_fit1_N1)
post2_N1 <- rstan::extract(rw_influenza_fit2_N1)
post3_N1 <- rstan::extract(rw_influenza_fit3_N1)

g_dis1 <- 1 / post1_N1$eta[1]
nb_dis1 <- 1 + 1 / post1_N1$phi
hist(nb_dis1 / g_dis1)
quantile(nb_dis1 / g_dis1, c(0.5, 0.025, 0.975))

g_dis2 <- 1 / post2_N1$eta[1]
nb_dis2 <- 1 + 1 / post2_N1$phi
hist(nb_dis2 / g_dis2)
quantile(nb_dis2 / g_dis2, c(0.5, 0.025, 0.975))

g_dis3 <- 1 / post3_N1$eta[1]
nb_dis3 <- 1 + 1 / post3_N1$phi
hist(nb_dis3 / g_dis3)
quantile(nb_dis3 / g_dis3, c(0.5, 0.025, 0.975))

###############################################################################
# Supplementary Figure 8
##############################################################################
# Get proportion of influenza-like illness positve for influenza
rw_mod_prop1 <- rw_proportion(rw_influenza_fit1,
                              num_days = nrow(df1),
                              time_labels = df1$week)
rw_mod_prop2 <- rw_proportion(rw_influenza_fit2,
                              num_days = nrow(df2),
                              time_labels = df2$week)
rw_mod_prop3 <- rw_proportion(rw_influenza_fit3,
                              num_days = nrow(df3),
                              time_labels = df3$week)

# Label countries and merge
rw_mod_prop1$country <- "Australia"
rw_mod_prop2$country <- "United States of America"
rw_mod_prop3$country <- "Singapore"

rw_mod_prop <- rbind(rw_mod_prop1, rw_mod_prop2, rw_mod_prop3)


cols <- c("navy", "red4", "gold3", "green4", "black")

# Plot figure and save
ggplot(rw_mod_prop[rw_mod_prop$pathogen == "Influenza", ]) +
  geom_line(aes(x = time, y = y, color = country)) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = country
  ), alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = country
  ), alpha = 0.2) +
  coord_cartesian(xlim = c(as.Date("2012-01-01") + 180, as.Date("2024-01-01") -
                             180)) +
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
  scale_color_brewer("", palette = "Dark2") +
  scale_fill_brewer("", palette = "Dark2") +
  xlab("Date") +
  ylab("Influenza positive proportion") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("Figures/InfluenzaPP1.png",
       width = 10,
       height = 8)
