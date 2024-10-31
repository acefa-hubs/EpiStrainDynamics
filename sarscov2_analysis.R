

# Loading required packages
library(ggplot2)
library(rstan)
library(RColorBrewer)
library(cowplot)
library(patchwork)
library(lubridate)
library(tidyverse)
library(stats)

# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')

###############################################################################
## Reading and organising the data

# Set date limits to consider
min_date <- as.Date("2020-06-01") # This will be updated later in script
max_date <- as.Date("2022-12-31")

# Load the UK COVID data from dashboard
df1a <- read.csv('Data/newCasesBySpecimenDate_nation_2020.csv')
df1b <- read.csv('Data/newCasesBySpecimenDate_nation_2021.csv')
df1c <- read.csv('Data/newCasesBySpecimenDate_nation_2022.csv')
df1 <- rbind(df1a, df1b, df1c)
df1 <- df1[df1$date >= min_date & df1$date <= max_date, ]

#Combine the data for each nation into one time series
df <- data.frame()
for (i in 1:length(unique(df1$date))) {
  tmp_df <- df1[df1$date %in% unique(df1$date)[i], ]

  row_df <- data.frame(date = as.Date(unique(df1$date)[i]),
                       cases = sum(tmp_df$value))
  df <- rbind(df, row_df)
}

# Load the variant data available from Lythgoe et al and format
df1S <- read.csv('Data/SampleInfo.csv')
df1S <- df1S[df1S$major_lineage != "Unassingned", ]

lineages_considered <- c(
  "B.1.177",
  "B.1.1.7",
  "B.1.617.2",
  "BA.1",
  "BA.2",
  "BA.2.75",
  "BA.4",
  "BA.5",
  "BQ.1",
  "XBB"
)
df1S$major_lineage <- factor(df1S$major_lineage, levels = lineages_considered)

df1S$collection_date <- as.Date(df1S$collection_date)

uniq_dates <- seq.Date(min(df1S$collection_date),
                       max(df1S$collection_date),
                       by = 1)
uniq_dates <- seq.Date(min_date, max_date, by = 1)

mat <- matrix(
  data = NA,
  nrow = length(uniq_dates),
  ncol = length(lineages_considered) + 1
)
rownames(mat) <- uniq_dates
colnames(mat) <- c(lineages_considered, "Other")

for (i in 1:length(uniq_dates)) {
  print(i)
  df_tmp <- df1S[df1S$collection_date %in% uniq_dates[i], ]
  tab <- table(df_tmp$major_lineage, useNA = "always")
  mat[i, ] <- tab
}

# Combine all data into one
df <- cbind(df, mat)
df$total <- rowSums(mat)
df$t <- as.numeric(df$date)
df$t <- df$t - min(df$t) + 1

# Include data after threshold for number of variants determined per day is reached
df <- df[df$date >= as.Date("2020-09-23"), ]

##############################################################################################################################################
# Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Loading stan models used for SARS-CoV-2 analysis
ps_single_mod <- stan_model('stan/ps_single_final.stan')
ps_mp_mod <- stan_model('stan/ps_mp_finalV2.stan')

######################################################################################################################################################
## Fitting to all data
# Calculate the locations of equally spaced knots
knots <- get_knots(df$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data <- list(
  num_data = nrow(df),
  num_path = 11,
  num_knots = length(knots),
  knots = knots,
  spline_degree = 3,
  Y = df$cases,
  X = df$t,
  P = t(df[3:13]),
  # Might have to fix this line
  week_effect = 7,
  DOW = (df$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data
)
saveRDS(ps_mp_fit, "FitStanModels/ps_mp_fit.rds")

######################################################################################################################################################
## Subsetting data for analysis of Delta emergence fitting up to multiple time-points

df3 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-05-04"), ]
df4 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-05-11"), ]
df5 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-05-18"), ]
df6 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-05-25"), ]
df7 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-06-01"), ]
df8 <- df[df$date >= as.Date("2020-09-23") &
            df$date <= as.Date("2021-06-08"), ]

######################################################################################################################################################
## Models fit up to first time point

knots3 <- get_knots(df3$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data3 <- list(
  num_data = nrow(df3),
  num_path = 4,
  num_knots = length(knots3),
  knots = knots3,
  spline_degree = 3,
  Y = df3$cases,
  X = df3$t,
  P = t(matrix(
    data = c(df3$B.1.177, df3$B.1.1.7, df3$Other, df3$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df3$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit3 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data3
)
saveRDS(ps_mp_fit3, "FitStanModels/ps_mp_fit3.rds")


ps_single_data3 <- list(
  num_data = nrow(df3),
  num_knots = length(knots3),
  knots = knots3,
  spline_degree = 3,
  Y = df3$cases,
  X = df3$t,
  week_effect = 7,
  DOW = (df3$t %% 7) + 1
)

ps_single_fit3 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data3
)
saveRDS(ps_single_fit3, "FitStanModels/ps_single_fit3.rds")

######################################################################################################################################################
## Models fit up to second time point

knots4 <- get_knots(df4$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data4 <- list(
  num_data = nrow(df4),
  num_path = 4,
  num_knots = length(knots4),
  knots = knots4,
  spline_degree = 3,
  Y = df4$cases,
  X = df4$t,
  P = t(matrix(
    data = c(df4$B.1.177, df4$B.1.1.7, df4$Other, df4$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df4$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit4 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data4
)
saveRDS(ps_mp_fit4, "FitStanModels/ps_mp_fit4.rds") # Didn't converge super well!

ps_single_data4 <- list(
  num_data = nrow(df4),
  num_knots = length(knots4),
  knots = knots4,
  spline_degree = 3,
  Y = df4$cases,
  X = df4$t,
  week_effect = 7,
  DOW = (df4$t %% 7) + 1
)

ps_single_fit4 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data4
)
saveRDS(ps_single_fit4, "FitStanModels/ps_single_fit4.rds")

######################################################################################################################################################
## Models fit up to third time point

knots5 <- get_knots(df5$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data5 <- list(
  num_data = nrow(df5),
  num_path = 4,
  num_knots = length(knots5),
  knots = knots5,
  spline_degree = 3,
  Y = df5$cases,
  X = df5$t,
  P = t(matrix(
    data = c(df5$B.1.177, df5$B.1.1.7, df5$Other, df5$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df5$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit5 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data5
)
saveRDS(ps_mp_fit5, "FitStanModels/ps_mp_fit5.rds")

ps_single_data5 <- list(
  num_data = nrow(df5),
  num_knots = length(knots5),
  knots = knots5,
  spline_degree = 3,
  Y = df5$cases,
  X = df5$t,
  week_effect = 7,
  DOW = (df5$t %% 7) + 1
)

ps_single_fit5 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data5
)
saveRDS(ps_single_fit5, "FitStanModels/ps_single_fit5.rds")

######################################################################################################################################################
## Models fit up to fourth timepoint
knots6 <- get_knots(df6$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data6 <- list(
  num_data = nrow(df6),
  num_path = 4,
  num_knots = length(knots6),
  knots = knots6,
  spline_degree = 3,
  Y = df6$cases,
  X = df6$t,
  P = t(matrix(
    data = c(df6$B.1.177, df6$B.1.1.7, df6$Other, df6$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df6$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit6 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data6
)
saveRDS(ps_mp_fit6, "FitStanModels/ps_mp_fit6.rds")

ps_single_data6 <- list(
  num_data = nrow(df6),
  num_knots = length(knots6),
  knots = knots6,
  spline_degree = 3,
  Y = df6$cases,
  X = df6$t,
  week_effect = 7,
  DOW = (df6$t %% 7) + 1
)

ps_single_fit6 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data6
)
saveRDS(ps_single_fit6, "FitStanModels/ps_single_fit6.rds")

###############################################################################
## Models fit up to fifth time point

knots7 <- get_knots(df7$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data7 <- list(
  num_data = nrow(df7),
  num_path = 4,
  num_knots = length(knots7),
  knots = knots7,
  spline_degree = 3,
  Y = df7$cases,
  X = df7$t,
  P = t(matrix(
    data = c(df7$B.1.177, df7$B.1.1.7, df7$Other, df7$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df7$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit7 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data7
)
saveRDS(ps_mp_fit7, "FitStanModels/ps_mp_fit7.rds")

ps_single_data7 <- list(
  num_data = nrow(df7),
  num_knots = length(knots7),
  knots = knots7,
  spline_degree = 3,
  Y = df7$cases,
  X = df7$t,
  week_effect = 7,
  DOW = (df7$t %% 7) + 1
)

ps_single_fit7 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data7
)
saveRDS(ps_single_fit7, "FitStanModels/ps_single_fit7.rds")

######################################################################################################################################################
## Models fit up to sixth time point

knots8 <- get_knots(df8$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data8 <- list(
  num_data = nrow(df8),
  num_path = 4,
  num_knots = length(knots8),
  knots = knots8,
  spline_degree = 3,
  Y = df8$cases,
  X = df8$t,
  P = t(matrix(
    data = c(df8$B.1.177, df8$B.1.1.7, df8$Other, df8$B.1.617.2),
    ncol = 4
  )),
  week_effect = 7,
  DOW = (df8$t %% 7) + 1,
  cov_structure = 0,
  noise_structure = 0
)

ps_mp_fit8 <- sampling(
  ps_mp_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_mp_data8
)
saveRDS(ps_mp_fit8, "FitStanModels/ps_mp_fit8.rds")

ps_single_data8 <- list(
  num_data = nrow(df8),
  num_knots = length(knots8),
  knots = knots8,
  spline_degree = 3,
  Y = df8$cases,
  X = df8$t,
  week_effect = 7,
  DOW = (df8$t %% 7) + 1
)

ps_single_fit8 <- sampling(
  ps_single_mod,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  data = ps_single_data8
)
saveRDS(ps_single_fit8, "FitStanModels/ps_single_fit8.rds")

###############################################################################
## Plotting Figures
###############################################################################

##############################################################
## Figure 2
##############################################################

# Read in model fit
ps_mp_fit <- readRDS("FitStanModels/ps_mp_fit.rds")

# Get incidences
mod_inc <- ps_incidence(
  ps_mp_fit,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 11,
  pathogen_names = c(lineages_considered, "Other")
)

# Get growth rates
mod_gr <- ps_growth_rate(
  ps_mp_fit,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 11,
  pathogen_names = c(lineages_considered, "Other")
)

# Set mask for values not within date ranges
mod_inc$mask <- 0
mod_gr$mask <- 0

pathogens <- c(lineages_considered, "Other")
for (i in 1:length(pathogens)) {
  index <- df[colnames(df) == pathogens[i]] > 0
  min_date <- min(df[index, ]$date)
  max_date <- max(df[index, ]$date)
  mod_gr[mod_gr$pathogen == pathogens[i] &
           (mod_gr$time >= min_date & mod_gr$time <= max_date), ]$mask <- 1
  mod_inc[mod_inc$pathogen == pathogens[i] &
            (mod_inc$time >= min_date & mod_inc$time <= max_date), ]$mask <- 1
  print(pathogens[i])
  #mod_gr<-mod_gr[ (!(mod_gr$pathogen==pathogens[i])) | (mod_gr$time>min_date & mod_gr$time<max_date) ,]
}
mod_gr[mod_gr$pathogen == "Other", ]$mask <- 0
mod_gr[mod_gr$pathogen == "Other" &
         (mod_gr$time <= as.Date("2021-01-19") |
            mod_gr$time >= as.Date("2022-03-05")), ]$mask <- 1

mod_gr$grp <- 0
mod_gr[mod_gr$pathogen == "Other" &
         mod_gr$time <= as.Date("2021-01-19"), ]$grp <- 1

# Set colours
cols <- c(brewer.pal(11, "Paired"), "grey")
cols <- brewer.pal(12, "Paired")
cols <- c(cols[1:10], cols[12])

# Plot sub figures
plt1 <- ggplot(mod_inc[mod_inc$pathogen == "Total", ]) +
  geom_line(aes(x = time, y = y)) +
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
  geom_point(data = df, aes(x = date, y = cases), size = 0.2) +
  #geom_line(data=df, aes(x=date, y=cases))+
  ylab("Cases") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30),
                  ylim = c(0, max(df$cases))) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000, 250000)) +
  theme_bw(base_size = 14)

plt2 <- ggplot(mod_inc[mod_inc$pathogen != "Total" &
                         mod_inc$mask == 1, ]) +
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
  coord_cartesian(xlim = c(min_date + 30, max_date - 30),
                  ylim = c(0, max(mod_inc$y))) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000, 250000)) +
  ylab("Modelled cases") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw(base_size = 14)

plt3 <- ggplot(mod_gr[mod_gr$pathogen != "Total" &
                        mod_gr$mask == 1, ]) +
  geom_line(aes(
    x = time,
    y = y,
    color = pathogen,
    group = interaction(pathogen, grp)
  )) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_50,
    ymax = ub_50,
    fill = pathogen,
    group = interaction(pathogen, grp)
  ),
  alpha = 0.2) +
  geom_ribbon(aes(
    x = time,
    y = y,
    ymin = lb_95,
    ymax = ub_95,
    fill = pathogen,
    group = interaction(pathogen, grp)
  ),
  alpha = 0.2) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  scale_color_manual("Variant", values = cols) +
  scale_fill_manual("Variant", values = cols) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Date") +
  scale_y_continuous("Growth rate") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

plt1i <- plt1 +
  coord_cartesian(xlim = c(as.Date("2022-05-01"), max_date),
                  ylim = c(0, 34000)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    plot.background = element_rect(color = "black", fill = alpha("grey", 0.2))
  ) +
  scale_y_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000, 30000))

plt2i <- plt2 +
  coord_cartesian(xlim = c(as.Date("2022-05-01"), max_date),
                  ylim = c(0, 8000)) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    plot.background = element_rect(color = "black", fill = alpha("grey", 0.2))
  ) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000,
                                6000, 7000, 8000, 10000, 15000))

plt1 <- plt1 +
  labs(tag = "A") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.tag.position = c(0.01, 0.96)
  ) +
  annotate(
    "rect",
    xmin = as.Date("2022-05-01"),
    xmax = max_date + 2,
    ymin = -5000,
    ymax = 35000,
    color = "black",
    fill = 'grey',
    alpha = 0.2,
    linetype = "dashed"
  )

plt2 <- plt2 +
  labs(tag = "B") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.tag.position = c(0.01, 0.96)
  ) +
  annotate(
    "rect",
    xmin = as.Date("2022-05-01"),
    xmax = max_date + 2,
    ymin = -5000,
    ymax = 9000,
    color = "black",
    fill = "grey",
    alpha = 0.2,
    linetype = "dashed"
  )

plt3 <- plt3 +
  labs(tag = "C") +
  theme(legend.position = "bottom",
        plot.tag.position = c(0.01, 0.96))

plt4 <- plt4 +
  labs(tag = "A") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.tag.position = c(0.01, 0.96)
  )

plt1 + plt1i + plt2 + plt2i + plt3 + plot_spacer() +
  plot_layout(ncol = 2,
              widths = c(3, 1),
              heights = c(2, 2, 1.5))

patch1 <- plt1 + plt2 + plt3 +
  plot_layout(nrow = 3, heights = c(2, 2, 1.5))
patch2 <- plt1i + plt2i +
  plot_layout(nrow = 2)

cowplot::plot_grid(patch1,
                   patch2,
                   rel_widths = c(3, 1),
                   rel_heights = c(1, 0.8))

ggsave("Figures/SARS1.png", width = 14, height = 10)

##############################################################
## Figure 3
##############################################################

# Read in model fits
ps_mp_fit3 <- readRDS("FitStanModels/ps_mp_fit3.rds")
ps_single_fit3 <- readRDS("FitStanModels/ps_single_fit3.rds")

ps_mp_fit4 <- readRDS("FitStanModels/ps_mp_fit4.rds")
ps_single_fit4 <- readRDS("FitStanModels/ps_single_fit4.rds")

ps_mp_fit5 <- readRDS("FitStanModels/ps_mp_fit5.rds")
ps_single_fit5 <- readRDS("FitStanModels/ps_single_fit5.rds")

ps_mp_fit6 <- readRDS("FitStanModels/ps_mp_fit6.rds")
ps_single_fit6 <- readRDS("FitStanModels/ps_single_fit6.rds")

ps_mp_fit7 <- readRDS("FitStanModels/ps_mp_fit7.rds")
ps_single_fit7 <- readRDS("FitStanModels/ps_single_fit7.rds")

ps_mp_fit8 <- readRDS("FitStanModels/ps_mp_fit8.rds")
ps_single_fit8 <- readRDS("FitStanModels/ps_single_fit8.rds")


# Get incidence curves
mod_inc3 <- ps_incidence(
  ps_mp_fit3,
  df3$t,
  num_days = nrow(df3),
  time_labels = df3$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
mod_inc4 <- ps_incidence(
  ps_mp_fit4,
  df4$t,
  num_days = nrow(df4),
  time_labels = df4$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
mod_inc5 <- ps_incidence(
  ps_mp_fit5,
  df5$t,
  num_days = nrow(df5),
  time_labels = df5$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
mod_inc6 <- ps_incidence(
  ps_mp_fit6,
  df6$t,
  num_days = nrow(df6),
  time_labels = df6$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
mod_inc7 <- ps_incidence(
  ps_mp_fit7,
  df7$t,
  num_days = nrow(df7),
  time_labels = df7$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
mod_inc8 <- ps_incidence(
  ps_mp_fit8,
  df8$t,
  num_days = nrow(df8),
  time_labels = df8$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)

mod_single_inc3 <- ps_single_incidence(ps_single_fit3,
                                       df3$t,
                                       num_days = nrow(df3),
                                       time_labels = df3$date)
mod_single_inc4 <- ps_single_incidence(ps_single_fit4,
                                       df4$t,
                                       num_days = nrow(df4),
                                       time_labels = df4$date)
mod_single_inc5 <- ps_single_incidence(ps_single_fit5,
                                       df5$t,
                                       num_days = nrow(df5),
                                       time_labels = df5$date)
mod_single_inc6 <- ps_single_incidence(ps_single_fit6,
                                       df6$t,
                                       num_days = nrow(df6),
                                       time_labels = df6$date)
mod_single_inc7 <- ps_single_incidence(ps_single_fit7,
                                       df7$t,
                                       num_days = nrow(df7),
                                       time_labels = df7$date)
mod_single_inc8 <- ps_single_incidence(ps_single_fit8,
                                       df8$t,
                                       num_days = nrow(df8),
                                       time_labels = df8$date)

# Get time-varying growth rates
gr_inc3 <- ps_growth_rate(
  ps_mp_fit3,
  df3$t,
  num_days = nrow(df3),
  time_labels = df3$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
gr_inc4 <- ps_growth_rate(
  ps_mp_fit4,
  df4$t,
  num_days = nrow(df4),
  time_labels = df4$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
gr_inc5 <- ps_growth_rate(
  ps_mp_fit5,
  df5$t,
  num_days = nrow(df5),
  time_labels = df5$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
gr_inc6 <- ps_growth_rate(
  ps_mp_fit6,
  df6$t,
  num_days = nrow(df6),
  time_labels = df6$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)
gr_inc7 <- ps_growth_rate(
  ps_mp_fit7,
  df7$t,
  num_days = nrow(df7),
  time_labels = df7$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)

gr_inc8 <- ps_growth_rate(
  ps_mp_fit8,
  df8$t,
  num_days = nrow(df8),
  time_labels = df8$date,
  num_path = 4,
  pathogen_names = c(
    "B.1.177",
    "B.1.1.7 (Alpha)",
    "Other (Wildtype)",
    "B.1.617.2 (Delta)"
  )
)

gr_single_inc3 <- ps_single_growth_rate(ps_single_fit3,
                                        df3$t,
                                        num_days = nrow(df3),
                                        time_labels = df3$date)
gr_single_inc4 <- ps_single_growth_rate(ps_single_fit4,
                                        df4$t,
                                        num_days = nrow(df4),
                                        time_labels = df4$date)
gr_single_inc5 <- ps_single_growth_rate(ps_single_fit5,
                                        df5$t,
                                        num_days = nrow(df5),
                                        time_labels = df5$date)
gr_single_inc6 <- ps_single_growth_rate(ps_single_fit6,
                                        df6$t,
                                        num_days = nrow(df6),
                                        time_labels = df6$date)
gr_single_inc7 <- ps_single_growth_rate(ps_single_fit7,
                                        df7$t,
                                        num_days = nrow(df7),
                                        time_labels = df7$date)
gr_single_inc8 <- ps_single_growth_rate(ps_single_fit8,
                                        df8$t,
                                        num_days = nrow(df8),
                                        time_labels = df8$date)

# Label and merge dataframes

mod_inc3$group <- 3
mod_inc4$group <- 4
mod_inc5$group <- 5
mod_inc6$group <- 6
mod_inc7$group <- 7
mod_inc8$group <- 8

mod_inc_comb1 <- rbind(mod_inc3, mod_inc4, mod_inc5, mod_inc6, mod_inc7, mod_inc8)

mod_single_inc3$group <- 3
mod_single_inc4$group <- 4
mod_single_inc5$group <- 5
mod_single_inc6$group <- 6
mod_single_inc7$group <- 7
mod_single_inc8$group <- 8

mod_single_inc_comb1 <- rbind(
  mod_single_inc3,
  mod_single_inc4,
  mod_single_inc5,
  mod_single_inc6,
  mod_single_inc7,
  mod_single_inc8
)

mod_single_inc_comb1$pathogen <- "Total"

gr_inc3$group <- 3
gr_inc4$group <- 4
gr_inc5$group <- 5
gr_inc6$group <- 6
gr_inc7$group <- 7
gr_inc8$group <- 8

gr_comb1 <- rbind(gr_inc3, gr_inc4, gr_inc5, gr_inc6, gr_inc7, gr_inc8)

gr_single_inc3$group <- 3
gr_single_inc4$group <- 4
gr_single_inc5$group <- 5
gr_single_inc6$group <- 6
gr_single_inc7$group <- 7
gr_single_inc8$group <- 8

gr_single_comb1 <- rbind(
  gr_single_inc3,
  gr_single_inc4,
  gr_single_inc5,
  gr_single_inc6,
  gr_single_inc7,
  gr_single_inc8
)

# Relabel and subset to include only Delta vs Alpha (i.e. no wildtype)
mod_inc_plt <- mod_inc_comb1[mod_inc_comb1$pathogen %in% c("B.1.1.7 (Alpha)", "B.1.617.2 (Delta)", "Total"), ]
mod_inc_plt <- mod_inc_plt[mod_inc_plt$group > 2, ]

mod_single_inc_plt <- mod_single_inc_comb1
mod_single_inc_plt$pathogen <- "Total\n(variant data not included)"
mod_single_inc_plt <- mod_single_inc_plt[mod_single_inc_plt$group > 2, ]

gr_plt <- gr_comb1[gr_comb1$pathogen %in%
                     c("B.1.1.7 (Alpha)", "B.1.617.2 (Delta)", "Total"), ]
gr_plt <- gr_plt[gr_plt$group > 2, ]

gr_single_comb1$pathogen <- "Total\n(variant data not included)"
gr_single_plt <- gr_single_comb1[gr_single_comb1$pathogen %in%
                                   c("B.1.1.7 (Alpha)",
                                     "B.1.617.2 (Delta)",
                                     "Total\n(variant data not included)"), ]
gr_single_plt <- gr_single_plt[gr_single_plt$group > 2, ]

# Set colours
cols <- c("navy", "red4", "gold3", "green4", "black")

# Make subplots
pltV1 <- ggplot(mod_inc_plt) +
  facet_wrap(. ~ group, nrow = 6) +
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
  geom_point(data = df, aes(x = date, y = cases), size = 0.4) +
  coord_cartesian(xlim = c(as.Date("2021-04-04"), as.Date("2021-06-28")),
                  ylim = c(0, 18000)) +
  #geom_line(data=df, aes(x=date, y=cases))+
  ylab("Cases") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = cols[c(1, 2, 4, 5)]) +
  scale_fill_manual(values = cols[c(1, 2, 4, 5)]) +
  xlab("Date") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank()
  )

pltV2 <- ggplot(mod_single_inc_plt) +
  facet_wrap(. ~ group, nrow = 6) +
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
  geom_point(data = df, aes(x = date, y = cases), size = 0.4) +
  coord_cartesian(xlim = c(as.Date("2021-04-04"), as.Date("2021-06-28")),
                  ylim = c(0, 18000)) +
  #geom_line(data=df, aes(x=date, y=cases))+
  ylab("Cases") +
  scale_color_manual(values = cols[c(5)]) +
  scale_fill_manual(values = cols[c(5)]) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  theme_bw(base_size = 14) +
  xlab("Date") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank()
  )

gr_plt <- rbind(gr_plt[gr_plt$pathogen == "B.1.617.2 (Delta)" &
                         gr_plt$time >= as.Date("2021-04-28"), ],
                gr_plt[gr_plt$pathogen %in% c("B.1.1.7 (Alpha)"), ])

pltV3 <- ggplot(gr_plt) +
  facet_wrap(. ~ group, nrow = 6) +
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
  geom_line(data = gr_single_plt, aes(x = time, y = y, color = pathogen)) +
  geom_ribbon(
    data = gr_single_plt,
    aes(
      x = time,
      y = y,
      ymin = lb_50,
      ymax = ub_50,
      fill = pathogen
    ),
    alpha = 0.2
  ) +
  geom_ribbon(
    data = gr_single_plt,
    aes(
      x = time,
      y = y,
      ymin = lb_95,
      ymax = ub_95,
      fill = pathogen
    ),
    alpha = 0.2
  ) +
  coord_cartesian(xlim = c(as.Date("2021-04-28"), as.Date("2021-06-08")),
                  ylim = c(-0.15, 0.5)) +
  #geom_line(data=df, aes(x=date, y=cases))+
  ylab("Growth rate") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  scale_color_manual(values = cols[c(1, 2, 5, 5)]) +
  scale_fill_manual(values = cols[c(1, 2, 5, 5)]) +
  theme_bw(base_size = 14) +
  xlab("Date") +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

pltV1 <- pltV1 + labs(tag = "A") +  theme(
  legend.position = c(0.25, 0.95),
  legend.background = element_rect(color = "black"),
  plot.tag.position = c(0.0, 0.99)
) +
  coord_cartesian(xlim = c(as.Date("2021-04-20"), as.Date("2021-06-18")),
                  ylim = c(0, 11000))

pltV2 <- pltV2 + labs(tag = "B") + theme(
  legend.position = c(0.35, 0.95),
  legend.background = element_rect(color =
                                     "black"),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  plot.tag.position = c(-0.01, 0.99)
) +
  coord_cartesian(xlim = c(as.Date("2021-04-20"), as.Date("2021-06-18")),
                  ylim = c(0, 11000))

pltV3 <- pltV3 + labs(tag = "C") +  theme(plot.tag.position = c(0.0, 0.99))

pltV1 + pltV2 + pltV3 + plot_layout(ncol = 3, widths = c(2, 2, 1))

ggsave("Figures/SARS_V2.png", width = 12, height = 12)

##########################################################################################
# Supplementary figure 4 (raw data)
##########################################################################################

df_new <- pivot_longer(df, cols = colnames(df)[3:13])

cols <- c(brewer.pal(11, "Paired"), "grey")
cols <- brewer.pal(12, "Paired")
cols <- c(cols[1:10], cols[12])

raw1 <- ggplot(data = df_new) +
  geom_col(aes(x = date, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual(values = cols) +
  scale_y_continuous("Variants sequenced") +
  guides(fill = guide_legend("Variant", nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend("Variant", nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom",
        legend.background = element_rect(color = "black")) +
  annotate(
    "rect",
    xmin = as.Date("2021-03-01"),
    xmax = as.Date("2021-06-28"),
    ymin = -10,
    ymax = 100,
    color = "black",
    fill = NA,
    alpha = 0.2,
    linetype = "dashed"
  ) +
  labs(tag = "A")

raw2 <- ggplot(data = df_new) +
  geom_col(aes(x = date, y = value, fill = name)) +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(as.Date("2021-03-01"),
                           as.Date("2021-06-28")),
                  ylim = c(0, 50)) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  xlab("Date") +
  scale_fill_manual(values = cols) +
  scale_y_continuous("Variants sequenced") +
  guides(fill = guide_legend("Variant", nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend("Variant", nrow = 2, byrow = TRUE)) +
  theme(legend.position = "none") +
  geom_vline(xintercept = as.Date("2021-05-04"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-05-11"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-05-18"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-05-25"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-06-01"), linetype = "dashed") +
  geom_vline(xintercept = as.Date("2021-06-08"), linetype = "dashed") +
  labs(tag = "B")

raw1 / raw2
ggsave("Figures/SARS_RAW.png", width = 10, height = 12)

###############################################################################
# Analysis for SFig 12 and SFig 13
###############################################################################
# Some functions required for this specific analysis

source('R/SFig12-13_functions.R')

###############################################################################
# SupFig 12
###############################################################################
# Get posterior smooth function
ps_a <- ps_get_a(
  ps_mp_fit,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 11
)

# Get growth rate advantage between specified pathogen pairs (path_i vs path_j)
mod_gr_adv <- ps_growth_rate_advantage(
  ps_a,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  num_path = 11,
  path_i = c(2, 3, 4, 5, 8),
  path_j = c(11, 2, 3, 4, 5),
  pathogen_names = c(lineages_considered, "Other")
)

## Set some generation time distributions (gamma distributions)
# parameters for each GI distribution (n common for all)
n <- 2.5

mu_w <- 4.95
mu_a <- 4.35
mu_d <- 3.9
mu_ba1 <- 3.3
mu_ba2 <- 2.9
mu_ba5 <- 2.3

gammaDist <- function(b, n, a)
  (b ** n) * (a ** (n - 1)) * exp(-b * a) / gamma(n)

gi_dist_i = c(function(x)
  gammaDist(b = n / mu_a, n = n, x), function(x)
    gammaDist(b = n / mu_d, n = n, x), function(x)
      gammaDist(b = n / mu_ba1, n = n, x), function(x)
        gammaDist(b = n / mu_ba2, n = n, x), function(x)
          gammaDist(b = n / mu_ba5, n = n, x))

gi_dist_j = c(function(x)
  gammaDist(b = n / mu_w, n = n, x), function(x)
    gammaDist(b = n / mu_a, n = n, x), function(x)
      gammaDist(b = n / mu_d, n = n, x), function(x)
        gammaDist(b = n / mu_ba1, n = n, x), function(x)
          gammaDist(b = n / mu_ba2, n = n, x))

# Get Rt advantagees assuming same generation interval distribution
mod_Rt_adv <- ps_Rt_advantage(
  ps_a,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  tau_max = 7,
  gi_dist = function(x)
    gammaDist(b = n / mu_w, n = n, x),
  num_path = 11,
  path_i = c(2, 3, 4, 5, 8),
  path_j = c(11, 2, 3, 4, 5),
  pathogen_names = c(lineages_considered, "Other")
)

# Get Rt advantagees assuming varying generation interval distribution
mod_Rt_advVar <- ps_Rt_advantage2(
  ps_a,
  df$t,
  num_days = nrow(df),
  time_labels = df$date,
  tau_max = 7,
  gi_dist_i = gi_dist_i,
  gi_dist_j = gi_dist_j,
  num_path = 11,
  path_i = c(2, 3, 4, 5, 8),
  path_j = c(11, 2, 3, 4, 5),
  pathogen_names = c(lineages_considered, "Other")
)
mod_gr_adv2 <- mod_gr_adv
mod_Rt_adv2 <- mod_Rt_adv
mod_Rt_advVar2 <- mod_Rt_advVar

# Mask so only periods in which both variants are circulating display
mod_gr_adv$mask <- 0
mod_Rt_adv$mask <- 0
mod_Rt_advVar$mask <- 0

for (i in 1:nrow(mod_gr_adv)) {
  print(i)
  msk1 <- mod_gr[mod_gr$time == mod_gr_adv$time[i] &
                   mod_gr$pathogen == mod_gr_adv$pathogen1[i], ]$mask

  msk2 <- mod_gr[mod_gr$time == mod_gr_adv$time[i] &
                   mod_gr$pathogen == mod_gr_adv$pathogen2[i], ]$mask


  if (msk1 == 1 & msk2 == 1) {
    mod_gr_adv$mask[i] <- 1
  }

  print(i)

  if (i <= nrow(mod_Rt_adv)) {
    msk1 <- mod_gr[mod_gr$time == mod_Rt_adv$time[i] &
                     mod_gr$pathogen == mod_Rt_adv$pathogen1[i], ]$mask

    msk2 <- mod_gr[mod_gr$time == mod_Rt_adv$time[i] &
                     mod_gr$pathogen == mod_Rt_adv$pathogen2[i], ]$mask

    if (msk1 == 1 & msk2 == 1) {
      mod_Rt_adv$mask[i] <- 1
      mod_Rt_advVar$mask[i] <- 1
    }
  }
}

#### Some formatting before plotting
mod_gr_adv$pathogen1 <- factor(mod_gr_adv$pathogen1)
mod_gr_adv$pathogen2 <- factor(mod_gr_adv$pathogen2)

levels(mod_gr_adv$pathogen1) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "BA.5 (Omicron)"
)
levels(mod_gr_adv$pathogen2) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "Wildtype"
)

mod_gr_adv$pathogen <- paste(mod_gr_adv$pathogen1,
                             " over ",
                             mod_gr_adv$pathogen2)

####
mod_Rt_adv$pathogen1 <- factor(mod_Rt_adv$pathogen1)
mod_Rt_adv$pathogen2 <- factor(mod_Rt_adv$pathogen2)

levels(mod_Rt_adv$pathogen1) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "BA.5 (Omicron)"
)
levels(mod_Rt_adv$pathogen2) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "Wildtype"
)

mod_Rt_adv$pathogen <- paste(mod_Rt_adv$pathogen1,
                             " over ",
                             mod_Rt_adv$pathogen2)

####
mod_Rt_advVar$pathogen1 <- factor(mod_Rt_advVar$pathogen1)
mod_Rt_advVar$pathogen2 <- factor(mod_Rt_advVar$pathogen2)

levels(mod_Rt_advVar$pathogen1) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "BA.5 (Omicron)"
)
levels(mod_Rt_advVar$pathogen2) <- c(
  "B.1.1.7 (Alpha)",
  "B.1.617.2 (Delta)",
  "BA.1 (Omicron)",
  "BA.2 (Omicron)",
  "Wildtype"
)

mod_Rt_advVar$pathogen <- paste(mod_Rt_advVar$pathogen1,
                                " over ",
                                mod_Rt_advVar$pathogen2)

####
# Plot subfigures

plt_gr1 <- ggplot(mod_gr_adv[mod_gr_adv$mask == 1, ]) +
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
  scale_color_brewer("Variants", palette = "Dark2") +
  scale_fill_brewer("Variants", palette = "Dark2") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Date") +
  scale_y_continuous("Growth rate advantage") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30))

plt_gr2 <- ggplot(mod_Rt_adv[mod_Rt_adv$mask == 1, ]) +
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
  scale_color_brewer("Variants", palette = "Dark2") +
  scale_fill_brewer("Variants", palette = "Dark2") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Date") +
  scale_y_continuous("Multiplicative R advantage") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30),
                  ylim = c(0, 4)) +
  geom_label(aes(
    label = "Constant generation interval distribution",
    x = as.Date("2021-03-01"),
    y = Inf
  ),
  vjust = 1.5)

plt_gr3 <- ggplot(mod_Rt_advVar[mod_Rt_advVar$mask == 1, ]) +
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
  scale_color_brewer("Variants", palette = "Dark2") +
  scale_fill_brewer("Variants", palette = "Dark2") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30)) +
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y") +
  ylab("Modelled cases") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Date") +
  scale_y_continuous("Multiplicative R advantage") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom") +
  coord_cartesian(xlim = c(min_date + 30, max_date - 30),
                  ylim = c(0, 3)) +
  geom_label(
    aes(
      label = "Generation interval distribution changes between variants",
      x = as.Date("2021-03-01"),
      y = Inf
    ),
    vjust = 1.5
  )

plt_gr1 <- plt_gr1 +
  labs(tag = "A") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

plt_gr2 <- plt_gr2 +
  labs(tag = "B") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none"
  )

plt_gr3 <- plt_gr3 +
  labs(tag = "C")

plt_gr1 + plt_gr2 + plt_gr3 + plot_layout(nrow = 3)

ggsave("Figures/SARS_GR1.png", width = 14, height = 10)

#################################################################################
# SupFig 13
##################################################################################

x <- seq(0, 14)

df_gi1 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_w, n = n, x),
                     Variant = "Wildtype")
df_gi2 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_a, n = n, x),
                     Variant = "B.1.1.7 (Alpha)")
df_gi3 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_d, n = n, x),
                     Variant = "Delta (B.1.617.2)")
df_gi4 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_ba1, n = n, x),
                     Variant = "BA.1 (Omicron)")
df_gi5 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_ba2, n = n, x),
                     Variant = "BA.2 (Omicron)")
df_gi6 <- data.frame(t = x,
                     p = gammaDist(b = n / mu_ba5, n = n, x),
                     Variant = "BA.5 (Omicron)")

df_gi <- rbind(df_gi1, df_gi2, df_gi3, df_gi4, df_gi5, df_gi6)

df_int <- df_gi[df_gi$t == 0, ]

Vars = unique(df_gi$Variant)
ggplot(df_gi) +
  geom_line(aes(x = t, y = p, color = Variant)) +
  geom_vline(data = df_int,
             aes(
               xintercept = c(mu_w, mu_a, mu_d, mu_ba1, mu_ba2, mu_ba5),
               color = Variant
             ),
             linetype = "dashed") +
  theme_bw(base_size = 14) +
  scale_color_brewer("Variants", palette = "Dark2") +
  ylab("Probability density") +
  xlab("Time since infection") +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_rect(color = 'black')
  ) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(0.0168, 0.35))

ggsave("Figures/SARS_GI.png", width = 7, height = 7)
