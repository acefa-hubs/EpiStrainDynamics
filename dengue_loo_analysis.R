setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/EpidemicAnalytics/GitHubRelease")

# Loading required packages
library(ggplot2)
library(rstan)
library(patchwork)
library(cowplot)


# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')

##############################################################################################################################################

# Read in Taiwan dengue data and format it into new data frame df
df1 <- read.csv('Data/Dengue_Daily.csv')

df1 <- df1[c(1,3,22)]
colnames(df1) <- c("Onset_date", "Notification_date", "Serotype")
df1$Onset_date <- as.Date(df1$Onset_date)
df1$Notification_date <- as.Date(df1$Notification_date)
table(df1$Serotype)
df1$Serotype <- factor(df1$Serotype)
levels(df1$Serotype) <- c("No data", "Serotype 1", "Serotype 3", "Serotype 2", "Serotype 4")

df <- data.frame()
onset_dates <- seq.Date(min(df1$Onset_date), max(df1$Onset_date), by=1)

for(i in 1:length(onset_dates) ){
  df_tmp <- df1[df1$Onset_date == onset_dates[i],]
  
  table(df_tmp$Serotype)
  
  row_df <- data.frame(date = onset_dates[i],
                       cases = nrow(df_tmp),
                       sero1 = nrow(df_tmp[df_tmp$Serotype=="Serotype 1",]),
                       sero2 = nrow(df_tmp[df_tmp$Serotype=="Serotype 2",]),
                       sero3 = nrow(df_tmp[df_tmp$Serotype=="Serotype 3",]),
                       sero4 = nrow(df_tmp[df_tmp$Serotype=="Serotype 4",]),
                       unsero = nrow(df_tmp[df_tmp$Serotype=="No data",]))
  
  df <- rbind(df, row_df)
  
}

##############################################################################################################################################
# Subset data into specific periods to fit model to
################################################################
# 2014 season only
min_date <- as.Date("2014-04-01")
max_date <- as.Date("2015-04-01")

dfS1 <- df[df$date<max_date & df$date>=min_date,]
dfS1$time <- seq(1, nrow(dfS1))

################################################################
# 2015 season only
min_date <- as.Date("2015-04-01")
max_date <- as.Date("2016-04-01")

dfS2 <- df[df$date<max_date & df$date>=min_date,]
dfS2$time <- seq(1, nrow(dfS2))

################################################################
# 2006-2015 seasons (model not fit here see dengue_analysis.R)
# Set limits on dates to consider
min_date <- as.Date("2006-04-01")
max_date <- as.Date("2016-04-01")

df2 <- df[df$date<max_date & df$date>=min_date,]
df2$time <- seq(1, nrow(df2))

################################################################
# 2006 - 2013 seasons only
min_date <- as.Date("2006-04-01")
max_date <- as.Date("2014-04-01")

df <- df[df$date<max_date & df$date>=min_date,]
df$time <- seq(1, nrow(df))


##############################################################################################################################################
# Set stan options
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Loading stan model used for LOO dengue analysis
ps_mp_modV2 <- stan_model('stan/ps_mp_finalV2.stan')
ps_mp_modS <- stan_model('stan/ps_mp_final - CopyWithPrior.stan')

######################################################################################################################################################
# Fitting model to 2006-2013 seasons
knots <- get_knots(df$time, days_per_knot = 5, spline_degree = 3)

ps_mp_data <- list(num_data = nrow(df),
                   num_path = 4,
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = df$cases,
                   X = df$time,
                   P = t(matrix(data= c(df$sero1, df$sero2, df$sero3, df$sero4), ncol=4)),
                   week_effect = 1,
                   DOW = (df$time %% 1)+1,
                   cov_structure = 1,
                   noise_structure = 0) 

ps_mp_fit <- sampling(ps_mp_modV2,
                      iter= 5000,
                      warmup = 1000,
                      chains=4,
                      data = ps_mp_data)
saveRDS(ps_mp_fit, "FitStanModels/dengue_fit_loo_prior.rds")


######################################################################################################################################################
# Extract posterior distributions from fit to 2006-2013 seasons
post<-rstan::extract(ps_mp_fit)

# Get mean and standard deviations for tau and phi
tau_mn<- apply(post$tau, 2, mean)
tau_sd<- apply(post$tau, 2, sd)
phi_mn <- c(mean(post$phi))
phi_sd <- c(sd(post$phi))


##############################################################################################################################
# Fitting model to 2014 season using informed priors

knotsS1 <- get_knots(dfS1$time, days_per_knot = 5, spline_degree = 3)

ps_mp_dataS1 <- list(num_data = nrow(dfS1),
                     num_path = 4,
                     num_knots = length(knotsS1),
                     knots = knotsS1,
                     spline_degree=3,
                     Y = dfS1$cases,
                     X = dfS1$time,
                     P = t(matrix(data= c(dfS1$sero1, dfS1$sero2, dfS1$sero3, dfS1$sero4), ncol=4)),
                     week_effect = 1,
                     DOW = (dfS1$time %% 1)+1,
                     cov_structure = 1,
                     tau_mn = tau_mn,
                     tau_sd = tau_sd,
                     phi_mn = phi_mn,
                     phi_sd = phi_sd) 

ps_mp_fitS1 <- sampling(ps_mp_modS,
                        iter= 5000,
                        warmup = 1000,
                        chains=4,
                        data = ps_mp_dataS1)
saveRDS(ps_mp_fitS1, "FitStanModels/dengue_fit_loo_post1.rds")


##############################################################################################################################
# Fitting model to 2015 season using informed priors

knotsS2 <- get_knots(dfS2$time, days_per_knot = 5, spline_degree = 3)

ps_mp_dataS2 <- list(num_data = nrow(dfS2),
                     num_path = 4,
                     num_knots = length(knotsS2),
                     knots = knotsS2,
                     spline_degree=3,
                     Y = dfS2$cases,
                     X = dfS2$time,
                     P = t(matrix(data= c(dfS2$sero1, dfS2$sero2, dfS2$sero3, dfS2$sero4), ncol=4)),
                     week_effect = 1,
                     DOW = (dfS2$time %% 1)+1,
                     cov_structure = 1,
                     tau_mn = tau_mn,
                     tau_sd = tau_sd,
                     phi_mn = phi_mn,
                     phi_sd = phi_sd) 

ps_mp_fitS2 <- sampling(ps_mp_modS,
                        iter= 5000,
                        warmup = 1000,
                        chains=4,
                        data = ps_mp_dataS2)
saveRDS(ps_mp_fitS2, "FitStanModels/dengue_fit_loo_post2.rds")

######################################################################################################################################################


######################################################################################################################################################
## Figures
######################################################################################################################################################
# Sup Fig 17
###############################################################
# Read in model fits
ps_mp_fitA <- readRDS("FitStanModels/dengue_fit_loo_prior.rds")
ps_mp_fitB <- readRDS("FitStanModels/dengue_fit4.rds")
ps_mp_fitS1 <- readRDS("FitStanModels/dengue_fit_loo_post1.rds")
ps_mp_fitS2 <- readRDS("FitStanModels/dengue_fit_loo_post2.rds")

# Get incidence curves
mod_incA <- ps_incidence(ps_mp_fitA, df$time, num_days=nrow(df), time_labels = df$date,
                         num_path = 4,
                         days_per_knot = 5,
                         pathogen_names = c("1", "2", "3", "4"))

mod_incB <- ps_incidence(ps_mp_fitB, df2$time, num_days=nrow(df2), time_labels = df2$date,
                         num_path = 4,
                         days_per_knot = 5,
                         pathogen_names = c("1", "2", "3", "4"))

mod_inc1 <- ps_incidence(ps_mp_fitS1, dfS1$time, num_days=nrow(dfS1), time_labels = dfS1$date,
                         num_path = 4,
                         days_per_knot = 5,
                         pathogen_names = c("1", "2", "3", "4"))

mod_inc2 <- ps_incidence(ps_mp_fitS2, dfS2$time, num_days=nrow(dfS2), time_labels = dfS2$date,
                         num_path = 4,
                         days_per_knot = 5,
                         pathogen_names = c("1", "2", "3", "4"))


# Label and merge
mod_incA$label <- "2006-2014"
mod_incB$label <- "2006-2016 seasons"
mod_inc1$label <- "2014 season"
mod_inc2$label <- "2015 season"

mod_inc <- rbind(mod_incB, mod_inc1, mod_inc2)


# Plot all the subfigures

cols <- c("black", "gold","green4")

loo0 <- ggplot(mod_inc[mod_inc$pathogen=="Total",])+
  geom_line(aes(x=time, y=y, colour=label))+
  geom_point(data=df2, aes(x=date, y=cases), size=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=label), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = as.Date("2014-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2015-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2016-04-01"), linetype="dashed")+
  scale_color_manual("Model fit to:",values=cols)+
  scale_fill_manual("Model fit to:",values=cols)+
  ylab("Total cases")+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  coord_cartesian(xlim=c(as.Date("2014-04-01"), as.Date("2016-04-01") )  )+
  theme(legend.position=c(0.5,0.5),
        legend.background=element_rect(color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

loo1 <- ggplot(mod_inc[mod_inc$pathogen=="1",])+
  geom_line(aes(x=time, y=y, colour=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=label), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = as.Date("2014-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2015-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2016-04-01"), linetype="dashed")+
  scale_color_manual("Model fit to:",values=cols)+
  scale_fill_manual("Model fit to:",values=cols)+
  ylab("Serotype 1 cases")+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  coord_cartesian(xlim=c(as.Date("2014-04-01"), as.Date("2016-04-01") )  )+
  theme(legend.position="none",
        legend.background=element_rect(color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

loo2 <- ggplot(mod_inc[mod_inc$pathogen=="2",])+
  geom_line(aes(x=time, y=y, colour=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=label), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = as.Date("2014-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2015-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2016-04-01"), linetype="dashed")+
  scale_color_manual("Model fit to:",values=cols)+
  scale_fill_manual("Model fit to:",values=cols)+
  ylab("Serotype 2 cases")+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  coord_cartesian(xlim=c(as.Date("2014-04-01"), as.Date("2016-04-01") )  )+
  theme(legend.position="none",
        legend.background=element_rect(color="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

loo3 <- ggplot(mod_inc[mod_inc$pathogen=="3",])+
  geom_line(aes(x=time, y=y, colour=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=label), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = as.Date("2014-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2015-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2016-04-01"), linetype="dashed")+
  scale_color_manual("Model fit to:",values=cols)+
  scale_fill_manual("Model fit to:",values=cols)+
  ylab("Serotype 3 cases")+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  coord_cartesian(xlim=c(as.Date("2014-04-01"), as.Date("2016-04-01") )  )+
  theme(legend.position="none",
        legend.background=element_rect(color="black"))


loo4 <- ggplot(mod_inc[mod_inc$pathogen=="4",])+
  geom_line(aes(x=time, y=y, colour=label))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=label), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=label), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_vline(xintercept = as.Date("2014-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2015-04-01"), linetype="dashed")+
  geom_vline(xintercept = as.Date("2016-04-01"), linetype="dashed")+
  scale_color_manual("Model fit to:",values=cols)+
  scale_fill_manual("Model fit to:",values=cols)+
  ylab("Serotype 4 cases")+
  xlab("Date")+
  scale_x_date(date_breaks = "2 month", date_labels =  "%b\n%Y")+
  coord_cartesian(xlim=c(as.Date("2014-04-01"), as.Date("2016-04-01") )  )+
  theme(legend.position="none",
        legend.background=element_rect(color="black"))




leg <- get_legend(loo0)
loo0 <- loo0 + theme(legend.position = "none")

# Combine plots for overall figure
loo0 + leg + loo1 + loo2 + loo3 + loo4 +plot_layout(nrow=3)

ggsave("Figures/DengueLOO.png", width=12, height=8)
