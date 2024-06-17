setwd("C:/Users/EALESO/OneDrive - The University of Melbourne/Projects/EpidemicAnalytics")

# Loading required packages
library(ggplot2)
library(rstan)


# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')

# Load some real data 
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

## Loading and fitting Stan models example code for a few representative models

#############################################################################################################################################
## Example 1
# Random-walk model, influenza A not-subtyped nuance, fit to Australian influenza data (obtained from WHO)
rw_influenza_mod <- stan_model('stan/rw_influenza.stan')

rw_influenza_data1 <- list(num_data = nrow(df1),
                           num_path = 4,
                           Y = df1$ili,
                           P1 = t(matrix(data= c(df1$inf_A, df1$inf_B, df1$num_spec-df1$inf_all), ncol=3) ),
                           P2 = t(matrix(data= c(df1$inf_H3N2, df1$inf_H1N1), ncol=2) ),
                           week_effect = 1,
                           DOW = (df1$time %% 1)+1,
                           cov_structure = 1)

rw_influenza_fit <- sampling(rw_influenza_mod,
                   iter=2000,
                   warmup =500,
                   chains=4,
                   data = rw_influenza_data1)

######################################################################################################################################################
## Example 2 
#Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated dataset to demonstrate utility during flu-season/pandemic

ps_mp_mod <- stan_model('stan/ps_mp.stan')

# Fitting the model 'mid-season', only include first 140 days
dfS2 <- dfS1[dfS1$t < 140,]

# Calculate the locations of equally spaced knots
knots <- get_knots(dfS2$t, days_per_knot = 5, spline_degree = 3)

ps_mp_data <- list(num_data = nrow(dfS2),
                   num_path = 3,
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = dfS2$y,
                   X = dfS2$t,
                   P = t(matrix(data= c(dfS2$H3N2, dfS2$H1N1, dfS2$B), ncol=3)) ) 

ps_mp_fit <- sampling(ps_mp_mod,
                      iter= 2000,
                      warmup = 500,
                      chains=4,
                      data = ps_mp_data)


######################################################################################################################################################
## Example 3
#Penalised-spline model, multiple-pathogens, fit up to day 140 of simulated dataset to demonstrate utility during flu-season/pandemic

ps_single_mod <- stan_model('stan/ps_single.stan')

# Fitting the model 'mid-season', only include first 140 days

# Calculate the locations of equally spaced knots
knots <- get_knots(dfC1$time, days_per_knot = 5, spline_degree = 3)

ps_single_data <- list(num_data = nrow(dfC1),
                   num_knots = length(knots),
                   knots = knots,
                   spline_degree=3,
                   Y = dfC1$metric_value,
                   X = dfC1$time,
                   week_effect=7,
                   DOW = (dfC1$time %% 7)+1) 

ps_single_fit <- sampling(ps_single_mod,
                      iter= 2000,
                      warmup = 500,
                      chains=4,
                      data = ps_single_data)



######################################################################################################################################################
## Example 1 Figures 

rw_mod_inc <- rw_incidence(rw_influenza_fit, num_days = nrow(df1), time_labels = df1$week)

ggplot(rw_mod_inc)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=df1, aes(x=week, y=ili))+
  theme_bw()

ggplot(rw_mod_inc[rw_mod_inc$pathogen %in% c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"),])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=df1, aes(x=week, y=ili*(inf_A/num_spec)*inf_H3N2/(inf_H3N2+inf_H1N1),  color="Influenza A H3N2"))+
  geom_point(data=df1, aes(x=week, y=ili*(inf_A/num_spec)*inf_H1N1/(inf_H3N2+inf_H1N1),  color="Influenza A H1N1"))+
  geom_point(data=df1, aes(x=week, y=ili*(inf_B/num_spec),  color="Influenza B"))+
  theme_bw()

rw_mod_prop <- rw_proportion(rw_influenza_fit, num_days = nrow(df1), time_labels = df1$week)

ggplot(rw_mod_prop[rw_mod_prop$pathogen%in%c("Influenza", "Influenza A", "Influenza B"),])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=df1, aes(x=week, y=inf_all/num_spec, color="Influenza"))+
  geom_point(data=df1, aes(x=week, y=inf_A/num_spec, color="Influenza A"))+
  geom_point(data=df1, aes(x=week, y=inf_B/num_spec, color="Influenza B"))

ggplot(rw_mod_prop[rw_mod_prop$pathogen%in%c("H3N2", "H1N1"),])+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=df1, aes(x=week, y=inf_H3N2/(inf_H3N2+inf_H1N1), color="H3N2"))+
  geom_point(data=df1, aes(x=week, y=inf_H1N1/(inf_H3N2+inf_H1N1), color="H1N1"))


rw_mod_gr <- rw_growth_rate(rw_influenza_fit, num_days = nrow(df1), time_labels = df1$week)

ggplot(rw_mod_gr)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2")+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(
    "Growth rate", 
    sec.axis = sec_axis(~., breaks=c(log(2)*7/3, log(2)*7/5, log(2)*7/7, log(2)*7/14,0, -log(2)*7/14, -log(2)*7/7, -log(2)*7/5, -log(2)*7/3), labels = c("3", "5", "7", "14", expression(infinity/-infinity), "-14", "-7", "-5", "-3"), name = "Doubling(+) / Halving(-) time (days)")
  )

gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)

rw_mod_Rt <- rw_Rt(rw_influenza_fit, num_days = nrow(df1), time_labels = df1$week,
               gi_dist = function(x) gammaDist(b=1, n=1, x))

ggplot(rw_mod_Rt)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 1, linetype="dashed")

######################################################################################################################################################
## Example 2 - Figures
# Note that the model was fit to dfS2 data which only goes up to day 140
# Some of the plots include the full data from dfS1 which goes up to day 365 (plot functions don't need this)

mod_inc <- ps_incidence(ps_mp_fit, dfS2$t, num_days=nrow(dfS2), time_labels = dfS2$t,
                        num_path = 3,
                        pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))


ggplot(mod_inc)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=dfS1, aes(x=t+0.2, y=y))
#geom_point(data=dfS1, aes(x=t, y=y*H3N2/tests, color="Influenza A H3N2"))+ 
#geom_point(data=dfS1, aes(x=t, y=y*H1N1/tests, color="Influenza A H1N1"))+
#geom_point(data=dfS1, aes(x=t, y=y*B/tests, color="Influenza B"))

mod_prop <- ps_proportion(ps_mp_fit, dfS2$t, num_days=nrow(dfS2), time_labels = dfS2$t,
                          num_path = 3,
                          comb_num=list(c(1), c(2), c(3)),
                          comb_den=list(c(1,2,3), c(1,2,3), c(1,2,3)),
                          comb_names=c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))

ggplot(mod_prop)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  geom_point(data=dfS1, aes(x=t, y=H3N2/tests, color="Influenza A H3N2"))+
  geom_point(data=dfS1, aes(x=t, y=H1N1/tests, color="Influenza A H1N1"))+
  geom_point(data=dfS1, aes(x=t, y=B/tests, color="Influenza B"))


mod_gr <- ps_growth_rate(ps_mp_fit, dfS2$t, num_days=nrow(dfS2), time_labels = dfS2$t,
                         num_path = 3,
                         pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))


ggplot(mod_gr)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(
    "Growth rate", 
    sec.axis = sec_axis(~., breaks=c(log(2)/3, log(2)/5, log(2)/7, log(2)/14,0, -log(2)/14, -log(2)/7, -log(2)/5, -log(2)/3), labels = c("3", "5", "7", "14", expression(infinity/-infinity), "-14", "-7", "-5", "-3"), name = "Doubling(+) / Halving(-) time (days)")
  )


gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)

mod_Rt <- ps_Rt(ps_mp_fit, dfS2$t, num_days=nrow(dfS2), time_labels = dfS2$t,
                num_path = 3,
                tau_max=7,
                gi_dist = function(x) gammaDist(b=2, n=2, x),
                pathogen_names = c("Influenza A H3N2", "Influenza A H1N1", "Influenza B"))

ggplot(mod_Rt)+
  geom_line(aes(x=time, y=y, color=pathogen))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=pathogen), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=pathogen), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 1, linetype="dashed")



######################################################################################################################################################
## Example 3 - Figures

cov_mod_inc <- ps_single_incidence(ps_single_fit, dfC1$time, num_days=nrow(dfC1), time_labels = as.Date(dfC1$date) )


ggplot(cov_mod_inc)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  geom_point(data=dfC1, aes(x=as.Date(date), y=metric_value))+
  geom_line(data=dfC1, aes(x=as.Date(date), y=metric_value))


cov_mod_inc_dow <- ps_single_incidence_dow(ps_single_fit, dfC1$time, num_days=nrow(dfC1), time_labels = as.Date(dfC1$date),
                                       week_effect = 7, DOW = (dfC1$time %% 7)+1)

ggplot(cov_mod_inc_dow)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  geom_point(data=dfC1, aes(x=as.Date(date), y=metric_value))+
  geom_line(data=dfC1, aes(x=as.Date(date), y=metric_value))


cov_mod_gr <- ps_single_growth_rate(ps_single_fit, dfC1$time, num_days=nrow(dfC1), time_labels = as.Date(dfC1$date) )


ggplot(cov_mod_gr)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_y_continuous(
    "Growth rate", 
    sec.axis = sec_axis(~., breaks=c(log(2)/3, log(2)/5, log(2)/7, log(2)/14,0, -log(2)/14, -log(2)/7, -log(2)/5, -log(2)/3), labels = c("3", "5", "7", "14", expression(infinity/-infinity), "-14", "-7", "-5", "-3"), name = "Doubling(+) / Halving(-) time (days)")
  )


gammaDist <- function(b, n, a) (b**n) * (a**(n-1)) * exp(-b*a) / gamma(n)

cov_mod_Rt <- ps_single_Rt(ps_single_fit, dfC1$time, num_days=nrow(dfC1), time_labels = as.Date(dfC1$date),
                tau_max=7,
                gi_dist = function(x) gammaDist(b=2, n=2, x))

ggplot(cov_mod_Rt)+
  geom_line(aes(x=time, y=y))+
  geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95), alpha=0.2)+
  theme_bw(base_size = 14)+
  geom_hline(yintercept = 1, linetype="dashed")
