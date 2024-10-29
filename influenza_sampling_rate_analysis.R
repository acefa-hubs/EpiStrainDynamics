
# Loading required packages
library(ggplot2)
library(rstan)
library(cowplot)
library(patchwork)
library(lubridate)
library(tidyverse)

# Loading required functions
source('R/rw_analysis_scripts.R')
source('R/ps_analysis_scripts.R')
source('R/rw_single_analysis_scripts.R')
source('R/ps_single_analysis_scripts.R')


##############################################################################################################################################
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

df1 <- df1[df1$week<max_date & df1$week>=min_date,]
df2 <- df2[df2$week<max_date & df2$week>=min_date,]
df3 <- df3[df3$week<max_date & df3$week>=min_date,]

df1 <- df1[order(df1$week),]
df2 <- df2[order(df2$week),]
df3 <- df3[order(df3$week),]

df1$time <- seq(1, nrow(df1))
df2$time <- seq(1, nrow(df2))
df3$time <- seq(1, nrow(df3))
#############################################################################################################################################
## Sample from the original data for each country for supplementary analyses

inf_sampler <- function(df1, sr=100, seed=12345){
  set.seed(seed)
  
  df_new <- df1
  
  for(i in 1:nrow(df1)){
    print(i)
    
    p1 <- c(rep("A", df1$inf_A[i]), rep("B", df1$inf_B[i]), rep("N", df1$num_spec[i]-df1$inf_all[i]))
    
    if(length(p1)>=sr){
      p1_sample <- sample(p1, sr)
    } else{
      p1_sample <- p1
    }
    
    df_new$inf_A[i] <- length(p1_sample[p1_sample=="A"])
    df_new$inf_B[i] <- length(p1_sample[p1_sample=="B"])
    df_new$inf_all[i] <- df_new$inf_A[i] + df_new$inf_B[i]
    df_new$num_spec[i] <- min(sr, length(p1))
    
    p2 <- c(rep("H3N2", df1$inf_H3N2[i]), rep("H1N1", df1$inf_H1N1[i]), rep("US", df1$inf_A[i]-df1$inf_H3N2[i]-df1$inf_H1N1[i] ) )
    
    p2_sample <- sample(p2, df_new$inf_A[i])
    
    df_new$inf_H3N2[i] <- length(p2_sample[p2_sample=="H3N2"])
    df_new$inf_H1N1[i] <- length(p2_sample[p2_sample=="H1N1"])
    
  }
  
  df_new
  
}

df1_SR100 <- inf_sampler(df1, sr=100)
df1_SR50 <- inf_sampler(df1, sr=50)
df1_SR20 <- inf_sampler(df1, sr=20)

df2_SR100 <- inf_sampler(df2, sr=100)
df2_SR50 <- inf_sampler(df2, sr=50)
df2_SR20 <- inf_sampler(df2, sr=20)

df3_SR100 <- inf_sampler(df3, sr=100)
df3_SR50 <- inf_sampler(df3, sr=50)
df3_SR20 <- inf_sampler(df3, sr=20)

########################################################################################################################################################
## Set some stan settings
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = 4)

## Loading stan model used for influenza analysis

rw_influenza_mod <- stan_model('stan/rw_influenza_finalV2.stan')


##########################################################################################################################################
# Fitting the stan models
#############################################################################################################################################
## Models fit to Australia
# SR100
rw_influenza_data1_SR100 <- list(num_data = nrow(df1_SR100),
                                 num_path = 4,
                                 Y = df1_SR100$ili,
                                 P1 = t(matrix(data= c(df1_SR100$inf_A, df1_SR100$inf_B, df1_SR100$num_spec-df1_SR100$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df1_SR100$inf_H3N2, df1_SR100$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df1$time %% 1)+1,
                                 cov_structure = 1,
                                 noise_structure = 0)

rw_influenza_fit1_SR100 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data1_SR100)

saveRDS(rw_influenza_fit1_SR100, "FitStanModels/rw_influenza_fit1_SR100.rds")

# SR50
rw_influenza_data1_SR50 <- list(num_data = nrow(df1_SR50),
                                 num_path = 4,
                                 Y = df1_SR50$ili,
                                 P1 = t(matrix(data= c(df1_SR50$inf_A, df1_SR50$inf_B, df1_SR50$num_spec-df1_SR50$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df1_SR50$inf_H3N2, df1_SR50$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df1$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit1_SR50 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data1_SR50)

saveRDS(rw_influenza_fit1_SR50, "FitStanModels/rw_influenza_fit1_SR50.rds")

# SR20
rw_influenza_data1_SR20 <- list(num_data = nrow(df1_SR20),
                                 num_path = 4,
                                 Y = df1_SR20$ili,
                                 P1 = t(matrix(data= c(df1_SR20$inf_A, df1_SR20$inf_B, df1_SR20$num_spec-df1_SR20$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df1_SR20$inf_H3N2, df1_SR20$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df1$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit1_SR20 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data1_SR20)

saveRDS(rw_influenza_fit1_SR20, "FitStanModels/rw_influenza_fit1_SR20.rds")

#############################################################################################################################################
## Models fit to United States of America
# SR100
rw_influenza_data2_SR100 <- list(num_data = nrow(df2_SR100),
                                 num_path = 4,
                                 Y = df2_SR100$ili,
                                 P1 = t(matrix(data= c(df2_SR100$inf_A, df2_SR100$inf_B, df2_SR100$num_spec-df2_SR100$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df2_SR100$inf_H3N2, df2_SR100$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df2_SR100$time %% 1)+1,
                                 cov_structure = 1,
                                 noise_structure = 0)

rw_influenza_fit2_SR100 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data2_SR100)

saveRDS(rw_influenza_fit2_SR100, "FitStanModels/rw_influenza_fit2_SR100.rds")

# SR50
rw_influenza_data2_SR50 <- list(num_data = nrow(df2_SR50),
                                 num_path = 4,
                                 Y = df2_SR50$ili,
                                 P1 = t(matrix(data= c(df2_SR50$inf_A, df2_SR50$inf_B, df2_SR50$num_spec-df2_SR50$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df2_SR50$inf_H3N2, df2_SR50$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df2_SR50$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit2_SR50 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data2_SR50)

saveRDS(rw_influenza_fit2_SR50, "FitStanModels/rw_influenza_fit2_SR50.rds")

# SR20
rw_influenza_data2_SR20 <- list(num_data = nrow(df2_SR20),
                                 num_path = 4,
                                 Y = df2_SR20$ili,
                                 P1 = t(matrix(data= c(df2_SR20$inf_A, df2_SR20$inf_B, df2_SR20$num_spec-df2_SR20$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df2_SR20$inf_H3N2, df2_SR20$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df2_SR20$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit2_SR20 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data2_SR20)

saveRDS(rw_influenza_fit2_SR20, "FitStanModels/rw_influenza_fit2_SR20.rds")

#############################################################################################################################################
## Model fit to Singapore
#SR100
rw_influenza_data3_SR100 <- list(num_data = nrow(df3_SR100),
                                 num_path = 4,
                                 Y = df3_SR100$ili,
                                 P1 = t(matrix(data= c(df3_SR100$inf_A, df3_SR100$inf_B, df3_SR100$num_spec-df3_SR100$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df3_SR100$inf_H3N2, df3_SR100$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df3_SR100$time %% 1)+1,
                                 cov_structure = 1,
                                 noise_structure = 0)

rw_influenza_fit3_SR100 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data3_SR100)

saveRDS(rw_influenza_fit3_SR100, "FitStanModels/rw_influenza_fit3_SR100.rds")


#SR50
rw_influenza_data3_SR50 <- list(num_data = nrow(df3_SR50),
                                 num_path = 4,
                                 Y = df3_SR50$ili,
                                 P1 = t(matrix(data= c(df3_SR50$inf_A, df3_SR50$inf_B, df3_SR50$num_spec-df3_SR50$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df3_SR50$inf_H3N2, df3_SR50$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df3_SR50$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit3_SR50 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data3_SR50)

saveRDS(rw_influenza_fit3_SR50, "FitStanModels/rw_influenza_fit3_SR50.rds")


#SR20
rw_influenza_data3_SR20 <- list(num_data = nrow(df3_SR20),
                                 num_path = 4,
                                 Y = df3_SR20$ili,
                                 P1 = t(matrix(data= c(df3_SR20$inf_A, df3_SR20$inf_B, df3_SR20$num_spec-df3_SR20$inf_all), ncol=3) ),
                                 P2 = t(matrix(data= c(df3_SR20$inf_H3N2, df3_SR20$inf_H1N1), ncol=2) ),
                                 week_effect = 1,
                                 DOW = (df3_SR20$time %% 1)+1,
                                 cov_structure = 1,
                                noise_structure = 0)

rw_influenza_fit3_SR20 <- sampling(rw_influenza_mod,
                                    iter=5000,
                                    warmup =1000,
                                    chains=4,
                                    data = rw_influenza_data3_SR20)

saveRDS(rw_influenza_fit3_SR20, "FitStanModels/rw_influenza_fit3_SR20.rds")

######################################################################################################################################################
# Plotting figures
######################################################################################################################################################

############################################################################################
## Supplementary Figure 9
############################################################################################

rw_influenza_fit1_SR20 <- readRDS("FitStanModels/rw_influenza_fit1_SR20.rds")
rw_influenza_fit1_SR50 <- readRDS("FitStanModels/rw_influenza_fit1_SR50.rds")
rw_influenza_fit1_SR100 <- readRDS("FitStanModels/rw_influenza_fit1_SR100.rds")
rw_influenza_fit1_SRX <- readRDS("FitStanModels/rw_influenza_fit1_N0.rds")


rw_mod_inc1_SR20 <- rw_incidence(rw_influenza_fit1_SR20, num_days = nrow(df1), time_labels = df1$week)
rw_mod_inc1_SR50 <- rw_incidence(rw_influenza_fit1_SR50, num_days = nrow(df1), time_labels = df1$week)
rw_mod_inc1_SR100 <- rw_incidence(rw_influenza_fit1_SR100, num_days = nrow(df1), time_labels = df1$week)
rw_mod_inc1_SRX <- rw_incidence(rw_influenza_fit1_SRX, num_days = nrow(df1), time_labels = df1$week)

rw_mod_inc1_SR20$sampling <- "20/week max"
rw_mod_inc1_SR50$sampling <- "50/week max"
rw_mod_inc1_SR100$sampling <- "100/week max"
rw_mod_inc1_SRX$sampling <- "Original"

rw_mod_inc1 <- rbind(rw_mod_inc1_SR20,
                     rw_mod_inc1_SR50,
                     rw_mod_inc1_SR100,
                     rw_mod_inc1_SRX)

rw_mod_inc1$sampling <- factor(rw_mod_inc1$sampling, levels = c( "20/week max", "50/week max", "100/week max", "Original") )
rw_mod_inc1$pathogen <- factor(rw_mod_inc1$pathogen, levels=c("Influenza A H1N1", "Influenza A H3N2", "Influenza B", "Unknown", "Total"))


plt1<-ggplot(rw_mod_inc1[rw_mod_inc1$pathogen=="Total",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  geom_point(data=df1, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Total influenza-like illness", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt2<-ggplot(rw_mod_inc1[rw_mod_inc1$pathogen=="Influenza A H3N2",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H3N2", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt3<-ggplot(rw_mod_inc1[rw_mod_inc1$pathogen=="Influenza A H1N1",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H1N1", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt4<-ggplot(rw_mod_inc1[rw_mod_inc1$pathogen=="Influenza B",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza B", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt5<-ggplot(rw_mod_inc1[rw_mod_inc1$pathogen=="Unknown",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df1, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer("Sampling rate", palette = "Dark2")+
  scale_fill_brewer("Sampling rate", palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Unknown", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")


plt2i <- plt2+coord_cartesian(xlim=c(as.Date("2016-04-01"), as.Date("2017-01-01") ),ylim = c(0, 125))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt3i <- plt3+coord_cartesian(xlim=c(as.Date("2016-01-01"), as.Date("2017-01-01") ),ylim = c(0, 15))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt4i <- plt4+coord_cartesian(xlim=c(as.Date("2016-01-01"), as.Date("2017-01-01") ),ylim = c(0, 20))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))


plt1 <- plt1+ 
  labs(tag="A")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
plt2 <- plt2+ 
  labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2016-04-01"), xmax = as.Date("2017-01-01"), ymin = -2, ymax = 125, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt3 <- plt3+ 
  labs(tag="C")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2016-01-01"), xmax = as.Date("2017-01-01"), ymin = -2, ymax = 15, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt4 <- plt4+ 
  labs(tag="D")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2016-01-01"), xmax = as.Date("2017-01-01"), ymin = -2, ymax = 20, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt5 <- plt5+ 
  labs(tag="E")+
  theme(plot.tag.position = c(0.01,0.97),
        legend.position = "bottom")


patch1 <- plt1 + plt2 + plt3 + plt4 + plt5 + plot_layout(nrow=5)

patch2 <- plt2i + plt3i + plt4i +
  plot_layout(nrow=3)

cowplot::plot_grid(patch1, patch2, rel_widths = c(3,1), rel_heights = c(1,0.8))


ggsave("Figures/Influenza_SR_Aus.png", width=10, height=14)


############################################################################################
## Supplementary Figure 11
############################################################################################

rw_influenza_fit2_SR20 <- readRDS("FitStanModels/rw_influenza_fit2_SR20.rds")
rw_influenza_fit2_SR50 <- readRDS("FitStanModels/rw_influenza_fit2_SR50.rds")
rw_influenza_fit2_SR100 <- readRDS("FitStanModels/rw_influenza_fit2_SR100.rds")
rw_influenza_fit2_SRX <- readRDS("FitStanModels/rw_influenza_fit2_N0.rds")


rw_mod_inc2_SR20 <- rw_incidence(rw_influenza_fit2_SR20, num_days = nrow(df2), time_labels = df2$week)
rw_mod_inc2_SR50 <- rw_incidence(rw_influenza_fit2_SR50, num_days = nrow(df2), time_labels = df2$week)
rw_mod_inc2_SR100 <- rw_incidence(rw_influenza_fit2_SR100, num_days = nrow(df2), time_labels = df2$week)
rw_mod_inc2_SRX <- rw_incidence(rw_influenza_fit2_SRX, num_days = nrow(df2), time_labels = df2$week)

rw_mod_inc2_SR20$sampling <- "20/week max"
rw_mod_inc2_SR50$sampling <- "50/week max"
rw_mod_inc2_SR100$sampling <- "100/week max"
rw_mod_inc2_SRX$sampling <- "Original"

rw_mod_inc2 <- rbind(rw_mod_inc2_SR20,
                     rw_mod_inc2_SR50,
                     rw_mod_inc2_SR100,
                     rw_mod_inc2_SRX)

rw_mod_inc2$sampling <- factor(rw_mod_inc2$sampling, levels = c( "20/week max", "50/week max", "100/week max", "Original") )


rw_mod_inc2$pathogen <- factor(rw_mod_inc2$pathogen, levels=c("Influenza A H1N1", "Influenza A H3N2", "Influenza B", "Unknown", "Total"))



plt1<-ggplot(rw_mod_inc2[rw_mod_inc2$pathogen=="Total",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  geom_point(data=df2, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Total influenza-like illness", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt2<-ggplot(rw_mod_inc2[rw_mod_inc2$pathogen=="Influenza A H3N2",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df2, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H3N2", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt3<-ggplot(rw_mod_inc2[rw_mod_inc2$pathogen=="Influenza A H1N1",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df2, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H1N1", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt4<-ggplot(rw_mod_inc2[rw_mod_inc2$pathogen=="Influenza B",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df2, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza B", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt5<-ggplot(rw_mod_inc2[rw_mod_inc2$pathogen=="Unknown",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df2, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer("Sampling rate",palette = "Dark2")+
  scale_fill_brewer("Sampling rate",palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Unknown", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")



plt2i <- plt2+coord_cartesian(xlim=c(as.Date("2018-11-01"), as.Date("2019-06-01") ),ylim = c(0, 17000))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt3i <- plt3+coord_cartesian(xlim=c(as.Date("2018-11-01"), as.Date("2019-05-01") ),ylim = c(0, 32000))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt4i <- plt4+coord_cartesian(xlim=c(as.Date("2018-07-01"), as.Date("2019-07-01") ),ylim = c(0, 2000))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))


plt1 <- plt1+ 
  labs(tag="A")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
plt2 <- plt2+ 
  labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2018-11-01"), xmax = as.Date("2019-06-01"), ymin = -200, ymax = 17000, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt3 <- plt3+ 
  labs(tag="C")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2018-11-01"), xmax = as.Date("2019-05-01"), ymin = -200, ymax = 32000, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt4 <- plt4+ 
  labs(tag="D")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2018-07-01"), xmax = as.Date("2019-07-01"), ymin = -200, ymax = 2000, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt5 <- plt5+ 
  labs(tag="E")+
  theme(plot.tag.position = c(0.01,0.97),
        legend.position = "bottom")


patch1 <- plt1 + plt2 + plt3 + plt4 + plt5 + plot_layout(nrow=5)

patch2 <- plt2i + plt3i + plt4i +
  plot_layout(nrow=3)

cowplot::plot_grid(patch1, patch2, rel_widths = c(3,1), rel_heights = c(1,0.8))




ggsave("Figures/Influenza_SR_USA.png", width=10, height=14)



############################################################################################
## Supplementary Figure 10
############################################################################################

rw_influenza_fit3_SR20 <- readRDS("FitStanModels/rw_influenza_fit3_SR20.rds")
rw_influenza_fit3_SR50 <- readRDS("FitStanModels/rw_influenza_fit3_SR50.rds")
rw_influenza_fit3_SR100 <- readRDS("FitStanModels/rw_influenza_fit3_SR100.rds")
rw_influenza_fit3_SRX <- readRDS("FitStanModels/rw_influenza_fit3_N0.rds")

rw_mod_inc3_SR20 <- rw_incidence(rw_influenza_fit3_SR20, num_days = nrow(df3), time_labels = df3$week)
rw_mod_inc3_SR50 <- rw_incidence(rw_influenza_fit3_SR50, num_days = nrow(df3), time_labels = df3$week)
rw_mod_inc3_SR100 <- rw_incidence(rw_influenza_fit3_SR100, num_days = nrow(df3), time_labels = df3$week)
rw_mod_inc3_SRX <- rw_incidence(rw_influenza_fit3_SRX, num_days = nrow(df3), time_labels = df3$week)

rw_mod_inc3_SR20$sampling <- "20/week max"
rw_mod_inc3_SR50$sampling <- "50/week max"
rw_mod_inc3_SR100$sampling <- "100/week max"
rw_mod_inc3_SRX$sampling <- "Original"

rw_mod_inc3 <- rbind(rw_mod_inc3_SR20,
                     rw_mod_inc3_SR50,
                     rw_mod_inc3_SR100,
                     rw_mod_inc3_SRX)


rw_mod_inc3$sampling <- factor(rw_mod_inc3$sampling, levels = c("20/week max", "50/week max", "100/week max", "Original") )


rw_mod_inc3$pathogen <- factor(rw_mod_inc3$pathogen, levels=c("Influenza A H1N1", "Influenza A H3N2", "Influenza B", "Unknown", "Total"))



plt1<-ggplot(rw_mod_inc3[rw_mod_inc3$pathogen=="Total",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  geom_point(data=df3, aes(x=week, y=ili), size=0.2)+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Total influenza-like illness", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt2<-ggplot(rw_mod_inc3[rw_mod_inc3$pathogen=="Influenza A H3N2",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df3, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H3N2", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt3<-ggplot(rw_mod_inc3[rw_mod_inc3$pathogen=="Influenza A H1N1",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df3, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza A H1N1", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt4<-ggplot(rw_mod_inc3[rw_mod_inc3$pathogen=="Influenza B",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df3, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Influenza B", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")

plt5<-ggplot(rw_mod_inc3[rw_mod_inc3$pathogen=="Unknown",])+
  geom_line(aes(x=time, y=y, color=sampling), size=0.2)+
  #geom_ribbon(aes(x=time, y=y, ymin=lb_50, ymax=ub_50, fill=sampling), alpha=0.2)+
  geom_ribbon(aes(x=time, y=y, ymin=lb_95, ymax=ub_95, fill=sampling), alpha=0.2)+
  #facet_wrap(.~pathogen)+
  #geom_point(data=df3, aes(x=week, y=ili))+
  coord_cartesian(xlim=c(as.Date("2012-01-01")+180, as.Date("2024-01-01")-180))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%Y")+
  scale_color_brewer("Sampling rate", palette = "Dark2")+
  scale_fill_brewer("Sampling rate",palette = "Dark2")+
  xlab("Date")+
  ylab("Modelled cases")+
  theme_bw()+
  geom_label(aes(label = "Unknown", x=as.Date("2018-01-01"), y=Inf ), vjust=1.2)+
  theme(legend.position = "none")


plt2i <- plt2+coord_cartesian(xlim=c(as.Date("2017-07-01"), as.Date("2018-07-01") ),ylim = c(0, 100))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt3i <- plt3+coord_cartesian(xlim=c(as.Date("2017-07-01"), as.Date("2018-07-01") ),ylim = c(0, 100))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))

plt4i <- plt4+coord_cartesian(xlim=c(as.Date("2017-07-01"), as.Date("2018-07-01") ),ylim = c(0, 300))+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(color="black", fill=alpha("grey",0.2)))


plt1 <- plt1+ 
  labs(tag="A")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
plt2 <- plt2+ 
  labs(tag="B")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2017-07-01"), xmax = as.Date("2018-07-01"), ymin = -5, ymax = 100, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt3 <- plt3+ 
  labs(tag="C")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2017-07-01"), xmax = as.Date("2018-07-01"), ymin = -5, ymax = 100, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt4 <- plt4+ 
  labs(tag="D")+
  theme(plot.tag.position = c(0.01,0.97),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  annotate("rect", xmin = as.Date("2017-07-01"), xmax = as.Date("2018-07-01"), ymin = -5, ymax = 300, color="black", 
           fill='grey', alpha = 0.2, linetype = "dashed")

plt5 <- plt5+ 
  labs(tag="E")+
  theme(plot.tag.position = c(0.01,0.97),
        legend.position = "bottom")


patch1 <- plt1 + plt2 + plt3 + plt4 + plt5 + plot_layout(nrow=5)

patch2 <- plt2i + plt3i + plt4i +
  plot_layout(nrow=3)

cowplot::plot_grid(patch1, patch2, rel_widths = c(3,1), rel_heights = c(1,0.8))





ggsave("Figures/Influenza_SR_SIN.png", width=10, height=14)


############################################################################################
## Supplementary Figures 1-3 (raw data in models)
############################################################################################

###################################################
df1$sampling <- "Original testing rate"
df1_SR100$sampling <- "100 tests/week max"
df1_SR50$sampling <- "50 tests/week max"
df1_SR20$sampling <- "20 tests/week max"

df1T <- rbind(df1, df1_SR100, df1_SR50, df1_SR20)

df1T$num_spec <- df1T$num_spec - df1T$inf_all
df1T$inf_A_ns <- df1T$inf_A - df1T$inf_H3N2 - df1T$inf_H1N1

df_new1 <- pivot_longer(df1T, cols=colnames(df1T)[c(7,13,10,11,12)])

df_new1$name <- factor(df_new1$name, labels = c("Influenza A (not subtyped)", "Influenza B", "Influenza A H1N1", "Influenza A H3N2", "Negative"))
df_new1$name <- factor(df_new1$name, levels = c("Negative", "Influenza B", "Influenza A H1N1", "Influenza A H3N2","Influenza A (not subtyped)"))

df_new1$sampling <- factor(df_new1$sampling, levels = c("Original testing rate","100 tests/week max","50 tests/week max","20 tests/week max"))



###################################################
df2$sampling <- "Original testing rate"
df2_SR100$sampling <- "100 tests/week max"
df2_SR50$sampling <- "50 tests/week max"
df2_SR20$sampling <- "20 tests/week max"

df2T <- rbind(df2, df2_SR100, df2_SR50, df2_SR20)

df2T$num_spec <- df2T$num_spec - df2T$inf_all
df2T$inf_A_ns <- df2T$inf_A - df2T$inf_H3N2 - df2T$inf_H1N1

df_new2 <- pivot_longer(df2T, cols=colnames(df2T)[c(7,13,10,11,12)])

df_new2$name <- factor(df_new2$name, labels = c("Influenza A (not subtyped)", "Influenza B", "Influenza A H1N1", "Influenza A H3N2", "Negative"))
df_new2$name <- factor(df_new2$name, levels = c("Negative", "Influenza B", "Influenza A H1N1", "Influenza A H3N2","Influenza A (not subtyped)"))

df_new2$sampling <- factor(df_new2$sampling, levels = c("Original testing rate","100 tests/week max","50 tests/week max","20 tests/week max"))


###################################################
df3$sampling <- "Original testing rate"
df3_SR100$sampling <- "100 tests/week max"
df3_SR50$sampling <- "50 tests/week max"
df3_SR20$sampling <- "20 tests/week max"

df3T <- rbind(df3, df3_SR100, df3_SR50, df3_SR20)

df3T$num_spec <- df3T$num_spec - df3T$inf_all
df3T$inf_A_ns <- df3T$inf_A - df3T$inf_H3N2 - df3T$inf_H1N1

df_new3 <- pivot_longer(df3T, cols=colnames(df3T)[c(7,13,10,11,12)])

df_new3$name <- factor(df_new3$name, labels = c("Influenza A (not subtyped)", "Influenza B", "Influenza A H1N1", "Influenza A H3N2", "Negative"))
df_new3$name <- factor(df_new3$name, levels = c("Negative", "Influenza B", "Influenza A H1N1", "Influenza A H3N2","Influenza A (not subtyped)"))

df_new3$sampling <- factor(df_new3$sampling, levels = c("Original testing rate","100 tests/week max","50 tests/week max","20 tests/week max"))



cols2 <- c("grey", "navy", "red4", "gold3","green4")

###################################################
ggplot(data=df_new1)+
  geom_col(aes(x=week, y=value, fill=name))+
  facet_wrap(.~sampling, nrow=4, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Number of tests")+
  scale_fill_manual("Test result",values=cols2)+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  #guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))


ggsave("Figures/InfluenzaRaw1.png", width=10, height=12)


###################################################
ggplot(data=df_new2)+
  geom_col(aes(x=week, y=value, fill=name))+
  facet_wrap(.~sampling, nrow=4, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Number of tests")+
  scale_fill_manual("Test result",values=cols2)+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  #guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))

ggsave("Figures/InfluenzaRaw2.png", width=10, height=12)

###################################################
ggplot(data=df_new3)+
  geom_col(aes(x=week, y=value, fill=name))+
  facet_wrap(.~sampling, nrow=4, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_cartesian(xlim=c(min_date+30, max_date-30))+
  scale_x_date(date_breaks = "1 year", date_labels =  "%b\n%Y")+
  ylab("Modelled cases")+
  xlab("Date")+
  scale_y_continuous("Number of tests")+
  scale_fill_manual("Test result",values=cols2)+
  #guides(fill=guide_legend(nrow=2,byrow=TRUE))+
  #guides(color=guide_legend(nrow=2,byrow=TRUE))+
  theme(legend.position = "bottom",
        legend.background = element_rect(color="black"),
        strip.background = element_rect(color="black", fill=NA))

ggsave("Figures/InfluenzaRaw3.png", width=10, height=12)









