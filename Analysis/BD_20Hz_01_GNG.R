#####
# Author: Julius Kricheldorff
# Analysis of the Go-NoGo task
# Date 11.01.2023
####

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(hypr)
library(tidybayes)
library(emmeans)

# set directory
wd <-"D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

# load data
GoNoGo<- read_csv(file = "GoNoGo.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

GoNoGo <- GoNoGo %>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & GoNoGo == "Go" ~ "S130Hz_Go",
    Stim_verb == "130Hz" & GoNoGo == "NoGo - Go" ~ "S130Hz_NoGo_Go",
    Stim_verb == "130Hz" & GoNoGo == "NoGo - Stop" ~ "S130Hz_NoGo_Stop",
    Stim_verb == "20Hz" & GoNoGo == "Go" ~ "S20Hz_Go",
    Stim_verb == "20Hz" & GoNoGo == "NoGo - Go" ~ "S20Hz_NoGo_Go",
    Stim_verb == "20Hz" & GoNoGo == "NoGo - Stop" ~ "S20Hz_NoGo_Stop",
    Stim_verb == "OFF" & GoNoGo == "Go" ~ "SOFF_Go",
    Stim_verb == "OFF" & GoNoGo == "NoGo - Go" ~ "SOFF_NoGo_Go",
    Stim_verb == "OFF" & GoNoGo == "NoGo - Stop" ~ "SOFF_NoGo_Stop",
  )))


# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- GoNoGo %>%
  filter(Correct_Response == 1,
         RT < 3,
         RT > 0.2,
         GoNoGo == "NoGo - Go" | GoNoGo == "Go") %>%
  mutate(RT_ms = RT*1000) %>%
  droplevels() %>%# drop Stop trials
  mutate(S130Hz = ifelse(Stim_verb == "130Hz", 1, 0),
         SOFF = ifelse(Stim_verb == "OFF", 1, 0),
         Go_diff =ifelse(GoNoGo == "NoGo - Go", 0.5, -0.5))


RT_data2 <- GoNoGo %>%
  filter(Correct_Response == 1,
        
         GoNoGo == "NoGo - Go" | GoNoGo == "Go") %>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)
# create a contrast matrix for our comparisons of interest we are interested 
# if we can observe more slowing on NoGo - Go trials in the 20Hz condition, relatively
# to the Go condition in the 20Hz condition, than the 130Hz and Off condition

prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = Go_diff), 
  prior(normal(0,  0.3), class = b, coef = Go_diff:S130Hz),
  prior(normal(0,  0.3), class = b, coef = Go_diff:SOFF), 
  prior(normal(0,  0.3), class = b, coef = S130Hz),
  prior(normal(0,  0.3), class = b, coef = SOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object RT analysis
m1_GoNoGo_shift <- bf(RT_ms ~ 1  + Go_diff + S130Hz + SOFF + S130Hz*Go_diff + SOFF*Go_diff + (1|Part_nr))


#fit the first model
fit_shifted_log_GNG <- brm(formula = m1_GoNoGo_shift,
                               family = shifted_lognormal(),
                               data = RT_data,
                               prior = prior_weakly_informed,
                               warmup = 2000,
                               iter = 12000,# 20000 is the limit necessary for bridge sampling
                               cores = 4, seed = 423,
                               control = list(adapt_delta = 0.9),
                               save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                               chains =4
)

save(fit_shifted_log_GNG, file = "../Modelle/shifted_log_GNG.rda")
load(file = "../Modelle/shifted_log_GNG.rda")

##### Error Analysis

# brms formula 
m1_GoNoGo_log <- bf(Error ~ 1  + Contrast_F + (1|Part_nr)) 

# Prior definition
prior_weakly_informed_logreg2<- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = Contrast_FGo_Stop_vs_Go),
  prior(normal(0,  1.5), class = b, coef = Contrast_FGo_Stop_vs_Stop),
  prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Stop_20Hz_vs_130Hz),
  prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Stop_20Hz_vs_OFF),
  prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Go_20Hz_vs_130Hz), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Go_20Hz_vs_OFF),
  prior(normal(0,  1.5), class = b, coef = Contrast_FGo_20Hz_vs_130Hz), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FGo_20Hz_vs_OFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)

# fit models
fit_log_flanker2 <- brm(formula = m1_GoNoGo_log,
                        family = bernoulli(link = logit),
                        data = GoNoGo,
                        prior = prior_weakly_informed_logreg2,
                        warmup = 2000,
                        iter = 12000,# 20000 is the limit necessary for bridge sampling
                        cores = 4, seed = 423,
                        save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                        chains =4)

save(fit_log_flanker2, file = "../Modelle/log_reg_GNG.rda")
load(file = "../Modelle/log_reg_GNG.rda")
