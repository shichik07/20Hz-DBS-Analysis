#####
# Author: Julius Kricheldorff
# Analysis of the simple reaction time data
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
# Date 11.01.2023
####

# Address git : C/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(hypr)
library(tiybayes)
library(bayestestR)

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted')

# load data
SimpleRT<- read_csv(file = "SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))
# create a contrast matrix for our comparisons of interest
# Contrasts only for the list-wide effect only
SRT_Contrast <- hypr(
  S130_S20 = S20Hz ~ S130Hz, 
  S20_SOFF = S20Hz ~ SOFF, 
  levels = c("S130Hz", "SOFF", "S20Hz")
)
SRT_Contrast

# assign the generated contrast matrix to the List Wide Factor
contrasts(SimpleRT$StimCon) <- contr.hypothesis(SRT_Contrast)
contrasts(SimpleRT$StimCon)   


# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- SimpleRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

RT_data2 <-  SimpleRT%>%
  filter(Correct_Response == 1)%>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)


# Prior informed weakly Item Specific
prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.1), class = b, coef = StimConS130_S20), 
  prior(normal(0,  0.1), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  0.1), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_SRT <- bf(RT_ms ~ 1  + StimCon + (1|Part_nr))

#get_prior(formula = m1_SRT, data = RT_data, family = shifted_lognormal())

#fit the first model
fit_shifted_log_SRT <- brm(formula = m1_SRT,
                  family = shifted_lognormal(),
                  data = RT_data,
                  prior = prior_weakly_informed,
                  warmup = 2000,
                  iter = 12000,# 20000 is the limit necessary for bridge sampling
                  cores = 4, seed = 423,
                  #control = list(adapt_delta = 0.9),
                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                  chains =4
)

save(fit_shifted_log_SRT, file = "E:/20Hz/Data/Modelle/shifted_log_SRT.rda")
load(file = "E:/20Hz/Data/Modelle/shifted_log_SRT.rda")

# posteriro predictive checks
pp_check(fit_shifted_log_SRT, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_shifted_log_SRT, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_shifted_log_SRT, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_shifted_log_SRT, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_shifted_log_SRT, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "StimCon")
# looks good

# Great, the model looks perfect. Lets see the effect sizes, 

# first get variable names
get_variables(fit_shifted_log_SRT)

post_eff <- fit_shifted_log_SRT %>%
  spread_draws(b_Intercept, b_StimConS20_SOFF, b_StimConS130_S20, sigma, ndt) %>%
  mutate(Est_S130Hz = exp(b_Intercept + (-(2/3))*b_StimConS130_S20 + (1/3)*b_StimConS20_SOFF) + sigma/2 + ndt,
         Est_S20Hz = exp(b_Intercept + (1/3)*b_StimConS130_S20 + (1/3)*b_StimConS20_SOFF) + sigma/2 + ndt,
         Est_SOFF = exp(b_Intercept + (1/3)*b_StimConS130_S20 + (-(2/3))*b_StimConS20_SOFF) + sigma/2 + ndt,
         S20_vs_S130 = Est_S20Hz - Est_S130Hz,
         S20_vs_SOFF = Est_S20Hz - Est_SOFF) %>%
  summarise_draws()

# Okay, we see a lot of variance (makes sense since the timing measurements were imprecise)
# but it appears, there is no difference in timing

## Next, let us look at the accuracy data

prior_weakly_informed_logreg <- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = StimConS130_S20), 
  prior(normal(0,  1.5), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)


#get_prior(formula = m1_SRT_log, data = RT_data, family = bernoulli(link = logit))

# brmsformula object Item Specific
m1_SRT_log <- bf(Error ~ 1  + StimCon + (1|Part_nr)) 

#### Fit Inducer Models ####


# we should consider varying non-decision times between the groups
fit_logReg_SRT <- brm(formula = m1_SRT_log,
                                  family = bernoulli(link = logit),
                                  data = SimpleRT,
                                  prior = prior_weakly_informed_logreg,
                                  warmup = 2000,
                                  iter = 12000,# 20000 is the limit necessary for bridge sampling
                                  cores = 4, seed = 423,
                                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                  chains =4)

save(fit_logReg_SRT, file = "E:/20Hz/Data/Modelle/log_reg_SRT.rda")
load(file = "E:/20Hz/Data/Modelle/log_reg_SRT.rda")
# posteriro predictive checks
pp_check(fit_logReg_SRT, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_logReg_SRT, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_logReg_SRT, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_logReg_SRT, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_logReg_SRT, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "StimCon")
# looks good

# first get variable names
get_variables(fit_logReg_SRT)

post_eff_log_SRT <- fit_logReg_SRT %>%
  spread_draws(b_Intercept, b_StimConS20_SOFF, b_StimConS130_S20) %>%
  mutate(Est_S130Hz = plogis(b_Intercept + (-(2/3))*b_StimConS130_S20 + (1/3)*b_StimConS20_SOFF),
         Est_S20Hz = plogis(b_Intercept + (1/3)*b_StimConS130_S20 + (1/3)*b_StimConS20_SOFF),
         Est_SOFF = plogis(b_Intercept + (1/3)*b_StimConS130_S20 + (-(2/3))*b_StimConS20_SOFF),
         S20_vs_S130 = Est_S20Hz - Est_S130Hz,
         S20_vs_SOFF = Est_S20Hz - Est_SOFF) %>%
  summarise_draws()


