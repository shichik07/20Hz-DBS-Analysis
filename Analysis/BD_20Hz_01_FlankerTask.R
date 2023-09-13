#####
# Author: Julius Kricheldorff
# Analysis of the flanker task
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
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

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted')

# load data
flankerRT<- read_csv(file = "flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

contrasts(flankerRT$StimCon)



# create a contrast matrix for our comparisons of interest
# Contrasts only for the list-wide effect only
flanker_Contrast <- hypr(
  Congruency = (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3 ~ (S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3,
  Stim_20v130 = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (S130Hz_incongruent + S130Hz_congruent)/2,
  Stim_20vOFF = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (SOFF_incongruent + SOFF_congruent)/2,
  Stroop_20v130 = (S20Hz_incongruent - S20Hz_congruent) ~ (S130Hz_incongruent - S130Hz_congruent),
  Stroop_20vOFF = (S20Hz_incongruent - S20Hz_congruent) ~ (SOFF_incongruent - SOFF_congruent),
  levels = c("S130Hz_incongruent", "S130Hz_congruent", "SOFF_congruent", 
             "SOFF_incongruent", "S20Hz_incongruent", "S20Hz_congruent")
)

flanker_Contrast

# create a variable for the contrasts

flankerRT <- flankerRT %>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & Congruency == "congruent" ~ "S130Hz_congruent",
    Stim_verb == "130Hz" & Congruency == "incongruent" ~ "S130Hz_incongruent",
    Stim_verb == "20Hz" & Congruency == "congruent" ~ "S20Hz_congruent",
    Stim_verb == "20Hz" & Congruency == "incongruent" ~ "S20Hz_incongruent",
    Stim_verb == "OFF" & Congruency == "congruent" ~ "SOFF_congruent",
    Stim_verb == "OFF" & Congruency == "incongruent" ~ "SOFF_incongruent",
  )))

# assign the generated contrast matrix to the List Wide Factor
contrasts(flankerRT$Contrast_F) <- contr.hypothesis(flanker_Contrast)
contrasts(flankerRT$Contrast_F)   


# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- flankerRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

RT_data2 <-  flankerRT %>%
  filter(Correct_Response == 1)%>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)

# Prior informed weakly Item Specific
prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = Contrast_FCongruency), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FStroop_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FStroop_20vOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_flanker <- bf(RT_ms ~ 1  + Contrast_F + (Contrast_F|Part_nr))

#get_prior(formula = m1_flanker, data = RT_data, family = shifted_lognormal())

#fit the first model
fit_shifted_log_flanker <- brm(formula = m1_flanker,
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

save(fit_shifted_log_flanker, file = "E:/20Hz/Data/Modelle/shifted_log_flanker.rda")

# posteriro predictive checks
pp_check(fit_shifted_log_flanker, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_shifted_log_flanker, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_shifted_log_flanker, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_shifted_log_flanker, ndraws = 1000, type = "stat", stat = "mean")
# looks really good
pp_check(fit_shifted_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_F")
# somewhat resonable
pp_check(fit_shifted_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")

# Great, the model looks good. Lets see the effect sizes, 

# first get variable names
get_variables(fit_shifted_log_flanker)

post_eff <- fit_shifted_log_flanker %>%
  spread_draws(b_Intercept, b_Contrast_FCongruency, b_Contrast_FStim_20v130, b_Contrast_FStim_20vOFF, 
               b_Contrast_FStroop_20v130, b_Contrast_FStroop_20vOFF, sigma, ndt) %>%
  mutate(Est_S130_I = exp(b_Intercept + (-0.5)*b_Contrast_FCongruency + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (-(1/3))*b_Contrast_FStroop_20v130 + (1/6)*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Est_S130_C = exp(b_Intercept + (0.5)*b_Contrast_FCongruency + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (1/3)*b_Contrast_FStroop_20v130 + (-(1/6))*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Est_S20_I = exp(b_Intercept + (-0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (1/6)*b_Contrast_FStroop_20v130 + (1/6)*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Est_S20_C = exp(b_Intercept + (0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (-(1/6))*b_Contrast_FStroop_20v130 + (-(1/6))*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Est_SOFF_I = exp(b_Intercept + (-0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (1/6)*b_Contrast_FStroop_20v130 + (-(1/3))*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Est_SOFF_C = exp(b_Intercept + (0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (-(1/6))*b_Contrast_FStroop_20v130 + (1/3)*b_Contrast_FStroop_20vOFF) + sigma/2 + ndt,
         Congruency = (Est_S20_I + Est_S130_I + Est_SOFF_I)/3 - (Est_S20_C + Est_S130_C + Est_SOFF_C)/3,
         S20_vs_S130 = (Est_S20_I + Est_S20_C)/2 - (Est_S130_I + Est_S130_C)/2,
         S20_vs_SOFF = (Est_S20_I + Est_S20_C)/2 - (Est_SOFF_I + Est_SOFF_C)/2,
         Stroop_S20_vs_S130 = (Est_S20_I - Est_S20_C) - (Est_S130_I - Est_S130_C),
         Stroop_S20_vs_SOFF =(Est_S20_I - Est_SOFF_C) - (Est_SOFF_I - Est_SOFF_C)) %>%
  select(Congruency, S20_vs_S130, S20_vs_SOFF, Stroop_S20_vs_S130, Stroop_S20_vs_SOFF) %>%
  summarise_draws()

# Okay, we see a lot of variance (makes sense since the timing measurements were imprecise)
# but it appears, there is no difference in RTs for any of the analyses. There was also no diffferences in Stroop effects.

## Next, let us look at the accuracy data

prior_weakly_informed_logreg<- c(
  prior(normal(-0.6, 0.6), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = Contrast_FCongruency), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20vOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)


#get_prior(formula = m1_SRT, data = RT_data, family = bernoulli(link = logit))

# brmsformula object Item Specific
m1_flanker_log <- bf(Error ~ 1  + Contrast_F + (Contrast_F|Part_nr)) 

#### Fit Inducer Models ####


# we should consider varying non-decision times between the groups
fit_log_flanker <- brm(formula = m1_flanker_log,
                      family = bernoulli(link = logit),
                      data = flankerRT,
                      prior = prior_weakly_informed_logreg,
                      warmup = 2000,
                      iter = 12000,# 20000 is the limit necessary for bridge sampling
                      cores = 4, seed = 423,
                      save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                      chains =4)

save(fit_log_flanker, file = "E:/20Hz/Data/Modelle/log_reg_flanker.rda")

# posteriro predictive checks
pp_check(fit_log_flanker, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_log_flanker, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_log_flanker, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_log_flanker, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_F")
# looks good
pp_check(fit_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")

# first get variable names
get_variables(fit_log_flanker)



post_eff_logreg <- fit_log_flanker %>%
  spread_draws(b_Intercept, b_Contrast_FCongruency, b_Contrast_FStim_20v130, b_Contrast_FStim_20vOFF, 
               b_Contrast_FStroop_20v130, b_Contrast_FStroop_20vOFF) %>%
  mutate(Est_S130_I = plogis(b_Intercept + (-0.5)*b_Contrast_FCongruency + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (-(1/3))*b_Contrast_FStroop_20v130 + (1/6)*b_Contrast_FStroop_20vOFF),
         Est_S130_C = plogis(b_Intercept + (0.5)*b_Contrast_FCongruency + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (1/3)*b_Contrast_FStroop_20v130 + (-(1/6))*b_Contrast_FStroop_20vOFF),
         Est_S20_I = plogis(b_Intercept + (-0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (1/6)*b_Contrast_FStroop_20v130 + (1/6)*b_Contrast_FStroop_20vOFF),
         Est_S20_C = plogis(b_Intercept + (0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (-(1/6))*b_Contrast_FStroop_20v130 + (-(1/6))*b_Contrast_FStroop_20vOFF),
         Est_SOFF_I = plogis(b_Intercept + (-0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (1/6)*b_Contrast_FStroop_20v130 + (-(1/3))*b_Contrast_FStroop_20vOFF),
         Est_SOFF_C = plogis(b_Intercept + (0.5)*b_Contrast_FCongruency + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (-(1/6))*b_Contrast_FStroop_20v130 + (1/3)*b_Contrast_FStroop_20vOFF),
         Congruency = (Est_S20_I + Est_S130_I + Est_SOFF_I)/3 - (Est_S20_C + Est_S130_C + Est_SOFF_C)/3,
         S20_vs_S130 = (Est_S20_I + Est_S20_C)/2 - (Est_S130_I + Est_S130_C)/2,
         S20_vs_SOFF = (Est_S20_I + Est_SOFF_C)/2 - (Est_SOFF_I + Est_SOFF_C)/2,
         Stroop_S20_vs_S130 = (Est_S20_I - Est_S20_C) - (Est_S130_I - Est_S130_C),
         Stroop_S20_vs_SOFF =(Est_S20_I - Est_S20_C) - (Est_SOFF_I - Est_SOFF_C)) %>%
  select(Congruency, S20_vs_S130, S20_vs_SOFF, Stroop_S20_vs_S130, Stroop_S20_vs_SOFF) %>%
  summarise_draws()
