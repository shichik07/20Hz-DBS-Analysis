#####
# Author: Julius Kricheldorff
# Analysis of the Go-NoGo task
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
library(tidybayes)
library(emmeans)

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted')

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
  mutate(RT_ms = RT*1000)

RT_data2 <- GoNoGo %>%
  filter(Correct_Response == 1,
        
         GoNoGo == "NoGo - Go" | GoNoGo == "Go") %>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)
# create a contrast matrix for our comparisons of interest we are interested 
# if we can observe more slowing on NoGo - Go trials in the 20Hz condition, relatively
# to the Go condition in the 20Hz condition, than the 130Hz and Off condition

GoNoGo_Contrast_RT <- hypr(
  Go_NoGo_Go = (S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3,
  Stim_20v130 = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (S130Hz_NoGo_Go + S130Hz_Go)/2,
  Stim_20vOFF = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (SOFF_NoGo_Go + SOFF_Go)/2,
  GoDiff_20v130 = (S20Hz_NoGo_Go - S20Hz_Go) ~ (S130Hz_NoGo_Go - S130Hz_Go),
  GoDiff_20vOFF = (S20Hz_NoGo_Go - S20Hz_Go) ~ (SOFF_NoGo_Go - SOFF_Go),
  #levels = c("S130Hz_NoGo_Go", "S130Hz_Go", "SOFF_Go", 
            # "SOFF_NoGo_Go", "S20Hz_NoGo_Go", "S20Hz_Go")
  levels = c("S130Hz_NoGo_Stop", "S130Hz_NoGo_Go", "S130Hz_Go",
             "SOFF_Go", "SOFF_NoGo_Stop", "SOFF_NoGo_Go", 
             "S20Hz_NoGo_Stop", "S20Hz_NoGo_Go", "S20Hz_Go")
)

GoNoGo_Contrast_RT

# assign the generated contrast matrix to the List Wide Factor
contrasts(RT_data$Contrast_F) <- contr.hypothesis(GoNoGo_Contrast_RT)
contrasts(RT_data$Contrast_F)   


# Prior informed weakly Item Specific. 
# Note to self for the ndt prior: If I do not model an effect on the ndt parameter it is
# specified in milliseconds (so Y_min). However, once I decide to put in a group 
# intercept for example, it will be modeled using a log link function. Thus,
# for the prior specification. Specify the prior on the actual scale if you
# are not interested in any group effects on this parameter, but specify it on the
# log scale if you want to parametrize ndt further.
prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = Contrast_FGo_NoGo_Go), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FGoDiff_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FGoDiff_20vOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)

prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = Contrast_FGo_NoGo_Go), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  0.3), class = b, coef = Contrast_FGoDiff_20v130),
  prior(normal(0,  0.3), class = b, coef = Contrast_FGoDiff_20vOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object RT analysis
m1_GoNoGo_shift <- bf(RT_ms ~ 1  + Contrast_F + (Contrast_F|Part_nr))

#get_prior(formula =  m1_GoNoGo_shift, data = RT_data, family = shifted_lognormal())

#fit the first model
fit_shifted_log_GoNoGo <- brm(formula = m1_GoNoGo_shift,
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

save(fit_shifted_log_GoNoGo, file = "E:/20Hz/Data/Modelle/shifted_log_GoNoGo.rda")
load(file = "E:/20Hz/Data/Modelle/shifted_log_GoNoGo.rda")

# posteriro predictive checks
pp_check(fit_shifted_log_GoNoGo, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_shifted_log_GoNoGo, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_shifted_log_GoNoGo, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_shifted_log_GoNoGo, ndraws = 1000, type = "stat", stat = "mean")
# looks really good
pp_check(fit_shifted_log_GoNoGo, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_F")
# somewhat resonable
pp_check(fit_shifted_log_GoNoGo, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")

# Great, the model looks good. Lets see the effect sizes, 

# first get variable names
get_variables(fit_shifted_log_GoNoGo)

post_eff <- fit_shifted_log_GoNoGo %>%
  spread_draws(b_Intercept, b_Contrast_FGo_NoGo_Go, b_Contrast_FStim_20v130, b_Contrast_FStim_20vOFF, 
               b_Contrast_FGoDiff_20v130, b_Contrast_FGoDiff_20vOFF, sigma, ndt) %>%
  mutate(Est_S130_NoGo = exp(b_Intercept + (-0.5)*b_Contrast_FGo_NoGo_Go + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (-(1/3))*b_Contrast_FGoDiff_20v130 + (1/6)*b_Contrast_FGoDiff_20vOFF) + sigma/2 + ndt,
         Est_S130_Go = exp(b_Intercept + (0.5)*b_Contrast_FGo_NoGo_Go + (-(2/3))*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                          + (1/3)*b_Contrast_FGoDiff_20v130 + (-(1/6))*b_Contrast_FGoDiff_20vOFF) + sigma/2 + ndt,
         Est_S20_NoGo = exp(b_Intercept + (-0.5)*b_Contrast_FGo_NoGo_Go + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (1/6)*b_Contrast_FGoDiff_20v130 + (1/6)*b_Contrast_FGoDiff_20vOFF)+ sigma/2 + ndt,
         Est_S20_Go = exp(b_Intercept + (0.5)*b_Contrast_FGo_NoGo_Go + (1/3)*b_Contrast_FStim_20v130 + (1/3)*b_Contrast_FStim_20vOFF
                         + (-(1/6))*b_Contrast_FGoDiff_20v130 + (-(1/6))*b_Contrast_FGoDiff_20vOFF) + sigma/2 + ndt,
         Est_SOFF_NoGo = exp(b_Intercept + (-0.5)*b_Contrast_FGo_NoGo_Go + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (1/6)*b_Contrast_FGoDiff_20v130 + (-(1/3))*b_Contrast_FGoDiff_20vOFF) + sigma/2 + ndt,
         Est_SOFF_Go = exp(b_Intercept + (0.5)*b_Contrast_FGo_NoGo_Go + (1/3)*b_Contrast_FStim_20v130 + (-(2/3))*b_Contrast_FStim_20vOFF
                          + (-(1/6))*b_Contrast_FGoDiff_20v130 + (1/3)*b_Contrast_FGoDiff_20vOFF) + sigma/2 + ndt,
         Go_NoGo_Go = (Est_S20_NoGo + Est_S130_NoGo + Est_SOFF_NoGo)/3 - (Est_S20_Go + Est_S130_Go + Est_SOFF_Go)/3,
         S20_vs_S130 = (Est_S20_NoGo + Est_S20_Go)/2 - (Est_S130_NoGo + Est_S130_Go)/2,
         S20_vs_SOFF = (Est_S20_NoGo + Est_SOFF_Go)/2 - (Est_SOFF_NoGo + Est_SOFF_Go)/2,
         GoDiff_S20_vs_S130 = (Est_S20_NoGo - Est_S20_Go) - (Est_S130_NoGo - Est_S130_Go),
         GoDiff_S20_vs_SOFF =(Est_S20_NoGo - Est_S20_Go) - (Est_SOFF_NoGo - Est_SOFF_Go)) %>%
  #select(Go_NoGo_Go, S20_vs_S130, S20_vs_SOFF, GoDiff_S20_vs_S130, GoDiff_S20_vs_SOFF) %>%
  summarise_draws(mean, median, sd, ~quantile(.x, probs = c(0.025, 0.975)))


# first calculate the estimated marginal means
emm_GNG_RT <- emmeans(fit_shifted_log_GoNoGo, specs = ~ Contrast_F, epred = TRUE)

# define new contrasts for emmeans - there should be an easier way to this
S130Hz_NoGo_Go <- c(1, 0, 0, 0, 0, 0)
S130Hz_Go <- c(0, 1, 0, 0, 0, 0)
SOFF_Go <- c(0, 0, 1, 0, 0, 0)
SOFF_NoGo_Go <- c(0, 0, 0, 1, 0, 0)
S20Hz_NoGo_Go <- c(0, 0, 0, 0, 1, 0)
S20Hz_Go <- c(0, 0, 0, 0, 0, 1)

# and contrasts for the differences in Go effects
GoDiff_S20_vs_S130 <- (S20Hz_NoGo_Go - S20Hz_Go) - (S130Hz_NoGo_Go - S130Hz_Go)
GoDiff_S20_vs_SOFF <- (S20Hz_NoGo_Go - S20Hz_Go) - (SOFF_NoGo_Go - SOFF_Go)

# Next we calculate the contrasts of interest from these marginal means
contrast(emm_GNG_RT, method = list("S130Hz_NoGo_Go - S130Hz_Go" = S130Hz_NoGo_Go - S130Hz_Go,
                                   "S20Hz_NoGo_Go - S20Hz_Go" = S20Hz_NoGo_Go - S20Hz_Go,
                                   "SOFFHz_NoGo_Go - SOFFHz_Go" = SOFF_NoGo_Go - SOFF_Go,
                                   "GoDiff_S20_vs_S130" = GoDiff_S20_vs_S130,
                                   "GoDiff_S20_vs_SOFF" = GoDiff_S20_vs_SOFF)) 


# alright, we see participants definitely slowed more in the 20Hz condition compared to the OFF state, there is
# some indication that this may also be the case for the parameter estimate of the 20 vs. 130 Hz stim condition
# but the actual parameter does include the null. Surprisingly our mean estimate does not? not sure what that is about
# Maybe try the median? Regardless a BF should provide better interpretability. Note: I understand
## Next, let us look at the accuracy data

#### Analysis of the accuracy data
# Things are a little trickier here. Now in addition we are interested in the error rates on stop trials.
# GoNoGo_Contrast_Acc <- hypr(
#   Go_NoGo_Go = (S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3,
#   NoGo_Stop_NoGo_Go = (S130Hz_NoGo_Stop + SOFF_NoGo_Stop + S20Hz_NoGo_Stop)/3 ~ 
#     (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3,
#   Stim_20v130 = (S20Hz_NoGo_Go + S20Hz_Go +  S20Hz_NoGo_Stop)/3 ~ 
#     (S130Hz_NoGo_Go + S130Hz_Go + S130Hz_NoGo_Stop)/3,
#   Stim_20vOFF = (S20Hz_NoGo_Go + S20Hz_Go +  S20Hz_NoGo_Stop)/3 ~ 
#     (SOFF_NoGo_Go + SOFF_Go +  SOFF_NoGo_Stop)/3,
#   NoGo_Go_Diff_20v130 = S20Hz_NoGo_Go ~ S130Hz_NoGo_Go,
#   NoGo_Stopp_Diff_20v130 = S20Hz_NoGo_Stop ~ S130Hz_NoGo_Stop,
#   NoGo_Go_Diff_20vOFF = S20Hz_NoGo_Go ~ SOFF_NoGo_Go,
#   NoGo_Stopp_Diff_20vOFF = S20Hz_NoGo_Stop ~ SOFF_NoGo_Stop,
#   levels = c("S130Hz_NoGo_Go", "S130Hz_Go", "SOFF_Go", 
#              "SOFF_NoGo_Go", "S20Hz_NoGo_Go", "S20Hz_Go",
#              "S130Hz_NoGo_Stop", "S20Hz_NoGo_Stop", "SOFF_NoGo_Stop")
# )
# 
# GoNoGo_Contrast_Acc
# 
# 
# # assign the generated contrast matrix to the List Wide Factor
# contrasts(GoNoGo$Contrast_F) <- contr.hypothesis(GoNoGo_Contrast_Acc)
# contrasts(GoNoGo$Contrast_F)   
# 
# 
# prior_weakly_informed_logreg<- c(
#   prior(normal(-0.6, 0.6), class = Intercept),
#   prior(normal(0,  1.5), class = b, coef = Contrast_FGo_NoGo_Go), 
#   prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Go_Diff_20v130),
#   prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Go_Diff_20vOFF), 
#   prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Stop_NoGo_Go),
#   prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Stopp_Diff_20v130),
#   prior(normal(0,  1.5), class = b, coef = Contrast_FNoGo_Stopp_Diff_20vOFF), 
#   prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20v130),
#   prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20vOFF),
#   prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
# )
# 
# 
# #get_prior(formula = m1_GoNoGo_log, data = GoNoGo, family = bernoulli(link = logit))
# 
# # brmsformula object Item Specific
# m1_GoNoGo_log <- bf(Error ~ 1  + Contrast_F + (Contrast_F|Part_nr)) 
# 
# #### Fit Inducer Models ####
# 
# 
# # we should consider varying non-decision times between the groups
# fit_log_flanker <- brm(formula = m1_GoNoGo_log,
#                        family = bernoulli(link = logit),
#                        data = GoNoGo,
#                        prior = prior_weakly_informed_logreg,
#                        warmup = 2000,
#                        iter = 12000,# 20000 is the limit necessary for bridge sampling
#                        cores = 4, seed = 423,
#                        save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
#                        chains =4)
# 
# save(fit_log_flanker, file = "E:/20Hz/Data/Modelle/log_reg_GoNoGo.rda")
# load(file = "E:/20Hz/Data/Modelle/log_reg_GoNoGo.rda")
# 
# warp_em <- emmeans(fit_log_flanker, ~ Contrast_F, type = "response")
# cont <- contrast(warp_em, type = "response")
# #get the posterior draws from the contrasts 
# cont_posterior <- gather_emmeans_draws(cont)
# 
# # posteriro predictive checks
# pp_check(fit_log_flanker, ndraws = 11, type = "hist")
# # Looks good for this model
# pp_check(fit_log_flanker, ndraws = 100, type = "dens_overlay")
# # Looks good for this model
# pp_check(fit_log_flanker, type = "boxplot", ndraws = 10)
# # A few observation outside, but our particpants had a deadline - cannot be modelled
# pp_check(fit_log_flanker, ndraws = 1000, type = "stat", stat = "mean")
# # looks good
# pp_check(fit_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_F")
# # looks good
# pp_check(fit_log_flanker, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")

#NOTE: The parametrization that I chose worked terribly. For the final model
# I will use the generic parametrization and use the emmeans for the predicted contrasts of interest.

## Second Analysis with the default parametrization

#get_prior(formula = m1_GoNoGo_log2, data = GoNoGo, family = bernoulli(link = logit))

prior_weakly_informed_logreg2<- c(
  prior(normal(-0.6, 0.6), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = GoNoGoNoGoMGo), 
  prior(normal(0,  1.5), class = b, coef = GoNoGoNoGoMStop),
  prior(normal(0,  1.5), class = b, coef = Stim_verb20Hz), 
  prior(normal(0,  1.5), class = b, coef = Stim_verb20Hz:GoNoGoNoGoMGo),
  prior(normal(0,  1.5), class = b, coef = Stim_verb20Hz:GoNoGoNoGoMStop),
  prior(normal(0,  1.5), class = b, coef = Stim_verbOFF), 
  prior(normal(0,  1.5), class = b, coef = Stim_verbOFF:GoNoGoNoGoMGo),
  prior(normal(0,  1.5), class = b, coef = Stim_verbOFF:GoNoGoNoGoMStop),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)

# brmsformula object Item Specific
m1_GoNoGo_log2 <- bf(Error ~ 1  + Stim_verb*GoNoGo + (Stim_verb*GoNoGo|Part_nr)) 


# we should consider varying non-decision times between the groups
fit_log_flanker2 <- brm(formula = m1_GoNoGo_log2,
                       family = bernoulli(link = logit),
                       data = GoNoGo,
                       prior = prior_weakly_informed_logreg2,
                       warmup = 2000,
                       iter = 12000,# 20000 is the limit necessary for bridge sampling
                       cores = 4, seed = 423,
                       save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                       chains =4)

save(fit_log_flanker2, file = "E:/20Hz/Data/Modelle/log_reg_GoNoGo2.rda")
load(file = "E:/20Hz/Data/Modelle/log_reg_GoNoGo2.rda")

# posteriro predictive checks
pp_check(fit_log_flanker2, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_log_flanker2, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_log_flanker2, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Stim_verb")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "GoNoGo")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")



warp_em2 <- emmeans(fit_log_flanker2, ~ Stim_verb|GoNoGo, epred = TRUE)
cont2 <- contrast(regrid(warp_em2, transform = "response"), interaction = "pairwise")

# we need to rework the brms object to a grid object
prior_mod <- unupdate(fit_log_flanker2)
prior_emmgrid <- emmeans(prior_mod, ~ Stim_verb|GoNoGo)

# then we can estimate the Bayesfactor
bayesfactor_parameters(cont2, prior = prior_emmgrid)



sum_acc <- summary(cont2, type = "response", point.est = mean)
#summary(cont2, type = "response", point.est = median)
#pairs(warp_em2, simple = "Stim_verb", point.est = mean)


## Look at the contrast matrix for the drift diffusion model

GoNoGo_Contrast_DDM <- hypr(
  DT_130Hz = (S130Hz_Go) ~ (S130Hz_NoGo_Go + S130Hz_NoGo_NoGo)/2,
  Stim_20v130 = (S20Hz_Go) ~ (S20Hz_NoGo_Go + S20Hz_NoGo_NoGo)/2,
  Stim_20vOFF = (SOFF_Go) ~ (SOFF_NoGo_Go + SOFF_NoGo_NoGo)/2,
  levels = c("S130Hz_NoGo_Go", "S130Hz_Go", "S130Hz_NoGo_NoGo", 
             "S20Hz_NoGo_Go", "S20Hz_Go", "S20Hz_NoGo_NoGo",
             "SOFF_NoGo_Go", "SOFF_Go", "SOFF_NoGo_NoGo")
)

GoNoGo_Contrast_DDM