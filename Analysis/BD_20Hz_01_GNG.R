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
m1_GoNoGo_shift <- bf(RT_ms ~ 1  + Go_diff + S130Hz + SOFF + S130Hz*Go_diff + SOFF*Go_diff + (1 + Go_diff + S130Hz + SOFF |Part_nr))

#get_prior(formula =  m1_GoNoGo_shift, data = RT_data, family = shifted_lognormal())

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

save(fit_shifted_log_GNG, file = "E:/20Hz/Data/Modelle/shifted_log_GNG.rda")
load(file = "E:/20Hz/Data/Modelle/shifted_log_GNG.rda")

# Perform posterior predictive checks 

# first calculate the posterior predictions
a <- brms::posterior_epred(fit_shifted_log_GNG)

# summarize over conditions
png(file="C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/RT_GNG_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 1, 3))

S130Hz_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "S130Hz_NoGo_Go"])
S130Hz_Go <- rowMeans(a[,RT_data$Contrast_F == "S130Hz_Go"])
pred_130 <- S130Hz_NoGo_Go - S130Hz_Go
mean_130 <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_Go"])
hist(pred_130, main = "Go(NoGo) - Go:  130Hz Stimulation", xlab = c("difference in ms"))
abline(v = mean_130, col = "red", lwd = 3)

SOFF_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "SOFF_NoGo_Go"])
SOFF_Go <- rowMeans(a[,RT_data$Contrast_F == "SOFF_Go"])
pred_OFF <- SOFF_NoGo_Go - SOFF_Go
mean_OFF <-  mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_Go"])
hist(pred_OFF, main = "Go(NoGo) - Go:  OFF Stimulation", xlab = c("difference in ms"))
abline(v = mean_OFF, col = "red", lwd = 3)

S20Hz_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "S20Hz_NoGo_Go"])
S20Hz_Go <- rowMeans(a[,RT_data$Contrast_F == "S20Hz_Go"])
pred_20 <- S20Hz_NoGo_Go - S20Hz_Go
mean_20 <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_Go"])
hist(pred_20, main = "Go(NoGo) - Go:  20Hz Stimulation", xlab = c("difference in ms"))
abline(v = mean_20, col = "red", lwd = 3)

# calculate marginal means
emm_GNG_RT <- emmeans(fit_shifted_log_GNG, specs = ~ SOFF*Go_diff +  S130Hz*Go_diff, epred = TRUE)
S130Hz_NoGo_Go <- c(0, 0, 0, 0, 0, 0, 1, 0)
S130Hz_Go <- c(0, 0, 0, 0, 1, 0, 0, 0)
SOFF_Go <- c(0, 1, 0, 0, 0, 0, 0, 0)
SOFF_NoGo_Go <- c(0, 0, 0, 1, 0, 0, 0, 0)
S20Hz_NoGo_Go <- c(0, 0, 1, 0, 0, 0, 0, 0)
S20Hz_Go <- c(1, 0, 0, 0, 0, 0, 0, 0)

GoDiff_S20_vs_S130 <- (S20Hz_NoGo_Go - S20Hz_Go) - (S130Hz_NoGo_Go - S130Hz_Go)
GoDiff_S20_vs_SOFF <- (S20Hz_NoGo_Go - S20Hz_Go) - (SOFF_NoGo_Go - SOFF_Go)

# Next we calculate the contrasts of interest from these marginal means
test <- contrast(emm_GNG_RT, method = list("S130Hz_NoGo_Go - S130Hz_Go" = S130Hz_NoGo_Go - S130Hz_Go,
                                            "S20Hz_NoGo_Go - S20Hz_Go" = S20Hz_NoGo_Go - S20Hz_Go,
                                            "SOFFHz_NoGo_Go - SOFFHz_Go" = SOFF_NoGo_Go - SOFF_Go,
                                            "GoDiff_S20_vs_S130" = GoDiff_S20_vs_S130,
                                            "GoDiff_S20_vs_SOFF" = GoDiff_S20_vs_SOFF)) 


# alright, we see participants definitely slowed more in the 20Hz condition compared to the OFF state, there is
# some indication that this may also be the case for the parameter estimate of the 20 vs. 130 Hz stim condition
# but the actual parameter does include the null. Surprisingly our mean estimate does not? not sure what that is about
# Maybe try the median? Regardless a BF should provide better interpretability. Note: I understand
## Next, let us look at the accuracy data


prior_weakly_informed_logreg2<- c(
  prior(normal(-2, 1), class = Intercept),
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


#Third Analysis - Here I formulate the contrasts which we are actually interested in - just 20Hz vs other for all response options

GoNoGo_Contrast_Acc <- hypr(
  Go_Stop_vs_Go = (S20Hz_NoGo_Go + S130Hz_NoGo_Go + SOFF_NoGo_Go)/3 ~ (S20Hz_Go + S130Hz_Go + SOFF_Go)/3,
  Go_Stop_vs_Stop = (S20Hz_NoGo_Go + S130Hz_NoGo_Go + SOFF_NoGo_Go)/3 ~ (S20Hz_NoGo_Stop + S130Hz_NoGo_Stop + SOFF_NoGo_Stop)/3,
  NoGo_Stop_20Hz_vs_130Hz = S20Hz_NoGo_Stop ~ S130Hz_NoGo_Stop,
  NoGo_Stop_20Hz_vs_OFF = S20Hz_NoGo_Stop ~ SOFF_NoGo_Stop,
  NoGo_Go_20Hz_vs_130Hz = S20Hz_NoGo_Go ~ S130Hz_NoGo_Go,
  NoGo_Go_20Hz_vs_OFF = S20Hz_NoGo_Go ~ SOFF_NoGo_Go,
  Go_20Hz_vs_130Hz = S20Hz_Go ~ S130Hz_Go,
  Go_20Hz_vs_OFF = S20Hz_Go ~ SOFF_Go,
  levels = c("S130Hz_NoGo_Stop", "S130Hz_NoGo_Go", "S130Hz_Go", 
             "SOFF_Go", "SOFF_NoGo_Stop", "SOFF_NoGo_Go", 
             "S20Hz_NoGo_Stop", "S20Hz_NoGo_Go", "S20Hz_Go")
)

GoNoGo_Contrast_Acc

contrasts(GoNoGo$Contrast_F) <- contr.hypothesis(GoNoGo_Contrast_Acc)
contrasts(GoNoGo$Contrast_F)

# brmsformula object Item Specific
m1_GoNoGo_log <- bf(Error ~ 1  + Contrast_F + (Contrast_F|Part_nr)) 
m1_GoNoGo_log <- bf(Error ~ 1  + Stim_verb * GoNoGo + (Stim_verb * GoNoGo|Part_nr)) 

#get_prior(formula =  m1_GoNoGo_log, data = GoNoGo, family = bernoulli(link = logit))

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



# we should consider varying non-decision times between the groups
fit_log_flanker2 <- brm(formula = m1_GoNoGo_log,
                       family = bernoulli(link = logit),
                       data = GoNoGo,
                       prior = prior_weakly_informed_logreg2,
                       warmup = 2000,
                       iter = 12000,# 20000 is the limit necessary for bridge sampling
                       cores = 4, seed = 423,
                       save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                       chains =4)

save(fit_log_flanker2, file = "E:/20Hz/Data/Modelle/log_reg_GNG.rda")
load(file = "E:/20Hz/Data/Modelle/log_reg_GNG.rda")

# posteriro predictive checks
pp_check(fit_log_flanker2, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_log_flanker2, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_log_flanker2, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Contrast_F")
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Stim_verb")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "GoNoGo")
# looks good
pp_check(fit_log_flanker2, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")



warp_em2 <- emmeans(fit_log_flanker2, ~ Stim_verb|GoNoGo, epred = TRUE)
warp_em3 <- emmeans(fit_log_flanker2, ~ Contrast_F, epred = TRUE)
cont2 <- contrast(regrid(warp_em3, transform = "response"), interaction = "pairwise")

# we need to rework the brms object to a grid object
prior_mod <- unupdate(fit_log_flanker2)
prior_emmgrid <- emmeans(prior_mod, ~ Contrast_F, epred = TRUE)

# then we can estimate the Bayesfactor
bayesfactor_parameters(cont2, prior = prior_emmgrid)
mod <- ref_grid(fit_log_flanker2)
prior_mod <- update(prior_weakly_informed_logreg2, prior_PD = TRUE)
bayesfactor_parameters(pairs(mod), prior =)



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