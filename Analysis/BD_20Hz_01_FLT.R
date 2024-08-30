#####
# Author: Julius Kricheldorff
# Analysis of the flanker task
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
wd <-"D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

# load data
FLTRT<- read_csv(file = "flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))


# create a contrast matrix for our comparisons of interest
FLT_Contrast <- hypr(
  Congruency = (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3 ~ (S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3,
  Stim_20v130 = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (S130Hz_incongruent + S130Hz_congruent)/2,
  Stim_20vOFF = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (SOFF_incongruent + SOFF_congruent)/2,
  Stroop_20v130 = (S20Hz_incongruent - S20Hz_congruent) ~ (S130Hz_incongruent - S130Hz_congruent),
  Stroop_20vOFF = (S20Hz_incongruent - S20Hz_congruent) ~ (SOFF_incongruent - SOFF_congruent),
  levels = c("S130Hz_incongruent", "S130Hz_congruent", "SOFF_congruent", 
             "SOFF_incongruent", "S20Hz_incongruent", "S20Hz_congruent")
)

FLT_Contrast

# create a variable for the contrasts

FLTRT <- FLTRT %>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & Congruency == "congruent" ~ "S130Hz_congruent",
    Stim_verb == "130Hz" & Congruency == "incongruent" ~ "S130Hz_incongruent",
    Stim_verb == "20Hz" & Congruency == "congruent" ~ "S20Hz_congruent",
    Stim_verb == "20Hz" & Congruency == "incongruent" ~ "S20Hz_incongruent",
    Stim_verb == "OFF" & Congruency == "congruent" ~ "SOFF_congruent",
    Stim_verb == "OFF" & Congruency == "incongruent" ~ "SOFF_incongruent",
  )))

# assign the generated contrast matrix to the List Wide Factor
contrasts(FLTRT$Contrast_F) <- contr.hypothesis(FLT_Contrast)
contrasts(FLTRT$Contrast_F)   


# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- FLTRT %>%
  filter(Correct_Response == 1,
         RT < 3,
         RT > 0.2) %>%
   mutate(RT_ms = RT*1000) 

RT_data2 <-  FLTRT %>%
  filter(Correct_Response == 1)%>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)

# Prior informed weakly 
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
m1_FLT <- bf(RT_ms ~ 1  + Contrast_F + (1 |Part_nr))

#fit the first model
fit_shifted_log_FLT <- brm(formula = m1_FLT,
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

save(fit_shifted_log_FLT, file = "../Modelle/shifted_log_FLT1.rda")
load(file = "../Modelle/shifted_log_FLT.rda")


## Next, let us look at the accuracy data

prior_weakly_informed_logreg<- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = Contrast_FCongruency), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20vOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)

# brmsformula object Item Specific
m1_FLT_log <- bf(Error ~ 1  + Contrast_F + (1|Part_nr)) 

#### Fit Accuracy Model ####


# we should consider varying non-decision times between the groups
fit_log_FLT <- brm(formula = m1_FLT_log,
                      family = bernoulli(link = logit),
                      data = FLTRT,
                      prior = prior_weakly_informed_logreg,
                      warmup = 2000,
                      iter = 12000,# 20000 is the limit necessary for bridge sampling
                      cores = 4, seed = 423,
                      save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                      chains =4)

save(fit_log_FLT, file = "../Modelle/log_reg_FLT.rda")
load(file = "../Modelle/log_reg_FLT.rda")

