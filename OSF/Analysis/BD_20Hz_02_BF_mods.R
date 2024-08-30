#####
# Author: Julius Kricheldorff
# Fit Models for Bayes Factor calculation in the Go-NoGo task
# Date 26.09.2023
####

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)
library(rstudioapi)
library(readr)

# Set a seed for sake of reproducibility
set.seed(32936)


wd <-"D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis"
setwd(wd)

###### Functions:

pass_brms = function(save_name, prior, data, model) {
  # performs a shifted log regression
  #
  # Args:
  #   save_name: string of save location where the partial model is saved
  #   prior: list of of formulas for the prior distributions
  #   data: tibble with the data to be fit
  #   model: formula object with the model definition
  # Returns:
  #   The specified (brmsfit) 
  fit_model <- brm(formula = model,
                   family = shifted_lognormal(),
                   data = data,
                   prior = prior,
                   warmup = 2000,
                   iter = 12000,# 20000 is the limit necessary for bridge sampling
                   cores = 4, seed = 423,
                   control = list(adapt_delta = 0.95),
                   save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                   chains =4)
  
  # save the model
  save(fit_model, file = save_name)
}

pass_brms_log = function(save_name, prior, data, model) {
  # performs a logistic regression
  #
  # Args:
  #   save_name: string of save location where the partial model is saved
  #   prior: list of of formulas for the prior distributions
  #   data: tibble with the data to be fit
  #   model: formula object with the model definition
  # Returns:
  #   The specified (brmsfit) 
  fit_model <- brm(formula = model,
                   family = bernoulli(link = logit),
                   data = data,
                   prior = prior,
                   warmup = 2000,
                   iter = 12000,# 20000 is the limit necessary for bridge sampling
                   cores = 4, seed = 423,
                   control = list(adapt_delta = 0.95),
                   save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                   chains =4)
  
  # save the model
  save(fit_model, file = save_name)
}


# load data
GoNoGo<- read_csv(file = "Data/Extracted/GoNoGo.csv") %>%
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

# brmsformula object RT analysis
m1_GoNoGo_shift <- bf(RT_ms ~ 1  + Go_diff + S130Hz + SOFF + S130Hz*Go_diff + SOFF*Go_diff + (1|Part_nr))


# models we assess
GNG_mods <- c("GNG_min_GoDiff_20vOFF",
             "GNG_min_GoDiff_20v130")

# Same priors as for the model before
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

##### now let us loop though our models and save the results 
for(mods in 1:length(GNG_mods)){
  # define prior 
  Prior_weakly <- prior_weakly_informed
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+4)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/GNG_RT_", GNG_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = Prior_weakly, data = RT_data, model = m1_GoNoGo_shift)
} 

#### Now Lets do the same for the accuracy data

# Define formulas so we can loop through them
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

contrasts(GoNoGo$Contrast_F) <- contr.hypothesis(GoNoGo_Contrast_Acc)
contrasts(GoNoGo$Contrast_F)


GNG_acc_mods <- c("GNG_min_NoGo_Stop_20Hz_vs_130Hz", 
              "GNG_min_NoGo_Stop_20Hz_vs_OFF",
              "GNG_min_NoGo_Go_20Hz_vs_130Hz",
              "GNG_min_NoGo_Go_20Hz_vs_OFF",
              "GNG_min_Go_20Hz_vs_130Hz",
              "GNG_min_Go_20Hz_vs_OFF")

# Same priors as for the model before
prior_weakly_informed_log<- c(
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

# brmsformula object RT analysis
m1_GoNoGo_log <- bf(Error ~ 1  + Contrast_F + (1|Part_nr))

##### now let us loop though our models and save the results 
for(mods in 1:length(GNG_acc_mods)){
  # define prior
  Prior_weakly <- prior_weakly_informed_log
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+3)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/GNG_acc_", GNG_acc_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_log(save_name = save_name, prior = Prior_weakly, data = GoNoGo, model = m1_GoNoGo_log)
} 

##### Now the same for the Flanker RT data

# load data
FLTRT<- read_csv(file = "Data/Extracted/flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
  )))%>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & Congruency == "congruent" ~ "S130Hz_congruent",
    Stim_verb == "130Hz" & Congruency == "incongruent" ~ "S130Hz_incongruent",
    Stim_verb == "20Hz" & Congruency == "congruent" ~ "S20Hz_congruent",
    Stim_verb == "20Hz" & Congruency == "incongruent" ~ "S20Hz_incongruent",
    Stim_verb == "OFF" & Congruency == "congruent" ~ "SOFF_congruent",
    Stim_verb == "OFF" & Congruency == "incongruent" ~ "SOFF_incongruent",
  )))

RT_data_FLT <- FLTRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

# Define formulas so we can loop through them
FLT_Contrast <- hypr(
  Congruency = (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3 ~ (S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3,
  Stim_20v130 = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (S130Hz_incongruent + S130Hz_congruent)/2,
  Stim_20vOFF = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (SOFF_incongruent + SOFF_congruent)/2,
  Stroop_20v130 = (S20Hz_incongruent - S20Hz_congruent) ~ (S130Hz_incongruent - S130Hz_congruent),
  Stroop_20vOFF = (S20Hz_incongruent - S20Hz_congruent) ~ (SOFF_incongruent - SOFF_congruent),
  levels = c("S130Hz_incongruent", "S130Hz_congruent", "SOFF_congruent", 
             "SOFF_incongruent", "S20Hz_incongruent", "S20Hz_congruent")
)

# assign the generated contrast matrix to the List Wide Factor
contrasts(RT_data_FLT$Contrast_F) <- contr.hypothesis(FLT_Contrast)
contrasts(RT_data_FLT$Contrast_F) 

FLT_RT_mods <- c("FLT_min_Congruency", 
                  "FLT_min_Stim_20v130",
                  "FLT_min_Stim_20vOFF",
                  "FLT_min_Stroop_20v130",
                  "FLT_min_Stroop_20vOFF")

# Prior informed weakly Item Specific
prior_weakly_informed_FLT<- c(
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
m1_FLT <- bf(RT_ms ~ 1  + Contrast_F + (1|Part_nr))

##### now let us loop though our models and save the results 
for(mods in 1:length(FLT_RT_mods)){
  # define prior
  Prior_weakly <- prior_weakly_informed_FLT
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+3)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/FLT_RT_", FLT_RT_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = Prior_weakly, data = RT_data_FLT, model = m1_FLT)
} 


##### Now the same for the Flanker Error data


# Define formulas so we can loop through them
FLT_Contrast <- hypr(
  Congruency = (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3 ~ (S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3,
  Stim_20v130 = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (S130Hz_incongruent + S130Hz_congruent)/2,
  Stim_20vOFF = (S20Hz_incongruent + S20Hz_congruent)/2 ~ (SOFF_incongruent + SOFF_congruent)/2,
  Stroop_20v130 = (S20Hz_incongruent - S20Hz_congruent) ~ (S130Hz_incongruent - S130Hz_congruent),
  Stroop_20vOFF = (S20Hz_incongruent - S20Hz_congruent) ~ (SOFF_incongruent - SOFF_congruent),
  levels = c("S130Hz_incongruent", "S130Hz_congruent", "SOFF_congruent", 
             "SOFF_incongruent", "S20Hz_incongruent", "S20Hz_congruent")
)

# assign the generated contrast matrix to the List Wide Factor
contrasts(FLTRT$Contrast_F) <- contr.hypothesis(FLT_Contrast)
contrasts(FLTRT$Contrast_F) 


FLT_acc_mods <- c("FLT_min_Congruency", 
                 "FLT_min_Stim_20v130",
                 "FLT_min_Stim_20vOFF",
                 "FLT_min_Stroop_20v130",
                 "FLT_min_Stroop_20vOFF")

# Prior informed weakly Item Specific
prior_weakly_informed_FLT_acc<- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = Contrast_FCongruency), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStim_20vOFF), 
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20v130),
  prior(normal(0,  1.5), class = b, coef = Contrast_FStroop_20vOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)

# brmsformula object List Wide
m1_FLT <- bf(Error ~ 1  + Contrast_F + (1|Part_nr))

##### now let us loop though our models and save the results 
for(mods in 1:length(FLT_acc_mods)){
  # define prior
  Prior_weakly <- prior_weakly_informed_FLT_acc
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+1)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/FLT_acc_", FLT_acc_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_log(save_name = save_name, prior = Prior_weakly, data = FLTRT, model = m1_FLT)
} 

##### And lastly lets do it for the SRT data again

# load data
SimpleRT<- read_csv(file = "Data/Extracted/SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data_SRT <- SimpleRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

# create a contrast matrix for our comparisons of interest
SRT_Contrast <- hypr(
  S130_S20 = S20Hz ~ S130Hz, 
  S20_SOFF = S20Hz ~ SOFF, 
  levels = c("S130Hz", "SOFF", "S20Hz")
)
SRT_Contrast

contrasts(RT_data_SRT$StimCon) <- contr.hypothesis(SRT_Contrast)
contrasts(RT_data_SRT$StimCon)

# Model names
SRT_RT_mods <- c("SRT_min_Stim_20v130",
                 "SRT_min_Stim_20vOFF"
                 )


# Prior informed weakly Item Specific
prior_weakly_informed_SRT<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = StimConS130_S20), 
  prior(normal(0,  0.3), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_SRT <- bf(RT_ms ~ 1  + StimCon + (1|Part_nr))
#get_prior(formula = m1_SRT, data = RT_data_SRT, family = shifted_lognormal())

##### now let us loop though our models and save the results 
for(mods in 1:length(SRT_RT_mods)){
  # define prior
  Prior_weakly <- prior_weakly_informed_SRT
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+3)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/SRT_RT_", SRT_RT_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = Prior_weakly, data = RT_data_SRT, model = m1_SRT)
} 

##### Now the same for the SRT Error data

# Define formulas so we can loop through them
SRT_Contrast <- hypr(
  S130_S20 = S20Hz ~ S130Hz, 
  S20_SOFF = S20Hz ~ SOFF, 
  levels = c("S130Hz", "SOFF", "S20Hz")
)
SRT_Contrast

contrasts(SimpleRT$StimCon) <- contr.hypothesis(SRT_Contrast)
contrasts(SimpleRT$StimCon)




SRT_acc_mods <- c("SRT_min_Stim_20v130",
                  "SRT_min_Stim_20vOFF"
)

# Prior informed weakly Item Specific
prior_weakly_informed_SRT_acc<- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = StimConS130_S20), 
  prior(normal(0,  1.5), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_SRT <- bf(Error ~ 1  + StimCon + (1|Part_nr))

##### now let us loop though our models and save the results 
for(mods in 1:length(SRT_acc_mods)){
  # define prior
  Prior_weakly <- prior_weakly_informed_SRT_acc
  # set the slope of our variable of interest to zero
  Prior_weakly$prior[(mods+1)] <- "constant(0)"
  # get save name for variable
  save_name <- paste("Data/Modelle/BF_mods/SRT_acc_", SRT_acc_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_log(save_name = save_name, prior = Prior_weakly, data = SimpleRT, model = m1_SRT)
} 

