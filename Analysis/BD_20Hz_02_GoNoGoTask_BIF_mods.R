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

# Set a seed for sake of reproducibility
set.seed(32936)

setwd('E:/20Hz/Data/Modelle/BF_mods')

# load data
GoNoGo<- read_csv(file = "C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted/GoNoGo.csv") %>%
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
  droplevels() # drop Stop trials
  

# Define formulas so we can loop through them
GNG_formulas <- c(
  formula((S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3), # main effect Go effects
  formula((S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (S130Hz_NoGo_Go + S130Hz_Go)/2), # Overall effect LFS vs HFS
  formula((S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (SOFF_NoGo_Go + SOFF_Go)/2), # Overall effect LFS vs OFF
  formula((S20Hz_NoGo_Go - S20Hz_Go) ~ (S130Hz_NoGo_Go - S130Hz_Go)), # Difference Go effects LFS vs HFS
  formula((S20Hz_NoGo_Go - S20Hz_Go) ~ (SOFF_NoGo_Go - SOFF_Go)) # Difference Go effects LFS vs OFF
)

GoNoGo_Contrast_RT <- hypr(
  Go_NoGo_Go = (S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3,
  Stim_20v130 = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (S130Hz_NoGo_Go + S130Hz_Go)/2,
  Stim_20vOFF = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (SOFF_NoGo_Go + SOFF_Go)/2,
  GoDiff_20v130 = (S20Hz_NoGo_Go - S20Hz_Go) ~ (S130Hz_NoGo_Go - S130Hz_Go),
  GoDiff_20vOFF = (S20Hz_NoGo_Go - S20Hz_Go) ~ (SOFF_NoGo_Go - SOFF_Go),
  levels = c("S130Hz_NoGo_Go", "S130Hz_Go", "SOFF_Go", 
             "SOFF_NoGo_Go", "S20Hz_NoGo_Go", "S20Hz_Go")
)

# contrast names separately
GNG_contrast_names <- c(
  "Go_NoGo_Go",
  "Stim_20v130",
  "Stim_20vOFF",
  "GoDiff_20v130",
  "GoDiff_20vOFF"
)

# GNG levels
GNG_levels <- c("S130Hz_NoGo_Go", "S130Hz_Go", "SOFF_Go", 
                "SOFF_NoGo_Go", "S20Hz_NoGo_Go", "S20Hz_Go")

GNG_mods <- c("GNG_min_Go_NoGo_Go", 
             "GNG_min_Stim_20v130",
             "GNG_min_Stim_20vOFF",
             "GNG_min_GoDiff_20v130",
             "GNG_min_GoDiff_20vOFF")

# First let us get the model contrasts
for(mods in 1:length(GNG_mods)){
  # first create the contrast matrix
  temp_mat <- hypr(GNG_formulas[-mods],
                   levels = GNG_levels)
  # next add the appropriate variable names
  names(temp_mat) <- GNG_contrast_names[-mods]
  
  # Lastly rename the contrast matrix
  assign(GNG_mods[mods], temp_mat)
}

# Same priors as for the model before
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

# okay create a function to pass model parameter
pass_brms = function(save_name, prior, data, model) {
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

##### now let us loop though our models and save the results 
for(mods in 1:length(GNG_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(RT_data$Contrast_F) <- contr.hypothesis(eval(parse(text = GNG_mods[mods])))
  # define prior
  Prior_weakly <- prior_weakly_informed[-(mods+3),]
  # get save name for variable
  save_name <- paste("GNG_RT_", GNG_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = Prior_weakly, data = RT_data, model = m1_GoNoGo_shift)
} 

#### Now Lets do the same for the accuracy data

# Define formulas so we can loop through them
GNG_acc_formulas <- c(
  formula((S20Hz_NoGo_Go + S130Hz_NoGo_Go + SOFF_NoGo_Go)/3 ~ (S20Hz_Go + S130Hz_Go + SOFF_Go)/3), # main effect Go effects
  formula((S20Hz_NoGo_Go + S130Hz_NoGo_Go + SOFF_NoGo_Go)/3 ~ (S20Hz_NoGo_Stop + S130Hz_NoGo_Stop + SOFF_NoGo_Stop)/3), # Overall effect LFS vs HFS
  formula(S20Hz_NoGo_Stop ~ S130Hz_NoGo_Stop), # Overall effect LFS vs OFF
  formula(S20Hz_NoGo_Stop ~ SOFF_NoGo_Stop), # Difference Go effects LFS vs HFS
  formula(S20Hz_NoGo_Go ~ S130Hz_NoGo_Go), # Difference Go effects LFS vs OFF
  formula(S20Hz_NoGo_Go ~ SOFF_NoGo_Go),
  formula(S20Hz_Go ~ S130Hz_Go),
  formula(S20Hz_Go ~ SOFF_Go)
)

# contrast names separately
GNG_acc_contrast_names <- c(
  "Go_Stop_vs_Go",
  "Go_Stop_vs_Stop",
  "NoGo_Stop_20Hz_vs_130Hz",
  "NoGo_Stop_20Hz_vs_OFF",
  "NoGo_Go_20Hz_vs_130Hz",
  "NoGo_Go_20Hz_vs_OFF",
  "Go_20Hz_vs_130Hz",
  "Go_20Hz_vs_OFF"
)


# GNG levels
GNG_acc_levels <- c("S130Hz_NoGo_Stop", "S130Hz_NoGo_Go", "S130Hz_Go", 
                "SOFF_Go", "SOFF_NoGo_Stop", "SOFF_NoGo_Go", 
                "S20Hz_NoGo_Stop", "S20Hz_NoGo_Go", "S20Hz_Go")

GNG_acc_mods <- c("GNG_min_NoGo_Stop_20Hz_vs_130Hz", 
              "GNG_min_NoGo_Stop_20Hz_vs_OFF",
              "GNG_min_NoGo_Go_20Hz_vs_130Hz",
              "GNG_min_NoGo_Go_20Hz_vs_OFF",
              "GNG_min_Go_20Hz_vs_130Hz",
              "GNG_min_Go_20Hz_vs_OFF")

# First let us get the model contrasts
for(mods in 1:length(GNG_acc_mods)){
  # first create the contrast matrix
  temp_mat <- hypr(GNG_acc_formulas[-(mods+2)],
                   levels = GNG_acc_levels)
  # next add the appropriate variable names
  names(temp_mat) <- GNG_acc_contrast_names[-(mods+2)]
  
  # Lastly rename the contrast matrix
  assign(GNG_acc_mods[mods], temp_mat)
}

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
m1_GoNoGo_log <- bf(Error ~ 1  + Contrast_F + (Contrast_F|Part_nr))

# okay create a function to pass model parameter
pass_brms_log = function(save_name, prior, data, model) {
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

##### now let us loop though our models and save the results 
for(mods in 1:length(GNG_acc_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(GoNoGo$Contrast_F) <- contr.hypothesis(eval(parse(text = GNG_acc_mods[mods])))
  # define prior
  Prior_weakly <- prior_weakly_informed_log[-(mods+3),]
  # get save name for variable
  save_name <- paste("GNG_acc_", GNG_acc_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms_log(save_name = save_name, prior = Prior_weakly, data = GoNoGo, model = m1_GoNoGo_log)
} 
