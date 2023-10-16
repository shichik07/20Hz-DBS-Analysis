#####
# Author: Julius Kricheldorff
# Load models and calculate Bayes Factors for model components in the Go-NoGo task
# Date 05.10.2023
####
setwd('E:/20Hz/Data/Modelle/BF_mods')

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)
library(rstudioapi)
library(xtable)
library(stringr)


# Set a seed for sake of reproducibility
set.seed(32936)

load_part_mod <- function(loc, ana, model_t, param){
  # loads the partial model
  #
  # Args:
  #   loc: string of file location where the partial model is saved
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   model_t: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   The specified loaded model (brmsfit) with generically renamed as fit_model
  string <- paste(ana, model_t, param, sep = "_")
  fit <- load(file.path(loc, paste(string, ".rda", sep ="")))
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}

load_full_mod <- function(loc, ana, model_t){
  # loads the partial model
  #
  # Args:
  #   loc: string of file location where the full model is saved
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   model_t: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  # Returns:
  #   The specified loaded model (brmsfit) with generically renamed as fit_model

  if (model_t == "RT"){ 
    string <- file.path(loc, paste("shifted_log", ana, sep = "_"))
  } else if (model_t == "Acc"){
    string <- file.path(loc, paste("log_reg", ana, sep = "_"))
  }
  fit <- load(paste(string, ".rda", sep =""))
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}

BF_calc <- function(fullmod_loc, part_mod_loc, parameter, ana, model_t){
  # calculates the Bayes for all models
  #
  # Args:
  #   fullmod_loc: string of file location where the full model is saved
  #   part_mod_loc: string of file location where the partial model is saved
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   model_t: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: vector with parameter names as strings
  # Returns:
  #   a tibble containing the BFs and parameter names
  
  # tibble to save the data
  Model_BF <- tibble(
    Model = character(),
    Parameter = character(),
    BF = numeric(),
    Analysis = character()
  )
  
  #print(paste("Analyzing the", eff, "effect for the", itm, "items of the", model_type, "data"))
  # tibble to save the likelihood values
  Likelihoods <- tibble(
    Model = character(),
    Log_lik = numeric()
  )
  
  # partial models that we are going to load
  submodel <- paste(ana, "min", parameter, sep = "_")
  
  # load the full model
  full_mod <- load_full_mod(loc = fullmod_loc,
                            model_t = model_t,
                            ana = ana)
  
  # now we load the partial models
  for (par in submodel){
    # load the partial model
    def_mod <- load_part_mod(loc = part_mod_loc,
                             ana = ana, 
                             model_t = model_t, 
                             param = par)
    
    #calculate the BF using the brms function for it
    temp_BF <- bayes_factor(full_mod, def_mod)
    
    #save results as a tibble
    temp <- tibble(
      Model = ana,
      Parameter = str_sub(par, start =9),
      BF = temp_BF$bf,
      Analysis = model_t
    )
    # join tibble as row
    Model_BF <- Model_BF %>%
      bind_rows(temp)
  }
  return(Model_BF)
}

conditional_effect_calc_shift <- function(model){
  # calculate effect estimates in ms or percent
  # Args:
  #   model: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   summary(model_effects): Table with summarizing the posterior distribution of effect estimates, median, upper and lower 95% hpd interval boundaries
  
  
  # get posterior samples to calculate conditional effects
  # first calculate the estimated marginal means
  emm_GNG_RT <- emmeans(model, specs = ~ SOFF*Go_diff +  S130Hz*Go_diff, epred = TRUE)
  S130Hz_NoGo_Go <- c(0, 0, 0, 0, 0, 0, 1, 0)
  S130Hz_Go <- c(0, 0, 0, 0, 1, 0, 0, 0)
  SOFF_Go <- c(0, 1, 0, 0, 0, 0, 0, 0)
  SOFF_NoGo_Go <- c(0, 0, 0, 1, 0, 0, 0, 0)
  S20Hz_NoGo_Go <- c(0, 0, 1, 0, 0, 0, 0, 0)
  S20Hz_Go <- c(1, 0, 0, 0, 0, 0, 0, 0)
  
  GoDiff_20v130 <- (S20Hz_NoGo_Go - S20Hz_Go) - (S130Hz_NoGo_Go - S130Hz_Go)
  GoDiff_20vOFF <- (S20Hz_NoGo_Go - S20Hz_Go) - (SOFF_NoGo_Go - SOFF_Go)
  
  # Next we calculate the contrasts of interest from these marginal means
  model_effects <- contrast(emm_GNG_RT, method = list("GoDiff_20v130" = GoDiff_20v130,
                                             "GoDiff_20vOFF" = GoDiff_20vOFF)) 
  return(summary(model_effects))
}

conditional_effect_calc_acc <- function(model){
  # calculate effect estimates in ms or percent
  # Args:
  #   model: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   summary(model_effects): Table with summarizing the posterior distribution of effect estimates, median, upper and lower 95% hpd interval boundaries
  
  
  # get posterior samples to calculate conditional effects
  # first calculate the estimated marginal means
  emm_GNG_Acc <- emmeans(model, ~ Contrast_F, epred = TRUE)
  S130Hz_NoGo_Go <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  S130Hz_NoGo <- c(1, 0, 0, 0, 0, 0, 0, 0, 0)
  SOFF_NoGo <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
  SOFF_NoGo_Go <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
  S20Hz_NoGo_Go <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)
  S20Hz_NoGo <- c(0, 0, 0, 0, 0, 0, 1, 0, 0)
 
  
  # Next we calculate the contrasts of interest from these marginal means
  model_effects <- contrast(emm_GNG_Acc , method = list("NoGo_Stop_20Hz_vs_130Hz" = S20Hz_NoGo - S130Hz_NoGo,
                                                      "NoGo_Stop_20Hz_vs_OFF" = S20Hz_NoGo - SOFF_NoGo,
                                                      "NoGo_Go_20Hz_vs_130Hz" = S20Hz_NoGo_Go - S130Hz_NoGo_Go,
                                                      "NoGo_Go_20Hz_vs_OFF"= S20Hz_NoGo_Go - SOFF_NoGo_Go)) 
  return(summary(model_effects))
}

### First we get all the variables we need to call our readily calculated models

# contrast names seperately

fullmod_loc <- r"{E:\20Hz\Data\Modelle}"
part_mod_loc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
parameter_RT <- c("GoDiff_20v130",
               "GoDiff_20vOFF")

parameter_Acc <- c(
  "NoGo_Stop_20Hz_vs_130Hz",
  "NoGo_Stop_20Hz_vs_OFF",
  "NoGo_Go_20Hz_vs_130Hz",
  "NoGo_Go_20Hz_vs_OFF"
)

ana <- "GNG"
model_RT <- "RT"
model_Acc <- "Acc"

Partial_models_saveloc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
Full_models_saveloc <- r"{E:\20Hz\Data\Modelle}"

# calculate BFs for the RT data
BF_Results <- BF_calc(fullmod_loc = Full_models_saveloc, 
                      part_mod_loc = Partial_models_saveloc, 
                      parameter = parameter_RT, 
                      ana = ana, 
                      model_t = model_RT)

#save data

# Now the same for the Acc data 
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                      part_mod_loc = Partial_models_saveloc, 
                      parameter = parameter_Acc, 
                      ana = ana, 
                      model_t = model_Acc))



#save data
write.table(BF_Results , file = "E:/20Hz/Data/Modelle/BF_Results_GNG.csv")

# load data again
BF_Results <- read.csv(file = "E:/20Hz/Data/Modelle/BF_Results_GNG.csv", header = TRUE, sep = "")
#### Now that we have the table with our BFs, let us load the parameter estimates

# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(estimate = NA , lower.HPD = NA, upper.HPD = NA) 

# Loop through analysis and save effects
mods = c("RT", "Acc")

for (md in mods) {
  #load model
  full_model <- load_full_mod(loc = fullmod_loc, 
                              ana = ana, 
                              model_t = md)
  
  # calculate posterior conditional effect samples
  if (md == "RT"){
    eff_post <- conditional_effect_calc_shift(model = full_model)
  } else if (md == "Acc"){
    eff_post <- conditional_effect_calc_acc(model = full_model)
  }
  
  #integrate summary stats into BF table
  for (elem in eff_post$contrast){
    # find correct row
    idx <- Full_Model_Info$Analysis == md & Full_Model_Info$Parameter == elem
    Full_Model_Info$estimate[idx] <- eff_post$estimate[eff_post$contrast == elem]
    Full_Model_Info$lower.HPD[idx] <- eff_post$lower.HPD[eff_post$contrast == elem]
    Full_Model_Info$upper.HPD[idx] <- eff_post$upper.HPD[eff_post$contrast == elem]
  }
}

#save data
write.table(Full_Model_Info , file = "E:/20Hz/Data/Modelle/Full_Results_GNG.csv")
