#####
# Author: Julius Kricheldorff
# Load models and calculate Bayes Factors for model components in the Go-NoGo task
# Date 05.10.2023
####
#wd ='E:/20Hz/Data/Modelle/BF_mods'
#wd = 'C:/Users/doex9445/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_mods'
wd = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_mods"
setwd(wd)
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
library(emmeans)
library(tidybayes)
library(readr)


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
  fit <- try(load(file.path(loc, paste(string, ".rda", sep =""))))
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
  fit <- try(load(paste(string, ".rda", sep ="")))
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
    temp_BF <- try(brms::bayes_factor(full_mod, def_mod))
    
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

conditional_effect_calc_shift_GNG <- function(model){
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

brms_epred_calc_GNG_RT <- function(model){
  preds <- brms::posterior_epred(model)
  # load FLT data
  #GNG<- read_csv(file = "C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted/GoNoGo.csv") %>%
  GNG<- read_csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted/GoNoGo.csv") %>%
    mutate(Error = 1 - Correct_Response) %>%
    mutate(StimCon = as_factor(case_when(
      Stim_verb == "130Hz" ~ "S130Hz",
      Stim_verb == "20Hz" ~ "S20Hz",
      Stim_verb == "OFF" ~ "SOFF"
    ))) %>%
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
    ))) %>%
    filter(Correct_Response == 1,
           RT <3,
           RT > 0.2) %>%
    mutate(RT_ms = RT*1000)%>%
    droplevels() # drop Stop trials

  # calculate contrasts
  S130Hz_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "S130Hz_NoGo_Go"])
  S130Hz_Go <- rowMeans(preds[,GNG$Contrast_F == "S130Hz_Go"])
  SOFF_Go <- rowMeans(preds[,GNG$Contrast_F == "SOFF_Go"])
  SOFF_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "SOFF_NoGo_Go"])
  S20Hz_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "S20Hz_NoGo_Go"])
  S20Hz_Go <- rowMeans(preds[,GNG$Contrast_F == "S20Hz_Go"])
  
  # calculate and save contrasts
  model_effects <- tibble(contrast = c("GoDiff_20v130", "GoDiff_20vOFF"))
  GoDiff_20v130 <- tidybayes::mean_hdi(((S20Hz_NoGo_Go - S20Hz_Go) - (S130Hz_NoGo_Go - S130Hz_Go)))
  GoDiff_20vOFF <- tidybayes::mean_hdi(((S20Hz_NoGo_Go - S20Hz_Go) - (SOFF_NoGo_Go - SOFF_Go)))
  
  model_effects<- model_effects %>%
    bind_cols(bind_rows(GoDiff_20v130, GoDiff_20vOFF))%>%
    rename(estimate = y, lower.HPD = ymin, upper.HPD = ymax)
  
  return(model_effects)
}

brms_epred_calc_GNG_acc <- function(model){
  preds <- brms::posterior_epred(model)
  # load GNG data
  #GNG<- read_csv(file = "C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted/GoNoGo.csv") %>%
  GNG<- read_csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted/GoNoGo.csv") %>%
    mutate(Error = 1 - Correct_Response) %>%
    mutate(StimCon = as_factor(case_when(
      Stim_verb == "130Hz" ~ "S130Hz",
      Stim_verb == "20Hz" ~ "S20Hz",
      Stim_verb == "OFF" ~ "SOFF"
    ))) %>%
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
  
  # calculate contrasts
  S130Hz_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "S130Hz_NoGo_Go"])
  S130Hz_NoGo <- rowMeans(preds[,GNG$Contrast_F == "S130Hz_NoGo_Stop"])
  SOFF_NoGo <- rowMeans(preds[,GNG$Contrast_F == "SOFF_NoGo_Stop"])
  SOFF_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "SOFF_NoGo_Go"])
  S20Hz_NoGo_Go <- rowMeans(preds[,GNG$Contrast_F == "S20Hz_NoGo_Go"])
  S20Hz_NoGo <- rowMeans(preds[,GNG$Contrast_F == "S20Hz_NoGo_Stop"])
  
  # calculate and save contrasts
  model_effects <- tibble(contrast = c("NoGo_Stop_20Hz_vs_130Hz", "NoGo_Stop_20Hz_vs_OFF", "NoGo_Go_20Hz_vs_130Hz", "NoGo_Go_20Hz_vs_OFF"))
  NoGo_Stop_20Hz_vs_130Hz <- tidybayes::mean_hdi((S20Hz_NoGo - S130Hz_NoGo))
  NoGo_Stop_20Hz_vs_OFF<- tidybayes::mean_hdi((S20Hz_NoGo - SOFF_NoGo))
  NoGo_Go_20Hz_vs_130Hz <- tidybayes::mean_hdi((S20Hz_NoGo_Go - S130Hz_NoGo_Go))
  NoGo_Go_20Hz_vs_OFF<- tidybayes::mean_hdi((S20Hz_NoGo_Go - SOFF_NoGo_Go))
  
  
  model_effects<- model_effects %>%
    bind_cols(bind_rows(NoGo_Stop_20Hz_vs_130Hz, NoGo_Stop_20Hz_vs_OFF, NoGo_Go_20Hz_vs_130Hz, NoGo_Go_20Hz_vs_OFF)) %>%
    rename(estimate = y, lower.HPD = ymin, upper.HPD = ymax)
  
  return(model_effects)
}


conditional_effect_calc_acc_GNG <- function(model){
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

brms_epred_calc_SRT <- function(model, ana){
  preds <- brms::posterior_epred(model)
  # load FLT data
  #SimpleRT<- read_csv(file = "C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted/SimpleRT.csv") %>%
  SimpleRT<- read_csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted/SimpleRT.csv") %>%
    
    mutate(Error = 1 - Correct_Response) %>%
    mutate(StimCon = as_factor(case_when(
      Stim_verb == "130Hz" ~ "S130Hz",
      Stim_verb == "20Hz" ~ "S20Hz",
      Stim_verb == "OFF" ~ "SOFF"
    ))) 
  
  # filter the data more for the RT analysis
  if (ana == "RT") {SimpleRT  <- SimpleRT %>%
      filter(Correct_Response == 1,
             RT <3,
             RT > 0.2) %>%
      mutate(RT_ms = RT*1000)}
  
  # calculate contrasts
  S130Hz <- rowMeans(preds[,SimpleRT$StimCon == "S130Hz"])
  S20Hz <- rowMeans(preds[,SimpleRT$StimCon == "S20Hz"])
  SOFF <- rowMeans(preds[,SimpleRT$StimCon == "SOFF"])

  # calculate and save contrasts
  model_effects <- tibble(contrast = c("Stim_20v130", "Stim_20vOFF"))
  Stim_20v130 <- tidybayes::mean_hdi((S20Hz - S130Hz))
  Stim_20vOFF <- tidybayes::mean_hdi((S20Hz - SOFF))
  
  model_effects<- model_effects %>%
    bind_cols(bind_rows(Stim_20v130, Stim_20vOFF))%>%
    rename(estimate = y, lower.HPD = ymin, upper.HPD = ymax)
  
  return(model_effects)
}

conditional_effect_calc_SRT <- function(model){
  # calculate effect estimates in ms or percent
  # Args:
  #   model: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   summary(model_effects): Table with summarizing the posterior distribution of effect estimates, median, upper and lower 95% hpd interval boundaries
  
  
  # get posterior samples to calculate conditional effects
  # first calculate the estimated marginal means
  emm_SRT <- emmeans(model, ~ StimCon, epred = TRUE)
  S130Hz <- c(1, 0, 0)
  SOFF <- c(0, 1, 0)
  S20Hz <- c(0, 0, 1)
  
  
  # Next we calculate the contrasts of interest from these marginal means
  model_effects <- contrast(emm_SRT , method = list("Stim_20v130" = S20Hz - S130Hz,
                                                    "Stim_20vOFF" = S20Hz - SOFF
  )) 
  return(summary(model_effects))
}

brms_epred_calc_FLT <- function(model, ana){
  preds <- brms::posterior_epred(model)
  # load FLT data
  #FLTRT<- read_csv(file = "C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted/flanker.csv") %>%
  FLTRT<- read_csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted/flanker.csv") %>%
    mutate(Error = 1 - Correct_Response) %>%
    mutate(StimCon = as_factor(case_when(
      Stim_verb == "130Hz" ~ "S130Hz",
      Stim_verb == "20Hz" ~ "S20Hz",
      Stim_verb == "OFF" ~ "SOFF"
    ))) 
  
  # filter the data more for the RT analysis
  if (ana == "RT") {FLTRT <- FLTRT %>%
    filter(Correct_Response == 1,
           RT <3,
           RT > 0.2) %>%
    mutate(RT_ms = RT*1000)}
  
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
  
  # calculate contrasts
  S130Hz_congruent <- rowMeans(preds[,FLTRT$Contrast_F == "S130Hz_congruent"])
  S130Hz_incongruent <- rowMeans(preds[,FLTRT$Contrast_F == "S130Hz_incongruent"])
  S20Hz_congruent <- rowMeans(preds[,FLTRT$Contrast_F == "S20Hz_congruent"])
  S20Hz_incongruent <- rowMeans(preds[,FLTRT$Contrast_F == "S20Hz_incongruent"])
  SOFF_congruent <- rowMeans(preds[,FLTRT$Contrast_F == "SOFF_congruent"])
  SOFF_incongruent <- rowMeans(preds[,FLTRT$Contrast_F == "SOFF_incongruent"])
  
  
  model_effects <- tibble(contrast = c("Congruency", "Stim_20v130", "Stim_20vOFF", "Stroop_20v130", "Stroop_20vOFF"))
  Congruency <- tidybayes::mean_hdi(((S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3 -
    (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3))
  Stim_20v130 <- tidybayes::mean_hdi(((S20Hz_incongruent + S20Hz_congruent)/2 - (S130Hz_incongruent + S130Hz_congruent)/2))
  Stim_20vOFF <- tidybayes::mean_hdi(((S20Hz_incongruent + S20Hz_congruent)/2 - (SOFF_incongruent + SOFF_congruent)/2))
  Stroop_20v130 <- tidybayes::mean_hdi(((S20Hz_incongruent - S20Hz_congruent) - (S130Hz_incongruent - S130Hz_congruent)))
  Stroop_20vOFF <- tidybayes::mean_hdi(((S20Hz_incongruent - S20Hz_congruent) - (SOFF_incongruent - SOFF_congruent)))
  
  model_effects<- model_effects %>%
    bind_cols(bind_rows(Congruency, Stim_20v130, Stim_20vOFF, Stroop_20v130, Stroop_20vOFF))%>%
    rename(estimate = y, lower.HPD = ymin, upper.HPD = ymax)
  
  return(model_effects)
}

conditional_effect_calc_FLT <- function(model){
  
  model <- full_model
  # calculate effect estimates in ms or percent
  # Args:
  #   model: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   summary(model_effects): Table with summarizing the posterior distribution of effect estimates, median, upper and lower 95% hpd interval boundaries
  
  
  # get posterior samples to calculate conditional effects
  # first calculate the estimated marginal means
  emm_FLT <- emmeans(model, ~ Contrast_F, epred = TRUE)
  S130Hz_incongruent <- c(1, 0, 0, 0, 0, 0)
  S130Hz_congruent <- c(0, 1, 0, 0, 0, 0)
  SOFF_congruent <- c(0, 0, 1, 0, 0, 0)
  SOFF_incongruent <- c(0, 0, 0, 1, 0, 0)
  S20Hz_incongruent <- c(0, 0, 0, 0, 1, 0)
  S20Hz_congruent <- c(0, 0, 0, 0, 0, 1)
  
  
  # Next we calculate the contrasts of interest from these marginal means
  model_effects <- contrast(emm_FLT , method = list("Congruency" = (S130Hz_incongruent + SOFF_incongruent + S20Hz_incongruent)/3 -
                                                      (S130Hz_congruent + SOFF_congruent + S20Hz_congruent)/3,
                                                    "Stim_20v130" = (S20Hz_incongruent + S20Hz_congruent)/2 - (S130Hz_incongruent + S130Hz_congruent)/2,
                                                    "Stim_20vOFF" = (S20Hz_incongruent + S20Hz_congruent)/2 - (SOFF_incongruent + SOFF_congruent)/2,
                                                    "Stroop_20v130"=(S20Hz_incongruent - S20Hz_congruent) - (S130Hz_incongruent - S130Hz_congruent),
                                                    "Stroop_20vOFF"=(S20Hz_incongruent - S20Hz_congruent) - (SOFF_incongruent - SOFF_congruent)))
  
  return(summary(model_effects))
}

### First we get all the variables we need to call our readily calculated models

# contrast names seperately

#fullmod_loc <- r"{E:\20Hz\Data\Modelle}"
#part_mod_loc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
parameter_RT_GNG <- c("GoDiff_20v130",
               "GoDiff_20vOFF")

parameter_Acc_GNG <- c(
  "NoGo_Stop_20Hz_vs_130Hz",
  "NoGo_Stop_20Hz_vs_OFF",
  "NoGo_Go_20Hz_vs_130Hz",
  "NoGo_Go_20Hz_vs_OFF"
)

parameter_SRT <- c(
  "Stim_20v130",
  "Stim_20vOFF"
)

parameter_FLT <- c(
  "Congruency",
  "Stim_20v130",
  "Stim_20vOFF",
  "Stroop_20v130",
  "Stroop_20vOFF"
)

#model_RT <- "RT"
#model_Acc <- "Acc"

#Partial_models_saveloc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
#Full_models_saveloc <- r"{E:\20Hz\Data\Modelle}"
#Partial_models_saveloc <-'C:/Users/doex9445/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_mods'
#Full_models_saveloc <- 'C:/Users/doex9445/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle'
Partial_models_saveloc <- "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_mods"
Full_models_saveloc <- "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/"

# calculate BFs for the RT data of the GNG task
BF_Results <- BF_calc(fullmod_loc = Full_models_saveloc, 
                      part_mod_loc = Partial_models_saveloc, 
                      parameter = parameter_RT_GNG, 
                      ana = "GNG", 
                      model_t = "RT")

#save data

# Now the same for the Acc data of the GNG task
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                      part_mod_loc = Partial_models_saveloc, 
                      parameter = parameter_Acc_GNG, 
                      ana = "GNG", 
                      model_t = "Acc"))

# Now the same for the SRT task for the RT analysis
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                    part_mod_loc = Partial_models_saveloc, 
                    parameter = parameter_SRT, 
                    ana = "SRT", 
                    model_t = "RT"))

# Now the same for the SRT task for the accuracy analysis
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                    part_mod_loc = Partial_models_saveloc, 
                    parameter = parameter_SRT, 
                    ana = "SRT", 
                    model_t = "Acc"))

# Now the same for the SRT task for the RT analysis
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                    part_mod_loc = Partial_models_saveloc, 
                    parameter = parameter_FLT, 
                    ana = "FLT", 
                    model_t = "RT"))


# Now the same for the FLT task for the accuracy analysis
BF_Results <- BF_Results %>%
  bind_rows(BF_calc(fullmod_loc = Full_models_saveloc, 
                    part_mod_loc = Partial_models_saveloc, 
                    parameter = parameter_FLT, 
                    ana = "FLT", 
                    model_t = "Acc"))

#save data
write.table(BF_Results , file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_Results.csv")

# load data again
BF_Results <- read.csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/BF_Results.csv", header = TRUE, sep = "")
#### Now that we have the table with our BFs, let us load the parameter estimates

# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(estimate = NA , lower.HPD = NA, upper.HPD = NA) 

# Loop through analysis and save effects
tasks = c("GNG", "SRT", "FLT")
mods = c("RT", "Acc")
for (tsk in tasks){
  for (md in mods) {
    #load model
    full_model <- load_full_mod(loc = Full_models_saveloc, 
                                ana = tsk, 
                                model_t = md)
    
    # calculate posterior conditional effect samples
    if (tsk == "GNG"){
      if (md == "RT"){
        #eff_post <- conditional_effect_calc_shift_GNG(model = full_model)
        eff_post <- brms_epred_calc_GNG_RT(model = full_model)
      } else if (md == "Acc"){
        #eff_post <- conditional_effect_calc_acc_GNG(model = full_model)
        eff_post <- brms_epred_calc_GNG_acc(model = full_model)
      }
    } else if (tsk == "SRT"){
      #eff_post <- conditional_effect_calc_SRT(model = full_model)
      eff_post <- brms_epred_calc_SRT(model = full_model, ana = md)
    }
    else if (tsk == "FLT"){
      #eff_post <- conditional_effect_calc_FLT(model = full_model)
      eff_post <- brms_epred_calc_FLT(model = full_model, ana = md)
    }
    #integrate summary stats into BF table
    for (elem in eff_post$contrast){
      # find correct row
      idx <- Full_Model_Info$Analysis == md & Full_Model_Info$Parameter == elem & Full_Model_Info$Model == tsk
      Full_Model_Info$estimate[idx] <- eff_post$estimate[eff_post$contrast == elem]
      Full_Model_Info$lower.HPD[idx] <- eff_post$lower.HPD[eff_post$contrast == elem]
      Full_Model_Info$upper.HPD[idx] <- eff_post$upper.HPD[eff_post$contrast == elem]
    }
  }
}

#save data
#write.table(Full_Model_Info , file = "E:/20Hz/Data/Modelle/Full_Results.csv")
write.table(Full_Model_Info , file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/Full_Results1.csv")
Full_Model_Info2 <- read.csv(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/Full_Results1.csv", header = TRUE, sep = "")
Full_Model_Info2 <- Full_Model_Info2 %>%
  mutate(BF_new = round(BF, 1),
         estimate_new = round(estimate, 3),
         lower.HPD_new = round(lower.HPD, 3),
         upper.HPD_new = round(upper.HPD, 3))