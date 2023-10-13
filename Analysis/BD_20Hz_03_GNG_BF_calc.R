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


# Set a seed for sake of reproducibility
set.seed(32936)

Bayes_factor_calc <- function(Likelihoods, param, ana){
  # function to calculate Bayes factor
  #
  # Args:
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   Likelihoods: tibble containing likelihoods for all models
  #   param: parameter we want to calculate BF for
  # Returns:
  #   Bayes Factor (numeric) for the specified parameter
  sub_param <- paste(ana, "min", param, sep = "_")
    mod_with <- Likelihoods %>%
      filter(Model == "full") %>%
      pull(Log_lik)
    mod_without <-  Likelihoods %>%
      filter(Model == sub_param) %>%
      pull(Log_lik) 
  
    BF <- exp(mod_with - mod_without) # we have to log transform again for division
    return(BF)
}

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
  # calculate the log likelihood
  full_mod_loglik <- bridge_sampler(full_mod, silent = TRUE)
  temp_t <- tibble(Model = "full",
                   Log_lik = full_mod_loglik$logml)
  Likelihoods <- bind_rows(Likelihoods,temp_t)
  
  # now we load the partial models
  for (par in submodel){
    # load the partial model
    def_mod <- load_part_mod(loc = part_mod_loc,
                             ana = ana, 
                             model_t = model_t, 
                             param = par)
    # calculate the log liklihood
    def_mod_loglik <- bridge_sampler(def_mod, silent = TRUE)
    temp_t <- tibble(Model = par,
                     Log_lik = def_mod_loglik$logml)
    Likelihoods <- bind_rows(Likelihoods,temp_t)
  }
  # now use the log-likelihoods to calculate the Bayes inclusion factor
  for (params in parameter){
    print(params)
    BF <- Bayes_factor_calc(Likelihoods = Likelihoods,
                            param = params,
                            ana = ana)
    #save results as a tibble
    temp <- tibble(
      Model = ana,
      Parameter = params,
      BF = BF,
    )
    # join tibble as row
    Model_BF <- Model_BF %>%
      bind_rows(temp)
  }
  return(Model_BF)
}

fullmod_loc <- r"{E:\20Hz\Data\Modelle}"
parameter <- c("Go_NoGo_Go", 
               "Stim_20v130",
               "Stim_20vOFF",
               "GoDiff_20v130",
               "GoDiff_20vOFF")

model<- load_full_mod(fullmod_loc, ana, model_t)

summary_table <- conditional_effect_calc_shift(model)


conditional_effect_calc_shift <- function(model){
  # calculate effect estimates in ms or percent
  # Args:
  #   model: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: string indicating which parameter the model should be lacking and depends on ana
  # Returns:
  #   summary(model_effects): Table with summarizing the posterior distribution of effect estimates, median, upper and lower 95% hpd interval boundaries
  
  
  # get posterior samples to calculate conditional effects
  # first calculate the estimated marginal means
  emm <- emmeans(model, specs = ~ Contrast_F, epred = TRUE)
  
  # define new contrasts for emmeans - there should be an easier way to this
  S130Hz_NoGo_Go <- c(1, 0, 0, 0, 0, 0)
  S130Hz_Go <- c(0, 1, 0, 0, 0, 0)
  SOFF_Go <- c(0, 0, 1, 0, 0, 0)
  SOFF_NoGo_Go <- c(0, 0, 0, 1, 0, 0)
  S20Hz_NoGo_Go <- c(0, 0, 0, 0, 1, 0)
  S20Hz_Go <- c(0, 0, 0, 0, 0, 1)
  
  # and contrasts for the differences in Go effects
  Go_NoGo_Go <- (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3 - (S130Hz_Go + SOFF_Go + S20Hz_Go)/3
  Stim_20v130 <- (S20Hz_NoGo_Go + S20Hz_Go)/2 - (S130Hz_NoGo_Go + S130Hz_Go)/2
  Stim_20vOFF <- (S20Hz_NoGo_Go + S20Hz_Go)/2 - (SOFF_NoGo_Go + SOFF_Go)/2
  GoDiff_20v130 <- (S20Hz_NoGo_Go - S20Hz_Go) - (S130Hz_NoGo_Go - S130Hz_Go)
  GoDiff_20vOFF <- (S20Hz_NoGo_Go - S20Hz_Go) - (SOFF_NoGo_Go - SOFF_Go)
  
  # Next we calculate the contrasts of interest from these marginal means
  model_effects <- contrast(emm, method = list("Go_NoGo_Go" = Go_NoGo_Go,
                                     "Stim_20v130" = Stim_20v130,
                                     "Stim_20vOFF" = Stim_20vOFF,
                                     "GoDiff_20v130" = GoDiff_20v130,
                                     "GoDiff_20vOFF" = GoDiff_20vOFF))
  return(summary(model_effects))
}

### First we get all the variables we need to call our readily calculated models

# contrast names seperately

fullmod_loc <- r"{E:\20Hz\Data\Modelle}"
part_mod_loc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
parameter <- c("Go_NoGo_Go", 
               "Stim_20v130",
               "Stim_20vOFF",
               "GoDiff_20v130",
               "GoDiff_20vOFF")

ana <- "GNG"
model_t <- "RT"

Partial_models_saveloc <- r"{E:\20Hz\Data\Modelle\BF_mods}"
Full_models_saveloc <- r"{E:\20Hz\Data\Modelle}"

# calculate BFs
BF_Results <- BF_calc(fullmod_loc = Full_models_saveloc, 
                      part_mod_loc = Partial_models_saveloc, 
                      parameter = parameter, 
                      ana = ana, 
                      model_t = model_t)
#save data
write.table(BF_Results , file = "r{E:\20Hz\Data\Modelle\BF_Results.csv}")

# load data again
BF_Results <- read.csv(file = "r{E:\20Hz\Data\Modelle\BF_Results.csv}", header = TRUE, sep = "")
#### Now that we have the table with our BFs, let us load the parameter estimates

# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(estimate = NA , lower.HPD = NA, upper.HPD = NA) 

# Load and save the data

md = c("RT")
ana <- "GNG"
model_t <- "RT"


for (md in mods) {
      #load model
      full_model <- load_full_mod(loc = fullmod_loc, 
                             ana = ana, 
                             model_t = model_t)
      
      # calculate posterior conditional effect samples
      if (md == "RT"){
        eff_post <- conditional_effect_calc_shift(model = full_model)
      } else if (md == "Acc"){
        eff_post <- conditional_effect_calc_acc(model = full_model)
      }
      #translate into summary stats
      sum_t <- post_sum_calc(eff_post)
      
      #integrate summary stats into BF table
      
      temp_t <- sum_t %>% filter(!grepl("Stroop", parameter)) # filter the stroop effects which we are not interested in to show
      for (vars in temp_t$parameter){
        group_in <- substr(vars, start = nchar(vars)-1, stop = nchar(vars))
        
        #NOTE: by accident when i created the models for the BF calculation, also the models
        # get a single row
        row_var <- sum_t %>% filter(parameter == vars)
        if (grepl("Congruency" ,vars)){
          par = "Congruency"
          # Write values
          Full_Model_Info$mean[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$mean
          Full_Model_Info$lower95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$lower95
          Full_Model_Info$upper95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$upper95
        } 
        if (grepl("Block" ,vars)){
          if (eff == "LW"){
            par = "Listwide" 
          } else {
            par = "Itemspecific"
          }
          
          # Write values
          Full_Model_Info$mean[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$mean
          Full_Model_Info$lower95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$lower95
          Full_Model_Info$upper95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$upper95
        } 
        if (grepl("Control" ,vars)) {
          if (eff == "LW"){
            par = "LW_Block" 
          } else {
            par = "IS_Block"
          }
          
          # Write values
          Full_Model_Info$mean[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$mean
          Full_Model_Info$lower95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$lower95
          Full_Model_Info$upper95[Full_Model_Info$Item_type == itm & Full_Model_Info$Model == md & Full_Model_Info$Effect == eff & Full_Model_Info$Group == group_in & Full_Model_Info$Parameter == par] <- row_var$upper95
        }
        if (grepl("CI" ,vars) | grepl("CC" ,vars) ){
          var_n <- substr(vars, start = 1, stop = 2)
          temp <- tibble(Model = md,
                         Item_type = itm,
                         Effect = eff,
                         Group = group_in,
                         Parameter = var_n,
                         BF = NA,
                         mean = row_var$mean,
                         lower95 = row_var$lower95,
                         upper95 = row_var$upper95)
          Full_Model_Info <- bind_rows(Full_Model_Info, temp)
          
        }
        if (grepl("Inc_diff" ,vars) | grepl("Con_diff" ,vars) ){
          var_n <- substr(vars, start = 1, stop = 3)
          temp <- tibble(Model = md,
                         Item_type = itm,
                         Effect = eff,
                         Group = group_in,
                         Parameter = var_n,
                         BF = NA,
                         mean = row_var$mean,
                         lower95 = row_var$lower95,
                         upper95 = row_var$upper95)
          Full_Model_Info <- bind_rows(Full_Model_Info, temp)
        }
        
      }
    }
  }
}
FM_old <- Full_Model_Info %>%
  mutate(mean = round(mean,1),
         lower95 = round(lower95,1),
         upper95 = round(upper95,1)) %>%
  filter(Effect == "IS")
#NOTE: by accident when i created the models for the BF calculation, also the models
#Note: When i created the models for the BF calculation, by accident, I named them the same 
# for the item specific effect as the list wise effect. So the itemspecific model actually is named listwise (variable names are correct though)
# same goes for the interaction. Thus in the above code only LW_Block and Listwise are used. Does not affect the results though.
#For better readability I correct that below
Full_Model_Info_fin <- Full_Model_Info #%>%
#mutate(Parameter = case_when(
#Effect == "IS" & Parameter == "Listwide" ~ "Itemspecific",
#Effect == "IS" & Parameter == "LW_Block" ~ "IS_Block",
#TRUE ~ Parameter))

# Export as latex table
print(xtable(Full_Model_Info_fin, type = "latex"), file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/Bheav_Summary.tex")

# Filter only the RT Analysis Results

Full_Model_Info_fin_Acc <- Full_Model_Info_fin %>%
  filter(Model == "Acc") %>%
  filter(Effect == "LW") %>%
  mutate(mean = round(mean,1)) %>%
  mutate(lower95 = round(lower95,1)) %>%
  mutate(upper95 = round(upper95,1)) %>%
  mutate(BF = BF)

print(xtable(Full_Model_Info_fin_Acc, type = "latex"), file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/Bheav_Summary_RT.tex")
