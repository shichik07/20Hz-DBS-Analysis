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



# Bridge_comp <- function(full_model, deficient_model){
#   #calculate marginal logliklihood using bridge sampling
#   margLogLik_full <- bridge_sampler(full_model, silent = TRUE)
#   margLogLik_null <- bridge_sampler(deficient_model, silent = TRUE)
#   # next calculate the BF
#   BF <- bayes_factor(margLogLik_full, margLogLik_null)
#   return(BF$bf)
# }


# para <- c("LW_min_CO_Intercept",
#           "LW_min_CO_Congruency", 
#           "LW_min_CO_BlocknItem",
#           "LW_min_CO_2Factor",
#           "LW_min_PD_Intercept",
#           "LW_min_PD_Congruency",
#           "LW_min_PD_BlocknItem",
#           "LW_min_PD_2Factor")
# Liklihoods <- tibble(
#   Model = character(),
#   Log_lik = numeric()
# )
# for(i in c(1:8)){
#   temp_t <- tibble(Model = para[i],
#                    Log_lik =i)
#   Liklihoods <- bind_rows(Liklihoods, temp_t)
# 
# }

# because we need to add very small probabilities (exponentiate the log),
# this is possible with this function. Source:
# https://stats.stackexchange.com/questions/379335/adding-very-small-probabilities-how-to-compute
logsum <- function(l1, l2){
  max(l1,l2) + log1p(exp(-abs(l1-l2)))
}

Bayes_factor_calc <- function(Liklihoods, param){
  # function to calculate bayes inclusion factor across matched models
  Group_var = substr(param,1,2)
  if(grepl("IS_Block|LW_Block", param)){
    mod_with <- Liklihoods %>%
      filter(Model == "full") %>%
      pull(Log_lik)
    mod_without <-  Liklihoods %>%
      filter(grepl(Group_var, Model)) %>%
      filter(grepl("2Factor", Model)) %>%
      pull(Log_lik) 
    BF <- exp(mod_with - mod_without)
  } 
  else{ 
    if (grepl("Congruency",param)){
      filter_with <- "2Factor|Congruency"
      filter_without <- "Intercept|BlocknItem"}
    else if (grepl("Itemspecific|Listwide", param)){
      filter_with <- "2Factor|BlocknItem"
      filter_without <- "Intercept|Congruency"}
    # get the models with the effect of interest
    mod_with_vals <- Liklihoods %>% 
      filter(grepl(Group_var, Model)) %>%
      filter(grepl(filter_with, Model)) %>%
      pull(Log_lik)
    mod_with <- logsum(mod_with_vals[1], mod_with_vals[2])
    # get the models without the effect of interest
    mod_without_vals <- Liklihoods %>% 
      filter(grepl(Group_var, Model)) %>%
      filter(grepl(filter_without, Model)) %>%
      pull(Log_lik)
    mod_without <- logsum(mod_without_vals[1], mod_without_vals[2])
    BF <- exp(mod_with - mod_without) # we have to log transform again for division
  }
  return(BF)
}


load_part_mod <- function(ana, model_t, parameter){
  # loads the partial model
  #
  # Args:
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   model_t: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  #   parameter: indicates which parameter the model should be lacking and depends on ana
  # Returns:
  #   The specified load model with generically renamed as fit_model
  string <- paste(ana, model_t, ana, parameter, sep = "_")
  fit <- load(paste(string, ".rda", sep =""))
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}



load_full_mod <- function(loc, ana, model_t){
  # loads the partial model
  #
  # Args:
  #   loc: file location where the full model is saved
  #   ana: string for the task we are interested in, can take on the values "GNG", "SRT", "SST", "FLT"
  #   model_t: string for the model we are interested in "RT" for shifted log-normal, "Acc" for the logistic regression
  # Returns:
  #   The specified load model with generically renamed as fit_model

  if (model_t == "RT"){ 
    string <- file.path(loc, paste("shifted_log", ana, sep = "_"))
  } else if (model_t == "Acc"){
    string <- file.path(loc, paste("log_reg", ana, sep = "_"))
  }
  fit <- load(paste(string, ".rda", sep =""))
  fit_model <- eval(parse(text = fit)) # rename the model
  return(fit_model)
}

# BF_calc <- function(fullmod_loc, part_mod_loc){
#   # tibble to save the data
#   Model_BF <- tibble(
#     Model = character(),
#     Item_type = character(),
#     Effect = character(),
#     Group = character(),
#     Parameter = character(),
#     BF = numeric(),
#   )
#   mods = c("RT", "Acc")
#   item_type = c("inducer","diagnostic")
#   effect = c("LW", "IS")
#   for (model_type in mods){
#     for (itm in item_type){
#       for(eff in effect){
#         if (eff == "LW"){
#           parameter <- c(
#             "CO_Congruency",
#             "CO_Listwide",
#             "CO_LW_Block",
#             "PD_Congruency",
#             "PD_Listwide",
#             "PD_LW_Block"
#           )
#         } else if (eff == "IS"){
#           parameter <- c(
#             "CO_Congruency",
#             "CO_Itemspecific",
#             "CO_IS_Block",
#             "PD_Congruency",
#             "PD_Itemspecific",
#             "PD_IS_Block"
#           )
#         }
#         # load the full model
#         full_mod <- load_full_mod(loc = fullmod_loc, 
#                                   model_t = model_type,
#                                   item_type = itm, 
#                                   effect = eff)
#         for (par in parameter){
#           # load the partial model
#           def_mod <- load_part_mod(loc = part_mod_loc, 
#                                    model_t = model_type, 
#                                    item_type = itm, 
#                                    effect = eff, 
#                                    parameter = par)
#           # claculate the BF
#           BF_inc <- Bridge_comp(full_model = full_mod,
#                                 deficient_model = def_mod)
#           #save results as a tibble
#           temp <- tibble(
#             Model = model_type,
#             Item_type = itm,
#             Effect = eff,
#             Group = substr(par, start = 1, stop = 2),
#             Parameter = substr(par, start = 4, stop = nchar(par)),
#             BF = BF_inc,
#           )
#           # join tibble as row 
#           Model_BF <- bind_rows(Model_BF, temp)
#         }
#       }
#     }
#   }
#   return(Model_BF)
# }

BF_calc <- function(fullmod_loc, part_mod_loc){
  # tibble to save the data
  Model_BF <- tibble(
    Model = character(),
    Item_type = character(),
    Effect = character(),
    Group = character(),
    Parameter = character(),
    BF = numeric(),
  )
  
  #mods = c("RT", "Acc") # Acc models have yet to be fit
  mods = "RT"
  item_type = c("inducer","diagnostic")
  effect = c("LW", "IS")
  for (model_type in mods){
    for (itm in item_type){
      for(eff in effect){
        print(paste("Analyzing the", eff, "effect for the", itm, "items of the", model_type, "data"))
        # tibble to save the liklihood values
        Liklihoods <- tibble(
          Model = character(),
          Log_lik = numeric()
        )
        if (eff == "LW"){
          # partial models that we are going to load
          submodel <- c("LW_min_CO_Intercept",
                        "LW_min_CO_Congruency", 
                        "LW_min_CO_BlocknItem",
                        "LW_min_CO_2Factor",
                        "LW_min_PD_Intercept",
                        "LW_min_PD_Congruency",
                        "LW_min_PD_BlocknItem",
                        "LW_min_PD_2Factor")
          # parameters that we need to calculate the BIF for
          parameter <- c(
            "CO_Congruency",
            "CO_Listwide",
            "CO_LW_Block",
            "PD_Congruency",
            "PD_Listwide",
            "PD_LW_Block"
          )
        } else if (eff == "IS"){
          # partial models that we are going to load
          submodel <- c("IS_min_CO_Intercept",
                        "IS_min_CO_Congruency", 
                        "IS_min_CO_BlocknItem",
                        "IS_min_CO_2Factor",
                        "IS_min_PD_Intercept",
                        "IS_min_PD_Congruency",
                        "IS_min_PD_BlocknItem",
                        "IS_min_PD_2Factor")
          # parameters that we need to calculate the BIF for
          parameter <- c(
            "CO_Congruency",
            "CO_Itemspecific",
            "CO_IS_Block",
            "PD_Congruency",
            "PD_Itemspecific",
            "PD_IS_Block"
          )
        }
        # load the full model
        full_mod <- load_full_mod(loc = fullmod_loc, 
                                  model_t = model_type,
                                  item_type = itm, 
                                  effect = eff)
        # calculate the log liklihood
        full_mod_loglik <- bridge_sampler(full_mod, silent = TRUE)
        temp_t <- tibble(Model = "full",
                         Log_lik = full_mod_loglik$logml)
        Liklihoods <- bind_rows(Liklihoods,temp_t)
        
        # now we load the partial models
        for (par in submodel){
          # load the partial model
          def_mod <- load_part_mod(loc = part_mod_loc, 
                                   model_t = model_type, 
                                   item_type = itm, 
                                   effect = eff, 
                                   parameter = par)
          # calculate the log liklihood
          def_mod_loglik <- bridge_sampler(def_mod, silent = TRUE)
          temp_t <- tibble(Model = par,
                           Log_lik = def_mod_loglik$logml)
          Liklihoods <- bind_rows(Liklihoods,temp_t)
        }
        # now use the logliklihoods to calculate the bayes inclusion factor
        for (params in parameter){
          print(params)
          BF <- Bayes_factor_calc(Liklihoods = Liklihoods, 
                                               param = params)
          #save results as a tibble
          temp <- tibble(
            Model = model_type,
            Item_type = itm,
            Effect = eff,
            Group = substr(params, start = 1, stop = 2),
            Parameter = substr(params, start = 4, stop = nchar(params)),
            BF = BF,
          )
          # join tibble as row 
          Model_BF <- Model_BF %>%
            bind_rows(temp)
        }
        
      }
    }
  }
  return(Model_BF)
}




# calculate effect estimates in ms
conditional_effect_calc_shift <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block) + m_post$sigma/2)
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      CC_CO = MC_C_CO - MI_C_CO,
      CC_PD = MC_C_PD - MI_C_PD,
      CI_CO = MI_I_CO - MC_I_CO,
      CI_PD = MI_I_PD - MC_I_PD,
      Inc_diff_CO_PD = CI_CO - CI_PD,
      Con_diff_CO_PD = CC_CO - CC_PD)
  } else if (effect == "IS"){
    MC_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_C_PD <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MI_I_PD <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block) + m_post$sigma/2)
    MC_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MC_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_C_CO <- exp((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    MI_I_CO <- exp((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block) + m_post$sigma/2)
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      CC_CO = MC_C_CO - MI_C_CO,
      CC_PD = MC_C_PD - MI_C_PD,
      CI_CO = MI_I_CO - MC_I_CO,
      CI_PD = MI_I_PD - MC_I_PD,
      Inc_diff_CO_PD = CI_CO - CI_PD,
      Con_diff_CO_PD = CC_CO - CC_PD)
  }
  return(model_effects)
}

conditional_effect_calc_acc <- function(effect, model){
  # get posterior samples to calculate conditional effects
  m_post <- posterior_samples(model)
  # calculate conditional effects depending on model
  if (effect == "LW"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency + 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  - 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWPD_Congruency - 0.5 * m_post$b_Contrast_LWPD_Listwide  + 0.5*m_post$b_Contrast_LWPD_LW_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency + 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  - 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_LWCO_Congruency - 0.5 * m_post$b_Contrast_LWCO_Listwide  + 0.5*m_post$b_Contrast_LWCO_LW_Block))*100
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      CC_CO = MC_C_CO - MI_C_CO,
      CC_PD = MC_C_PD - MI_C_PD,
      CI_CO = MI_I_CO - MC_I_CO,
      CI_PD = MI_I_PD - MC_I_PD,
      Inc_diff_CO_PD = CI_CO - CI_PD,
      Con_diff_CO_PD = CC_CO - CC_PD)
  } else if (effect == "IS"){
    MC_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency + 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_C_PD <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  - 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MI_I_PD <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISPD_Congruency - 0.5 * m_post$b_Contrast_ISPD_Itemspecific  + 0.5*m_post$b_Contrast_ISPD_IS_Block))*100
    MC_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MC_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency + 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_C_CO <- plogis((m_post$b_Intercept + 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  - 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    MI_I_CO <- plogis((m_post$b_Intercept - 0.5* m_post$b_Contrast_ISCO_Congruency - 0.5 * m_post$b_Contrast_ISCO_Itemspecific  + 0.5*m_post$b_Contrast_ISCO_IS_Block))*100
    model_effects <- tibble(
      Congruency_PD = (MC_I_PD + MI_I_PD)/2 - (MC_C_PD + MI_C_PD)/2,
      Congruency_CO = (MC_I_CO + MI_I_CO)/2 - (MC_C_CO + MI_C_CO)/2,
      Block_CO = (MC_I_CO + MC_C_CO)/2 - (MI_I_CO + MI_C_CO)/2,
      Block_PD = (MC_I_PD + MC_C_PD)/2 - (MI_I_PD + MI_C_PD)/2,
      Stroop_MC_PD = MC_I_PD - MC_C_PD,
      Stroop_MI_PD = MI_I_PD - MI_C_PD,
      Stroop_MC_CO = MC_I_CO - MC_C_CO,
      Stroop_MI_CO = MI_I_CO - MI_C_CO,
      Control_PD = Stroop_MC_PD - Stroop_MI_PD,
      Control_CO = Stroop_MC_CO - Stroop_MI_CO,
      CC_CO = MC_C_CO - MI_C_CO,
      CC_PD = MC_C_PD - MI_C_PD,
      CI_CO = MI_I_CO - MC_I_CO,
      CI_PD = MI_I_PD - MC_I_PD,
      Inc_diff_CO_PD = CI_CO - CI_PD,
      Con_diff_CO_PD = CC_CO - CC_PD)
  }
  return(model_effects)
}

# function calculate summary statistics from posterior distributions
post_sum_calc <- function(summary_table){
  col_n <- colnames(summary_table)
  model_summary <- tibble(
    parameter = character(),
    mean = numeric(),
    lower95 = numeric(),
    upper95 = numeric()
  )
  for (var_n in col_n){
    var_t <- eval(parse(text = paste('summary_table$', var_n)))
    temp <- tibble(
      parameter = var_n,
      mean = mean(var_t),
      lower95 = quantile(var_t, probs = 0.025),
      upper95 = quantile(var_t, probs = 0.975)
    )
    # join new row with summary table
    model_summary <- bind_rows(model_summary, temp)
  }
  return(model_summary)
}

### First we get all the variables we need to call our readily calculated models

# contrast names seperately

Partial_models_saveloc <- "E:/AdaptiveControl/Data/BehaviorResults"
#Full_models_saveloc <- "E:/AdaptiveControl/Data/BehaviorResults"
Full_models_saveloc <- "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data"

# calculate BFs
BF_Results <- BF_calc( 
  fullmod_loc = Full_models_saveloc, 
  part_mod_loc = Partial_models_saveloc)

#save data
write.table(BF_Results , file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/BF_Results.csv")

# load data again
BF_Results <- read.csv(file = "C:/Users/doex9445/Dateien/Julius/AdaptiveControl/Data/BF_Results.csv", header = TRUE, sep = "")
#### Now that we have the table with our BFs, let us load the parameter estimates

# New Tibble
Full_Model_Info <- BF_Results %>%
  mutate(mean = NA , lower95 = NA, upper95 = NA) 

# Load and save the data

mods = c("RT")
item_type = c("inducer","diagnostic")
effect = c("LW", "IS")

for (md in mods) {
  for (eff in effect){
    for (itm in item_type){
      #load model
      model <- load_full_mod(loc = Full_models_saveloc,
                             model_t = md,
                             item_type = itm,
                             effect = eff)
      
      # calculate posterior conditional effect samples
      if (md == "RT"){
        eff_post <- conditional_effect_calc_shift(effect = eff,
                                                  model = model)
      } else if (md == "Acc"){
        eff_post <- conditional_effect_calc_acc(effect = eff,
                                                model = model)
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
