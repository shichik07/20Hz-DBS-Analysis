#####
# Author: Julius Kricheldorff
# Bayesian analysis of the reaction time data
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
# Date 11.01.2023
####

# First we determine our contrasts of interest. When we are not looking for an 
# interaction this is easy. We simply have two contrast coded variables for 20Hz
# versus 130Hz and 20Hz versus OFF condition. For the ineractions with congruency
# it gets a little more tricky. We again have two effect coded contrasts for the main
# effect of Stimulation setting. One for the main effect of congruency and two to
# compare the interaction effects 20Hz/130Hz and 20Hz/OFF. 

setwd('C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted')

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)

#library(chkptstanr) #package to interupt and restart sampling with stanr/brmslibrary(posterior)

# Set a seed for sake of reproducibility
set.seed(32946)

# create a hypothetical contrast matrix for the flanker task
TreatmentContrast <- hypr(
  S130_S20 = (S130Hz_incongruent + S130Hz_congruent)/2 ~ (S20Hz_incongruent + S20Hz_congruent)/2, # main effect Congruence healthy controls
  SOFF_S20 = (SOFFHz_incongruent + SOFFHz_congruent)/2 ~ (S20Hz_incongruent + S20Hz_congruent)/2, # main effect Block Listwide control
  Congruency = (S130Hz_incongruent + S20Hz_incongruent +SOFFHz_incongruent)/3 ~ 
    (S130Hz_congruent + S20Hz_congruent +SOFFHz_congruent)/3, # interaction listwide effect
  Int_S20_S130 = (S130Hz_incongruent - S130Hz_congruent) ~ (S20Hz_incongruent - S20Hz_congruent),
  Int_S20_S130 = (SOFFHz_incongruent - SOFFHz_congruent) ~ (S20Hz_incongruent - S20Hz_congruent),
  levels = c("S130Hz_incongruent", "S130Hz_congruent", "S20Hz_incongruent", "S20Hz_congruent",
             "SOFFHz_incongruent", "SOFFHz_congruent") 
  
)
TreatmentContrast

