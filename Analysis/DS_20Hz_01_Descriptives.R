#####
# Author: Julius Kricheldorff
# Descriptive Stats Demopgraphics
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
# Date 11.01.2023
####

# Address git : C/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted')

# load data
Demographics<- read_csv2(file = "PartChars.csv") 

Summary_quant <-  Demographics %>% 
  select(Alter, DauerPD, HnY, MMST, LEDD) %>%
  summarise(round(across(everything(), list(mean, sd)),1))



