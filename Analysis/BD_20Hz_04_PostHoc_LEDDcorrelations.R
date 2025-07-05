#####
# Author: Julius Kricheldorff
# Post-Hoc Analysis - LEDD and effect sizes
# Includes task data the Flanker-, Stop-Signal-, Go-NoGo and Simple RT task
# Date 25.11.2024
####

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# set directory
setwd('D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted')

# load data
Demographics<- read_csv2(file = "PartChars.csv") |>
  select(ProbandNEU,ProbandALT, LEDD)

# Next we need to extract the effects of interest for our analysis - in this case it makes the most sense to only do this for the GoNoGo analys