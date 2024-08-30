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
library(ggplot2)

# set directory
setwd('D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted')

# load data
Demographics<- read_csv2(file = "PartChars.csv") 

Summary_quant <-  Demographics %>% 
  select(Alter, DauerPD, DauerDBS, HnY, MMST, LEDD) %>%
  summarise(round(across(everything(), list(mean, sd)),1))

# get infos for gender, dominant side and disease subtype
summary_gender <- Demographics %>%
  group_by( Geschlecht) %>%
  summarise(number = n())

summary_dominant <- Demographics %>%
  group_by(Dominanzseite) %>%
  summarise(number = n())

summary_dst <- Demographics %>%
  group_by(Krankheitstyp) %>%
  summarise(number = n())

# Now the same for the UPDRS data

# load data
path <- file.path("D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Original/Daten final/Ergebnisse/UPDRS.sav")
UPDRS_data <- read_sav(file = path) %>%
  mutate(Stim_verb = as.factor(case_when(
    Stim == 1 ~ "130Hz",
    Stim == 2 ~ "20Hz",
    Stim == 3 ~ "OFF"
  ))) 

# get mean values and sd

UPDRS_data %>%
  group_by(Stim_verb) %>%
  summarize(mean_UPDRS = mean(UPDRS), sd_UPDRS = sd(UPDRS))

# Plot the data just for a visual check and confirm that only 130Hz significantly reduced UPDRS
ggplot(data =UPDRS_data, aes( y = UPDRS, x = Stim_verb), colour = Stim_verb) +
  geom_boxplot()

# Perform t-test for inferential statistics
UPD_20Hz <- UPDRS_data %>%
  filter(Stim_verb == "20Hz") %>%
  arrange(Proband)
UPD_130Hz <- UPDRS_data %>%
  filter(Stim_verb == "130Hz")%>%
  arrange(Proband)
UPD_OFF <- UPDRS_data %>%
  filter(Stim_verb == "OFF")%>%
  arrange(Proband)

# Bayesian Version
BF_UPD_20_130 <- ttestBF(x = UPD_130Hz$UPDRS - UPD_20Hz$UPDRS)
BF_UPD_20_OFF <- ttestBF(x = UPD_OFF$UPDRS - UPD_20Hz$UPDRS)
BF_UPD_OFF_130 <- ttestBF(x = UPD_130Hz$UPDRS - UPD_OFF$UPDRS)


# Frequentist Version
t.test(UPD_130Hz$UPDRS - UPD_20Hz$UPDRS, alternative = "two.sided")
t.test(UPD_OFF$UPDRS - UPD_20Hz$UPDRS, alternative = "two.sided")
t.test(UPD_130Hz$UPDRS - UPD_OFF$UPDRS, alternative = "two.sided")





