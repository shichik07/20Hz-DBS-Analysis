# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(hypr)
library(tidybayes)
library(emmeans)

# set directory
wd <-"D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

#### Posterior predictve Checks GoNoGo Reaction Times

# Load the data
# load data
GoNoGo<- read_csv(file = "GoNoGo.csv") %>%
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

# Load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/shifted_log_GNG.rda")

# first calculate the posterior predictions
a <- brms::posterior_epred(fit_shifted_log_GNG)

# summarize over conditions
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/RT_GNG_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 1, 3))

S130Hz_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "S130Hz_NoGo_Go"])
S130Hz_Go <- rowMeans(a[,RT_data$Contrast_F == "S130Hz_Go"])
pred_130 <- S130Hz_NoGo_Go - S130Hz_Go
mean_130 <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_Go"])
hist(pred_130, main = "Go(NoGo) - Go:  130Hz Stimulation", xlab = c("difference in ms"))
abline(v = mean_130, col = "red", lwd = 3)

SOFF_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "SOFF_NoGo_Go"])
SOFF_Go <- rowMeans(a[,RT_data$Contrast_F == "SOFF_Go"])
pred_OFF <- SOFF_NoGo_Go - SOFF_Go
mean_OFF <-  mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_Go"])
hist(pred_OFF, main = "Go(NoGo) - Go:  OFF Stimulation", xlab = c("difference in ms"))
abline(v = mean_OFF, col = "red", lwd = 3)

S20Hz_NoGo_Go <- rowMeans(a[,RT_data$Contrast_F == "S20Hz_NoGo_Go"])
S20Hz_Go <- rowMeans(a[,RT_data$Contrast_F == "S20Hz_Go"])
pred_20 <- S20Hz_NoGo_Go - S20Hz_Go
mean_20 <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_NoGo_Go"]) - mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_Go"])
hist(pred_20, main = "Go(NoGo) - Go:  20Hz Stimulation", xlab = c("difference in ms"))
abline(v = mean_20, col = "red", lwd = 3)

#save
dev.off()

# clean up workspace a little
rm(a, fit_shifted_log_GNG)

# Now the same for the Accuracy Data

#load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/log_reg_GNG.rda")

# get predictions 
acc_GNG_pred <- brms::posterior_epred(fit_log_flanker2)

# set up the figure
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/RT_GNG_acc_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 3, 3))

# get data for 130Hz
S130Hz_NoGo_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S130Hz_NoGo_Go"])
S130Hz_NoGo_Stop <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S130Hz_NoGo_Stop"])
S130Hz_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S130Hz_Go"])
mean_130_NoGo_Go <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "S130Hz_NoGo_Go"])
mean_130_NoGo_Stop <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "S130Hz_NoGo_Stop"]) 
mean_130_Go <-mean(GoNoGo$Error[GoNoGo$Contrast_F == "S130Hz_Go"])

# Plot 130 Figures
hist(S130Hz_NoGo_Go, main = "NoGo (Go):  130Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_130_NoGo_Go, col = "red", lwd = 3)
hist(S130Hz_NoGo_Stop, main = "NoGo (Stop):  130Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_130_NoGo_Stop, col = "red", lwd = 3)
hist(S130Hz_Go, main = "Go:  130Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_130_Go, col = "red", lwd = 3)

# get data for 20Hz
S20Hz_NoGo_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S20Hz_NoGo_Go"])
S20Hz_NoGo_Stop <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S20Hz_NoGo_Stop"])
S20Hz_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "S20Hz_Go"])
mean_20_NoGo_Go <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "S20Hz_NoGo_Go"])
mean_20_NoGo_Stop <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "S20Hz_NoGo_Stop"]) 
mean_20_Go <-mean(GoNoGo$Error[GoNoGo$Contrast_F == "S20Hz_Go"])

# Plot 20 Figures
hist(S20Hz_NoGo_Go, main = "NoGo (Go):  20Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_20_NoGo_Go, col = "red", lwd = 3)
hist(S20Hz_NoGo_Stop, main = "NoGo (Stop):  20Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_20_NoGo_Stop, col = "red", lwd = 3)
hist(S20Hz_Go, main = "Go:  20Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_20_Go, col = "red", lwd = 3)

# get data for OFF
SOFF_NoGo_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "SOFF_NoGo_Go"])
SOFF_NoGo_Stop <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "SOFF_NoGo_Stop"])
SOFF_Go <- rowMeans(acc_GNG_pred[,GoNoGo$Contrast_F == "SOFF_Go"])
mean_OFF_NoGo_Go <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "SOFF_NoGo_Go"])
mean_OFF_NoGo_Stop <- mean(GoNoGo$Error[GoNoGo$Contrast_F == "SOFF_NoGo_Stop"]) 
mean_OFF_Go <-mean(GoNoGo$Error[GoNoGo$Contrast_F == "SOFF_Go"])

# Plot OFF Figures
hist(SOFF_NoGo_Go, main = "NoGo (Go):  OFF Stimulation", xlab = c("Error in %"))
abline(v = mean_OFF_NoGo_Go, col = "red", lwd = 3)
hist(SOFF_NoGo_Stop, main = "NoGo (Stop):  OFF Stimulation", xlab = c("Error in %"))
abline(v = mean_OFF_NoGo_Stop, col = "red", lwd = 3)
hist(SOFF_Go, main = "Go:  OFF Stimulation", xlab = c("Error in %"))
abline(v = mean_OFF_Go, col = "red", lwd = 3)

# save 
dev.off()

# clean up workspace
rm(acc_GNG_pred, fit_log_flanker2)

#### Posterior Predictive Checks Response Selection Task Reaction Times

#load the data
# load data
SimpleRT<- read_csv(file = "SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

RT_data <- SimpleRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

# Load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/shifted_log_SRT.rda")

# first calculate the posterior predictions
RT_SRT_pred <- brms::posterior_epred(fit_shifted_log_SRT)

# summarize over conditions
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/RT_SRT_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 1, 2))

# get data 
S130Hz <- rowMeans(RT_SRT_pred[,RT_data$StimCon == "S130Hz"])
S20Hz <- rowMeans(RT_SRT_pred[,RT_data$StimCon == "S20Hz"])
SOFF <- rowMeans(RT_SRT_pred[,RT_data$StimCon == "SOFF"])
pred_20_130 <- S20Hz - S130Hz
pred_20_OFF <- S20Hz - SOFF
mean_20_130 <- mean(RT_data$RT_ms[RT_data$StimCon == "S20Hz"]) - mean(RT_data$RT_ms[RT_data$StimCon == "S130Hz"])
mean_20_OFF <- mean(RT_data$RT_ms[RT_data$StimCon == "S20Hz"]) - mean(RT_data$RT_ms[RT_data$StimCon == "SOFF"])

# Plot Figures
hist(pred_20_130, main = "RT 20Hz versus 130Hz", xlab = c("difference in ms"))
abline(v = mean_20_130, col = "red", lwd = 3)
hist(pred_20_OFF, main = "RT 20Hz versus OFF", xlab = c("difference in ms"))
abline(v = mean_20_OFF, col = "red", lwd = 3)

#save
dev.off()

# clean up workspace a little
rm(RT_SRT_pred, fit_shifted_log_SRT)

# Now the same for the Accuracy Data

# Load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/log_reg_SRT.rda")

# first calculate the posterior predictions
acc_SRT_pred <- brms::posterior_epred(fit_logReg_SRT)

# summarize over conditions
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/acc_SRT_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 1, 2))

# get data 
S130Hz <- rowMeans(acc_SRT_pred[,SimpleRT$StimCon == "S130Hz"])
S20Hz <- rowMeans(acc_SRT_pred[,SimpleRT$StimCon == "S20Hz"])
SOFF <- rowMeans(acc_SRT_pred[,SimpleRT$StimCon == "SOFF"])
pred_20_130 <- S20Hz - S130Hz
pred_20_OFF <- S20Hz - SOFF
mean_20_130 <- mean(SimpleRT$Error[SimpleRT$StimCon == "S20Hz"]) - mean(SimpleRT$Error[SimpleRT$StimCon == "S130Hz"])
mean_20_OFF <- mean(SimpleRT$Error[SimpleRT$StimCon == "S20Hz"]) - mean(SimpleRT$Error[SimpleRT$StimCon == "SOFF"])

# Plot Figures
hist(pred_20_130, main = "RT 20Hz versus 130Hz", xlab = c("difference in ms"))
abline(v = mean_20_130, col = "red", lwd = 3)
hist(pred_20_OFF, main = "RT 20Hz versus OFF", xlab = c("difference in ms"))
abline(v = mean_20_OFF, col = "red", lwd = 3)

#save
dev.off()

# clean up workspace a little
rm(acc_SRT_pred, fit_logReg_SRT)

#### Posterior Predictive Checks Flanker Task Reaction Times

# load the data
# load data
FLTRT<- read_csv(file = "flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

FLTRT <- FLTRT %>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & Congruency == "congruent" ~ "S130Hz_congruent",
    Stim_verb == "130Hz" & Congruency == "incongruent" ~ "S130Hz_incongruent",
    Stim_verb == "20Hz" & Congruency == "congruent" ~ "S20Hz_congruent",
    Stim_verb == "20Hz" & Congruency == "incongruent" ~ "S20Hz_incongruent",
    Stim_verb == "OFF" & Congruency == "congruent" ~ "SOFF_congruent",
    Stim_verb == "OFF" & Congruency == "incongruent" ~ "SOFF_incongruent",
  )))

RT_data <- FLTRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

# load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/shifted_log_FLT.rda")

# calculate predictions
RT_FLT_pred <- brms::posterior_epred(fit_shifted_log_FLT)

# summarize over conditions
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/RT_FLT_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 3, 2))

# get data for 130Hz
S130Hz_Congruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "S130Hz_congruent"])
S130Hz_Incongruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "S130Hz_incongruent"])
mean_130Hz_Congruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_congruent"])
mean_130Hz_Incongruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S130Hz_incongruent"]) 

# Plot 130 Figures
hist(S130Hz_Congruent, main = "Congruent:  130Hz Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_130Hz_Congruent, col = "red", lwd = 3)
hist(S130Hz_Incongruent, main = "Incongruent:  130Hz Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_130Hz_Incongruent, col = "red", lwd = 3)

# get data for 20Hz
S20Hz_Congruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "S20Hz_congruent"])
S20Hz_Incongruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "S20Hz_incongruent"])
mean_20Hz_Congruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_congruent"])
mean_20Hz_Incongruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "S20Hz_incongruent"]) 

# Plot 20 Figures
hist(S20Hz_Congruent, main = "Congruent:  20Hz Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_20Hz_Congruent, col = "red", lwd = 3)
hist(S20Hz_Incongruent, main = "Incongruent:  20Hz Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_20Hz_Incongruent, col = "red", lwd = 3)

# get data for OFF
SOFF_Congruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "SOFF_congruent"])
SOFF_Incongruent <- rowMeans(RT_FLT_pred[,RT_data$Contrast_F == "SOFF_incongruent"])
mean_OFF_Congruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_congruent"])
mean_OFF_Incongruent <- mean(RT_data$RT_ms[RT_data$Contrast_F == "SOFF_incongruent"]) 

# Plot OFF Figures
hist(SOFF_Congruent, main = "Congruent:  OFF Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_OFF_Congruent, col = "red", lwd = 3)
hist(SOFF_Incongruent, main = "Incongruent:  OFF Stimulation", xlab = c("RT in ms"), breaks = 30)
abline(v = mean_OFF_Incongruent, col = "red", lwd = 3)


# save 
dev.off()

# clean up workspace
rm(RT_FLT_pred, fit_shifted_log_FLT)

# And lastly for the accuracy data

# load the model
load(file = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Modelle/log_reg_FLT.rda")

# calculate predictions
acc_FLT_pred <- brms::posterior_epred(fit_log_FLT)

# summarize over conditions
png(file="D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Posterior_Checks/acc_FLT_posterior_pred.png",
    width=18, height=12, units="cm",res = 400)
par(mfrow = c( 3, 2))

# get data for 130Hz
S130Hz_Congruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "S130Hz_congruent"])
S130Hz_Incongruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "S130Hz_incongruent"])
mean_130Hz_Congruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "S130Hz_congruent"])
mean_130Hz_Incongruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "S130Hz_incongruent"]) 

# Plot 130 Figures
hist(S130Hz_Congruent, main = "Congruent:  130Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_130Hz_Congruent, col = "red", lwd = 3)
hist(S130Hz_Incongruent, main = "Incongruent:  130Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_130Hz_Incongruent, col = "red", lwd = 3)

# get data for 20Hz
S20Hz_Congruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "S20Hz_congruent"])
S20Hz_Incongruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "S20Hz_incongruent"])
mean_20Hz_Congruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "S20Hz_congruent"])
mean_20Hz_Incongruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "S20Hz_incongruent"]) 

# Plot 20 Figures
hist(S20Hz_Congruent, main = "Congruent:  20Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_20Hz_Congruent, col = "red", lwd = 3)
hist(S20Hz_Incongruent, main = "Incongruent:  20Hz Stimulation", xlab = c("Error in %"))
abline(v = mean_20Hz_Incongruent, col = "red", lwd = 3)

# get data for OFF
SOFF_Congruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "SOFF_congruent"])
SOFF_Incongruent <- rowMeans(acc_FLT_pred[,FLTRT$Contrast_F == "SOFF_incongruent"])
mean_OFF_Congruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "SOFF_congruent"])
mean_OFF_Incongruent <- mean(FLTRT$Error[FLTRT$Contrast_F == "SOFF_incongruent"]) 

# Plot OFF Figures
hist(SOFF_Congruent, main = "Congruent:  OFF Stimulation", xlab = c("Error in %"))
abline(v = mean_OFF_Congruent, col = "red", lwd = 3)
hist(SOFF_Incongruent, main = "Incongruent:  OFF Stimulation", xlab = c("Error in %"))
abline(v = mean_OFF_Incongruent, col = "red", lwd = 3)

# save 
dev.off()

# clean up workspace
rm(acc_FLT_pred, fit_log_FLT)