#####
# Author: Julius Kricheldorff
# Analysis of the Stop-Signal task
# Date 26.01.2023
####

# Address git : C/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(hypr)
library(tidybayes)
library(emmeans)
library(purrr)
library(BayesFactor)

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted')

# load data
SST<- read_csv(file = "StoppSignal.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

SST <- SST %>%
  mutate(Contrast_F = as_factor(case_when(
    Stim_verb == "130Hz" & StopTrl == "Go" ~ "S130Hz_Go",
    Stim_verb == "130Hz" & StopTrl == "Stop" ~ "S130Hz_Stop",
    Stim_verb == "20Hz" & StopTrl == "Go" ~ "S20Hz_Go",
    Stim_verb == "20Hz" & StopTrl == "Stop" ~ "S20Hz_Stop",
    Stim_verb == "OFF" & StopTrl == "Go" ~ "SOFF_Go",
    Stim_verb == "OFF" & StopTrl == "Stop" ~ "SOFF_Stop",
  )))

SST_Error <- SST %>%
  filter(StopTrl == "Go")

# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- SST %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2,
         StopTrl == "Go") %>%
  mutate(RT_ms = RT*1000)

RT_data2 <-  SST %>%
  filter(Correct_Response == 1,
         StopTrl == "Go") %>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)




# Prior informed weakly Item Specific
prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.1), class = b, coef = StimConS20Hz), 
  prior(normal(0,  0.1), class = b, coef = StimConSOFF),
  prior(normal(0,  0.1), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_SST <- bf(RT_ms ~ 1  + StimCon + (StimCon|Part_nr))

#get_prior(formula = m1_SST, data = RT_data, family = shifted_lognormal())

#fit the first model
fit_shifted_log_SST <- brm(formula = m1_SST,
                           family = shifted_lognormal(),
                           data = RT_data,
                           prior = prior_weakly_informed,
                           warmup = 200,
                           iter = 1200,# 20000 is the limit necessary for bridge sampling
                           cores = 4, seed = 423,
                           #control = list(adapt_delta = 0.9),
                           save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                           chains =4
)

save(fit_shifted_log_SST, file = "E:/20Hz/Data/Modelle/shifted_log_SST.rda")
load(file = "E:/20Hz/Data/Modelle/shifted_log_SST.rda")

# Next check how well the posteriors fit the actual data

# posteriro predictive checks
pp_check(fit_shifted_log_SST, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_shifted_log_SST, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_shifted_log_SST, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_shifted_log_SST, ndraws = 1000, type = "stat", stat = "mean")
# looks really good
pp_check(fit_shifted_log_SST, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "StimCon")
# somewhat resonable
pp_check(fit_shifted_log_SST, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "Part_nr")

# Lastly use emmeans to get the contrasts
prior_sum <- prior_summary(fit_shifted_log_SST)

warp_em2 <- pairs(emmeans(fit_shifted_log_SST, ~ StimCon))

# we need to rework the brms object to a grid object
prior_mod <- unupdate(fit_shifted_log_SST)
prior_emmgrid <- emmeans(prior_mod, ~ StimCon)

# then we can estimate the Bayesfactor
bayesfactor_parameters(cont2, prior = prior_emmgrid)


## Next, let us look at the accuracy data

prior_weakly_informed_logreg <- c(
  prior(normal(-0.6, 0.6), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = StimConS20Hz), 
  prior(normal(0,  1.5), class = b, coef = StimConSOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)


#get_prior(formula = m1_SRT_log, data = RT_data, family = bernoulli(link = logit))

# brmsformula object Item Specific
m1_SRT_log <- bf(Error ~ 1  + StimCon + (StimCon|Part_nr)) 

#### Fit Inducer Models ####


# we should consider varying non-decision times between the groups
fit_logReg_SST <- brm(formula = m1_SRT_log,
                      family = bernoulli(link = logit),
                      data = SST_Error,
                      prior = prior_weakly_informed_logreg,
                      warmup = 2000,
                      iter = 12000,# 20000 is the limit necessary for bridge sampling
                      cores = 4, seed = 423,
                      save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                      chains =4)

save(fit_logReg_SST, file = "E:/20Hz/Data/Modelle/logreg_SST.rda")
load(file = "E:/20Hz/Data/Modelle/logreg_SST.rda")
# posteriro predictive checks
pp_check(fit_logReg_SST, ndraws = 11, type = "hist")
# Looks good for this model
pp_check(fit_logReg_SST, ndraws = 100, type = "dens_overlay")
# Looks good for this model
pp_check(fit_logReg_SST, type = "boxplot", ndraws = 10)
# A few observation outside, but our particpants had a deadline - cannot be modelled
pp_check(fit_logReg_SST, ndraws = 1000, type = "stat", stat = "mean")
# looks good
pp_check(fit_logReg_SST, ndraws = 1000, type = "stat_grouped", stat = "mean", group = "StimCon")
# looks good

warp_em2 <- emmeans(fit_logReg_SST, ~ StimCon, epred = TRUE)
cont2 <- contrast(regrid(warp_em2, transform = "response"), interaction = "pairwise")
summary(cont2, type = "response", point.est = mean)


# In summary, we do not see any difference inr eaction times on Go trials, which is not really surprising

# Next, let us go for the accuracy rates on Go trials (Stop trials should be at a solid 50%)


#next let us try the stop signal integration method
library(splithalfr)
data("ds_sst", package = "splithalfr")
?ds_sst

# get data of the first participant in the 130Hz condition

fn_scoretidy <- function(.data) {
  # Mean SSD
  mean_ssd <- .data %>% filter(StopTrl == "Stop") %>%
    summarise(meanssd = mean(delay)) %>%
    pull(meanssd)
  # Proportion of failed nogos
  p_failed_nogo <- .data %>% filter(StopTrl == "Stop") %>%
    summarise(meanerr = mean(Error)) %>%
    pull(meanerr)
  # Go RTs
  go_rts <- .data %>% filter(StopTrl != "Stop",
                                      RT > 0)
  # n-th percentile of Go RTs
  rt_quantile <- go_rts %>%
    summarise(rt_quantile = quantile(RT, p_failed_nogo, names = FALSE)*1000) %>%
    pull(rt_quantile)
  # SSRTi
  return(round(rt_quantile - mean_ssd))
}

SST_filtered <- SST %>%
  filter(RT < 3)

res <- SST_filtered %>% 
  group_by(Part_nr, Stim_verb) %>%
  group_modify(~ tibble::enframe(fn_scoretidy(.), name = NULL)) %>%
  ungroup()

mean_res <- res %>%
  group_by(Stim_verb) %>%
  summarise(SCT = round(mean(value), 2))

SCT_20Hz <- res %>%
  filter(Stim_verb == "20Hz") %>%
  arrange(Part_nr)
SCT_130Hz <- res %>%
  filter(Stim_verb == "130Hz")%>%
  arrange(Part_nr)
SCT_Off <- res %>%
  filter(Stim_verb == "OFF")%>%
  arrange(Part_nr)
SCT_20Hz$value - SCT_130Hz$value

# post-hoc t-test
BF_SCT_20_130 <- ttestBF(SCT_20Hz$value - SCT_130Hz$value)
BF_SCT_20_OFF <- ttestBF(SCT_20Hz$value - SCT_Off$value)



# Get the corresponding p-value
t.test(SCT_20Hz$value - SCT_130Hz$value, alternative = "two.sided")
t.test(SCT_20Hz$value - SCT_Off$value, alternative = "two.sided")




