#####
# Author: Julius Kricheldorff
# Analysis of the simple reaction time data
# Date 11.01.2023
####

# load packages
library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(brms)
library(hypr)
library(tidybayes)
library(bayestestR)

# set directory
wd <-"D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

# load data
SimpleRT<- read_csv(file = "SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))
# create a contrast matrix for our comparisons of interest
SRT_Contrast <- hypr(
  S130_S20 = S20Hz ~ S130Hz, 
  S20_SOFF = S20Hz ~ SOFF, 
  levels = c("S130Hz", "SOFF", "S20Hz")
)
SRT_Contrast

# assign the generated contrast matrix to the List Wide Factor
contrasts(SimpleRT$StimCon) <- contr.hypothesis(SRT_Contrast)
contrasts(SimpleRT$StimCon)


# for the RT analysis we filter only correct trials, and exclude trials shorter 
# than 200ms a longer than 2s
RT_data <- SimpleRT %>%
  filter(Correct_Response == 1,
         RT <3,
         RT > 0.2) %>%
  mutate(RT_ms = RT*1000)

RT_data2 <-  SimpleRT%>%
  filter(Correct_Response == 1)%>%
  mutate(RT_ms = RT*1000)

data_excluded <- 1-nrow(RT_data)/nrow(RT_data2)


# Prior informed weakly Item Specific
prior_weakly_informed<- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0 ,0.5), class = sigma, lb = 0),
  prior(uniform(0, min_Y), class = ndt),
  prior(normal(0,  0.3), class = b, coef = StimConS130_S20), 
  prior(normal(0,  0.3), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  0.3), class = sd, coef = Intercept, group = Part_nr)
)


# brmsformula object List Wide
m1_SRT <- bf(RT_ms ~ 1  + StimCon + (1| Part_nr))


#fit the first model
fit_shifted_log_SRT <- brm(formula = m1_SRT,
                  family = shifted_lognormal(),
                  data = RT_data,
                  prior = prior_weakly_informed,
                  warmup = 2000,
                  iter = 12000,# 20000 is the limit necessary for bridge sampling
                  cores = 4, seed = 423,
                  #control = list(adapt_delta = 0.9),
                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                  chains =4
)

save(fit_shifted_log_SRT, file = "../Modelle/shifted_log_SRT.rda")
load(file = "../Modelle/shifted_log_SRT.rda")

## Next, let us look at the accuracy data
prior_weakly_informed_logreg <- c(
  prior(normal(-2, 1), class = Intercept),
  prior(normal(0,  1.5), class = b, coef = StimConS130_S20), 
  prior(normal(0,  1.5), class = b, coef = StimConS20_SOFF),
  prior(normal(0,  1.5), class = sd, coef = Intercept, group = Part_nr)
)

# brmsformula object Item Specific
m1_SRT_log <- bf(Error ~ 1  + StimCon + (1|Part_nr)) 

#### Fit Error Models ####



fit_logReg_SRT <- brm(formula = m1_SRT_log,
                                  family = bernoulli(link = logit),
                                  data = SimpleRT,
                                  prior = prior_weakly_informed_logreg,
                                  warmup = 4000,
                                  iter = 14000,# 20000 is the limit necessary for bridge sampling
                                  cores = 4, seed = 423,
                                  save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                                  chains =4)

save(fit_logReg_SRT, file = "../Modelle/log_reg_SRT.rda")
load(file = "../Modelle/log_reg_SRT.rda")
