#####
# Author: Julius Kricheldorff
# Calculate Bayes Factors for contrasts of interest in the Go-NoGo task
# Date 26.09.2023
####

# Load packages
library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)
library(rstudioapi)

# Set a seed for sake of reproducibility
set.seed(32936)

setwd('C:/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Data/Extracted')

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
  mutate(RT_ms = RT*1000)

# Define formulas so we can loop through them
GNG_formulas <- c(
  formula((S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3), # main effect Go effects
  formula((S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (S130Hz_NoGo_Go + S130Hz_Go)/2), # Overall effect LFS vs HFS
  formula((S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (SOFF_NoGo_Go + SOFF_Go)/2), # Overall effect LFS vs OFF
  formula((S20Hz_NoGo_Go - S20Hz_Go) ~ (S130Hz_NoGo_Go - S130Hz_Go)), # Difference Go effects LFS vs HFS
  formula((S20Hz_NoGo_Go - S20Hz_Go) ~ (SOFF_NoGo_Go - SOFF_Go)) # Difference Go effects LFS vs OFF
)


GoNoGo_Contrast_RT <- hypr(
  Go_NoGo_Go = (S130Hz_Go + SOFF_Go + S20Hz_Go)/3 ~ (S130Hz_NoGo_Go + SOFF_NoGo_Go + S20Hz_NoGo_Go)/3,
  Stim_20v130 = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (S130Hz_NoGo_Go + S130Hz_Go)/2,
  Stim_20vOFF = (S20Hz_NoGo_Go + S20Hz_Go)/2 ~ (SOFF_NoGo_Go + SOFF_Go)/2,
  GoDiff_20v130 = (S20Hz_NoGo_Go - S20Hz_Go) ~ (S130Hz_NoGo_Go - S130Hz_Go),
  GoDiff_20vOFF = (S20Hz_NoGo_Go - S20Hz_Go) ~ (SOFF_NoGo_Go - SOFF_Go),
  levels = c("S130Hz_NoGo_Stop", "S130Hz_NoGo_Go", "S130Hz_Go",
             "SOFF_Go", "SOFF_NoGo_Stop", "SOFF_NoGo_Go", 
             "S20Hz_NoGo_Stop", "S20Hz_NoGo_Go", "S20Hz_Go")
)


# contrast names seperately
GNG_contrast_names <- c(
  "Go_NoGo_Go",
  "Stim_20v130",
  "Stim_20vOFF",
  "GoDiff_20v130",
  "GoDiff_20vOFF"
)

# LW levels
GNG_levels <- c("S130Hz_NoGo_Stop", "S130Hz_NoGo_Go", "S130Hz_Go",
                "SOFF_Go", "SOFF_NoGo_Stop", "SOFF_NoGo_Go", 
                "S20Hz_NoGo_Stop", "S20Hz_NoGo_Go", "S20Hz_Go")



BF_mods <- c("GNG_min_Go_NoGo_Go", 
             "GNG_min_Stim_20v130",
             "GNG_min_Stim_20vOFF",
             "GNG_min_GoDiff_20v130",
             "GNG_min_GoDiff_20vOFF")


# First let us get the listwide model contrasts
for(mods in 1:length(BF_mods)){
  # first create the contrast matrix
  temp_mat <- hypr(GNG_formulas[-mods],
                   levels = GNG_levels)
  # next add the appropriate variable names
  names(temp_mat) <- GNG_contrast_names[-mods]
  
  # Lastly rename the contrast matrix
  assign(BF_mods[mods], temp_mat)
}


# Now we create an extra column for the contrasts for the List wide effect 
Data <- Data %>%
  mutate(Contrast_LW = case_when(
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "CO"  ~ "CO_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_LW_incon_I",
    Analysis_type == "MC" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_con_C",
    Analysis_type == "MC" & Congruency == "incongruent" & Exp_group == "PD"  ~ "PD_LW_con_I",
    Analysis_type == "MI" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_LW_incon_C",
    Analysis_type == "MI" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_LW_incon_I",
  )) %>%
  mutate(Contrast_LW  = as_factor(Contrast_LW)) 

# Now we create an extra column for the contrasts for the Item Specific effect
Data <- Data %>%
  mutate(Contrast_IS = case_when(
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "CO" ~ "CO_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "CO" ~ "CO_IS_incon_I",
    Analysis_type == "main_con" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_con_C",
    Analysis_type == "main_con" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_con_I",
    Analysis_type == "main_incon" & Congruency == "congruent"  & Exp_group == "PD" ~ "PD_IS_incon_C",
    Analysis_type == "main_incon" & Congruency == "incongruent"  & Exp_group == "PD" ~ "PD_IS_incon_I",
  )) %>%
  mutate(Contrast_IS  = as_factor(Contrast_IS))


# Filter only young participants, remove practice trials and only keep correct trials
Data_y_inducer_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type != "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group)) 

Data_y_diagnostic_LWPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "List_wide") %>% # only Listwide items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_inducer_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # only Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type != "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

Data_y_diagnostic_ISPC <- Data %>%
  dplyr::filter(Exp_group != "CY") %>% # Control and old
  dplyr::filter(Item_specific == "Item_spec") %>% # only Item_specific items
  dplyr::filter(Block != 'practice') %>% # no practice trials
  dplyr::filter(Trl_type == "Diagnostic")%>%
  dplyr::filter(Correct_A == 1) %>% # only correct responses
  mutate(Exp_group  = as_factor(Exp_group))

# Same priors as for the model before

# Prior informed weakly List wide
prior_weakly_informed_LW <- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0), # group effect non-decision time CO is 0 in the contrast we assume on average 200ms
  prior(normal(0,0.5), class = sigma, lb = 0), 
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_Listwide),
  prior(normal(0,0.3), class = b, coef = Contrast_LWCO_LW_Block),
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_Listwide),
  prior(normal(0,0.3), class = b, coef = Contrast_LWPD_LW_Block),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0, 0.3), class = sd, coef = Intercept, group = Item)
)


# Prior informed weakly Item Specific
prior_weakly_informed_IS <- c(
  prior(normal(6.5, 0.5), class = Intercept, lb = 0),
  prior(normal(0,0.5), class = sigma, lb = 0), 
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupCO, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(5.3, 0.5), class = b, coef = Exp_groupPD, dpar = ndt), # group effect non-decision time CO is 5.3 in the contrast we assume on average 200ms
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_Congruency), # Priors of Contrasts for CO
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_Itemspecific),
  prior(normal(0,0.3), class = b, coef = Contrast_ISCO_IS_Block),
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_Congruency), # Priors of Contrasts for PD patients
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_Itemspecific),
  prior(normal(0,0.3), class = b, coef = Contrast_ISPD_IS_Block),
  prior(normal(0,0.3), class = sd, coef = Intercept, group = Subject),
  prior(normal(0,0.3), class = sd, coef = Intercept, group = Item)
)

# brmsformula object List Wide
m1_LW <- bf(RT ~ 1  + Contrast_LW + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group)

# brmsformula object Item Specific
m1_IS <- bf(RT ~ 1  + Contrast_IS + (1|Subject) + (1|Item),
            ndt ~  0 + Exp_group) 

# okay create a function to pass model parameter
pass_brms = function(save_name, prior, data, model) {
  fit_model <- brm(formula = model,
                   family = shifted_lognormal(),
                   data = data,
                   prior = prior,
                   warmup = 2000,
                   iter = 12000,# 20000 is the limit necessary for bridge sampling
                   cores = 4, seed = 423,
                   control = list(adapt_delta = 0.95),
                   save_pars = save_pars(all = TRUE), # must be set to true for bridgesampling
                   chains =4)
  
  # save the model
  save(fit_model, file = save_name)
}


##### now let us loop though our models and save the results for the listwide models

for(mods in 1:length(LW_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_LWPC$Contrast_LW) <- contr.hypothesis(eval(parse(text = LW_mods[mods])))
  # define prior
  LW_weakly <- prior_weakly_informed_LW[-(mods+4),]
  # get save name for variable
  save_name <- paste("fit_shifted_inducer_", LW_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = LW_weakly, data = Data_y_inducer_LWPC, model = m1_LW)
} 

# Now the same with the item specific models
for(mods in 1:length(IS_mods)){
  # first assign appropriate contrasts to our dataset 
  contrasts(Data_y_inducer_ISPC$Contrast_IS) <- contr.hypothesis(eval(parse(text = IS_mods[mods])))
  # define prior
  IS_weakly <- prior_weakly_informed_IS[-(mods+4),]
  # get save name for variable
  save_name <- paste("fit_shifted_inducer_", IS_mods[mods], ".rda", sep = "")
  # fit the model
  pass_brms(save_name = save_name, prior = IS_weakly, data = Data_y_inducer_ISPC, model = m1_IS)
  
}