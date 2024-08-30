#####
# Author: Julius Kricheldorff
# Install packagaes necessary for recreating the analysis using renv
# Date 25.07.2024
####
wd = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis"
setwd(wd)

library(brms)
library(dplyr)
library(tidyr)
library(haven)
library(hypr)
library(rstan)# package to generate contrasts
library(StanHeaders)
library(rstudioapi)
library(xtable)
library(stringr)
library(emmeans)
library(tidybayes)
library(ggtext)
library(colorspace)
library(ggplot2)
library(forcats)# so we can simply reorder the variables with fct_inorder
library(ggsci)
library(plotrix) # to calculate standard error
library(ggpubr)
library(stringr)
library(ggdist)
library(cowplot)
library(purrr)

