library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(ggtext)
library(colorspace)
library(ggplot2)
library(forcats)# so we can simply reorder the variables with fct_inorder
library(ggsci)
library(plotrix) # to calculate standard error
library(ggpubr)
library(stringr)
library(ggdist)

# set directory
wd = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

# load data
GoNoGo<- read_csv(file = "GoNoGo.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  ))) %>%
  mutate(Part_nr = as.factor(as.numeric(Part_nr)))


GoNoGo_cleanRT <- GoNoGo %>% filter(Correct_Response == 1,
                                    GoNoGo != "NoGo - Stop",
                                    RT < 3,
                                    RT > 0.2)%>%
  mutate(RT = RT*1000)

# load clinical data
Clinic<- read_csv2(file = "PartChars.csv") %>%
  select(ProbandALT, Alter, DauerPD, Krankheitstyp, Dominanzseite, HnY, DauerDBS, LEDD) %>%
  rename(Part_nr = ProbandALT,
         Age = Alter,
         Dur_DBS = DauerDBS,
         DiseaseSubtype = Krankheitstyp,
         DomSide = Dominanzseite,
         Dur_PD = DauerPD) %>%
  mutate(Part_nr = as.factor(Part_nr),
         DiseaseSubtype = as.factor(DiseaseSubtype),
         DomSide = as.factor(DomSide),
         HnY = as.factor(HnY)) %>%
  mutate(new_DiseaseSubType = recode_factor(DiseaseSubtype, "Ä"="ET", "T"="TD", "HR"="AR"),
         new_DomSide = recode_factor(DomSide, 'L' = "left", 'R'="right"))




GNG_individuals_RT <- GoNoGo_cleanRT %>%
  group_by(Stim_verb, GoNoGo, Part_nr) %>%
  summarise(RT = mean(RT)) %>%
  mutate(id = as_factor(case_when(
    GoNoGo == "Go" ~ as.numeric(Part_nr)+21,
    GoNoGo == "NoGo - Stop" ~ as.numeric(Part_nr),
    GoNoGo == "NoGo - Go" ~ as.numeric(Part_nr) +41
  ))) %>%
  ungroup()

# now join the clinical data
Outcome_RT <- GNG_individuals_RT%>%
  group_by(Stim_verb, Part_nr) %>%
  reframe(diff_RT = RT[GoNoGo == "NoGo - Go"] - RT[GoNoGo == "Go"])

# now calculate the difference for OFF vs 20Hz and 130Hz vs 20Hz
Outcome_RT_20_v_OFF <- Outcome_RT %>% 
  group_by(Part_nr) %>%
  reframe(diff_RT = diff_RT[Stim_verb == "20Hz"] - diff_RT[Stim_verb == "OFF"]) %>%
  mutate(comparison = "20_v_OFF")

Outcome_RT_all <- Outcome_RT %>% 
  group_by(Part_nr) %>%
  reframe(diff_RT = diff_RT[Stim_verb == "20Hz"] - diff_RT[Stim_verb == "130Hz"]) %>%
  mutate(comparison = "20_v_130") %>%
  bind_rows(Outcome_RT_20_v_OFF) %>%
  full_join(Clinic, by="Part_nr")

# Plot by LEDD
# now let us plot as a test vs LEDD
sp <- ggscatter(Outcome_RT_all, x = "diff_RT", y = "LEDD",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  labs(y = "LEDD in mg", x= "RT slowing in ms") +
  ggtitle("Association LEED and RT Slowing in ms")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 20, label.y = 25)

# now let us plot as a test vs Disease Duration
sp <- ggscatter(Outcome_RT_all, x = "diff_RT", y = "Dur_PD",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  labs(y = "Disease Duration in Years", x= "RT slowing in ms") +
  ggtitle("Association Disease Duration and RT Slowing in ms")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 20, label.y = 25)

# now let us plot as a test vs Duration DBS
sp <- ggscatter(Outcome_RT_all, x = "diff_RT", y = "Dur_DBS",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  labs(y = "Time since implantation in Years", x= "RT slowing in ms") +
  ggtitle("Association Time since implantation and RT Slowing in ms")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 20, label.y = 25)

# now let us plot as a test vs Age
sp <- ggscatter(Outcome_RT_all, x = "diff_RT", y = "Age",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  labs(y = "Age in Years", x= "RT slowing in ms") +
  ggtitle("Association Age and RT Slowing in ms")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 20, label.y = 25)

## Next see if there is a difference by affected sides
sides <- ggplot(data = Outcome_RT_all, 
                aes(x = fct_inorder(new_DomSide), y = diff_RT, color = fct_inorder(new_DomSide)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  labs(x = "Affected Side", y = "RT slowing in ms in ms") +
  ggtitle("Affected Side and RT Slowing in ms") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Affected Side",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

## Next see if there is a difference by disease subtype
DST <- ggplot(data = Outcome_RT_all, 
                aes(x = fct_inorder(new_DiseaseSubType), y = diff_RT, color = fct_inorder(new_DiseaseSubType)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  labs(x = "Disease Subtype", y = "RT slowing in ms in ms") +
  ggtitle("Disease Subtype and RT Slowing in ms") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Disease Subtype",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

## Next see if there is a difference by HNY
HNY <- ggplot(data = Outcome_RT_all, 
              aes(x = fct_inseq(HnY), y = diff_RT, color = fct_inorder(HnY)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  labs(x = "Hoehn and Yahr", y = "RT slowing in ms in ms") +
  ggtitle("Hoehn and Yahr and RT Slowing in ms") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Hoehn and Yahr",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

GNG_individual_err <- GoNoGo %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) 

# now calculate the difference for OFF vs 20Hz and 130Hz vs 20Hz
Outcome_err_20_v_OFF <- GNG_individual_err %>% 
  group_by(Part_nr) %>%
  reframe(diff_err = prop_error[Stim_verb == "20Hz"] - prop_error[Stim_verb == "OFF"]) %>%
  mutate(comparison = "20_v_OFF")

Outcome_err_all <- GNG_individual_err %>% 
  group_by(Part_nr) %>%
  reframe(diff_err = prop_error[Stim_verb == "20Hz"] - prop_error[Stim_verb == "130Hz"]) %>%
  mutate(comparison = "20_v_130") %>%
  bind_rows(Outcome_err_20_v_OFF) %>%
  full_join(Clinic, by="Part_nr")

# Plot by LEDD
# now let us plot as a test vs LEDD
sp <- ggscatter(Outcome_err_all, x = "diff_err", y = "LEDD",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.5, .1, by = .2),
    expand = c(.001, .001)
  ) +
  labs(y = "LEDD in mg", x= "Error Reduction in %") +
  ggtitle("Association LEED and Error Reduction")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = -0.3, label.y = 1500)

# now let us plot as a test vs Disease Duration
sp <- ggscatter(Outcome_err_all, x = "diff_err", y = "Dur_PD",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.5, .1, by = .2),
    expand = c(.001, .001)
  ) +
  labs(y = "Disease Duration in Years", x= "Error Reduction in %") +
  ggtitle("Association Disease Duration and Error Reduction in %")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = -0.3, label.y = 2)

# now let us plot as a test vs Duration DBS
sp <- ggscatter(Outcome_err_all, x = "diff_err", y = "Dur_DBS",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.5, .1, by = .2),
    expand = c(.001, .001)
  ) +
  labs(y = "Time since implantation in Years", x= "Error Reduction in %") +
  ggtitle("Association Time since implantation and Error Reduction in %")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = -0.3, label.y = 2)

# now let us plot as a test vs Age
sp <- ggscatter(Outcome_err_all, x = "diff_err", y = "Age",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval,
) + theme_bw(base_size = 14) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  scale_x_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.5, .1, by = .2),
    expand = c(.001, .001)
  ) +
  labs(y = "Age in Years", x= "Error Reduction in %") +
  ggtitle("Association Age and Error Reduction in %")
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = -0.3, label.y = 45)

## Next see if there is a difference by affected sides
sides <- ggplot(data = Outcome_err_all, 
                aes(x = fct_inorder(new_DomSide), y = diff_err, color = fct_inorder(new_DomSide)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.4, .2, by = .2),
    expand = c(.001, .001)
  ) +
  labs(x = "Affected Side", y = "Error Reduction in %") +
  ggtitle("Affected Side and Error Reduction in %") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Affected Side",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

## Next see if there is a difference by disease subtype
DST <- ggplot(data = Outcome_err_all, 
              aes(x = fct_inorder(new_DiseaseSubType), y = diff_err, color = fct_inorder(new_DiseaseSubType)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.4, .2, by = .2),
    expand = c(.001, .001)
  ) +
  labs(x = "Disease Subtype", y = "Error Reduction in %") +
  ggtitle("Disease Subtype and Error Reduction in %") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Disease Subtype",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

## Next see if there is a difference by HNY
HNY <- ggplot(data = Outcome_err_all, 
              aes(x = fct_inseq(HnY), y = diff_err, color = fct_inorder(HnY)))  +
  geom_point(size=0.8, show.legend =FALSE) +
  facet_grid(cols = vars(fct_inorder(comparison))) +
  geom_boxplot(alpha = 0.4, show.legend =FALSE) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(-0.5, .1),
    breaks = seq(-0.4, .2, by = .2),
    expand = c(.001, .001)
  ) +
  labs(x = "Hoehn and Yahr", y = "Error Reduction in %") +
  ggtitle("Hoehn and Yahr and Error Reduction in %") + 
  theme_minimal(base_family = "Zilla Slab", base_size = 8)+
  guides(color=guide_legend(title = "Hoehn and Yahr",
                            override.aes = list(shape = 15,
                                                size = 3))) +
  theme(
    panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.text = element_text(fac
                               e = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face="bold"),
    plot.title = element_text(face = "bold"),
    plot.title.position  =  "panel"
  )

