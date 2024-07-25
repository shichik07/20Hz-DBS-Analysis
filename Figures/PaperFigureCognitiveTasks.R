library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
#(tidybayes)
library(ggtext)
library(colorspace)
#library(ragg)
library(ggplot2)
library(forcats)# so we can simply reorder the variables with fct_inorder
library(ggsci)
library(plotrix) # to calculate standard error
library(ggpubr)
library(stringr)
library(ggdist)


# set directory

wd = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
#wd = "/home/jules/Dropbox/PhD_Thesis/UniOL/Julius/20Hz/Data/Extracted"
#wd = "C:/Users/doex9445/Dropbox/PhD_Thesis/UniOL/Julius/20Hz/Data/Extracted"
setwd(wd)

sequential_hcl(4, palette = "YlGnBu")

#vars
dot_size <- 0.5
text_size <- 8
guide_size <- 2
move_sum_by <- -5.4


# load data
SimpleRT<- read_csv(file = "SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

SimpleRT_cleanRT <- SimpleRT %>% filter(Correct_Response == 1,
                                        RT < 3,
                                        RT > 0.2) %>%
  mutate(RT = RT*1000)

SimpleRT_individuals <- SimpleRT_cleanRT %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(RT = median(RT))


pal_futurama("planetexpress")(4)

# Colors
pal1 <- c("#FF6F00FF", "#008EA0FF", "#8A4198FF")
#pal <- c( "#4DBBD5FF", "#E64B35FF")
pal <- c("#253494", "#56B4E9")


SRT_RT_plt <- SimpleRT_cleanRT %>% 
  group_by(Stim_verb, Part_nr) %>% 
  mutate(Median_RT = mean(RT)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = RT)) + 
  ggdist::stat_halfeye(
    aes(fill = Stim_verb,
        fill = after_scale(lighten(fill, .5))),
    adjust = .8, 
    width = .5, 
    .width = 0,
    justification = -.5, 
    point_color = NA,
    position = position_dodge(width = 1)
  )+
  geom_boxplot(#aes(color = Stim_verb,
                  # ),
    width = .35, 
    outlier.shape = NA,
    show_guide = FALSE
  ) +
  geom_point(data = SimpleRT_individuals,
             aes(fill = Stim_verb,
                 fill = after_scale(darken(fill, .1, space = "HLS"))),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = dot_size,
             alpha = .6,
             position = position_dodge(width = 1),
             show_guide = FALSE
  ) +
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(..y..),
        color = Stim_verb,
        fill = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by-3,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  #coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    limits = c(200, 2000),
    breaks = seq(200, 2000, by = 400),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal1, guide = "none") +
  scale_fill_manual(values = pal1, guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0, 
                                                 label.position = "bottom"))
  ) +
  labs(
    x = NULL,
    y = "RT in ms",
    title = "Response Selection Task",
    subtitle = "RT"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(family = "Roboto Mono",size = text_size),
    axis.text.x = element_text(
      color = rev(darken(pal1, .1, space = "HLS")), 
      size = text_size,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = text_size),
   # plot.title = element_markdown(face = "bold", size = text_size),
    plot.title = element_text(face = "bold", size = text_size,margin = margin(0, 0, 3, 0)),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.caption = element_markdown(
      color = "grey40", lineheight = 1.2,
      margin = margin(20, 0, 0, 0)),
    plot.margin = margin(15, 15, 10, 15)
  )

## Now Next Let us plot the Error Rates

SimpleError_individual <- SimpleRT %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1))

SRT_Error_plt <- SimpleRT %>%
  group_by(Stim_verb) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error 
             #color = Stim_verb,
             #fill = after_scale(lighten(color, .5))
             ),
         width = .25) + 
  geom_bar(aes(color = Stim_verb,
               fill = after_scale(lighten(color, .5))), stat = "identity", position = position_dodge(), width = .5) +
  geom_point(data = SimpleError_individual,
             aes(fill = Stim_verb),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = dot_size,
             alpha = 0.6,
             position = position_dodge(width = 0.8),
             show_guide = FALSE
  )  +
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = str_c(format(round(..y.., 3)*100 , nsmall =1)),
        y = stage(prop_error, after_stat = 0.3),
        color = Stim_verb,
        color = after_scale(darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by+2,
    #hjust = -0.7,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) + 
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.63),
    breaks = seq(0.1, 1.0, by = .2),
    expand = c(.001, .001)
  )  +
  scale_fill_manual(values = pal1,  guide = "none") +
  scale_color_manual(values = pal1, name="Stimulation\nCondition")+
  guides(color = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0))
  ) +
  labs(
    x = NULL,
    y = "Error in %",
    #title = "Choice Selection Task",
    subtitle = "Error",
    fill = "Stimulation"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    legend.title = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    legend.key.size = unit(0.3, 'cm'),
    #legend.position="bottom",
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(family = "Roboto Mono", size = text_size),
    axis.text.x = element_text(
      color = rev(darken(pal1, .1, space = "HLS")), 
      size = text_size,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = text_size),
    plot.title = element_markdown(face = "bold", size = text_size),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.caption = element_markdown(
      color = "grey40", lineheight = 1.2,
      margin = margin(20, 0, 0, 0)),
    plot.margin = margin(15, 15, 10, 15)
  ) +
  geom_errorbar(aes(ymin = prop_error - se, ymax = prop_error + se), 
                alpha = 0.5, position = position_dodge(width = 0.8),
                width = 0.5)

### Now the Flanker data
# load data
Flanker<- read_csv(file = "flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

#vars
dot_size <- 0.5
text_size <- 8
guide_size <- 2
move_sum_by <- -6.4

Flanker_cleanRT <- Flanker %>% filter(Correct_Response == 1,
                                      RT < 3,
                                      RT > 0.2)%>%
  mutate(RT = RT*1000)

Flanker_individuals <- Flanker_cleanRT %>%
  group_by(Stim_verb, Congruency, Part_nr) %>%
  summarise(RT = median(RT)) %>%
  mutate(id = as_factor(case_when(
    Congruency == "incongruent" ~ as.numeric(Part_nr)+21,
    Congruency == "congruent" ~ as.numeric(Part_nr)
  ))) 

FLT_RT_plt <- Flanker_cleanRT %>% 
  group_by(Stim_verb, Part_nr) %>% 
  ggplot(aes(x = fct_rev(Stim_verb), y = RT, fill = Congruency)) + 
  ggdist::stat_halfeye(
    aes(fill = Congruency,
        fill = after_scale(lighten(fill, .5))),
    adjust = .8, 
    width = .5, 
    .width = 0,
    justification = -.5, 
    point_color = NA,
    position = position_dodge(width = 1)
  ) + 
  geom_boxplot(
    aes(fill = Congruency,
        fill = after_scale(lighten(fill, 1, space = "HLS"))),
    width = .3,
    outlier.shape = NA,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  geom_point(data = Flanker_individuals,
             aes(fill = Congruency,
                 fill = after_scale(darken(fill, .1, space = "HLS"))),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = dot_size,
             alpha = .6,
             position = position_dodge(width = 1),
             show_guide = FALSE
  ) +
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(..y..),
        color = Congruency,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by-3,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  #coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    limits = c(200, 2000),
    breaks = seq(200, 2000, by = 400),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, labels = c("Congruent", "Incongruent" )) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0, 
                                                 label.position = "bottom"))
  ) +
  labs(
    x = NULL,
    y = "RT in ms",
    title = "Flanker Task",
    subtitle = "RT"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(family = "Roboto Mono",size = text_size),
    axis.text.x = element_text(
      color = rev(darken(pal1, .1, space = "HLS")), 
      size = text_size,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = text_size),
    plot.title = element_text(face = "bold", size = text_size,margin = margin(0, 0, 3, 0)),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.caption = element_markdown(
      color = "grey40", lineheight = 1.2,
      margin = margin(20, 0, 0, 0)),
    plot.margin = margin(15, 15, 10, 15)
  )

## Now Next Let us plot the Error Rates

SimpleError_individual <- Flanker %>%
  group_by(Stim_verb, Congruency, Part_nr) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  mutate(id = as_factor(case_when(
    Congruency == "incongruent" ~ as.numeric(Part_nr)+21,
    Congruency == "congruent" ~ as.numeric(Part_nr)
  ))) 

FLT_Error_plt <- Flanker %>%
  group_by(Stim_verb, Congruency) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error, fill = Congruency),
         #fill = Stim_verb,
         #fill = after_scale(desaturate(lighten(fill, .8), .4))),
         #color = Stim_verb,
         #color = after_scale(darken(color, .1, space = "HLS")),
         width = .25) + 
  geom_bar(aes(fill = Congruency, 
               fill = after_scale(lighten(fill, .5))),
           stat = "identity", 
           position = position_dodge(width = 0.8),
           width = 0.5
  ) +
  geom_point(data = SimpleError_individual,
             aes(fill = Congruency,
                 fill = after_scale(darken(fill, .1, space = "HLS"))),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = dot_size,
             alpha = 0.6,
             position = position_dodge(width = 0.8),
             show_guide = FALSE
  ) + 
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = str_c(format(round(..y.., 3)*100 , nsmall =1)),
        y = stage(prop_error, after_stat = 0.3),
        color = Congruency,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by+2,
    #hjust = -0.7,
    position = position_dodge(width = 1),
    show_guide = FALSE
  )  +
  #coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.62),
    breaks = seq(0.1, 1.0, by = .2),
    expand = c(.001, .001)
  )  +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, labels = c("Congruent", "Incongruent" )) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0))
  ) +
  labs(
    x = NULL,
    y = "Error in %",
    #title = "Flanker Task",
    subtitle = "Error"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    legend.title = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    legend.key.size = unit(0.3, 'cm'),
    #legend.position="bottom",
    panel.grid.minor = element_blank(),
    #panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(family = "Roboto Mono", size = text_size),
    axis.text.x = element_text(
      color = rev(darken(pal1, .1, space = "HLS")), 
      size = text_size,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = text_size),
    plot.title = element_markdown(face = "bold", size = text_size),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ),
    plot.title.position = "plot",
    plot.caption = element_markdown(
      color = "grey40", lineheight = 1.2,
      margin = margin(20, 0, 0, 0)),
    plot.margin = margin(15, 15, 10, 15)
  ) +
  geom_errorbar(aes(ymin = prop_error - se, ymax = prop_error + se), 
                alpha = 0.5, position = position_dodge(width = 0.8),
                width = 0.5)

# combine the plots

SRT_comb <- ggarrange(SRT_RT_plt + rremove("legend"), SRT_Error_plt, FLT_RT_plt + rremove("legend"),
                      FLT_Error_plt,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, align = "h", widths = c(3,4),
                      font.label=list(color="black",size=10))

save_n <- "Flanker_fin.tiff"
#save_path <- "C:/Users/doex9445/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Paper"
#save_path <- "/home/jules/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Paper"
save_path <- "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Paper"
ggsave(path = save_path, filename = save_n,  dpi=600,  units = "mm", height =  110, width = 160)
