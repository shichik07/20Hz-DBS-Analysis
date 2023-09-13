library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(tiybayes)
library(ggtext)
library(colorspace)
library(ragg)
library(ggplot2)
library(forcats)# so we can simply reorder the variables with fct_inorder
library(ggsci)
library(plotrix) # to calculate standard error
library(ggpubr)

#colorspace::choose_palette()
sequential_hcl(4, palette = "YlGnBu")

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted')

# load data
SimpleRT<- read_csv(file = "SimpleRT.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

SimpleRT_cleanRT <- SimpleRT %>% filter(Correct_Response == 1,
                    RT < 2,
                    RT > 0.2) %>%
  mutate(RT = RT*1000)

SimpleRT_individuals <- SimpleRT_cleanRT %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(RT = median(RT))


pal_futurama("planetexpress")(4)


pal <- c("#FF8C00", "#A034F0", "#159090")
pal <- c("#26185F", "#0095AF", "#9ADCBB")
pal <- c("#FF6F00FF", "#008EA0FF", "#8A4198FF")

add_sample <- function(x){
  return(c(y = max(x) + .025, 
           label = length(x)))
}

jitter <- position_jitter(width = 0.12, seed = 123)

RT_plt <- SimpleRT_cleanRT %>% 
  group_by(Stim_verb, Part_nr) %>% 
  mutate(Median_RT = mean(RT)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = RT)) + 
  ggdist::stat_halfeye(
    aes(color = Stim_verb,
        fill = after_scale(lighten(color, .5))),
    adjust = .5, 
    width = .6, 
    .width = 0,
    justification = -.4, 
    point_color = NA) + 
  geom_boxplot(
    aes(color = Stim_verb,
        color = after_scale(darken(color, .1, space = "HLS"))),
        #fill = after_scale(desaturate(lighten(color, .8), .4))),
    width = .35, 
    outlier.shape = NA
  ) +
  geom_point(data = SimpleRT_individuals,
    aes(color = Stim_verb,
        color = after_scale(darken(color, .1, space = "HLS"))),
    fill = "white",
    stat = "identity",
    shape = 21,
    stroke = .4,
    size = 2,
    position = jitter
  ) + 
  geom_point(data = SimpleRT_individuals,
    aes(fill = Stim_verb),
    color = "transparent",
    stat = "identity",
    shape = 21,
    stroke = .4,
    size = 2,
    alpha = .3,
    position = jitter
  ) + 
  geom_path(data = SimpleRT_individuals, 
            aes(group = Part_nr),
            stat = "identity",
            color = "black",
            alpha = .15,
            position = jitter
            ) +
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = round(..y..),
        color = Stim_verb,
        color = after_scale(darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size =9,
    vjust = -2
  ) +
  coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    limits = c(200, 2000),
    breaks = seq(200, 2000, by = 400),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(
    x = NULL,
    y = "RT in milliseconds"#,
    #title = "RTs - Simple Reaction Time Task"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono", face = "bold", size = 20),
    axis.text.y = element_text(
      color = rev(darken(pal, .1, space = "HLS")), 
      size = 28,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = 25),
    plot.title = element_markdown(face = "bold", size = 21),
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

Error_plt <- SimpleRT %>%
  group_by(Stim_verb) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error, 
             color = Stim_verb,
             fill = after_scale(lighten(color, .5))),
         width = .25) + 
  geom_bar(stat = "identity", position = position_dodge(), width = .4) +
  geom_point(data = SimpleError_individual,
             aes(color = Stim_verb,
                 color = after_scale(darken(color, .1, space = "HLS"))),
             fill = "white",
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2
  ) + 
  geom_point(data = SimpleError_individual,
             aes(fill = Stim_verb),
             color = "transparent",
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2,
             alpha = .3
  ) +
  geom_path(data = SimpleError_individual, 
            aes(group = Part_nr),
            stat = "identity",
            color = "black",
            alpha = .15
  ) +
  stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = str_c(format(round(..y.., 3)*100 , nsmall =1)),
        y = stage(prop_error, after_stat = 0.3),
        color = Stim_verb,
        color = after_scale(darken(color, .1, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 9,
    vjust = -0.3,
    show_guide = FALSE
  ) +
  coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.7),
    breaks = seq(0.1, 1.0, by = .2),
    expand = c(.001, .001)
  )  +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(
    x = NULL,
    y = "Error in Percent"#,
    #title = "Error Rates - Simple Reaction Time Task"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono",face = "bold", size = 20),
    axis.text.y = element_text(
      color = rev(darken(pal, .1, space = "HLS")), 
      size = 28,
      face = "bold"
    ),
    axis.title.x = element_text(margin = margin(t = 10),
                                size = 24),
    plot.title = element_markdown(face = "bold", size = 21),
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
                alpha = 0.6, position = position_dodge(width = 0.9),
                width = 0.4)

# combine the plots

SRT_comb <- ggarrange(RT_plt, Error_plt + rremove("y.text"),
                        labels = c("RT", "Error"),
                        ncol = 2, nrow = 1, align = "h", widths = c(1,1),
                      font.label=list(color="black",size=28))

annotate_figure(SRT_comb,
                top =text_grob( "Selection Task", face = "bold", size = 32))



save_n <- "SRT.png"
save_path <- "C:/Users/doex9445/Dateien/Julius/20Hz/Figures/Poster"

ggsave(path = save_path, filename = save_n,  dpi=300,  units = "mm", height =  160, width = 380)
