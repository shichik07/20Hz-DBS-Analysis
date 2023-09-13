library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(tidybayes)
library(ggtext)
library(colorspace)
library(ragg)
library(ggplot2)
library(forcats)# so we can simply reorder the variables with fct_inorder
library(ggsci)
library(plotrix) # to calculate standard error

pal_npg("nrc")(2)
hcl_palettes(palette = "reds", plot = TRUE)
sequential_hcl(8,"reds")

# set directory
setwd('C:/Users/doex9445/Dateien/Julius/20Hz/Data/Extracted')

# load data
GoNoGo<- read_csv(file = "GoNoGo.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

GoNoGo_cleanRT <- GoNoGo %>% filter(Correct_Response == 1,
                                    GoNoGo != "NoGo - Stop",
                                      RT < 2,
                                      RT > 0.2)%>%
  mutate(RT = RT*1000)

GoNoGo_individuals <- GoNoGo_cleanRT %>%
  group_by(Stim_verb, GoNoGo, Part_nr) %>%
  summarise(RT = median(RT)) %>%
  mutate(id = as_factor(case_when(
    GoNoGo == "Go" ~ as.numeric(Part_nr)+21,
    GoNoGo == "NoGo - Stop" ~ as.numeric(Part_nr),
    GoNoGo == "NoGo - Go" ~ as.numeric(Part_nr) +41
  ))) 

# Colors
pal1 <- c("#FF6F00FF", "#008EA0FF", "#8A4198FF")
pal <- c( "#0072B5FF", "#EE4C97FF")
pal2 <- c( "#0072B5FF", "#EE4C97FF", "#BC3C29FF")
#pal2 <- c("#0072B5FF", "#A30034", "#6D0026")
pal2 <- c( "#0072B5FF", "#EE4C97FF", "#A30034")


RT_plt <- GoNoGo_cleanRT %>% 
  group_by(Stim_verb, Part_nr) %>% 
  ggplot(aes(x = fct_rev(Stim_verb), y = RT, fill = GoNoGo)) + 
  ggdist::stat_halfeye(
    aes(fill = GoNoGo,
        fill = after_scale(lighten(fill, .5))),
    adjust = .8, 
    width = .5, 
    .width = 0,
    justification = -.5, 
    point_color = NA,
    position = position_dodge(width = 1)
  ) + 
  geom_boxplot(
    aes(fill = GoNoGo,
        fill = after_scale(lighten(fill, 1, space = "HLS"))),
    width = .3, 
    outlier.shape = NA,
    position = position_dodge(width = 1), 
    show_guide = FALSE
  ) +
  geom_point(data = GoNoGo_individuals %>%
               filter(GoNoGo != "NoGo - Stop"),
             aes(fill = GoNoGo,
                 fill = after_scale(darken(fill, .1, space = "HLS"))),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2,
             position = position_dodge(width = 1),
             show_guide = FALSE
  ) + 
  geom_point(data = GoNoGo_individuals %>%
               filter(GoNoGo != "NoGo - Stop"),
             aes(fill = GoNoGo),
             color = "transparent",
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2,
             alpha = .3,
             position = position_dodge(width = 1),
             show_guide = FALSE
  ) + 
  stat_summary(
    geom = "text",
    fun = "median",
    aes(label = round(..y..),
        color = GoNoGo,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 9,
    vjust = -0.8,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    limits = c(200, 2000),
    breaks = seq(200, 2000, by = 400),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, labels = c("Go", "Go (NoGo)" )) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = 6,
                                                 linetype = 0))
  ) +
  labs(
    x = NULL,
    y = "RT in milliseconds"#,
    #title = "RTs - Flanker Task"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono", face = "bold", size = 20),
    axis.text.y = element_text(
      color = rev(darken(pal1, .1, space = "HLS")), 
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
  )

## Now Next Let us plot the Error Rates

SimpleError_individual <- GoNoGo %>%
  group_by(Stim_verb, GoNoGo, Part_nr) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  mutate(id = as_factor(case_when(
    GoNoGo == "Go" ~ as.numeric(Part_nr)+21,
    GoNoGo == "NoGo - Stop" ~ as.numeric(Part_nr),
    GoNoGo == "NoGo - Go" ~ as.numeric(Part_nr) +41
  ))) 

# width of dodge
wd <- 0.8


Error_plt <- GoNoGo %>%
  group_by(Stim_verb, GoNoGo) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error, fill = GoNoGo),
         #fill = Stim_verb,
         #fill = after_scale(desaturate(lighten(fill, .8), .4))),
         #color = Stim_verb,
         #color = after_scale(darken(color, .1, space = "HLS")),
         width = .25) + 
  geom_bar(aes(fill = GoNoGo, 
               fill = after_scale(lighten(fill, .5))),
           stat = "identity", 
           position = position_dodge(width = wd),
           width = 0.5
  ) +
  geom_point(data = SimpleError_individual,
             aes(fill = GoNoGo,
                 fill = after_scale(darken(fill, .1, space = "HLS"))),
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2,
             position = position_dodge(width = wd),
             show_guide = FALSE
  ) + 
  geom_point(data = SimpleError_individual,
             aes(fill = GoNoGo),
             color = "transparent",
             stat = "identity",
             shape = 21,
             stroke = .4,
             size = 2,
             alpha = .3,
             position = position_dodge(width = wd),
             show_guide = FALSE
  ) + 
   stat_summary(
    geom = "text",
    fun = "mean",
    aes(label = str_c(format(round(..y.., 3)*100 , nsmall =1)),
        y = stage(prop_error, after_stat = 0.3),
        color = GoNoGo,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 9,
    vjust = -0.3,
    #hjust = -0.7,
    position = position_dodge(width = wd),
    show_guide = FALSE
  )  +
  coord_flip(xlim = c(1.2, NA), clip = "on") +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.7),
    breaks = seq(0.1, 1.0, by = .2),
    expand = c(.001, .001)
  )  +
  scale_color_manual(values = pal2, guide = "none") +
  scale_fill_manual(values = pal2, 
                    labels = c("Go", "Go (NoGo)", "NoGo (NoGo)"), 
                    name = "Trial Type") +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = 6,
                                                 linetype = 0
                                                 ))
  ) +
  labs(
    x = NULL,
    y = "Error in Percent"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = 15) +
  theme(
    legend.title = element_text(size=24),
    legend.text = element_text(size=21),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono", face = "bold", size = 20),
    axis.text.y = element_text(
      color = darken(pal1, .1, space = "HLS"), 
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
                alpha = 0.5, position = position_dodge(width = wd),
                width = 0.5)
# combine the plots

SRT_comb <- ggarrange(RT_plt + rremove("legend"), Error_plt + rremove("y.text"),
                      labels = c("RT", "Error"),
                      ncol = 2, nrow = 1, align = "h", widths = c(1,1),
                      font.label=list(color="black",size=28))

annotate_figure(SRT_comb,
                top =text_grob( "Go-NoGo Task", face = "bold", size = 32))



save_n <- "GoNoGo.png"
save_path <- "C:/Users/doex9445/Dateien/Julius/20Hz/Figures/Poster"

ggsave(path = save_path, filename = save_n,  dpi=300,  units = "mm", height =  160, width = 380)

