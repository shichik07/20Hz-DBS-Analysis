library(haven) # import SPSS files
library(dplyr)
library(tidyr)
library(readr)
library(ggtext)
library(colorspace)
library(ggplot2)
library(forcats)
library(ggsci)
library(plotrix) # to calculate standard error
library(ggpubr)
library(stringr)
library(ggdist)

# set directory
wd = "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Data/Extracted"
setwd(wd)

# vars
dot_size <- 0.5
text_size <- 8
guide_size <- 2
move_sum_by <- -5.4

# load data
Flanker <- read_csv(file = "Flanker.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
  )))

# Clean RT data
Flanker_cleanRT <- Flanker %>% 
  filter(Correct_Response == 1,
         RT < 2,
         RT > 0.2) %>%
  mutate(RT = RT*1000)

# Calculate congruency effect (incongruent - congruent RT) for each participant and stimulation condition
Flanker_congruency_effect <- Flanker_cleanRT %>%
  group_by(Stim_verb, Congruency, Part_nr) %>%
  summarise(mean_RT = mean(RT), .groups = "drop") %>%
  pivot_wider(names_from = Congruency, values_from = mean_RT) %>%
  mutate(congruency_effect = incongruent - congruent) %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(congruency_effect = mean(congruency_effect), .groups = "drop")

# Calculate summary statistics for plotting RT congruency effect
Flanker_congruency_summary <- Flanker_congruency_effect %>%
  group_by(Stim_verb) %>%
  summarise(
    mean_effect = mean(congruency_effect),
    se = std.error(congruency_effect),
    .groups = "drop"
  )

# Calculate congruency effect for error rates (incongruent - congruent error rates)
Flanker_error_congruency_effect <- Flanker %>%
  group_by(Stim_verb, Congruency, Part_nr) %>%
  summarise(error_rate = mean(Error), .groups = "drop") %>%
  pivot_wider(names_from = Congruency, values_from = error_rate) %>%
  mutate(error_congruency_effect = incongruent - congruent) %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(error_congruency_effect = mean(error_congruency_effect), .groups = "drop")

# Calculate summary statistics for plotting error congruency effect
Flanker_error_congruency_summary <- Flanker_error_congruency_effect %>%
  group_by(Stim_verb) %>%
  summarise(
    mean_effect = mean(error_congruency_effect),
    se = std.error(error_congruency_effect),
    .groups = "drop"
  )

# Color palette
pal1 <- c("#FF6F00FF", "#008EA0FF", "#8A4198FF")

# Create the RT congruency effect plot
FLT_Congruency_RT_plt <- ggplot(Flanker_congruency_effect, aes(x = fct_rev(Stim_verb), y = congruency_effect)) + 
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
  geom_boxplot(
    width = .35, 
    outlier.shape = NA,
    show_guide = FALSE
  ) +
  geom_point(
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
  scale_y_continuous(
    limits = c(0, 300),
    breaks = seq(0, 300, by = 50),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal1, guide = "none") +
  scale_fill_manual(values = pal1, guide = "none") +
  labs(
    x = NULL,
    y = "Congruency Effect (ms)",
    subtitle = "RT Congruency Effect"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    legend.title = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid.minor = element_blank(),
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
    plot.margin = margin(15, 15, 10, 15),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_errorbar(
    data = Flanker_congruency_summary,
    aes(y = mean_effect, ymin = mean_effect - se, ymax = mean_effect + se), 
    alpha = 0.5, width = 0.5
  )

# Create the Error congruency effect plot
FLT_Congruency_Error_plt <- ggplot(Flanker_error_congruency_effect, aes(x = fct_rev(Stim_verb), y = error_congruency_effect)) + 
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
  geom_boxplot(
    width = .35, 
    outlier.shape = NA,
    show_guide = FALSE
  ) +
  geom_point(
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
    aes(label = round(..y..*100, 3),
        color = Stim_verb,
        fill = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by-3,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal1, guide = "none") +
  scale_fill_manual(values = pal1, guide = "none") +
  labs(
    x = NULL,
    y = "Error Congruency Effect",
    subtitle = "Error Congruency Effect"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    legend.title = element_text(size=text_size),
    legend.text = element_text(size=text_size),
    legend.key.size = unit(0.3, 'cm'),
    panel.grid.minor = element_blank(),
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
    plot.margin = margin(15, 15, 10, 15),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  geom_errorbar(
    data = Flanker_error_congruency_summary,
    aes(y = mean_effect, ymin = mean_effect - se, ymax = mean_effect + se), 
    alpha = 0.5, width = 0.5
  )

# Combine the plots
FLT_Congruency_plt <- ggarrange(
  FLT_Congruency_RT_plt, FLT_Congruency_Error_plt,
  ncol = 2,
  labels = c("A", "B"),
  font.label = list(size = text_size, family = "Zilla Slab", face = "bold")
) %>%
  annotate_figure(top = text_grob("Flanker Task Congruency Effect", 
                                 face = "bold", 
                                 size = text_size, 
                                 family = "Zilla Slab"))

# Create a function to add white background to any plot
add_white_bg <- function(plot) {
  plot_with_bg <- ggplot() + 
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)) +
    annotation_custom(ggplotGrob(plot))
  return(plot_with_bg)
}

# Apply white background to the combined figure
FLT_Congruency_plt <- add_white_bg(FLT_Congruency_plt)

# Save the figure
save_path <- "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Paper"
ggsave(path = save_path, filename = "Flanker_Congruency_Effect.tiff", plot = FLT_Congruency_plt, dpi = 600, units = "mm", height = 110, width = 160)
ggsave(path = save_path, filename = "Flanker_Congruency_Effect.png", plot = FLT_Congruency_plt, dpi = 600, units = "mm", height = 110, width = 160)
