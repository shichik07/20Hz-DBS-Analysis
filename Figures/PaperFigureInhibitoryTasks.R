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

#vars
dot_size <- 0.5
text_size <- 8
guide_size <- 2
move_sum_by <- -5.4


# load data
StoppSignal<- read_csv(file = "StoppSignal.csv") %>%
  mutate(Error = 1 - Correct_Response) %>%
  mutate(StimCon = as_factor(case_when(
    Stim_verb == "130Hz" ~ "S130Hz",
    Stim_verb == "20Hz" ~ "S20Hz",
    Stim_verb == "OFF" ~ "SOFF",
    
  )))

StoppSignal_cleanRT <- StoppSignal %>% filter(Correct_Response == 1,
                                              RT < 2,
                                              RT > 0.2)%>%
  mutate(RT = RT*1000)

# Integration Medthod to get Stopsignal Reaction times

StoppSignal_individuals <- StoppSignal_cleanRT %>%
  filter(StopTrl != "Stop") %>%
  group_by(Stim_verb, StopTrl, Part_nr) %>%
  summarise(RT = median(RT))

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

SSRT <- StoppSignal %>% 
  group_by(Part_nr, Stim_verb) %>%
  group_modify(~ tibble::enframe(fn_scoretidy(.), name = NULL))


pal <- c("#FF8C00", "#A034F0", "#159090")
pal <- c("#26185F", "#0095AF", "#9ADCBB")
pal1 <- c("#FF6F00FF", "#008EA0FF", "#8A4198FF")


SSRT_txt <- SSRT %>%
  group_by(Stim_verb) %>%
  summarise(medSSRT = median(value))


Stop_RT_plt <- SSRT %>%
  group_by(Stim_verb) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = value)) + 
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
    title = "Stop-Change Task",
    subtitle = "SCRT"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(panel.grid.minor = element_blank(),
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
        plot.margin = margin(15, 15, 10, 15),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
  )

## Now Next Let us plot the Error Rates

SimpleError_individual <- StoppSignal %>%
  filter(StopTrl != "Stop") %>%
  group_by(Stim_verb, Part_nr) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1))


Stop_Error_plt <- StoppSignal %>%
  filter(StopTrl != "Stop") %>%
  group_by(Stim_verb) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error 
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
    subtitle = "Error",
    fill = "Stimulation"
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
  )

### Now the GoNoGo data
# set directory

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
                                    RT < 3,
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


GNG_RT_plt <-  GoNoGo_cleanRT %>% 
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
  geom_point(data = GoNoGo_individuals,
             aes(fill = GoNoGo,
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
        color = GoNoGo,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by-3,
    position = position_dodge(width = 1),
    show_guide = FALSE
  ) +
  # Add significance annotation for 20Hz certain Go vs uncertain Go trials (asterisk)
  geom_path(data = data.frame(x = c(2.2, 2.2, 2.8, 2.8), 
                             y = c(1850, 1900, 1900, 1850), 
                             GoNoGo = c("Go", "Go", "uncertain Go", "uncertain Go")),
           aes(x = x, y = y), 
           inherit.aes = FALSE,
           color = "gray30",
           linewidth = 0.8) +
  # Add asterisk
  annotate("text", x = 2.5, y = 1950, label = "*", 
           size = 2.5, fontface = "bold") +
  scale_y_continuous(
    limits = c(200, 2000),
    breaks = seq(200, 2000, by = 400),
    expand = c(.001, .001)
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, labels = c("certain Go", "uncertain Go" )) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0, 
                                                 label.position = "bottom"))
  ) +
  labs(
    x = NULL,
    y = "RT in ms",
    title = "Go-NoGo Task",
    subtitle = "RT"
  ) +
  theme_minimal(base_family = "Zilla Slab", base_size = text_size) +
  theme(
    panel.grid.minor = element_blank(),
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
    plot.margin = margin(15, 15, 10, 15),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
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

GNG_Error_plt <- GoNoGo %>%
  group_by(Stim_verb, GoNoGo) %>%
  summarise(prop_error = abs(round(mean(Correct_Response-1), 4)), se = std.error(Correct_Response-1)) %>%
  ggplot(aes(x = fct_rev(Stim_verb), y = prop_error, fill = GoNoGo),
         width = .25) + 
  geom_bar(aes(fill = GoNoGo, 
               fill = after_scale(lighten(fill, .5))),
           stat = "identity", 
           position = position_dodge(width = 0.8),
           width = 0.5
  ) +
  geom_point(data = SimpleError_individual,
             aes(fill = GoNoGo,
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
        color = GoNoGo,
        color = after_scale(darken(color, .3, space = "HLS"))),
    family = "Roboto Mono",
    fontface = "bold",
    size = 2.5,
    vjust = move_sum_by+2,
    position = position_dodge(width = 1),
    show_guide = FALSE
  )  +
  # Add significance annotation for 20Hz uncertain Go trials vs 130Hz uncertain Go trials
  geom_path(data = data.frame(x = c(2-0.3, 2-0.3, 3-0.3, 3-0.3), 
                             y = c(0.25, 0.35, 0.35, 0.25), 
                             GoNoGo = rep("uncertain Go", 4)),
           aes(x = x, y = y), 
           inherit.aes = FALSE,
           color = "gray30",
           linewidth = 0.8) +
  # Add asterisk
  annotate("text", x = 2.5-0.3, y = 0.37, label = "*", 
           size = 2.5, fontface = "bold") +
           
  # Add significance annotation for 20Hz uncertain Go trials vs OFF uncertain Go trials
  geom_path(data = data.frame(x = c(1-0.3, 1-0.3, 2-0.3, 2-0.3), 
                             y = c(0.25-0.05, 0.42-0.05, 0.42-0.05, 0.25-0.05), 
                             GoNoGo = rep("uncertain Go", 4)),
           aes(x = x, y = y), 
           inherit.aes = FALSE,
           color = "gray30",
           linewidth = 0.8) +
  # Add asterisk
  annotate("text", x = 1.5-0.3, y = 0.44-0.05, label = "*", 
           size = 2.5, fontface = "bold") +
           
  # Add significance annotation for 20Hz NoGo trials vs 130Hz NoGo trials
  geom_path(data = data.frame(x = c(2+0.3, 2+0.3, 3+0.3, 3+0.3), 
                             y = c(0.55, 0.57, 0.57, 0.55), 
                             GoNoGo = rep("uncertain Go", 4)),
           aes(x = x, y = y), 
           inherit.aes = FALSE,
           color = "gray30",
           linewidth = 0.8) +
  # Add asterisk
  annotate("text", x = 2.5+0.3, y = 0.6, label = "*", 
           size = 2.5, fontface = "bold") +
           
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 0.62),
    breaks = seq(0.1, 1.0, by = .2),
    expand = c(.001, .001)
  )  +
  scale_color_manual(values = pal2, guide = "none") +
  scale_fill_manual(values = pal2, labels = c("certain Go", "uncertain Go", "NoGo")) +
  guides(fill = guide_legend(override.aes = list(shape = NA, 
                                                 size = guide_size,
                                                 linetype = 0))
  ) +
  labs(
    x = NULL,
    y = "Error in %",
    subtitle = "Error"
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
  )

# combine the plots

# Make sure all plots have proper error bars
Stop_Error_plt <- Stop_Error_plt + 
  geom_errorbar(aes(ymin = prop_error - se, ymax = prop_error + se), 
                alpha = 0.5, position = position_dodge(width = 0.8),
                width = 0.5)

GNG_Error_plt <- GNG_Error_plt + 
  geom_errorbar(aes(ymin = prop_error - se, ymax = prop_error + se), 
                alpha = 0.5, position = position_dodge(width = 0.8),
                width = 0.5)

SRT_comb <- ggarrange(Stop_RT_plt + rremove("legend"), Stop_Error_plt, GNG_RT_plt + rremove("legend"),
                      GNG_Error_plt,
                      labels = c("A", "B", "C", "D"),
                      ncol = 2, nrow = 2, align = "h", widths = c(3,4),
                      font.label=list(color="black",size=10))

# Save as TIFF
save_n <- "Inhibition_fin.tiff"
save_path <- "D:/Data/Dropbox/PhD_Thesis/UniOL/Julius/20Hz-DBS-Analysis/Figures/Paper"
ggsave(path = save_path, filename = save_n, dpi=600, units = "mm", height = 110, width = 160)

# Save as PNG
save_n_png <- "Inhibition_fin.png"
ggsave(path = save_path, filename = save_n_png, dpi=600, units = "mm", height = 110, width = 160)
