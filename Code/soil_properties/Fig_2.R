# Set the working directory
setwd("G:/R_script")

path_properties <- file.path(getwd(), "Data/Soil_properties")


# 1. Loading package and data----


# Package
library(readxl)
library(tidyverse)
library(lme4)
library(performance)
library(emmeans)
library(ggforce)
library(ggpubr)


# Data
data.enzymes <- read_excel(file.path(path_properties,
                     "soil.properties.xlsx"), sheet = "Data") %>%
  arrange(Topography, Pairs, Status) %>%
  mutate(Pairs = factor(Pairs)) %>%
  mutate(
    mass_normalized_activity = rowMeans(
      as.matrix(across(c(AG, BG, CBH, BX, NAG, PER, PHO), ~ as.numeric(scale(.x)))),
      na.rm = TRUE),
    SOC_normalized_activity = rowMeans(
      as.matrix(across(c(AG, BG, CBH, BX, NAG, PER, PHO), ~ as.numeric(scale(.x / SOC)))),
      na.rm = TRUE)
  ) %>%
  select(Sample.ID, Pairs, Status, Topography,
         mass_normalized_activity, SOC_normalized_activity) %>%
  data.frame()





# 2. Mass_normalized_activity----


## 2.1) Model construction----
model.mass.normalized <- lmer(mass_normalized_activity ~ Status * Topography + (1 | Pairs),
                              data = data.enzymes)


# Normality test
check_normality(model.mass.normalized)

# Heteroscedasticity test
check_heteroscedasticity(model.mass.normalized)


# Variance analysis
emmeans(model.mass.normalized, ~ Status | Topography) %>%
  pairs(adjust = "tukey")



## 2.2) Visualization----

# Color
mapping.fill <- c("Ridge_Alive" = "#FDD49E", "Ridge_Dead" = "#FC8D59",
                  "Valley_Alive" = "#C6DBEF", "Valley_Dead" = "#168ED4")

mapping.color <- c("Ridge_Alive"  = "#FCAE91",
                   "Ridge_Dead"   = "#A63603",
                   "Valley_Alive" = "#9ECAE1",
                   "Valley_Dead"  = "#08519C")

errorbar.color <- c("Ridge_Alive"  = "#F16913",
                    "Ridge_Dead"   = "#D94801",
                    "Valley_Alive" = "#6BAED6",
                    "Valley_Dead"  = "#3182BD")   

# Fig_2a----
### 2.2.1 Ridge----

p1 <- data.enzymes %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = mass_normalized_activity)) +
  geom_violin(aes(fill = Group, colour = Group),
              linewidth = 1, show.legend = F, bw = 0.5,
              trim = F, width = 0.5, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.enzymes[, c("mass_normalized_activity", "Pairs")], by = c("y" = "mass_normalized_activity")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })




p.Ridge <- p1 + geom_bezier0(aes(x, y, group = Pairs),
                             data = Ridge.line.data, colour = "grey70") +
  labs(x = NULL, y = "mass_normalized_activity"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 3.4,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 3.6,
           label = expression(italic(p) == 0.0015),
           size = 6) +
  stat_summary(geom = "errorbar", fun.data = mean_se,
               color = "#D94801", linewidth = 1,
               width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean,
               size = 3, color = "#D94801",
               show.legend = F) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


p.Ridge


### 2.2.2 Valley----

p2 <- data.enzymes %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = mass_normalized_activity)) +
  geom_violin(aes(fill = Group, colour = Group),
              linewidth = 1, show.legend = F, bw = 0.3,
              trim = F, width = 0.5, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.enzymes[, c("mass_normalized_activity", "Pairs")], by = c("y" = "mass_normalized_activity")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })




p.Valley <- p2 + geom_bezier0(aes(x, y, group = Pairs),
                              data = Valley.line.data, colour = "grey70") +
  labs(x = NULL, y = "mass_normalized_activity"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 2.2,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 2.4,
           label = expression(italic(p) == 0.6112),
           size = 6) +
  stat_summary(geom = "errorbar", fun.data = mean_se,
               color = "blue", linewidth = 1,
               width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean,
               size = 3, color = "blue",
               show.legend = F) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


p.Valley





# 3. SOC_normalized_activity----


## 3.1) Model construction----
model.SOC.normalized <- lmer(SOC_normalized_activity ~ Status * Topography + (1 | Pairs),
                             data = data.enzymes)


# Normality test
check_normality(model.SOC.normalized)

# Heteroscedasticity test
check_heteroscedasticity(model.SOC.normalized)


# variance analysis
emmeans(model.SOC.normalized, ~ Status | Topography) %>%
  pairs(adjust = "tukey")



## 3.2) Visualization----

# Fig_2b----
### 3.2.1 Ridge----

p1 <- data.enzymes %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = SOC_normalized_activity)) +
  geom_violin(aes(fill = Group, colour = Group),
              linewidth = 1, show.legend = F, bw = 0.3,
              trim = F, width = 0.5, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.enzymes[, c("SOC_normalized_activity", "Pairs")], by = c("y" = "SOC_normalized_activity")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })




p.Ridge <- p1 + geom_bezier0(aes(x, y, group = Pairs),
                             data = Ridge.line.data, colour = "grey70") +
  labs(x = NULL, y = "SOC_normalized_activity"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 2.7,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 2.9,
           label = expression(italic(p) == 0.4403),
           size = 6) +
  stat_summary(geom = "errorbar", fun.data = mean_se,
               color = "#D94801", linewidth = 1,
               width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean,
               size = 3, color = "#D94801",
               show.legend = F) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


p.Ridge


### 3.2.2 Valley----

p2 <- data.enzymes %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = SOC_normalized_activity)) +
  geom_violin(aes(fill = Group, colour = Group),
              linewidth = 1, show.legend = F, bw = 0.5,
              trim = F, width = 0.5, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 1))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.enzymes[, c("SOC_normalized_activity", "Pairs")], by = c("y" = "SOC_normalized_activity")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })




p.Valley <- p2 + geom_bezier0(aes(x, y, group = Pairs),
                              data = Valley.line.data, colour = "grey70") +
  labs(x = NULL, y = "SOC_normalized_activity"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 2.9,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 3.1,
           label = expression(italic(p) == 0.0315),
           size = 6) +
  stat_summary(geom = "errorbar", fun.data = mean_se,
               color = "blue", linewidth = 1,
               width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean,
               size = 3, color = "blue",
               show.legend = F) +
  theme_bw(base_size = 18) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


p.Valley



