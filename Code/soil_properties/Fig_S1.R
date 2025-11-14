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
library(clubSandwich)


# Data
data.SOC1 <- read_excel(file.path(path_properties,
                                  "soil.properties.xlsx"), sheet = "Data") %>%
  arrange(Topography, Pairs, Status) %>%
  mutate(Pairs = factor(Pairs)) %>%
  select(Sample.ID, Pairs, Status, Topography,
         POC, MAOC, Moisture) %>%
  data.frame()


# Fig_S1a----

# 2. POC----

## 2.1) Model construction----
model.POC <- lmer(POC ~ Status * Topography + (1 | Pairs),
                  data = data.SOC1)


# Normality test
check_normality(model.POC)

# Heteroscedasticity test
check_heteroscedasticity(model.POC)



# Variance analysis

vcovCR.POC <- vcovCR(model.POC, type = "CR2", cluster = data.SOC1$Pairs)

emmeans(model.POC, ~ Status | Topography,
        vcov. = vcovCR.POC) %>%
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


### 2.2.1 Ridge----

p1 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = POC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 1, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 18), breaks = seq(0, 18, 3))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.SOC1[, c("POC", "Pairs")], by = c("y" = "POC")) %>%
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
  labs(x = NULL, y = "POC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 15.5,
               label = "",
               color = "black",
               size = 0.7,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 16.2,
           label = "ns",
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

p2 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = POC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 1.5, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 6))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.SOC1[, c("POC", "Pairs")], by = c("y" = "POC")) %>%
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
  labs(x = NULL, y = "POC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 27,
               label = "",
               color = "black",
               size = 0.7,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 28,
           label = expression(italic(p) == 0.0421),
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




# Fig_S1b----
# 3. MAOC----


## 3.1) Model construction----
model.MAOC <- lmer(MAOC ~ Status * Topography + (1 | Pairs),
                   data = data.SOC1)


# Normality test
check_normality(model.MAOC)

# Heteroscedasticity test
check_heteroscedasticity(model.MAOC)


# Variance analysis
emmeans(model.MAOC, ~ Status | Topography) %>%
  pairs(adjust = "tukey")



## 3.2) Visualization----

### 3.2.1 Ridge----

p1 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = MAOC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 2, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.SOC1[, c("MAOC", "Pairs")], by = c("y" = "MAOC")) %>%
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
  labs(x = NULL, y = "MAOC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 42,
               label = "",
               color = "black",
               size = 0.7,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 44,
           label = "ns",
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

p2 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = MAOC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 2, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.SOC1[, c("MAOC", "Pairs")], by = c("y" = "MAOC")) %>%
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
  labs(x = NULL, y = "MAOC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 44,
               label = "",
               color = "black",
               size = 0.7,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 46,
           label = "ns",
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


# Fig_S1c----
# 4. Moisture----


## 4.1) Model construction----
model.Moisture <- lmer(Moisture ~ Status * Topography + (1 | Pairs),
                       data = data.SOC1)


# Normality test
check_normality(model.Moisture)

# Heteroscedasticity test
check_heteroscedasticity(model.Moisture)


# Variance analysis
emmeans(model.Moisture, ~ Status | Topography) %>%
  pairs(adjust = "tukey")



## 4.2) Visualization----


### 4.2.1 Ridge----

p1 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = Moisture)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 1.2, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(5, 30), breaks = seq(5, 30, 5))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.SOC1[, c("Moisture", "Pairs")], by = c("y" = "Moisture")) %>%
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
  labs(x = NULL, y = "Moisture"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 26,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.08) +
  annotate("text",
           x = 1.5, y = 27,
           label = expression(italic(p) < 0.001),
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


### 4.2.2 Valley----

p2 <- data.SOC1 %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = Moisture)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 2, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(5, 37), breaks = seq(5, 37, 5))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.SOC1[, c("Moisture", "Pairs")], by = c("y" = "Moisture")) %>%
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
  labs(x = NULL, y = "Moisture"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 35,
               label = "",
               color = "black",
               size = 0.5,
               tip.length = 0.03) +
  annotate("text",
           x = 1.5, y = 36,
           label = "ns",
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


