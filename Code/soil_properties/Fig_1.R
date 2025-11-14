
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
data.SOC <- read_excel(file.path(path_properties,
                                 "soil.properties.xlsx"), sheet = "Data") %>%
  select(Sample.ID, Pairs, Species, Status, Topography, SOC) %>%
  arrange(Topography, Species, Pairs, Status) %>%
  mutate(Pairs = factor(Pairs)) %>%
  data.frame()




# 2. SOC----

## 2.1) Model construction----

model.SOC <- lmer(SOC ~ Status * Topography + (1 | Pairs),
                  data = data.SOC)


# Normality test
check_normality(model.SOC)

# Heteroscedasticity test
check_heteroscedasticity(model.SOC)


# Variance analysis
emmeans(model.SOC, ~ Status | Topography) %>%
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

p1 <- data.SOC %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Ridge") %>%
  ggplot(aes(x = Status, y = SOC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 3, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10))


Ridge.line.data <- ggplot_build(p1)$data[[2]] %>%
  left_join(data.SOC[, c("SOC", "Pairs")], by = c("y" = "SOC")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })



### Fig_1a----
p.Ridge <- p1 + geom_bezier0(aes(x, y, group = Pairs),
                             data = Ridge.line.data, colour = "grey70") +
  labs(x = NULL, y = "SOC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 49,
               label = "",
               color = "black",
               size = 0.75,
               tip.length = 0.05) +
  annotate("text",
           x = 1.5, y = 51,
           label = expression(italic(p) == 0.0137),
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

p2 <- data.SOC %>%
  mutate(Group = paste(Topography, Status, sep = "_")) %>%
  subset(Topography == "Valley") %>%
  ggplot(aes(x = Status, y = SOC)) +
  geom_violin(aes(fill = Group, colour = Group),
              adjust = 1.5, bw = 3, alpha = 1,
              linewidth = 1, show.legend = F,
              trim = F, width = 0.6, scale = "width") +
  geom_point(aes(group = Status, fill = Group, colour = Group),
             shape = 21, size = 6, color = "black",
             position = position_jitterdodge(jitter.width = 0.4,
                                             seed = 1234),
             show.legend = F) +
  scale_color_manual(values = mapping.color) +
  scale_fill_manual(values = mapping.fill) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, 10))


Valley.line.data <- ggplot_build(p2)$data[[2]] %>%
  left_join(data.SOC[, c("SOC", "Pairs")], by = c("y" = "SOC")) %>%
  select(x, y, group, Pairs) %>%
  group_by(Pairs) %>%
  group_modify(function (data, group) {
    data.frame(
      x = c(data$x[1], data$x[1] + 0.3, data$x[2] - 0.3, data$x[2]),
      y = c(data$y[1], data$y[1], data$y[2], data$y[2])
    ) 
  })



### Fig_1c----
p.Valley <- p2 + geom_bezier0(aes(x, y, group = Pairs),
                              data = Valley.line.data, colour = "grey70") +
  labs(x = NULL, y = "SOC"~(g~kg^-1~soil)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 60,
               label = "",
               color = "black",
               size = 0.7,
               tip.length = 0.03) +
  annotate("text",
           x = 1.5, y = 62.5,
           label = expression(italic(p) == 0.0128),
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




# 3. ΔSOC vs Δaggregate----


### Fig_1b----
## 4.1) Ridge----
delta.aggregate.ridge <- read_excel(file.path(path_properties,
                                              "soil.properties.xlsx"),
                                    sheet = "Data") %>%
  subset(Topography == "Ridge") %>%
  select(Pairs, Species, Status, Topography, SOC, POC, MAOC) %>%
  arrange(Pairs, Status) %>% 
  data.frame() %>%
  group_by(Species, Topography, Pairs) %>%
  group_modify(function(data, Status){
    data.frame(SOC = data$SOC[2] - data$SOC[1],
               POC = data$POC[2] - data$POC[1],
               MAOC = data$MAOC[2] - data$MAOC[1])
  }) %>% data.frame()



### 4.1.1 Model construction----
model.delta.aggregate.ridge <- delta.aggregate.ridge %>%
  lm(SOC ~ POC * MAOC, data = .)

check_normality(model.delta.aggregate.ridge)
check_heteroscedasticity(model.delta.aggregate.ridge)




### 4.1.2 POC marginal dependence----

newdata.delta.POC.ridge <-
  expand.grid(POC = seq(min(delta.aggregate.ridge$POC),
                        max(delta.aggregate.ridge$POC),
                        length.out = 100),
              MAOC = mean(delta.aggregate.ridge$MAOC)) %>%
  data.frame()


predict.delta.POC.ridge <- predict.lm(model.delta.aggregate.ridge,
                                      newdata = newdata.delta.POC.ridge,
                                      se.fit = T)


plot.data.delta.POC.ridge <- newdata.delta.POC.ridge %>%
  mutate(SOC.fit = predict.delta.POC.ridge$fit,
         SE = predict.delta.POC.ridge$se.fit,
         lower.fit = SOC.fit - 1.96 * SE,
         upper.fit = SOC.fit + 1.96 * SE)




# Marginal slope
emtrends(model.delta.aggregate.ridge, ~ 1, var = "POC",
         infer = c(TRUE, TRUE)) %>%
  summary()


# Visualization
delta.aggregate.ridge %>%
  ggplot(aes(x = POC, y = SOC)) +
  geom_point(size = 6, fill = "#FCD29E", shape = 21) + 
  geom_ribbon(data = plot.data.delta.POC.ridge,
              aes(x = POC, y = SOC.fit,
                  ymin = lower.fit, ymax = upper.fit),
              fill = "#BFA27D", inherit.aes = F,
              alpha = 0.5) +
  geom_line(data = plot.data.delta.POC.ridge,
            aes(x = POC, y = SOC.fit), colour = "#E60012",
            inherit.aes = F, linewidth = 1, linetype = 2) +
  geom_text(x = -0.4, y = 7,
            size = 7,
            label = "Slope = 0.541; p = 0.4228") +
  labs(x = "Changes in POC on ridges",
       y = "Changes in SOC on ridges") +
  scale_x_continuous(limits = c(-2, 2.1), breaks = seq(-2, 2, 1)) +
  scale_y_continuous(limits = c(-12, 8), breaks = seq(-12, 8, 4)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),  # x轴刻度线朝里（负值朝内）
        axis.ticks.length.y = unit(-0.25, "cm"))



### 4.1.3 MAOC marginal dependence----

newdata.delta.MAOC.ridge <-
  expand.grid(MAOC = seq(min(delta.aggregate.ridge$MAOC),
                         max(delta.aggregate.ridge$MAOC),
                         length.out = 100),
              POC = mean(delta.aggregate.ridge$POC)) %>%
  data.frame()


predict.delta.MAOC.ridge <- predict.lm(model.delta.aggregate.ridge,
                                       newdata = newdata.delta.MAOC.ridge,
                                       se.fit = T)


plot.data.delta.MAOC.ridge <- newdata.delta.MAOC.ridge %>%
  mutate(SOC.fit = predict.delta.MAOC.ridge$fit,
         SE = predict.delta.MAOC.ridge$se.fit,
         lower.fit = SOC.fit - 1.96 * SE,
         upper.fit = SOC.fit + 1.96 * SE)



# Marginal slope
emtrends(model.delta.aggregate.ridge, ~ 1, var = "MAOC",
         infer = c(TRUE, TRUE)) %>%
  summary()



# Visualization
delta.aggregate.ridge %>%
  ggplot(aes(x = MAOC, y = SOC)) +
  geom_point(size = 6, fill = "#FCD29E", shape = 21) + 
  geom_ribbon(data = plot.data.delta.MAOC.ridge,
              aes(x = MAOC, y = SOC.fit,
                  ymin = lower.fit, ymax = upper.fit),
              fill = "#BFA27D", inherit.aes = F, alpha = 0.5) +
  geom_line(data = plot.data.delta.MAOC.ridge,
            aes(x = MAOC, y = SOC.fit), colour = "#E60012",
            inherit.aes = F, linewidth = 1, linetype = 1) +
  geom_text(x = -3, y = 7,
            size = 7,
            label = "Slope = 0.621; p = 0.0080") +
  labs(x = "Changes in MAOC on ridges",
       y = "Changes in SOC on ridges") +
  scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, 4)) +
  scale_y_continuous(limits = c(-12, 8), breaks = seq(-12, 8, 4)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),  # x轴刻度线朝里（负值朝内）
        axis.ticks.length.y = unit(-0.25, "cm"))


### Fig_1d----

## 4.2) Valley----
delta.aggregate.valley <- read_excel(file.path(path_properties,
                                               "soil.properties.xlsx"),
                                     sheet = "Data") %>%
  subset(Topography == "Valley") %>%
  select(Pairs, Species, Status, Topography, SOC, POC, MAOC) %>%
  arrange(Pairs, Status) %>% 
  data.frame() %>%
  group_by(Species, Topography, Pairs) %>%
  group_modify(function(data, Status){
    data.frame(SOC = data$SOC[2] - data$SOC[1],
               POC = data$POC[2] - data$POC[1],
               MAOC = data$MAOC[2] - data$MAOC[1])
  }) %>% data.frame()



### 4.2.1 Model construction----
model.delta.aggregate.valley <- delta.aggregate.valley %>%
  lm(SOC ~ POC * MAOC, data = .)

check_normality(model.delta.aggregate.valley)
check_heteroscedasticity(model.delta.aggregate.valley)


### 4.2.2 POC marginal dependence----

newdata.delta.POC.valley <-
  expand.grid(POC = seq(min(delta.aggregate.valley$POC),
                        max(delta.aggregate.valley$POC),
                        length.out = 100),
              MAOC = mean(delta.aggregate.valley$MAOC)) %>%
  data.frame()


predict.delta.POC.valley <- predict.lm(model.delta.aggregate.valley,
                                       newdata = newdata.delta.POC.valley,
                                       se.fit = T)


plot.data.delta.POC.valley <- newdata.delta.POC.valley %>%
  mutate(SOC.fit = predict.delta.POC.valley$fit,
         SE = predict.delta.POC.valley$se.fit,
         lower.fit = SOC.fit - 1.96 * SE,
         upper.fit = SOC.fit + 1.96 * SE)



# Marginal slope
emtrends(model.delta.aggregate.valley, ~ 1, var = "POC",
         infer = c(TRUE, TRUE)) %>%
  summary()



# Visualization
delta.aggregate.valley %>%
  ggplot(aes(x = POC, y = SOC)) +
  geom_point(size = 6, fill = "#C6DAED", shape = 21) + 
  geom_ribbon(data = plot.data.delta.POC.valley,
              aes(x = POC, y = SOC.fit,
                  ymin = lower.fit, ymax = upper.fit),
              fill = "#768CA0", inherit.aes = F, alpha = 0.5) +
  geom_line(data = plot.data.delta.POC.valley,
            aes(x = POC, y = SOC.fit), colour = "#3666FD",
            inherit.aes = F, linewidth = 1, linetype = 1) +
  geom_text(x = 0, y = 11,
            size = 7,
            label = "Slope = 0.529; p = 0.0166") +
  labs(x = "Changes in POC on valleys",
       y = "Changes in SOC on valleys") +
  scale_x_continuous(limits = c(-6, 9), breaks = seq(-6, 9, 3)) +
  scale_y_continuous(limits = c(-6, 12), breaks = seq(-6, 12, 3)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"), 
        axis.ticks.length.y = unit(-0.25, "cm"))



### 4.2.3 MAOC marginal dependence----

newdata.delta.MAOC.valley <-
  expand.grid(MAOC = seq(min(delta.aggregate.valley$MAOC),
                         max(delta.aggregate.valley$MAOC),
                         length.out = 100),
              POC = mean(delta.aggregate.valley$POC)) %>%
  data.frame()


predict.delta.MAOC.valley <- predict.lm(model.delta.aggregate.valley,
                                        newdata = newdata.delta.MAOC.valley,
                                        se.fit = T)


plot.data.delta.MAOC.valley <- newdata.delta.MAOC.valley %>%
  mutate(SOC.fit = predict.delta.MAOC.valley$fit,
         SE = predict.delta.MAOC.valley$se.fit,
         lower.fit = SOC.fit - 1.96 * SE,
         upper.fit = SOC.fit + 1.96 * SE)



# Marginal slope
emtrends(model.delta.aggregate.valley, ~ 0, var = "MAOC",
         infer = c(TRUE, TRUE)) %>%
  summary()


# Visualization
delta.aggregate.valley %>%
  ggplot(aes(x = MAOC, y = SOC)) +
  geom_point(size = 6, fill = "#C6DAED", shape = 21) + 
  geom_ribbon(data = plot.data.delta.MAOC.valley,
              aes(x = MAOC, y = SOC.fit,
                  ymin = lower.fit, ymax = upper.fit),
              fill = "#768CA0", inherit.aes = F, alpha = 0.5) +
  geom_line(data = plot.data.delta.MAOC.valley,
            aes(x = MAOC, y = SOC.fit), colour = "#3666FD",
            inherit.aes = F, linewidth = 1, linetype = 1) +
  geom_text(x = -3, y = 11,
            size = 7,
            label = "Slope = 0.200; p = 0.3508") +
  labs(x = "Changes in MAOC on valleys",
       y = "Changes in SOC on valleys") +
  scale_x_continuous(limits = c(-8, 6), breaks = seq(-8, 6, 2)) +
  scale_y_continuous(limits = c(-6, 12), breaks = seq(-6, 12, 3)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"), 
        axis.ticks.length.y = unit(-0.25, "cm"))

