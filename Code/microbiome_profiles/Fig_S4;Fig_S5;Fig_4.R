# Set the working directory
setwd("G:/R_script")

path_microbiome <- file.path(getwd(), "Data/microbiome_profiles")
path_properties <- file.path(getwd(), "Data/Soil_properties")
path_output <- file.path(getwd(), "Output")


# load packages
library(tidyverse)
library(caret)
library(xgboost)
library(SHAPforxgboost)
library(ggpubr)
library(patchwork)
library(edgeR)
library(readxl)
library(pheatmap)
library(RColorBrewer)


# 1. Loading data----

## 1.1) Sample information----
sample.info <- read_excel(file.path(path_properties,
                                    "soil.properties.xlsx"), sheet = "Data") %>%
  arrange(Topography, Pairs, Status) %>%
  mutate(Pairs = factor(Pairs))


sample.ridge <- sample.info %>%
  subset(Topography == "Ridge") %>%
  pull(Sample.ID)


sample.valley <- sample.info %>%
  subset(Topography == "Valley") %>%
  pull(Sample.ID)



## 1.2) microbial_composition----

load(file.path(path_output,
               "microbial_composition_results.Rdata"))




## 1.3) Target gene----
bacteria.target.gene <- c("Exo_amylase", "Endo_amylase", "Exo_cellulase",
                          "Endo_cellulase", "Exo_hemicellulose",
                          "Endo_hemicellulose", "Exo_chitosanase", "Endo_chitosanase",
                          "TEA_Oxygen", "TEA_Nitrate")


fungi.target.gene <- c("Exo_amylase", "Endo_amylase", "Exo_cellulase",
                       "Endo_cellulase", "Exo_hemicellulose", "Endo_hemicellulose",
                       "Exo_chitosanase", "Endo_chitosanase", "Ligninase")


## 1.4) Gene copy number matrix----

### 1.4.1 Bacteria----
bacteria.data <-  read.csv(file.path(path_microbiome,
                                     "bacterial.taxa.genome.gene.copy.number.matrix.csv"),
                           header = TRUE) %>%
  group_by(Silva.Species) %>%
  summarise(across(all_of(bacteria.target.gene), ~mean(.x)), .groups = "drop") %>%
  left_join(bacteria.wilcox.result_df[, c("taxa", "log2FC.ridge", "log2FC.valley")],
             join_by("Silva.Species" == "taxa")) %>%
  column_to_rownames("Silva.Species") %>%
  select(all_of(c("log2FC.ridge", "log2FC.valley", bacteria.target.gene)))



### 1.4.1 Fungi----
fungi.data <-  read.csv(file.path(path_microbiome,
                                     "fungal.taxa.genome.gene.copy.number.matrix.csv"),
                           header = TRUE) %>%
  group_by(Silva.Genus) %>%
  summarise(across(all_of(fungi.target.gene), ~mean(.x)), .groups = "drop") %>%
  left_join(fungi.wilcox.result_df[, c("Genus", "Species", "log2FC.ridge", "log2FC.valley")],
            join_by("Silva.Genus" == "Genus")) %>%
  column_to_rownames("Species") %>%
  select(all_of(c("log2FC.ridge", "log2FC.valley", fungi.target.gene)))



## 1.5) Hyperparameter Matrix----


param_grid <- expand.grid(
  eta = c(0.02, 0.05, 0.1),
  max_depth = c(2, 3, 4),
  min_child_weight = c(3, 4, 5),
  gamma = c(0.1, 0.2, 0.3),
  subsample = c(0.6, 0.8),
  colsample_bytree = c(0.6, 0.8),
  lambda = c(0.3, 1, 3),
  alpha  = c(0.1, 0.25, 0.4)
)




nrow(param_grid)



# Train the XGBoost model----

# 2. Bacteria.ridge----

# Data preparation
bacteria.ridge.total.data <- bacteria.data %>%
  select(all_of(c("log2FC.ridge", bacteria.target.gene))) %>% 
  rename(Abundance_shift = log2FC.ridge)



## 2.1) Training set and testing set----
set.seed(20250508)
bacteria.ridge.train_idx <- createDataPartition(bacteria.ridge.total.data$Abundance_shift,
                                                p = 0.8, list = FALSE)

# Training set
bacteria.ridge.train_data <- bacteria.ridge.total.data[bacteria.ridge.train_idx, ]

# testing set
bacteria.ridge.test_data <- bacteria.ridge.total.data[-bacteria.ridge.train_idx, ]



## 2.2) Convert to DMatrix----
bacteria.ridge.dtrain <- xgb.DMatrix(
  data = as.matrix(bacteria.ridge.train_data[, bacteria.target.gene]),
  label = bacteria.ridge.train_data$Abundance_shift
)



## 2.3) Hyperparameter tuning----


bacteria.ridge.cv_results <- param_grid %>%
  mutate(best_rmse_min = NA_real_,
         best_rmse_tail10 = NA_real_,
         best_nrounds = NA_integer_)


set.seed(20250508)
for (i in 1:nrow(param_grid)) {
  cat("The", i, "round of cycles,", "Total", nrow(param_grid),"rounds", "\n")
  params <- list(booster = "gbtree",
                 objective = "reg:squarederror",
                 eval_metric="rmse",
                 eta = param_grid$eta[i],
                 max_depth = param_grid$max_depth[i],
                 subsample = param_grid$subsample[i],
                 colsample_bytree = param_grid$colsample_bytree[i],
                 lambda = param_grid$lambda[i],
                 alpha = param_grid$alpha[i],
                 nthread = 4)
  
  cv <- xgb.cv(params = params,
               data = bacteria.ridge.dtrain,
               nrounds = 500,
               nfold = 10,
               early_stopping_rounds = 50,
               verbose = 0)
  
  eval_log <- cv$evaluation_log$test_rmse_mean
  bacteria.ridge.cv_results$best_rmse_min[i] <- min(eval_log)
  bacteria.ridge.cv_results$best_rmse_tail10[i] <- mean(tail(eval_log, 10))
  bacteria.ridge.cv_results$best_nrounds[i] <- cv$best_iteration
}




## 2.4) Extract the optimal parameter----
bacteria.ridge.best_idx <- which.min(bacteria.ridge.cv_results$best_rmse_tail10)
bacteria.ridge.best_params <- param_grid[bacteria.ridge.best_idx, ]
bacteria.ridge.best_nrounds <- bacteria.ridge.cv_results$best_nrounds[bacteria.ridge.best_idx]


## 2.5) Train the optimal model----

set.seed(20250508)
bacteria.ridge.final_model <- xgb.train(
  params = as.list(bacteria.ridge.best_params) %>% c(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse"),
  data = bacteria.ridge.dtrain,
  nrounds = bacteria.ridge.best_nrounds,
  verbose = 1
)



## 2.6) Testing model performance----
bacteria.ridge.pred_test <- predict(bacteria.ridge.final_model,
                                    as.matrix(bacteria.ridge.test_data[, bacteria.target.gene]))

bacteria.ridge.test.result <- data.frame(pred = bacteria.ridge.pred_test,
                                         obs = bacteria.ridge.test_data$Abundance_shift)


## 2.6.1 Fig_S4a----

postResample(pred = bacteria.ridge.test.result$pred, obs = bacteria.ridge.test.result$obs)


# Visualization
bacteria.ridge.test.result %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 7, shape = 21, fill = "#FEA6BF") +
  stat_cor(size = 10) +
  stat_smooth(geom = "smooth", formula = "y ~ x", alpha = 0.5,
              se = T, method = "lm", fill = "#7E243D",
              colour = "#FC215F", linewidth = 1.5) +
  labs(x = "Observed value", y = "Predicted value") +
  scale_x_continuous(limits = c(-5, 1), breaks = seq(-5, 1, 1)) +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 0.5)) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.20, "cm"),
        axis.ticks.length.y = unit(-0.20, "cm"))




## 2.7) SHAP Analysis----

bacteria.ridge.shap_long <- shap.prep(xgb_model = bacteria.ridge.final_model,
                                      X_train = as.matrix(bacteria.ridge.train_data[, bacteria.target.gene])) %>%
  data.frame() %>%
  group_by(variable) %>%
  mutate(across(stdfvalue, ~ (.-min(.)) / (max(.) - min(.))))


# Feature Importance
bacteria.ridge_importance <- bacteria.ridge.shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value))) %>%
  arrange(desc(mean_abs_shap))

head(bacteria.ridge.shap_long)


## 2.8) Fig_S5a----


### 2.8.1 swarm diagram----
bacteria.ridge_SHAP_Swarm_Plot <- bacteria.ridge.shap_long %>%
  ggplot(aes(x = value, y = reorder(variable, mean_value))) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(aes(fill = stdfvalue),
             shape = 21, size = 3.5,
             stroke = 0, alpha = 0.4, color = "black",
             position = position_jitter(height = 0.06, seed = 1234)) +
  scale_fill_gradientn(colours = c("blue", "red"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       name = "Feature value",
                       breaks = c(0, 1),
                       labels = c("Low", "High"),
                       guide = guide_colourbar(
                         title.position = "bottom",
                         title.hjust = 0.5,
                         title.vjust = 4,
                         barwidth = unit(8, "cm"),
                         barheight = unit(0.2, "cm"))) +
  labs(x = "SHAP value (Impact on model output)",
       y = NULL, fill = NULL,
       title = "Bacteria ridge") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x.bottom = element_text(vjust = -1),
        legend.box.margin = margin(r = 10),
        plot.margin = margin(10, 0, 5, 5))


### 2.8.2 Variable Importance----
bacteria.ridge_Feature_Importance_Plot <- bacteria.ridge_importance %>%
  ggplot(aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap))) +
  geom_col(fill = "grey80", width = 0.6, color = "black") +
  scale_x_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05),
                     expand =  expansion(0, 0), position = "top") +
  labs(x = "Mean |SHAP|", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),   # x 方向刻度线
        plot.margin = margin(10, 10, 5, 0),
        axis.title.x.top = element_text(vjust = 2))


### 2.8.3 Merge Graphics----
bacteria.ridge_SHAP_Swarm_Plot + 
  bacteria.ridge_Feature_Importance_Plot + 
  plot_layout(widths = c(3.5, 1.5))


bacteria.ridge.important.variable.top6 <- bacteria.ridge_importance$variable[1:6]


## Fig_4a----
bacteria.ridge.shap_long %>%
  subset(variable %in%  bacteria.ridge.important.variable.top6) %>%
  ggplot(aes(x = rfvalue, y = value)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1) +
  geom_point(size = 7, shape = 21, fill = "gray90") +
  geom_smooth(method = "glm", span = 1, formula = y ~ log(x + 1e-6),
              linewidth = 2, alpha = 1) + 
  scale_x_continuous(expand = expansion(mult = 0.1, add = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0.1)) +
  facet_wrap(~ variable, scale = "free", nrow = 1) +
  labs(x = "Gene copy number", y = "SHAP value") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))



## 2.9) Save cross-validation results----
save(bacteria.ridge.cv_results,
     file = file.path(path_microbiome, 
                      "xgboost.parameter.cv_results.Rdata"))



# 3. Bacteria.valley----

# Data preparation
bacteria.valley.total.data <- bacteria.data %>%
  select(all_of(c("log2FC.valley", bacteria.target.gene))) %>% 
  rename(Abundance_shift = log2FC.valley)



## 3.1) Training set and testing set----
set.seed(20250508)
bacteria.valley.train_idx <- createDataPartition(bacteria.valley.total.data$Abundance_shift,
                                                 p = 0.8, list = FALSE)

# Training set
bacteria.valley.train_data <- bacteria.valley.total.data[bacteria.valley.train_idx, ]

# testing set
bacteria.valley.test_data <- bacteria.valley.total.data[-bacteria.valley.train_idx, ]



## 3.2) Convert to DMatrix----
bacteria.valley.dtrain <- xgb.DMatrix(
  data = as.matrix(bacteria.valley.train_data[, bacteria.target.gene]),
  label = bacteria.valley.train_data$Abundance_shift
)



## 3.3) Hyperparameter tuning----


bacteria.valley.cv_results <- param_grid %>%
  mutate(best_rmse_min = NA_real_,
         best_rmse_tail10 = NA_real_,
         best_nrounds = NA_integer_)


set.seed(20250508)
for (i in 1:nrow(param_grid)) {
  cat("The", i, "round of cycles,", "Total", nrow(param_grid),"rounds", "\n")
  params <- list(booster = "gbtree",
                 objective = "reg:squarederror",
                 eval_metric="rmse",
                 eta = param_grid$eta[i],
                 max_depth = param_grid$max_depth[i],
                 subsample = param_grid$subsample[i],
                 colsample_bytree = param_grid$colsample_bytree[i],
                 lambda = param_grid$lambda[i],
                 alpha = param_grid$alpha[i],
                 nthread = 4)
  
  cv <- xgb.cv(params = params,
               data = bacteria.valley.dtrain,
               nrounds = 500,
               nfold = 10,
               early_stopping_rounds = 50,
               verbose = 0)
  
  eval_log <- cv$evaluation_log$test_rmse_mean
  bacteria.valley.cv_results$best_rmse_min[i] <- min(eval_log)
  bacteria.valley.cv_results$best_rmse_tail10[i] <- mean(tail(eval_log, 10))
  bacteria.valley.cv_results$best_nrounds[i] <- cv$best_iteration
}




## 3.4) Extract the optimal parameter----
bacteria.valley.best_idx <- which.min(bacteria.valley.cv_results$best_rmse_tail10)
bacteria.valley.best_params <- param_grid[bacteria.valley.best_idx, ]
bacteria.valley.best_nrounds <- bacteria.valley.cv_results$best_nrounds[bacteria.valley.best_idx]


## 3.5) Train the optimal model----

set.seed(20250508)
bacteria.valley.final_model <- xgb.train(
  params = as.list(bacteria.valley.best_params) %>% c(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse"),
  data = bacteria.valley.dtrain,
  nrounds = bacteria.valley.best_nrounds,
  verbose = 1
)



## 3.6) Testing model performance----
bacteria.valley.pred_test <- predict(bacteria.valley.final_model,
                                     as.matrix(bacteria.valley.test_data[, bacteria.target.gene]))

bacteria.valley.test.result <- data.frame(pred = bacteria.valley.pred_test,
                                          obs = bacteria.valley.test_data$Abundance_shift)


### 3.6.1 Fig_S4b----

postResample(pred = bacteria.valley.test.result$pred, obs = bacteria.valley.test.result$obs)


# Visualization
bacteria.valley.test.result %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 7, shape = 21, fill = "#FEA6BF") +
  stat_cor(size = 10) +
  stat_smooth(geom = "smooth", formula = "y ~ x", alpha = 0.5,
              se = T, method = "lm", fill = "#7E243D",
              colour = "#FC215F", linewidth = 1.5) +
  labs(x = "Observed value", y = "Predicted value") +
  scale_x_continuous(limits = c(-1, 3), breaks = seq(-1, 3, 1)) +
  scale_y_continuous(limits = c(-0.5, 1), breaks = seq(-0.5, 1, 0.25)) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.20, "cm"),
        axis.ticks.length.y = unit(-0.20, "cm"))




## 3.7) SHAP Analysis----

bacteria.valley.shap_long <- shap.prep(xgb_model = bacteria.valley.final_model,
                                       X_train = as.matrix(bacteria.valley.train_data[, bacteria.target.gene])) %>%
  data.frame() %>%
  group_by(variable) %>%
  mutate(across(stdfvalue, ~ (.-min(.)) / (max(.) - min(.))))


# Feature Importance
bacteria.valley_importance <- bacteria.valley.shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value))) %>%
  arrange(desc(mean_abs_shap))

head(bacteria.valley.shap_long)


## 3.8) Fig_S5b----


### 3.8.1 swarm diagram----
bacteria.valley_SHAP_Swarm_Plot <- bacteria.valley.shap_long %>%
  ggplot(aes(x = value, y = reorder(variable, mean_value))) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(aes(fill = stdfvalue),
             shape = 21, size = 3.5,
             stroke = 0, alpha = 0.4, color = "black",
             position = position_jitter(height = 0.06, seed = 1234)) +
  scale_fill_gradientn(colours = c("blue", "red"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       name = "Feature value",
                       breaks = c(0, 1),
                       labels = c("Low", "High"),
                       guide = guide_colourbar(
                         title.position = "bottom",
                         title.hjust = 0.5,
                         title.vjust = 4,
                         barwidth = unit(8, "cm"),
                         barheight = unit(0.2, "cm"))) +
  labs(x = "SHAP value (Impact on model output)",
       y = NULL, fill = NULL,
       title = "Bacteria valley") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x.bottom = element_text(vjust = -1),
        legend.box.margin = margin(r = 10),
        plot.margin = margin(10, 0, 5, 5))


### 3.8.2 Variable Importance----
bacteria.valley_Feature_Importance_Plot <- bacteria.valley_importance %>%
  ggplot(aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap))) +
  geom_col(fill = "grey80", width = 0.6, color = "black") +
  scale_x_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, 0.05),
                     expand =  expansion(0, 0), position = "top") +
  labs(x = "Mean |SHAP|", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),   # x 方向刻度线
        plot.margin = margin(10, 10, 5, 0),
        axis.title.x.top = element_text(vjust = 2))


### 3.8.3 Merge Graphics----
bacteria.valley_SHAP_Swarm_Plot + 
  bacteria.valley_Feature_Importance_Plot + 
  plot_layout(widths = c(3.5, 1.5))


### Fig_4b----
bacteria.valley.important.variable.top6 <- bacteria.valley_importance$variable[1:6]

bacteria.valley.shap_long %>%
  subset(variable %in%  bacteria.valley.important.variable.top6) %>%
  ggplot(aes(x = rfvalue, y = value)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1) +
  geom_point(size = 7, shape = 21, fill = "gray90") +
  geom_smooth(method = "glm", span = 1, formula = y ~ log(x + 1e-6),
              linewidth = 2, alpha = 1) + 
  scale_x_continuous(expand = expansion(mult = 0.1, add = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1, add = 0.1)) +
  facet_wrap(~ variable, scale = "free", nrow = 1) +
  labs(x = "Gene copy number", y = "SHAP value") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length = unit(-0.15, "cm"))


## 3.9) Save cross-validation results----
save(bacteria.ridge.cv_results,
     bacteria.valley.cv_results,
     file = file.path(path_microbiome, 
                      "xgboost.parameter.cv_results.Rdata"))




# 4. Fungi.ridge----

# Data preparation
fungi.ridge.total.data <- fungi.data %>%
  select(all_of(c("log2FC.ridge", fungi.target.gene))) %>% 
  rename(Abundance_shift = log2FC.ridge)



## 4.1) Training set and testing set----
set.seed(20250508)
fungi.ridge.train_idx <- createDataPartition(fungi.ridge.total.data$Abundance_shift,
                                             p = 0.8, list = FALSE)

# Training set
fungi.ridge.train_data <- fungi.ridge.total.data[fungi.ridge.train_idx, ]

# testing set
fungi.ridge.test_data <- fungi.ridge.total.data[-fungi.ridge.train_idx, ]



## 4.2) Convert to DMatrix----
fungi.ridge.dtrain <- xgb.DMatrix(
  data = as.matrix(fungi.ridge.train_data[, fungi.target.gene]),
  label = fungi.ridge.train_data$Abundance_shift
)



## 4.3) Hyperparameter tuning----



fungi.ridge.cv_results <- param_grid %>%
  mutate(best_rmse_min = NA_real_,
         best_rmse_tail10 = NA_real_,
         best_nrounds = NA_integer_)


set.seed(20250508)
for (i in 1:nrow(param_grid)) {
  cat("The", i, "round of cycles,", "Total", nrow(param_grid),"rounds", "\n")
  params <- list(booster = "gbtree",
                 objective = "reg:squarederror",
                 eval_metric="rmse",
                 eta = param_grid$eta[i],
                 max_depth = param_grid$max_depth[i],
                 subsample = param_grid$subsample[i],
                 colsample_bytree = param_grid$colsample_bytree[i],
                 lambda = param_grid$lambda[i],
                 alpha = param_grid$alpha[i],
                 nthread = 4)
  
  cv <- xgb.cv(params = params,
               data = fungi.ridge.dtrain,
               nrounds = 500,
               nfold = 10,
               early_stopping_rounds = 50,
               verbose = 0)
  
  eval_log <- cv$evaluation_log$test_rmse_mean
  fungi.ridge.cv_results$best_rmse_min[i] <- min(eval_log)
  fungi.ridge.cv_results$best_rmse_tail10[i] <- mean(tail(eval_log, 10))
  fungi.ridge.cv_results$best_nrounds[i] <- cv$best_iteration
}




## 4.4) Extract the optimal parameter----
fungi.ridge.best_idx <- which.min(fungi.ridge.cv_results$best_rmse_tail10)
fungi.ridge.best_params <- param_grid[fungi.ridge.best_idx, ]
fungi.ridge.best_nrounds <- fungi.ridge.cv_results$best_nrounds[fungi.ridge.best_idx]


## 4.5) Train the optimal model----

set.seed(20250508)
fungi.ridge.final_model <- xgb.train(
  params = as.list(fungi.ridge.best_params) %>% c(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse"),
  data = fungi.ridge.dtrain,
  nrounds = fungi.ridge.best_nrounds,
  verbose = 1
)



## 4.6) Testing model performance----
fungi.ridge.pred_test <- predict(fungi.ridge.final_model,
                                 as.matrix(fungi.ridge.test_data[, fungi.target.gene]))

fungi.ridge.test.result <- data.frame(pred = fungi.ridge.pred_test,
                                      obs = fungi.ridge.test_data$Abundance_shift)




### 4.6.1 Fig_S4c----
postResample(pred = fungi.ridge.test.result$pred, obs = fungi.ridge.test.result$obs)


# Visualization
fungi.ridge.test.result %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 7, shape = 21, fill = "#F3BAF2") +
  stat_cor(size = 10) +
  stat_smooth(geom = "smooth", formula = "y ~ x", alpha = 0.5,
              se = T, method = "lm", fill = "#792D78",
              colour = "#E053DF", linewidth = 1.5) +
  labs(x = "Observed value", y = "Predicted value") +
  scale_x_continuous(limits = c(-10, 8), breaks = seq(-10, 8, 2)) +
  scale_y_continuous(limits = c(-2, 1), breaks = seq(-2, 1, 0.5)) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.20, "cm"),
        axis.ticks.length.y = unit(-0.20, "cm"))





## 4.7) SHAP Analysis----

fungi.ridge.shap_long <- shap.prep(xgb_model = fungi.ridge.final_model,
                                   X_train = as.matrix(fungi.ridge.train_data[, fungi.target.gene])) %>%
  data.frame() %>%
  group_by(variable) %>%
  mutate(across(stdfvalue, ~ (.-min(.)) / (max(.) - min(.))))


# Feature Importance
fungi.ridge_importance <- fungi.ridge.shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value))) %>%
  arrange(desc(mean_abs_shap))

head(fungi.ridge.shap_long)


## 4.8) SHAP Visualization----


### 4.8.1 swarm diagram----
fungi.ridge_SHAP_Swarm_Plot <- fungi.ridge.shap_long %>%
  ggplot(aes(x = value, y = reorder(variable, mean_value))) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(aes(fill = stdfvalue),
             shape = 21, size = 3.5,
             stroke = 0, alpha = 0.4, color = "black",
             position = position_jitter(height = 0.06, seed = 1234)) +
  scale_fill_gradientn(colours = c("blue", "red"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       name = "Feature value",
                       breaks = c(0, 1),
                       labels = c("Low", "High"),
                       guide = guide_colourbar(
                         title.position = "bottom",
                         title.hjust = 0.5,
                         title.vjust = 4,
                         barwidth = unit(8, "cm"),
                         barheight = unit(0.2, "cm"))) +
  labs(x = "SHAP value (Impact on model output)",
       y = NULL, fill = NULL,
       title = "Fungi ridge") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x.bottom = element_text(vjust = -1),
        legend.box.margin = margin(r = 10),
        plot.margin = margin(10, 0, 5, 5))


### 4.8.2 Variable Importance----
fungi.ridge_Feature_Importance_Plot <- fungi.ridge_importance %>%
  ggplot(aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap))) +
  geom_col(fill = "grey80", width = 0.6, color = "black") +
  scale_x_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.05),
                     expand =  expansion(0, 0), position = "top") +
  labs(x = "Mean |SHAP|", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),   # x 方向刻度线
        plot.margin = margin(10, 10, 5, 0),
        axis.title.x.top = element_text(vjust = 2))


### 4.8.3 Merge Graphics----
fungi.ridge_SHAP_Swarm_Plot + 
  fungi.ridge_Feature_Importance_Plot + 
  plot_layout(widths = c(3.5, 1.5))


## 4.9) Save cross-validation results----
save(bacteria.ridge.cv_results,
     bacteria.valley.cv_results,
     fungi.ridge.cv_results,
     file = file.path(path_microbiome, 
                      "xgboost.parameter.cv_results.Rdata"))



# 5. Fungi.valley----

# Data preparation
fungi.valley.total.data <- fungi.data %>%
  select(all_of(c("log2FC.valley", fungi.target.gene))) %>% 
  rename(Abundance_shift = log2FC.valley)



## 5.1) Training set and testing set----
set.seed(20250508)
fungi.valley.train_idx <- createDataPartition(fungi.valley.total.data$Abundance_shift,
                                              p = 0.8, list = FALSE)

# Training set
fungi.valley.train_data <- fungi.valley.total.data[fungi.valley.train_idx, ]

# testing set
fungi.valley.test_data <- fungi.valley.total.data[-fungi.valley.train_idx, ]



## 5.2) Convert to DMatrix----
fungi.valley.dtrain <- xgb.DMatrix(
  data = as.matrix(fungi.valley.train_data[, fungi.target.gene]),
  label = fungi.valley.train_data$Abundance_shift
)



## 5.3) Hyperparameter tuning----


fungi.valley.cv_results <- param_grid %>%
  mutate(best_rmse_min = NA_real_,
         best_rmse_tail10 = NA_real_,
         best_nrounds = NA_integer_)


set.seed(20250508)
for (i in 1:nrow(param_grid)) {
  cat("The", i, "round of cycles,", "Total", nrow(param_grid),"rounds", "\n")
  params <- list(booster = "gbtree",
                 objective = "reg:squarederror",
                 eval_metric="rmse",
                 eta = param_grid$eta[i],
                 max_depth = param_grid$max_depth[i],
                 subsample = param_grid$subsample[i],
                 colsample_bytree = param_grid$colsample_bytree[i],
                 lambda = param_grid$lambda[i],
                 alpha = param_grid$alpha[i],
                 nthread = 4)
  
  cv <- xgb.cv(params = params,
               data = fungi.valley.dtrain,
               nrounds = 500,
               nfold = 10,
               early_stopping_rounds = 50,
               verbose = 0)
  
  eval_log <- cv$evaluation_log$test_rmse_mean
  fungi.valley.cv_results$best_rmse_min[i] <- min(eval_log)
  fungi.valley.cv_results$best_rmse_tail10[i] <- mean(tail(eval_log, 10))
  fungi.valley.cv_results$best_nrounds[i] <- cv$best_iteration
}




## 5.4) Extract the optimal parameter----
fungi.valley.best_idx <- which.min(fungi.valley.cv_results$best_rmse_tail10)
fungi.valley.best_params <- param_grid[fungi.valley.best_idx, ]
fungi.valley.best_nrounds <- fungi.valley.cv_results$best_nrounds[fungi.valley.best_idx]


## 5.5) Train the optimal model----

set.seed(20250508)
fungi.valley.final_model <- xgb.train(
  params = as.list(fungi.valley.best_params) %>% c(
    booster = "gbtree",
    objective = "reg:squarederror",
    eval_metric = "rmse"),
  data = fungi.valley.dtrain,
  nrounds = fungi.valley.best_nrounds,
  verbose = 1
)



## 5.6) Testing model performance----
fungi.valley.pred_test <- predict(fungi.valley.final_model,
                                  as.matrix(fungi.valley.test_data[, fungi.target.gene]))

fungi.valley.test.result <- data.frame(pred = fungi.valley.pred_test,
                                       obs = fungi.valley.test_data$Abundance_shift)



### 5.6.1 Fig_S4d----

postResample(pred = fungi.valley.test.result$pred, obs = fungi.valley.test.result$obs)


# Visualization
fungi.valley.test.result %>%
  ggplot(aes(x = obs, y = pred)) +
  geom_point(size = 7, shape = 21, fill = "#F3BAF2") +
  stat_cor(size = 10) +
  stat_smooth(geom = "smooth", formula = "y ~ x", alpha = 0.5,
              se = T, method = "lm", fill = "#792D78",
              colour = "#E053DF", linewidth = 1.5) +
  labs(x = "Observed value", y = "Predicted value") +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, 2)) +
  scale_y_continuous(limits = c(-0.8, 1.2), breaks = seq(-0.8, 1.2, 0.4)) +
  theme_bw(base_size = 16) +
  theme(panel.border = element_rect(colour = "black",
                                    linewidth = 1,
                                    fill = NA),
        panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.20, "cm"),
        axis.ticks.length.y = unit(-0.20, "cm"))





## 5.7) SHAP Analysis----

fungi.valley.shap_long <- shap.prep(xgb_model = fungi.valley.final_model,
                                    X_train = as.matrix(fungi.valley.train_data[, fungi.target.gene])) %>%
  data.frame() %>%
  group_by(variable) %>%
  mutate(across(stdfvalue, ~ (.-min(.)) / (max(.) - min(.))))


# Feature Importance
fungi.valley_importance <- fungi.valley.shap_long %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = mean(abs(value))) %>%
  arrange(desc(mean_abs_shap))

head(fungi.valley.shap_long)


## 5.8) SHAP Visualization----


### 5.8.1 swarm diagram----
fungi.valley_SHAP_Swarm_Plot <- fungi.valley.shap_long %>%
  ggplot(aes(x = value, y = reorder(variable, mean_value))) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  geom_point(aes(fill = stdfvalue),
             shape = 21, size = 3.5,
             stroke = 0, alpha = 0.4, color = "black",
             position = position_jitter(height = 0.06, seed = 1234)) +
  scale_fill_gradientn(colours = c("blue", "red"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       name = "Feature value",
                       breaks = c(0, 1),
                       labels = c("Low", "High"),
                       guide = guide_colourbar(
                         title.position = "bottom",
                         title.hjust = 0.5,
                         title.vjust = 4,
                         barwidth = unit(8, "cm"),
                         barheight = unit(0.2, "cm"))) +
  labs(x = "SHAP value (Impact on model output)",
       y = NULL, fill = NULL,
       title = "Fungi valley") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x.bottom = element_text(vjust = -1),
        legend.box.margin = margin(r = 10),
        plot.margin = margin(10, 0, 5, 5))


### 5.8.2 Variable Importance----
fungi.valley_Feature_Importance_Plot <- fungi.valley_importance %>%
  ggplot(aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap))) +
  geom_col(fill = "grey80", width = 0.6, color = "black") +
  scale_x_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01),
                     expand =  expansion(0, 0), position = "top") +
  labs(x = "Mean |SHAP|", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),   # x 方向刻度线
        plot.margin = margin(10, 10, 5, 0),
        axis.title.x.top = element_text(vjust = 2))


### 5.8.3 Merge Graphics----
fungi.valley_SHAP_Swarm_Plot + 
  fungi.valley_Feature_Importance_Plot + 
  plot_layout(widths = c(3.5, 1.5))


## 5.9) Save cross-validation results----
save(bacteria.data,
     fungi.data,
     bacteria.ridge.cv_results,
     bacteria.valley.cv_results,
     fungi.ridge.cv_results,
     fungi.valley.cv_results,
     file = file.path(path_microbiome, 
                      "xgboost.parameter.cv_results.Rdata"))


# 6. Difference analysis of gene abundance----

## 6.1) gene.catalog----

gene.catalog <- 
  read_excel(file.path(path_microbiome,
                       "carbon_acquisition_TEA_utilization_gene_catalog.xlsx")) %>%
  select(KO, Type) %>%
  data.frame()


gene.list <- c("Exo_amylase", "Exo_cellulase", "Exo_hemicellulose", "Exo_chitosanase",
               "Endo_amylase", "Endo_cellulase", "Endo_hemicellulose", "Endo_chitosanase",
               "Ligninase", "Oxygen", "Nitrate", "Sulfite")



## 6.2) KEGG annotation list----
KEGG.annotation.list <- read.table(file.path(path_microbiome, "gene.ko.annotation.txt"),
                                   header = T)


## 6.3) Genes count table----

genes.count.cpm.table <- read.table(
  file.path(path_microbiome, "genes_count_table.txt"), header = TRUE) %>%
  column_to_rownames("SeqID") %>%
  cpm() %>% data.frame() %>%
  mutate(SeqID = row.names(.))


carbon.degrading.gene.count <- genes.count.cpm.table %>%
  subset(SeqID %in% KEGG.annotation.list$SeqID) %>%
  left_join(KEGG.annotation.list, by = "SeqID") %>%
  subset(KO %in% gene.catalog$KO) %>%
  left_join(gene.catalog, by = "KO") %>%
  group_by(Type) %>%
  summarise(across(sample.info$Sample.ID, ~sum(.x)), .groups = "drop") %>%
  data.frame() %>%
  column_to_rownames("Type")





## 6.4) Ridge----

group.ridge <- sample.info %>%
  subset(Topography == "Ridge") %>%
  pull(Status)

paired.ridge <- sample.info %>%
  subset(Topography == "Ridge") %>%
  pull(Pairs) %>%
  factor() %>%
  as.character()


design.ridge <- model.matrix(~ group.ridge + paired.ridge)


count.table.ridge <- carbon.degrading.gene.count %>%
  select(all_of(sample.ridge))


### 6.4.1 Result----
genes.count.diff.ridge <- DGEList(counts = count.table.ridge,
                                  samples = sample.ridge,
                                  group = group.ridge) %>%
  estimateDisp(design = design.ridge) %>%
  glmFit() %>%
  glmLRT(coef = "group.ridgeDead") %>%
  topTags(sort.by = "none", n = nrow(count.table.ridge),
          adjust.method = "BH") %>%
  data.frame() %>%
  dplyr::slice(match(gene.list, row.names(.)))


genes.count.diff.ridge


## 6.5) Valley----

group.valley <- sample.info %>%
  subset(Topography == "Valley") %>%
  pull(Status)

paired.valley <- sample.info %>%
  subset(Topography == "Valley") %>%
  pull(Pairs) %>%
  factor() %>%
  as.character()


design.valley <- model.matrix(~ group.valley + paired.valley)


count.table.valley <- carbon.degrading.gene.count %>%
  select(all_of(sample.valley))


### 6.5.1 Results----

genes.count.diff.valley <- DGEList(counts = count.table.valley,
                                   samples = sample.valley,
                                   group = group.valley) %>%
  estimateDisp(design = design.valley) %>%
  glmFit() %>%
  glmLRT(coef = "group.valleyDead") %>%
  topTags(sort.by = "none", n = nrow(count.table.valley),
          adjust.method = "BH") %>%
  data.frame() %>%
  dplyr::slice(match(gene.list, row.names(.)))


genes.count.diff.ridge
genes.count.diff.valley



# 7 Fig_4c----

## 7.1) Heatmap----


plot.data <- data.frame(Ridge = genes.count.diff.ridge$logFC,
                        Valley = genes.count.diff.valley$logFC,
                        row.names = row.names(genes.count.diff.ridge)) %>%
  mutate(across(c("Ridge", "Valley"), ~ round(.x, 3))) %>%
  subset(row.names(.) %in% gene.list) %>%
  dplyr::slice(match(gene.list, row.names(.)))




plot.signif <- data.frame(Ridge = genes.count.diff.ridge$FDR,
                          Valley = genes.count.diff.valley$FDR,
                          row.names = row.names(genes.count.diff.ridge)) %>%
  mutate(across(c("Ridge", "Valley"), ~ ifelse(.x < 0.05, "*", ""))) %>%
  subset(row.names(.) %in% gene.list) %>%
  dplyr::slice(match(gene.list, row.names(.)))



pheatmap(t(plot.data), 
         display_numbers = t(plot.signif),  
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")[-c(4, 5, 7, 8)]))(99),  # 自定义颜色
         fontsize_number = 30,
         border_color = "#606060",
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-0.6, 0.6,length.out = 100),
         legend_breaks = seq(-0.6, 0.6, 0.2))


pheatmap(t(plot.data), 
         display_numbers = t(plot.data),  # 在每个格子中显示显著性标记
         color = colorRampPalette(rev(brewer.pal(11, "PiYG")[-c(4, 5, 7, 8)]))(99),  # 自定义颜色
         fontsize_number = 16,
         border_color = "#606060",
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-0.6, 0.6,length.out = 100),
         legend_breaks = seq(-0.6, 0.6, 0.2))



# 8. save results----

save(carbon.degrading.gene.count,
     genes.count.diff.ridge, genes.count.diff.valley,
     file = file.path(path_output, "carbon.degrading.gene.count.Rdata"))


