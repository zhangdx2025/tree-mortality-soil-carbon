
# Set the working directory
setwd("G:/R_script")

path_properties <- file.path(getwd(), "Data/Soil_properties")


# Loading package

library(readxl)
library(tidyverse)
library(lme4)
library(performance)
library(car)
library(emmeans)
library(clubSandwich)



# 1. Loading Data----

soil.properties <- read_excel(file.path(path_properties,
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
  select(Sample.ID, Pairs, Status, Topography, SOC, POC, MAOC, Moisture,
         mass_normalized_activity, SOC_normalized_activity) %>%
  data.frame() 



## Specify variable name

variable_names <- c("SOC", "POC", "MAOC", "Moisture",
                    "mass_normalized_activity", "SOC_normalized_activity")




# 2. Normality and Heteroscedasticity test----

for (i in variable_names) {
  
  df <- soil.properties %>%
    select(all_of(c("Sample.ID", "Pairs", "Status", "Topography", i))) %>%
    rename("y" = i)
  

  print(paste("Normality test:", i))
  df %>%
    lmer(y ~ Status * Topography + (1 | Pairs), data = .) %>%
    check_normality() %>%
    print()
  
  
  print(paste("Heteroscedasticity test:", i))
  df %>%
    lmer(y ~ Status * Topography + (1 | Pairs), data = .) %>%
    check_heteroscedasticity() %>%
    print()
  }


# Table_S1----
# 2. Main effect analysis----

for (i in variable_names) {
  
  df <- soil.properties %>%
    select(all_of(c("Sample.ID", "Pairs", "Status", "Topography", i))) %>%
    rename("y" = i)

  print(i)
  
  if(i == "POC") {
    
    vcovCR <- df %>%
      lmer(y ~ Status * Topography + (1 | Pairs), data = .) %>%
      vcovCR(type = "CR2", cluster = df$Pairs)
    
    df %>%
      lmer(y ~ Status * Topography + (1 | Pairs), data = .) %>%
      Anova(vcov. = vcovCR) %>%
      data.frame() %>%
      add_column(Term = row.names(.), .before = 1) %>%
      {row.names(.) = NULL;.} %>%
      setNames(c("Term", "Chi_squared", "Degrees of freedom", "Pr_Chisq")) %>%
      mutate(Pr_Chisq = ifelse(Pr_Chisq < 0.001, "p < 0.001", round(Pr_Chisq, 3)),
             Chi_squared  = round(Chi_squared, 4))  %>%
      data.frame() %>%
      add_column(Enzymes = i, .before = 1) %>%
      print()
  } else {
    
    print(i)
    df %>%
      lmer(y ~ Status * Topography + (1 | Pairs), data = .) %>%
      Anova() %>%
      data.frame() %>%
      add_column(Term = row.names(.), .before = 1) %>%
      {row.names(.) = NULL;.} %>%
      setNames(c("Term", "Chi_squared", "Degrees of freedom", "Pr_Chisq")) %>%
      mutate(Pr_Chisq = ifelse(Pr_Chisq < 0.001, "p < 0.001", round(Pr_Chisq, 3)),
             Chi_squared  = round(Chi_squared, 4))  %>%
      data.frame() %>%
      add_column(Enzymes = i, .before = 1) %>%
      print()
    }}
 
 

