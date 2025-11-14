# Set the working directory
setwd("G:/R_script")


path_microbiome <- file.path(getwd(), "Data/microbiome_profiles")


# load packages
library(tidyverse)
library(readxl)




target.gene <- c("Exo_amylase", "Endo_amylase", "Exo_cellulase", "Endo_cellulase",
                 "Exo_hemicellulose", "Endo_hemicellulose", "Exo_chitosanase",
                 "Endo_chitosanase", "Ligninase", "TEA_Oxygen", "TEA_Nitrate",
                 "TEA_Sulfite")


# 1. Baterial intraspecific variation----


bacteria.data <- read.csv(file.path(path_microbiome,
                                    "Bacteria_for_gene_copy_number_variation.csv"),
                          header = TRUE)


nrow(bacteria.data)
bacteria.data$Ncbi.Species %>% unique() %>% length()




Species_results <- list()
for (i in seq_along(target.gene)) {
  
  df <- bacteria.data %>%
    select(all_of(c("Ncbi.Species", "Strain", target.gene[i]))) %>%
    rename("Gene" = target.gene[i],
           "Species" = "Ncbi.Species") %>%
    data.frame()
  
  
  formula <- as.formula(paste("Gene", "~ Species"))
  model <- aov(formula, data = df)
  
  between_ss <- summary(model)[[1]][[1, "Mean Sq"]]
  aov.residual <- summary(model)[[1]][[2, "Mean Sq"]]
  
  Species_results[[i]] <- data.frame(Between_ss = between_ss,
                                     Residual = aov.residual)
  
  names(Species_results)[i] <- target.gene[i]
  
}



do.call(rbind, Species_results) %>% data.frame() %>%
  mutate(Winthin_Total = Residual / (Between_ss + Residual),
         Between_ss = round(Between_ss, 2),
         Residual = round(Residual, 2),
         Winthin_Total = round(Winthin_Total, 4))



# 2. Fungal intraspecific variation----

fungi.data <- read.csv(file.path(path_microbiome,
                                 "Fungi_for_gene_copy_number_variation.csv"),
                       header = TRUE)


nrow(fungi.data)
fungi.data$Ncbi.Genus %>% unique() %>% length()
fungi.data$Ncbi.Species %>% unique() %>% length()



Genus_results <- list()
for (i in seq_along(target.gene)) {
  
  df <- fungi.data %>%
    select(all_of(c("Ncbi.Genus", "Ncbi.Species", "Strain", target.gene[i]))) %>%
    rename("Gene" = target.gene[i],
           "Genus" = "Ncbi.Genus",
           "Species" = "Ncbi.Species") %>%
    group_by(Genus, Species) %>%
    summarise(Gene = mean(Gene), .groups = "drop") %>%
    data.frame()
  
  
  formula <- as.formula(paste("Gene", "~ Genus"))
  fit <- aov(formula, data = df)
  
  between_genus <- summary(fit)[[1]][[1, "Mean Sq"]]
  within_genus <- summary(fit)[[1]][[2, "Mean Sq"]]
  
  Genus_results[[i]] <- data.frame(between_genus = between_genus,
                                   within_genus = within_genus)
  
  names(Genus_results)[i] <- target.gene[i]
  
}




do.call(rbind, Genus_results) %>% 
  data.frame(check.names = FALSE) %>%
  mutate(Within_genus_total = (within_genus / (within_genus + between_genus)),
         between_genus = round(between_genus, 2),
         within_genus = round(within_genus, 2),
         Within_genus_total = round(Within_genus_total, 4)
  )


