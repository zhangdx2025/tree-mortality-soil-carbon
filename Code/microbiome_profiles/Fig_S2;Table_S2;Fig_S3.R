
# Set the working directory
setwd("G:/R_script")


path_properties <- file.path(getwd(), "Data/Soil_properties")
path_microbiome <- file.path(getwd(), "Data/microbiome_profiles")
path_output <- file.path(getwd(), "Output")


# load packages
library(tidyverse)
library(readxl)
library(vegan)
library(smatr)
library(ggpmisc)

rename_unclassified <- function(vec, colname) {
  is_unclassified <- vec == "Unclassified"
  vec[is_unclassified] <- paste0("Unc.", colname,
                                 seq_len(sum(is_unclassified)))
  return(vec)
}



# 1. Loading Data----

## 1.1) Taxonomic.level----
Taxonomic.level <- c("Kingdom", "Phylum", "Class", "Order",
                     "Family", "Genus", "Species")


## 1.2) soil.properties----
soil.properties <- read_excel(file.path(path_properties,
                                        "soil.properties.xlsx"),
                              sheet = "Data") %>%
  arrange(Topography, Pairs, Status) %>%
  column_to_rownames("Sample.ID")


## 1.3) Bacteria----

### 1.3.1 Relative.abundance----

Relative.abu.Bac <- read.csv(file.path(path_microbiome, "FeatureTable_16S.txt"),
                             header = TRUE, sep = "\t") %>%
  separate(taxonomy, into = Taxonomic.level, sep = ";\\s*", fill = "right") %>%
  mutate(across(Kingdom:Species, ~str_remove(., "^[a-z]__"))) %>%
  filter(Kingdom == "Bacteria") %>%
  mutate(
    Phylum  = rename_unclassified(Phylum, "Phylum"),
    Class   = rename_unclassified(Class, "Class"),
    Order   = rename_unclassified(Order, "Order"),
    Family  = rename_unclassified(Family, "Family"),
    Genus   = rename_unclassified(Genus, "Genus"),
    Species = rename_unclassified(Species, "Species"),
    Species = if_else(Species == "Rhodoplanes_piscinae", "Rhodoplanes_serenus", Species)
  ) %>%
  group_by(across(all_of(Taxonomic.level))) %>%
  summarise(across(all_of(row.names(soil.properties)), sum, na.rm = TRUE), .groups = "drop") %>%
  mutate(ASV.ID = paste("ASV_B", row_number(), sep = "")) %>%
  relocate(ASV.ID, .before = 1) %>%
  mutate(across(where(is.numeric),
                ~ decostand(., method = "total", MARGIN = 2))) %>%
  data.frame()


Taxonomic.levels <- c("Phylum","Class", "Order", "Family","Genus","Species")


### 1.3.1 Abundance filtering----
bacteria.table <- Relative.abu.Bac %>%
  filter(rowMeans(select(., where(is.numeric))) > 0.0001,
         rowSums(across(where(is.numeric), ~ .x > 0)) > 3) %>%
  data.frame()


## 1.4) Fungi----

### 1.4.1 Relative.abundance----
Relative.abu.Fungi <- read.csv(file.path(path_microbiome,
                                         "FeatureTable_ITS.txt"),
                               header = TRUE, sep = "\t") %>%
  separate(taxonomy,
           into = Taxonomic.level,
           sep = ";\\s*", fill = "right") %>%
  mutate(across(Kingdom:Species, ~str_remove(., "^[a-z]__"))) %>%
  subset(Kingdom == "Fungi") %>%
  mutate(Phylum = rename_unclassified(Phylum, "Phylum"),
         Class = rename_unclassified(Class, "Class"),
         Order = rename_unclassified(Order, "Order"),
         Family = rename_unclassified(Family, "Family"),
         Genus = rename_unclassified(Genus, "Genus"),
         Species = rename_unclassified(Species, "Species")) %>%
  add_column(ASV.ID = paste("ASV_F", 1:nrow(.), sep = ""),
             .before = 1) %>%
  select(all_of(c("ASV.ID", Taxonomic.level, row.names(soil.properties)))) %>%
  mutate(across(where(is.numeric),
                ~ decostand(., method = "total", MARGIN = 2)))


### 1.4.2 Abundance filtering----
fungi.table <- Relative.abu.Fungi %>%
  filter(rowMeans(select(., where(is.numeric))) > 0.0001,
         rowSums(across(where(is.numeric), ~ .x > 0)) > 3) %>%
  data.frame()


# 2. Fig_S2----

## 2.1) Bateria----

### 2.1.1 Phylum----


# Taxonomic ranking
bacteria.phylum.list <- Relative.abu.Bac %>%
  data.frame() %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Phylum_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Phylum_Abundance) %>%
  dplyr::slice(which(Phylum_Abundance > 0.02)) %>%
  arrange(desc(Phylum_Abundance)) %>%
  pull(Phylum)


# Abundance table
bacteria.Phylum.table <- Relative.abu.Bac %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% bacteria.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(bacteria.phylum.list, "Others1"))) %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Phylum_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Phylum_Abundance) %>%
  arrange(Phylum) %>%
  add_column(Level = "Phylum", .before = 1) %>%
  select(Level, Taxa = Phylum, Phylum = Phylum,
         Abundance = Phylum_Abundance)



### 2.1.2 Class----

# Taxonomic ranking
bacteria.Class.list <- Relative.abu.Bac %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% bacteria.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(bacteria.phylum.list, "Others1"))) %>%
  group_by(Phylum, Class) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Class_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Class_Abundance) %>%
  arrange(Phylum, desc(Class_Abundance)) %>% 
  data.frame() %>%
  subset(Phylum != "Others1") %>%
  pull(Class)

bacteria.Class.list


# Abundance table
bacteria.Class.table <- Relative.abu.Bac %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% bacteria.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(bacteria.phylum.list, "Others1")),
         Class = ifelse(Class %in% bacteria.Class.list, Class, "Others2")) %>%
  group_by(Phylum, Class) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Class_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Class_Abundance) %>%
  arrange(Phylum, desc(Class_Abundance)) %>%
  add_column(Level = "Class", .before = 1) %>%
  select(Level, Taxa = Class, Phylum = Phylum,
         Abundance = Class_Abundance)



### 2.1.3 Order----

# Taxonomic ranking
bacteria.Order.list <- Relative.abu.Bac %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% bacteria.phylum.list, Phylum, "Others1"),
         Class = ifelse(Class %in% bacteria.Class.list, Class, "Others2"),
         Phylum = factor(Phylum, levels = c(bacteria.phylum.list, "Others1")),
         Class = factor(Class, levels = c(bacteria.Class.list, "Others2")),
         Order = ifelse(Phylum %in% bacteria.phylum.list, Order, "Others3")) %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Order_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Order, Order_Abundance) %>%
  arrange(Phylum, Class, desc(Order_Abundance)) %>%
  data.frame() %>%
  pull(Order)

bacteria.Order.list

# Abundance table
bacteria.Order.table <- Relative.abu.Bac %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% bacteria.phylum.list, Phylum, "Others1"),
         Class = ifelse(Class %in% bacteria.Class.list, Class, "Others2"),
         Phylum = factor(Phylum, levels = c(bacteria.phylum.list, "Others1")),
         Class = factor(Class, levels = c(bacteria.Class.list, "Others2")),
         Order = ifelse(Phylum %in% bacteria.phylum.list, Order, "Others3")) %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Order_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Order, Order_Abundance) %>%
  arrange(Phylum, Class, desc(Order_Abundance)) %>%
  add_column(Level = "Order", .before = 1) %>%
  select(Level, Taxa = Order, Phylum = Phylum,
         Abundance = Order_Abundance)




plot.bacteria.data  <- rbind(bacteria.Phylum.table,
                   bacteria.Class.table,
                   bacteria.Order.table) %>%
  data.frame() %>%
  mutate(fill = case_when(Phylum == "Proteobacteria" ~ "#FF505C",
                          Phylum == "Acidobacteriota" ~ "#00FFFF",
                          Phylum == "Planctomycetota" ~ "#38AFCE",
                          Phylum == "Actinobacteriota" ~ "#9B7FFA",
                          Phylum == "Verrucomicrobiota" ~ "#FBE13F",
                          Phylum == "Desulfobacterota" ~ "#FB9C6E",
                          Phylum == "Chloroflexi" ~ "#AB3E09",
                          Phylum == "Firmicutes" ~ "#5C1A17",
                          .default = "#2E210C"),
         Level = factor(Level, levels = c("Phylum", "Class", "Order")))



### 2.1.4 Stacked bar chart----
plot.bacteria.data  %>%
  ggplot(aes(x = Level, y = Abundance, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill", show.legend = F,
           color = "grey70", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.025))) +
  scale_fill_manual(values = plot.bacteria.data $fill) +
  labs(x = "Taxonomic Level",
       y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


## 2.2) Fungi----

### 2.2.1 Phylum----


# Taxonomic ranking
fungi.phylum.list <- Relative.abu.Fungi %>%
  data.frame() %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Phylum_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Phylum_Abundance) %>%
  dplyr::slice(which(Phylum_Abundance > 0.02)) %>%
  arrange(desc(Phylum_Abundance)) %>%
  pull(Phylum)


# Abundance table
fungi.Phylum.table <- Relative.abu.Fungi %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% fungi.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(fungi.phylum.list, "Others1"))) %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Phylum_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Phylum_Abundance) %>%
  arrange(Phylum) %>%
  add_column(Level = "Phylum", .before = 1) %>%
  select(Level, Taxa = Phylum, Phylum = Phylum,
         Abundance = Phylum_Abundance)



### 2.2.2 Class----

# Taxonomic ranking
fungi.Class.list <- Relative.abu.Fungi %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% fungi.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(fungi.phylum.list, "Others1"))) %>%
  group_by(Phylum, Class) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Class_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Class_Abundance) %>%
  arrange(Phylum, desc(Class_Abundance)) %>% 
  data.frame() %>%
  subset(Phylum != "Others1") %>%
  pull(Class)

fungi.Class.list


# Abundance table
fungi.Class.table <- Relative.abu.Fungi %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% fungi.phylum.list, Phylum, "Others1"),
         Phylum = factor(Phylum, levels = c(fungi.phylum.list, "Others1")),
         Class = ifelse(Class %in% fungi.Class.list, Class, "Others2")) %>%
  group_by(Phylum, Class) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Class_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Class_Abundance) %>%
  arrange(Phylum, desc(Class_Abundance)) %>%
  add_column(Level = "Class", .before = 1) %>%
  select(Level, Taxa = Class, Phylum = Phylum,
         Abundance = Class_Abundance)



### 2.2.3 Order----

# Taxonomic ranking
fungi.Order.list <- Relative.abu.Fungi %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% fungi.phylum.list, Phylum, "Others1"),
         Class = ifelse(Class %in% fungi.Class.list, Class, "Others2"),
         Phylum = factor(Phylum, levels = c(fungi.phylum.list, "Others1")),
         Class = factor(Class, levels = c(fungi.Class.list, "Others2")),
         Order = ifelse(Phylum %in% fungi.phylum.list, Order, "Others3")) %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Order_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Order, Order_Abundance) %>%
  arrange(Phylum, Class, desc(Order_Abundance)) %>%
  data.frame() %>%
  pull(Order)

fungi.Order.list

# Abundance table
fungi.Order.table <- Relative.abu.Fungi %>%
  data.frame() %>%
  mutate(Phylum = ifelse(Phylum %in% fungi.phylum.list, Phylum, "Others1"),
         Class = ifelse(Class %in% fungi.Class.list, Class, "Others2"),
         Phylum = factor(Phylum, levels = c(fungi.phylum.list, "Others1")),
         Class = factor(Class, levels = c(fungi.Class.list, "Others2")),
         Order = ifelse(Phylum %in% fungi.phylum.list, Order, "Others3")) %>%
  group_by(Phylum, Class, Order) %>%
  summarise(across(where(is.numeric), sum, .names = "{col}"), .groups = "drop") %>%
  mutate(Order_Abundance = rowMeans(select(., where(is.numeric)))) %>%
  select(Phylum, Class, Order, Order_Abundance) %>%
  arrange(Phylum, Class, desc(Order_Abundance)) %>%
  add_column(Level = "Order", .before = 1) %>%
  select(Level, Taxa = Order, Phylum = Phylum,
         Abundance = Order_Abundance)




plot.fungi.data  <- rbind(fungi.Phylum.table,
                          fungi.Class.table,
                          fungi.Order.table) %>%
  data.frame() %>%
  mutate(fill = case_when(Phylum == "Basidiomycota" ~ "#FF505C",
                          Phylum == "Ascomycota" ~ "#00FFFF",
                          Phylum == "Mortierellomycota" ~ "#38AFCE",
                          Phylum == "Rozellomycota" ~ "#9B7FFA",
                          Phylum == "Mucoromycota" ~ "#FBE13F",
                          Phylum == "Glomeromycota" ~ "#FB9C6E",
                          Phylum == "Chytridiomycota" ~ "#AB3E09",
                          .default = "#2E210C"),
         Level = factor(Level, levels = c("Phylum", "Class", "Order")))



### 2.2.4 Stacked bar chart----
plot.fungi.data  %>%
  ggplot(aes(x = Level, y = Abundance, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill", show.legend = F,
           color = "grey70", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.025))) +
  scale_fill_manual(values = plot.fungi.data $fill) +
  labs(x = "Taxonomic Level",
       y = "Relative Abundance (%)",
       fill = "Phylum") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


# 3. Table_S2----

## 3.1) Sample grouping----
sample.ridge.alive <- soil.properties %>%
  subset(Topography == "Ridge" & Status == "Alive") %>%
  arrange(Pairs) %>%
  row.names()


sample.ridge.dead <- soil.properties %>%
  subset(Topography == "Ridge" & Status == "Dead") %>%
  arrange(Pairs) %>%
  row.names() 


sample.valley.alive <- soil.properties %>%
  subset(Topography == "Valley" & Status == "Alive") %>%
  arrange(Pairs) %>%
  row.names()

sample.valley.dead <- soil.properties %>%
  subset(Topography == "Valley" & Status == "Dead") %>%
  arrange(Pairs) %>%
  row.names()

Group <- soil.properties %>%
  select(Topography, Status)


## 3.2) Bacteria----
bacteria.abundance.matrix <- bacteria.table %>%
  data.frame() %>%
  select(where(is.numeric)) %>%
  t() %>% data.frame()


row.names(bacteria.abundance.matrix) == row.names(soil.properties)


#### 3.2.1 PERMANOVA----

set.seed(20250508)
adonis2(bacteria.abundance.matrix ~ Topography * Status,
        data = Group, permutations = 9999,
        method = "bray", by = "terms") %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  rename(`Sum of Squares` = SumOfSqs,
         `R²`              = R2,
         `F value`         = F,
         `P value`         = `Pr..F.`)


## 3.3) Fungi----
fungi.abundance.matrix <- fungi.table %>%
  data.frame() %>%
  select(where(is.numeric)) %>%
  t() %>% data.frame()


row.names(fungi.abundance.matrix) == row.names(soil.properties)


#### 3.3.1 PERMANOVA----

set.seed(20250508)
adonis2(fungi.abundance.matrix ~ Topography * Status,
        data = Group, permutations = 9999,
        method = "bray", by = "terms") %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Term") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
  rename(`Sum of Squares` = SumOfSqs,
         `R²`              = R2,
         `F value`         = F,
         `P value`         = `Pr..F.`)



# 4 Fig_S3----

## 4.1) Bacteria----

### 4.1.1 Wilcox Difference Analysis----
bacteria.wilcox.results <- apply(bacteria.table, 1, function(row) {
  
  ridge.alive <- as.numeric(row[sample.ridge.alive])
  ridge.dead  <- as.numeric(row[sample.ridge.dead])
  valley.alive <- as.numeric(row[sample.valley.alive])
  valley.dead  <- as.numeric(row[sample.valley.dead])
  
  
  
  test.ridge <- wilcox.test(ridge.alive, ridge.dead, paired = TRUE)
  test.valley <- wilcox.test(valley.alive, valley.dead, paired = TRUE)
  
  log2FC.ridge <- log2(mean(ridge.dead) / mean(ridge.alive))
  log2FC.valley <- log2(mean(valley.dead) / mean(valley.alive))
  
  data.frame(
    statistic.ridge = test.ridge$statistic,
    statistic.valley = test.valley$statistic,
    log2FC.ridge = log2FC.ridge, 
    log2FC.valley = log2FC.valley,
    p.value.ridge = test.ridge$p.value,
    p.value.valley = test.valley$p.value
  )
})



bacteria.wilcox.result_df <- do.call(rbind, bacteria.wilcox.results) %>%
  data.frame() %>%
  mutate(taxa = bacteria.table$Species,
         FDR.ridge = p.adjust(p.value.ridge, method = "BH"),
         FDR.valley = p.adjust(p.value.valley, method = "BH")) %>%
  select(taxa, everything())



### 4.1.2 Cross-topographic relationships----

bacteria.wilcox.result_df %>%
  sma(log2FC.valley ~ log2FC.ridge, data = ., method="SMA") %>%
  summary()



bacteria.wilcox.result_df %>%
  ggplot(aes(x = log2FC.ridge, y = log2FC.valley)) +
  geom_point(size = 7, shape = 21, fill = "#FEA6BF") +
  stat_ma_line(method = "SMA", se = TRUE, level = 0.95, linewidth = 1,
               fill = "#7E243D", colour = "#FC215F", alpha = 1,
               fm.values = T, nperm = 500) +
  labs(x = expression("Log"[2]*" fold change of abundance (Ridge)"),
       y = expression("Log"[2]*" fold change of abundance (Valley)")) +
  scale_y_continuous(limits = c(-2.5, 5), breaks = seq(-3, 5, 1),
                     expand = expansion(mult = c(0.05, 0))) +
  scale_x_continuous(limits = c(-5, 3), breaks = seq(-4, 3, 1),
                     expand = expansion(mult = c(0.05, 0))) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black"),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))




## 4.2) Fungi----

### 4.2.1 Wilcox Difference Analysis----
fungi.wilcox.results <- apply(fungi.table, 1, function(row) {
  
  ridge.alive <- as.numeric(row[sample.ridge.alive])
  ridge.dead  <- as.numeric(row[sample.ridge.dead])
  valley.alive <- as.numeric(row[sample.valley.alive])
  valley.dead  <- as.numeric(row[sample.valley.dead])
  
  
  test.ridge <- wilcox.test(ridge.alive, ridge.dead, paired = TRUE)
  test.valley <- wilcox.test(valley.alive, valley.dead, paired = TRUE)
  
  eps_r <- ifelse(min(c(ridge.alive, ridge.dead)) == 0, 
                  1e-6, min(c(ridge.alive, ridge.dead)) * 0.5)
  eps_v <- ifelse(min(c(valley.alive, valley.dead)) == 0, 
                  1e-6, min(c(valley.alive, valley.dead)) * 0.5)
  
  
  log2FC.ridge <- log2((mean(ridge.dead) + eps_r) / (mean(ridge.alive) + eps_r))
  log2FC.valley <- log2((mean(valley.dead) + eps_v) / (mean(valley.alive) + eps_v))
  
  data.frame(
    statistic.ridge = test.ridge$statistic,
    statistic.valley = test.valley$statistic,
    log2FC.ridge = log2FC.ridge, 
    log2FC.valley = log2FC.valley,
    p.value.ridge = test.ridge$p.value,
    p.value.valley = test.valley$p.value
  )
})



fungi.wilcox.result_df <- do.call(rbind, fungi.wilcox.results) %>%
  data.frame() %>%
  mutate(Genus = fungi.table$Genus,
         Species = fungi.table$Species,
         FDR.ridge = p.adjust(p.value.ridge, method = "BH"),
         FDR.valley = p.adjust(p.value.valley, method = "BH")) %>%
  select(Genus, Species, everything())



### 4.1.2 Cross-topographic relationships----
fungi.wilcox.result_df %>% 
  sma(log2FC.valley ~ log2FC.ridge, data = ., method= "SMA") %>%
  summary()


fungi.wilcox.result_df %>%
  ggplot(aes(x = log2FC.ridge, y = log2FC.valley)) +
  geom_point(size = 7, shape = 21, fill = "#F3BAF2") +
  stat_ma_line(method = "SMA", se = TRUE, level = 0.95, linewidth = 1,
               fill = "#792D78", colour = "#E053DF",alpha = 1,
               fm.values = T, nperm = 500) +
  labs(x = expression("Log"[2]*" fold change of abundance (Ridge)"),
       y = expression("Log"[2]*" fold change of abundance (Valley)")) +
  scale_y_continuous(limits = c(-15, 15), breaks = seq(-15, 15, 5),
                     expand = expansion(mult = c(0.05, 0))) +
  scale_x_continuous(limits = c(-15, 10), breaks = seq(-15, 10, 5),
                     expand = expansion(mult = c(0.05, 0))) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.y.left = element_line(color = "black"),
        axis.line.x.bottom = element_line(color = "black"),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))




# 5. Save results----

save(Relative.abu.Bac, Relative.abu.Fungi,
     bacteria.table, fungi.table, 
     bacteria.wilcox.result_df, fungi.wilcox.result_df,
     file = file.path(path_output,
                      "microbial_composition_results.Rdata"))

