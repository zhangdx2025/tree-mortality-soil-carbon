# Set the working directory
setwd("G:/R_script")


path_properties <- file.path(getwd(), "Data/Soil_properties")
path_microbiome <- file.path(getwd(), "Data/microbiome_profiles")



# load packages
library(tidyverse)
library(readxl)
library(metacoder)


# load funtion
cell_fun = function(j, i, x, y, width, height, fill) {
  grid.rect(x = x, y = y,
            width = width, height = height,
            gp = gpar(col = "grey80", lwd = 0.5, fill = NA))  
}


# 1. Loading Data----

## 1.1) microbial composition----
load(file.path(path_microbiome, "microbial_composition_results.Rdata"))


## 1.2) gene copy number in the genome----


# gene copy number in the genome of bacteria amplicon taxa

### 1.2.1 Bacterial genomes----
Bacterial.taxa.genome.gene.copy.number.matrix <- 
  read.csv(file.path(path_microbiome,
                     "bacterial.taxa.genome.gene.copy.number.matrix.txt"),
           header = TRUE)


### 1.2.1 Fungal genomes----
Fungal.taxa.genome.gene.copy.number.matrix <- 
  read.csv(file.path(path_microbiome,
                     "fungal.taxa.genome.gene.copy.number.matrix.txt"),
           header = TRUE)





## 1.3) soil.properties----
soil.properties <- read_excel(file.path(path_properties,
                                        "soil.properties.xlsx"),
                              sheet = "Data") %>%
  arrange(Topography, Pairs, Status) %>%
  column_to_rownames("Sample.ID")

sample.ridge <- soil.properties %>%
  subset(Topography == "Ridge") %>%
  row.names()

sample.valley <- soil.properties %>%
  subset(Topography == "Valley") %>%
  row.names()



# 2. Table_S3----

## 2.1) Bacterial genome profiles----

### 2.1.1 Genome number----
Bacterial.taxa.genome.gene.copy.number.matrix %>% nrow()


### 2.1.2 Genus number----

# Total genus number
bacteria.table$Genus %>% unique() %>% length()

# Number of genus with genomes
Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Genus %>%
  unique() %>% length()


### 2.1.3 Species number----

# Total species number
bacteria.table$Species %>% unique() %>% length()

# Number of species with genomes
Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Species %>%
  unique() %>% length()


### 2.1.4 Abundance of covered Genus----

# Mean
bacteria.table %>%
  subset(Genus %in% Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Genus) %>%
  select(where(is.numeric)) %>%
  colSums() %>% mean()

# SD
bacteria.table %>%
  subset(Genus %in% Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Genus) %>%
  select(where(is.numeric)) %>%
  colSums() %>% sd()


### 2.1.5 Abundance of covered Species----

# Mean
bacteria.table %>%
  subset(Species %in% Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Species) %>%
  select(where(is.numeric)) %>%
  colSums() %>% mean()

# SD
bacteria.table %>%
  subset(Species %in% Bacterial.taxa.genome.gene.copy.number.matrix$Silva.Species) %>%
  select(where(is.numeric)) %>%
  rowMeans() %>% sd()



## 2.2) Fungal genome profiles----

### 2.2.1 Genome number----
Fungal.taxa.genome.gene.copy.number.matrix %>% nrow()


### 2.2.2 Genus number----

# Total genus number
fungi.table$Genus %>% unique() %>% length()


# Number of genus with genomes
Fungal.taxa.genome.gene.copy.number.matrix$Silva.Genus %>%
  unique() %>% length()


### 2.2.3 Species number----

# Total species number
fungi.table$Species %>% unique() %>% length()

# Number of species with genomes
Fungal.taxa.genome.gene.copy.number.matrix$Silva.Species %>%
  unique() %>% length()


### 2.2.4 Abundance of covered Genus----

# Mean
fungi.table %>%
  subset(Genus %in% Fungal.taxa.genome.gene.copy.number.matrix$Silva.Genus) %>%
  select(where(is.numeric)) %>%
  colSums() %>% mean()

# SD
fungi.table %>%
  subset(Genus %in% Fungal.taxa.genome.gene.copy.number.matrix$Silva.Genus) %>%
  select(where(is.numeric)) %>%
  colSums() %>% sd()


### 2.2.5 Abundance of covered Species----

# Mean
fungi.table %>%
  subset(Species %in% Fungal.taxa.genome.gene.copy.number.matrix$Silva.Species) %>%
  select(where(is.numeric)) %>%
  colSums() %>% mean()

# SD
fungi.table %>%
  subset(Species %in% Fungal.taxa.genome.gene.copy.number.matrix$Silva.Species) %>%
  select(where(is.numeric)) %>%
  rowMeans() %>% sd()


## 2.3) Fig_3a----

bacteria.object <- Bacterial.taxa.genome.gene.copy.number.matrix %>%
  select(Phylum = Silva.Phylum, Class = Silva.Class, Order = Silva.Order) %>%
  add_column(Bacteria = "Bacteria", .before = 1) %>%
  parse_tax_data(class_cols = c("Bacteria", "Phylum",
                                "Class", "Order"),
                 class_key = "taxon_name",
                 named_by_rank = T)


### 2.3.1 Taxonomic distribution----

heat_tree(
  bacteria.object,
  node_label = taxon_names,
  node_size = n_obs,  
  node_color = n_obs,
  node_color_range = c("#FEE8EE", "#BF2652"),
  node_size_trans = "log10 area",
  edge_size_trans = "log10 area",
  make_node_legend = T,
  node_color_interval = c(1, 164),
  node_size_interval = c(1, 164),
  node_size_range = c(0.005, 0.06),
  node_label_size_range = c(0.02, 0.06),
  edge_size_range = c(0.003, 0.02),
  node_size_axis_label = "Number of species",
  node_color_axis_label = "Number of species"
) 


## 2.4) Fig_3b----

fungi.object <- Fungal.taxa.genome.gene.copy.number.matrix %>%
  select(Phylum = Silva.Phylum, Class = Silva.Class, Order = Silva.Order) %>%
  add_column(Fungi = "Fungi", .before = 1) %>%
  parse_tax_data(class_cols = c("Fungi", "Phylum",
                                "Class", "Order"),
                 class_key = "taxon_name",
                 named_by_rank = T)


### 2.4.1 Taxonomic distribution----

heat_tree(
  fungi.object,
  node_label = taxon_names, 
  node_size = n_obs,
  node_color = n_obs,
  node_color_range = c("#F7E9F7", "#8D348B"),
  node_size_trans = "log10 area",
  edge_size_trans = "log10 area",
  make_node_legend = T,
  node_color_interval = c(1, 223),
  node_size_interval = c(1, 223),
  node_size_range = c(0.005, 0.06),
  node_label_size_range = c(0.02, 0.06),
  edge_size_range = c(0.003, 0.02),
  node_size_axis_label = "Number of speceis",
  node_color_axis_label = "Number of speceis"
) 



# 3. Functional gene profiles----

target.gene <- c("Exo_amylase", "Endo_amylase", "Exo_cellulase", "Endo_cellulase",
                 "Exo_hemicellulose", "Endo_hemicellulose", "Exo_chitosanase",
                 "Endo_chitosanase", "Ligninase", "TEA_Oxygen", "TEA_Nitrate",
                 "TEA_Sulfite")


## 3.1) Bacteria----
mat.bacteria <- Bacterial.taxa.genome.gene.copy.number.matrix %>%
  group_by(Silva.Species) %>%
  summarise(across(any_of(target.gene), ~mean(.x)), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(.x >= 1, 1, 0))) %>%
  column_to_rownames("Silva.Species") %>%
  arrange(rowSums(.)) %>%
  as.matrix()


### 3.1.1 Fig_3c----


Heatmap(t(mat.bacteria), cluster_columns = F,
        cluster_rows = F,
        col = colorRamp2(c(0, 1), c("white", "#FC477A")),
        show_column_names = FALSE,
        show_row_names = TRUE,
        show_heatmap_legend = FALSE,
        cell_fun = cell_fun)



# 3.2) Fungi----
mat.fungi <- Fungal.taxa.genome.gene.copy.number.matrix %>%
  group_by(Ncbi.Species) %>%
  summarise(across(any_of(target.gene), ~mean(.x)), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(.x >= 1, 1, 0))) %>%
  column_to_rownames("Ncbi.Species") %>%
  arrange(rowSums(.)) %>%
  as.matrix()


### 3.2.1 Fig_3d----
Heatmap(t(mat.fungi), cluster_columns = F,
        cluster_rows = F,
        col = colorRamp2(c(0, 1), c("white", "#AD5EFC")),
        show_column_names = FALSE,
        show_row_names = TRUE,
        show_heatmap_legend = FALSE,
        cell_fun = cell_fun)


# 4. Fig_3e----

## 4.1) Bacteria----
mat.bacteria %>%
  rowSums() %>%
  data.frame(Gene.Count = .) %>%
  ggplot(aes(x = Gene.Count)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, fill = "#FC477A", color = "grey30", boundary = 0) +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "Number of functional genes per species",
       y = "Proportion of species") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


## 4.2) Fungi----
mat.fungi %>%
  rowSums() %>%
  data.frame(Gene.Count = .) %>%
  ggplot(aes(x = Gene.Count)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 1, fill = "#AD5EFC", color = "grey30", boundary = 0) +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.001, 0.05))) +
  labs(x = "Number of functional genes per species",
       y = "Proportion of species") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))


# 5. Kernel density plot----

Jaccard.distance.bacteria <- mat.bacteria %>%
  as.matrix() %>%
  vegdist(method = "jaccard", binary = TRUE, diag = FALSE) %>%
  as.vector()

Jaccard.distance.fungi <- mat.fungi %>%
  as.matrix() %>%
  vegdist(method = "jaccard", binary = TRUE, diag = FALSE) %>%
  as.vector()


## 5.1) Fig_3f----
data.frame(Group = c(rep("Bacteria", length(Jaccard.distance.bacteria)),
                                            rep("Fungi", length(Jaccard.distance.fungi))),
                                  Jaccard.distance = c(Jaccard.distance.bacteria, Jaccard.distance.fungi)) %>%
  subset(!is.na(Jaccard.distance)) %>%
  ggplot(aes(x = Jaccard.distance, fill = Group)) +
  geom_density(alpha = 1, show.legend = F) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0.001, 0.05))) +
  scale_fill_manual(values = c("Bacteria" = "#FC477A", "Fungi" = "#AD5EFC")) +
  labs(x = "Jaccard fistance",
       y = "Relative frequency") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        axis.ticks.length.x = unit(-0.25, "cm"),
        axis.ticks.length.y = unit(-0.25, "cm"))

