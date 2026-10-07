library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(cowplot)
library(maaslin3)
library(randomForest)
library(pROC)
library(tibble)
library(stringr)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/AP.csv")

# ===================== All Cohort ===================== #

### stool 

stool_data <- pdac %>%
  filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

### species 

species_matrix <- species %>%
  column_to_rownames(var = "ID")
species_raw <- as.data.frame(t(species_matrix))
species_raw <- species_raw[stool_data$ID, ]
species_raw <- species_raw / 100

abundance_threshold  <- 2e-4   # 0.0002
prevalence_threshold <- 0.10   # 10%

filter_group_species <- function(group_name, stool_df, species_mat,
                                 abundance_threshold = 2e-4,
                                 prevalence_threshold = 0.10) {
  group_df <- stool_df %>%
    filter(Group_new == group_name)
  
  species_group <- species_mat[group_df$ID, , drop = FALSE]
  
  taxa_list <- lapply(split(group_df$ID, group_df$Data), function(sample_ids) {
    species_sub <- species_group[sample_ids, , drop = FALSE]
    apply(species_sub, 2, function(x) {
      mean(x > abundance_threshold) >= prevalence_threshold
    })
  })
  
  taxa_keep <- Reduce("|", taxa_list)
  
  species_group_filtered <- species_group[, taxa_keep, drop = FALSE]
  species_group_filtered
}

species_control <- filter_group_species("Control", stool_data, species_raw,
                                        abundance_threshold, prevalence_threshold)

species_pdac <- filter_group_species("AP", stool_data, species_raw,
                                     abundance_threshold, prevalence_threshold)

species_union <- union(colnames(species_control), colnames(species_pdac))

species_stool <- species_raw[, species_union, drop = FALSE]
species_stool <- as.data.frame(species_stool)

meta_stool <- stool_data %>%
  filter(ID %in% rownames(species_stool)) %>%
  arrange(match(ID, rownames(species_stool))) %>%
  mutate(
    Group_new  = factor(Group_new, levels = c("Control", "AP")),
    Sex        = factor(Sex),
    smoker     = factor(smoker),
    alcohol    = factor(ifelse(alcohol == 1, "Yes", "No"), levels = c("No", "Yes")),
    diabetes   = factor(ifelse(diabetes == 1, "Yes", "No"), levels = c("No", "Yes")),
    antibiotic = factor(ifelse(antibiotic == 1, "Yes", "No"), levels = c("No", "Yes"))
  )

meta_stool <- as.data.frame(meta_stool)
rownames(meta_stool) <- meta_stool$ID 

all(rownames(species_stool) == rownames(meta_stool))

fit_stool <- maaslin3(
  input_data = species_stool,
  input_metadata = meta_stool,
  output = "result/Maaslin3_stool_ap",
  fixed_effects = c("Group_new"),
  reference = c("Group_new,Control"), 
  random_effects = c("Data"), 
  normalization = "NONE",
  transform = "LOG",
  correction = "BH",
  min_prevalence = 0.10
)

# ===================== Single Cohort ===================== #

dat5 <- pdac %>%
  filter(Data == "Dat5")

### stool 

stool_data <- dat5 %>%
  filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

### species 

species_matrix <- species %>%
  column_to_rownames(var = "ID")
species_raw <- as.data.frame(t(species_matrix))
species_raw <- species_raw[stool_data$ID, ]
species_raw <- species_raw / 100

abundance_threshold  <- 2e-4   # 0.0002
prevalence_threshold <- 0.10   # 10%

filter_group_species <- function(group_name, stool_df, species_mat,
                                 abundance_threshold = 2e-4,
                                 prevalence_threshold = 0.10) {
  group_df <- stool_df %>%
    filter(Group_new == group_name)
  
  species_group <- species_mat[group_df$ID, , drop = FALSE]
  
  taxa_list <- lapply(split(group_df$ID, group_df$Data), function(sample_ids) {
    species_sub <- species_group[sample_ids, , drop = FALSE]
    apply(species_sub, 2, function(x) {
      mean(x > abundance_threshold) >= prevalence_threshold
    })
  })
  
  taxa_keep <- Reduce("|", taxa_list)
  
  species_group_filtered <- species_group[, taxa_keep, drop = FALSE]
  species_group_filtered
}

species_control <- filter_group_species("Control", stool_data, species_raw,
                                        abundance_threshold, prevalence_threshold)

species_pdac <- filter_group_species("AP", stool_data, species_raw,
                                     abundance_threshold, prevalence_threshold)

species_union <- union(colnames(species_control), colnames(species_pdac))

species_stool <- species_raw[, species_union, drop = FALSE]
species_stool <- as.data.frame(species_stool)

meta_stool <- stool_data %>%
  filter(ID %in% rownames(species_stool)) %>%
  arrange(match(ID, rownames(species_stool))) %>%
  mutate(
    Group_new  = factor(Group_new, levels = c("Control", "AP")),
    Sex        = factor(Sex),
    smoker     = factor(smoker),
    alcohol    = factor(ifelse(alcohol == 1, "Yes", "No"), levels = c("No", "Yes")),
    diabetes   = factor(ifelse(diabetes == 1, "Yes", "No"), levels = c("No", "Yes")),
    antibiotic = factor(ifelse(antibiotic == 1, "Yes", "No"), levels = c("No", "Yes"))
  )

meta_stool <- as.data.frame(meta_stool)
rownames(meta_stool) <- meta_stool$ID 

all(rownames(species_stool) == rownames(meta_stool))

fit_stool <- maaslin3(
  input_data = species_stool,
  input_metadata = meta_stool,
  output = "result/Maaslin3_stool_dat4_ap",
  fixed_effects = c("Group_new"),
  reference = c("Group_new,Control"), 
  normalization = "NONE",
  transform = "LOG",
  correction = "BH",
  min_prevalence = 0.10
)

fit_stool <- maaslin3(
  input_data = species_stool,
  input_metadata = meta_stool,
  output = "result/Maaslin3_stool_dat5_ap",
  fixed_effects = c("Group_new", "Age", "Sex", "smoker", "alcohol", "diabetes", "antibiotic"),
  reference = c("Group_new,Control", "Sex,male","smoker,non-smoker","alcohol,No","diabetes,No","antibiotic,No"), 
  normalization = "NONE",
  transform = "LOG",
  correction = "BH",
  min_prevalence = 0.10
)
