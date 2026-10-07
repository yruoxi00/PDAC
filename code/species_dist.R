library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/PDAC.csv")

### stool 

stool_data <- pdac %>%
  filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

# species data 

species_matrix <- species %>%
  column_to_rownames(var = "ID")
species <- as.data.frame(t(species_matrix))
species <- species[stool_data$ID, ]
species <- species / 100

abundance_threshold <- 2e-4   # 0.0002
prevalence_threshold <- 0.10  # 10%

taxa <- lapply(split(stool_data$ID, stool_data$Data), function(sample_ids) {
  species_sub <- species[sample_ids, , drop = FALSE]  
  apply(species_sub, 2, function(x) {
    mean(x > abundance_threshold) >= prevalence_threshold
  })
})

taxa <- Reduce("|", taxa)
cat("Original number of taxa:", ncol(species), "\n")
species <- species[, taxa]
cat("Filtered number of taxa:", ncol(species), "\n")

plain <- function(x, ...) format(x, ..., scientific = FALSE, drop0trailing = TRUE)

studycol <- c(
  "Vilà-Quintana_2024" = "#7DABD1",
  "Kartal_2022"        = "#FFAC59",
  "Meng_2025"          = "#F6D963",
  "Matsukawa_2021"     = "#FAADD5",
  "Nagata_2022"        = "#BC8CC3"
)

map_study_names <- function(df, col = "study") {
  df %>%
    mutate(
      !!col := case_when(
        .data[[col]] == "_relabDat2" ~ "Vilà-Quintana_2024",
        .data[[col]] == "_relabDat5" ~ "Kartal_2022",
        .data[[col]] == "_relabDat6" ~ "Meng_2025",
        .data[[col]] == "_relabDat8" ~ "Matsukawa_2021",
        .data[[col]] == "_relabDat9" ~ "Nagata_2022",
        TRUE ~ .data[[col]]
      )
    ) %>%
    mutate(!!col := factor(.data[[col]],
                           levels = c("Vilà-Quintana_2024","Kartal_2022","Meng_2025","Matsukawa_2021","Nagata_2022")))
}

make_species_dist_plot <- function(stool_data_sub, species_all, title_text, top_n = 25) {
  species_sub <- species_all[stool_data_sub$ID, , drop = FALSE]
  
  species_presence <- species_sub > 0
  spe_all <- data.frame(species = colnames(species_sub))
  spe_all$present <- apply(species_presence, 2, function(x) {
    length(unique(stool_data_sub$Data[x]))
  })
  
  species_long <- species_sub %>%
    mutate(Data = stool_data_sub$Data) %>%
    group_by(Data) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(-Data, names_to = "species", values_to = "relab") %>%
    pivot_wider(names_from = Data, values_from = relab, names_prefix = "_relab")
  
  spe_all <- merge(spe_all, species_long, by = "species")
  
  remove_pattern <- "(SGB|GGB|unknown|uncultured|unclassified|metagenome)"
  spe_all <- spe_all %>%
    filter(!grepl(remove_pattern, species, ignore.case = TRUE))
  
  n_studies <- length(grep("^_relab", names(spe_all)))
  
  universal <- spe_all %>%
    filter(present == n_studies) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Universal")
  
  overlap <- spe_all %>%
    filter(present != n_studies & present != 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Overlapping")
  
  solo <- spe_all %>%
    filter(present == 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Singular")
  
  spe_plot <- bind_rows(universal, overlap, solo)
  
  spe_plot2 <- spe_plot %>%
    select(-mean_abund) %>%
    pivot_longer(cols = contains("_relab"), names_to = "study", values_to = "relab") %>%
    filter(relab > 0) %>%
    select(species, group, study, relab)
  
  spe_plot2 <- spe_plot2 %>%
    mutate(species = gsub("^s__", "", species),
           species = gsub("_", " ", species))
  
  spe_order <- spe_plot2 %>%
    group_by(group, species) %>%
    summarise(mean = mean(relab), .groups = "drop_last") %>%
    arrange(group, desc(mean))
  spe_plot2$species <- factor(spe_plot2$species, levels = spe_order$species)
  spe_plot2$group <- factor(spe_plot2$group, levels = c("Universal","Overlapping","Singular"))
  
  spe_plot2 <- map_study_names(spe_plot2, col = "study")
  
  g <- ggplot(spe_plot2, aes(x = species, y = relab)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
    geom_point(aes(fill = study), shape = 21, size = 2.6, alpha = 0.85,
               position = position_jitter(width = 0.08, height = 0)) +
    scale_fill_manual(values = studycol, drop = FALSE) +
    ylab("Relative abundance (%)") + xlab("") +
    scale_y_log10(breaks = c(10, 1, 0.1, 0.01), labels = plain) +
    facet_wrap(~ group, scales = "free_x") +
    theme(strip.text.x = element_text(hjust = 0.5, color = "black", size = 14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 11, face = "italic", color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(legend.position = c(.9, .85),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.background = element_blank(),
          legend.key.height = unit(.3, 'cm'),
          legend.key.width  = unit(.3, 'cm'),
          legend.key = element_rect(fill = "white")) +
    guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
    ggtitle(title_text)
  
  g
}

# ---- PC subset: Pancreatic Cancer ----
stool_PC <- stool_data %>%
  filter(Group_new == "Pancreatic Cancer")

spe_dist_PC1 <- make_species_dist_plot(
  stool_data_sub = stool_PC,
  species_all    = species,
  title_text     = "Stool species distribution: PC (Pancreatic Cancer)",
  top_n          = 25
)

# ---- Control subset: Control ----
stool_con <- stool_data %>%
  filter(Group_new == "Control")

spe_dist_con1 <- make_species_dist_plot(
  stool_data_sub = stool_con,
  species_all    = species,
  title_text     = "Stool species distributiob: Control",
  top_n          = 25
)

p1 <- spe_dist_PC1
p2 <- spe_dist_con1

### saliva 

saliva_data <- pdac %>%
  filter(Type == "Saliva")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

# species data 

species_matrix <- species %>%
  column_to_rownames(var = "ID")
species <- as.data.frame(t(species_matrix))
species <- species[saliva_data$ID, ]
species <- species / 100

abundance_threshold <- 2e-4   # 0.0002
prevalence_threshold <- 0.10  # 10%

taxa <- lapply(split(saliva_data$ID, saliva_data$Data), function(sample_ids) {
  species_sub <- species[sample_ids, , drop = FALSE]  
  apply(species_sub, 2, function(x) {
    mean(x > abundance_threshold) >= prevalence_threshold
  })
})

taxa <- Reduce("|", taxa)
cat("Original number of taxa:", ncol(species), "\n")
species <- species[, taxa]
cat("Filtered number of taxa:", ncol(species), "\n")

plain <- function(x, ...) format(x, ..., scientific = FALSE, drop0trailing = TRUE)

studycol <- c(
  "Vilà-Quintana_2024" = "#7DABD1",
  "Kartal_2022"        = "#FFAC59",
  "Matsukawa_2021"     = "#FAADD5",
  "Nagata_2022"        = "#BC8CC3"
)

map_study_names <- function(df, col = "study") {
  df %>%
    mutate(
      !!col := case_when(
        .data[[col]] == "_relabDat2" ~ "Vilà-Quintana_2024",
        .data[[col]] == "_relabDat5" ~ "Kartal_2022",
        .data[[col]] == "_relabDat8" ~ "Matsukawa_2021",
        .data[[col]] == "_relabDat9" ~ "Nagata_2022",
        TRUE ~ .data[[col]]
      )
    ) %>%
    mutate(!!col := factor(.data[[col]],
                           levels = c("Vilà-Quintana_2024","Kartal_2022","Matsukawa_2021","Nagata_2022")))
}

make_species_dist_plot <- function(saliva_data_sub, species_all, title_text, top_n = 25) {
  species_sub <- species_all[saliva_data_sub$ID, , drop = FALSE]
  
  species_presence <- species_sub > 0
  spe_all <- data.frame(species = colnames(species_sub))
  spe_all$present <- apply(species_presence, 2, function(x) {
    length(unique(saliva_data_sub$Data[x]))
  })
  
  species_long <- species_sub %>%
    mutate(Data = saliva_data_sub$Data) %>%
    group_by(Data) %>%
    summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(-Data, names_to = "species", values_to = "relab") %>%
    pivot_wider(names_from = Data, values_from = relab, names_prefix = "_relab")
  
  spe_all <- merge(spe_all, species_long, by = "species")
  
  n_studies <- length(grep("^_relab", names(spe_all)))
  
  universal <- spe_all %>%
    filter(present == n_studies) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Universal")
  
  overlap <- spe_all %>%
    filter(present != n_studies & present != 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Overlapping")
  
  solo <- spe_all %>%
    filter(present == 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>% slice_head(n = top_n) %>%
    mutate(group = "Singular")
  
  spe_plot <- bind_rows(universal, overlap, solo)
  
  spe_plot2 <- spe_plot %>%
    select(-mean_abund) %>%
    pivot_longer(cols = contains("_relab"), names_to = "study", values_to = "relab") %>%
    filter(relab > 0) %>%
    select(species, group, study, relab)
  
  spe_plot2 <- spe_plot2 %>%
    mutate(
      species_id = species,
      species_label = gsub("_", " ", gsub("^s__", "", species))
    )
  
  spe_order <- spe_plot2 %>%
    group_by(group, species_id) %>%
    summarise(mean = mean(relab), .groups = "drop_last") %>%
    arrange(group, desc(mean))
  
  spe_plot2$species_id <- factor(spe_plot2$species_id, levels = unique(spe_order$species_id))
  spe_plot2$group <- factor(spe_plot2$group, levels = c("Universal","Overlapping","Singular"))
  
  spe_plot2 <- map_study_names(spe_plot2, col = "study")
  
  label_map <- setNames(spe_plot2$species_label, spe_plot2$species_id)
  
  g <- ggplot(spe_plot2, aes(x = species_id, y = relab)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
    geom_point(aes(fill = study), shape = 21, size = 2.6, alpha = 0.85,
               position = position_jitter(width = 0.08, height = 0)) +
    scale_fill_manual(values = studycol, drop = FALSE) +
    scale_x_discrete(labels = label_map) +
    ylab("Relative abundance (%)") + xlab("") +
    scale_y_log10(breaks = c(10, 1, 0.1, 0.01), labels = plain) +
    facet_wrap(~ group, scales = "free_x") +
    theme(strip.text.x = element_text(hjust = 0.5, color = "black", size = 14),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 11, face = "italic", color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(legend.position = c(.9, .85),
          legend.title = element_blank(),
          legend.text = element_text(size = 10, color = "black"),
          legend.background = element_blank(),
          legend.key.height = unit(.3, 'cm'),
          legend.key.width  = unit(.3, 'cm'),
          legend.key = element_rect(fill = "white")) +
    guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
    ggtitle(title_text)
  
  g
}

# ---- PC subset ----
saliva_PC <- saliva_data %>%
  filter(Group_new == "Pancreatic Cancer")

spe_dist_PC2 <- make_species_dist_plot(
  saliva_data_sub = saliva_PC,
  species_all    = species,
  title_text     = "Saliva species distribution: PC (Pancreatic Cancer)",
  top_n          = 25
)

# ---- Control subset: Control ----
saliva_con <- saliva_data %>%
  filter(Group_new == "Control")

spe_dist_con2 <- make_species_dist_plot(
  saliva_data_sub = saliva_con,
  species_all    = species,
  title_text     = "Saliva species distributiob: Control",
  top_n          = 25
)

p3 <- spe_dist_PC2
p4 <- spe_dist_con2

p5 <- plot_grid(
  p1, p2, p3, p4, 
  labels = c("A", "B", "C", "D"),     
  label_size = 20,          
  ncol = 2,
  align = "hv"
)
print(p5)

p12 <- plot_grid(
  p1, p2,
  ncol = 2,
  align = "hv"
)
print(p12)

ggsave(p12, filename="fig/stool_distribution.pdf", width = 25, height = 8)
ggsave(p5, filename="fig/distribution.pdf", width = 25, height = 20)
