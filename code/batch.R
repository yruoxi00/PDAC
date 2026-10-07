library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(MMUPHin)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

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

## batch effect correction 

meta_data_stool <- stool_data %>%
  column_to_rownames("ID") %>%
  select(-Group) %>%
  mutate(
    Data = as.factor(Data),
    Group_new = factor(Group_new)
  )

# confounding check

batch_var <- "Data"
covariates <- c("Group_new")

for (cov in covariates) {
  tab <- table(meta_data_stool[[batch_var]], meta_data_stool[[cov]])
  unique_counts <- tapply(meta_data_stool[[cov]], meta_data_stool[[batch_var]], function(x) length(unique(x)))
  
  if (all(unique_counts <= 1)) {
    cat(cov, "is completely confounded with", batch_var, "\n")
  } else {
    cat(cov, "has variation within", batch_var, "\n")
  }
}

# cross-study correction 

species_adj <- adjust_batch(
  feature_abd = t(species),
  batch       = "Data",                         
  covariates  = c("Group_new"), 
  data        = meta_data_stool,
  control     = list(verbose = TRUE, diagnostics = TRUE)
)
species_adj_matrix <- species_adj$feature_abd_adj

species_adj_matrix <- sweep(species_adj_matrix, 2, colSums(species_adj_matrix), FUN = "/")
summary(colSums(species_adj_matrix))  # should all be 1

## validation (PERMANOVA and PCoA)

dist_before <- vegdist(species, method = "bray")
dist_after <- vegdist(t(species_adj_matrix), method = "bray")

adonis_before <- adonis2(dist_before ~ Data, data = meta_data_stool)
adonis_after <- adonis2(dist_after ~ Data, data = meta_data_stool)
print(adonis_before)
print(adonis_after)

# PCoA before 
pcoa_before <- cmdscale(dist_before, k = 2, eig = TRUE)
pcoa_before_df <- data.frame(
  Axis1 = pcoa_before$points[, 1],
  Axis2 = pcoa_before$points[, 2],
  Batch = meta_data_stool$Study
)
variance <- round(pcoa_before$eig / sum(pcoa_before$eig) * 100, 2)
xlab_txt <- paste0("PCoA 1 (", variance[1], "%)")
ylab_txt <- paste0("PCoA 2 (", variance[2], "%)")

p1 <- ggplot(pcoa_before_df, aes(x = Axis1, y = Axis2, fill = Batch)) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.3) +
  stat_ellipse(color = NA, fill = NA, size = 0) +
  labs(
    x = xlab_txt,
    y = ylab_txt,
    fill = "Batch"
  ) +
  scale_fill_manual(values = c(
    "Vilà-Quintana_2024" = "#7DABD1",
    "Kartal-1_2022"      = "#C7D3E5", 
    "Kartal-2_2022"      = "#F6D963",
    "Matsukawa_2021"     = "#FAADD5",
    "Nagata_2022"        = "#BC8CC3"
  )) +
  theme_classic() +
  theme(
    legend.position      = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    # larger legend text
    legend.text          = element_text(size = 14),
    legend.key.size      = unit(0.9, "lines"),
    plot.title           = element_text(hjust = 0.5),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    # larger axis text and labels
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = "Before~correction*\",\"~R^2==\"6.2\"*'%'",
    hjust = -0.05, vjust = 1.5,
    parse = TRUE,
    size = 7.5,        # larger annotation text
    fontface = "bold"
  )
print(p1)

# PCoA after 
pcoa_after <- cmdscale(dist_after, k = 2, eig = TRUE)
pcoa_after_df <- data.frame(
  Axis1 = pcoa_after$points[, 1],
  Axis2 = pcoa_after$points[, 2],
  Batch = meta_data_stool$Study
)
variance <- round(pcoa_after$eig / sum(pcoa_after$eig) * 100, 2)
xlab_txt <- paste0("PCoA 1 (", variance[1], "%)")
ylab_txt <- paste0("PCoA 2 (", variance[2], "%)")

p2 <- ggplot(pcoa_after_df, aes(x = Axis1, y = Axis2, fill = Batch)) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.3) +
  stat_ellipse(color = NA, fill = NA, size = 0) + 
  labs(
    x = xlab_txt,
    y = ylab_txt,
    fill = "Batch"
  ) +
  scale_fill_manual(values = c(
    "Vilà-Quintana_2024" = "#7DABD1",
    "Kartal-1_2022"      = "#C7D3E5", 
    "Kartal-2_2022"      = "#F6D963",
    "Matsukawa_2021"     = "#FAADD5",
    "Nagata_2022"        = "#BC8CC3"
  )) +
  theme_classic() +
  theme(
    legend.position      = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 14),
    legend.key.size      = unit(0.9, "lines"),
    plot.title           = element_text(hjust = 0.5),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = "After~correction*\",\"~R^2==\"2.7\"*'%'",
    hjust = -0.05, vjust = 1.5,
    parse = TRUE,
    size = 7.5,
    fontface = "bold"
  )
print(p2)

p1 <- p1 + stat_ellipse(level = 0.95, show.legend = F,linetype = 2)
p2 <- p2 + stat_ellipse(level = 0.95, show.legend = F,linetype = 2)
print(p1)
print(p2)

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

## batch effect correction 

meta_data_saliva <- saliva_data %>%
  tibble::column_to_rownames("ID") %>%
  select(-Group) %>%
  mutate(
    Data = as.factor(Data),
    Group_new = factor(Group_new)
  )

# confounding check

batch_var <- "Data"
covariates <- c("Group_new")

for (cov in covariates) {
  tab <- table(meta_data_saliva[[batch_var]], meta_data_saliva[[cov]])
  unique_counts <- tapply(meta_data_saliva[[cov]], meta_data_saliva[[batch_var]], function(x) length(unique(x)))
  
  if (all(unique_counts <= 1)) {
    cat(cov, "is completely confounded with", batch_var, "\n")
  } else {
    cat(cov, "has variation within", batch_var, "\n")
  }
}

# cross-study correction 

species_adj <- adjust_batch(
  feature_abd = t(species),
  batch       = "Data",                         
  covariates  = c("Group_new"), 
  data        = meta_data_saliva,
  control     = list(verbose = TRUE, diagnostics = TRUE)
)
species_adj_matrix <- species_adj$feature_abd_adj

species_adj_matrix <- sweep(species_adj_matrix, 2, colSums(species_adj_matrix), FUN = "/")
summary(colSums(species_adj_matrix))  # should all be 1

## validation (PERMANOVA and PCoA)

dist_before <- vegdist(species, method = "bray")
dist_after <- vegdist(t(species_adj_matrix), method = "bray")

adonis_before <- adonis2(dist_before ~ Data, data = meta_data_saliva)
adonis_after <- adonis2(dist_after ~ Data, data = meta_data_saliva)
print(adonis_before)
print(adonis_after)

# PCoA before 
pcoa_before <- cmdscale(dist_before, k = 2, eig = TRUE)
pcoa_before_df <- data.frame(
  Axis1 = pcoa_before$points[, 1],
  Axis2 = pcoa_before$points[, 2],
  Batch = meta_data_saliva$Study
)
variance <- round(pcoa_before$eig / sum(pcoa_before$eig) * 100, 2)
xlab_txt <- paste0("PCoA 1 (", variance[1], "%)")
ylab_txt <- paste0("PCoA 2 (", variance[2], "%)")

p3 <- ggplot(pcoa_before_df, aes(x = Axis1, y = Axis2, fill = Batch)) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.3) +
  stat_ellipse(color = NA, fill = NA, size = 0) +
  labs(
    x = xlab_txt,
    y = ylab_txt,
    fill = "Batch"
  ) +
  scale_fill_manual(values = c(
    "Vilà-Quintana_2024" = "#7DABD1",
    "Kartal-1_2022"      = "#C7D3E5", 
    "Matsukawa_2021"     = "#FAADD5",
    "Nagata_2022"        = "#BC8CC3"
  )) +
  theme_classic() +
  theme(
    legend.position      = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    # larger legend text
    legend.text          = element_text(size = 14),
    legend.key.size      = unit(0.9, "lines"),
    plot.title           = element_text(hjust = 0.5),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    # larger axis text and labels
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = "Before~correction*\",\"~R^2==\"11.0\"*'%'",
    hjust = -0.05, vjust = 1.5,
    parse = TRUE,
    size = 7.5,        # larger annotation text
    fontface = "bold"
  )
print(p3)

# PCoA after 
pcoa_after <- cmdscale(dist_after, k = 2, eig = TRUE)
pcoa_after_df <- data.frame(
  Axis1 = pcoa_after$points[, 1],
  Axis2 = pcoa_after$points[, 2],
  Batch = meta_data_saliva$Study
)
variance <- round(pcoa_after$eig / sum(pcoa_after$eig) * 100, 2)
xlab_txt <- paste0("PCoA 1 (", variance[1], "%)")
ylab_txt <- paste0("PCoA 2 (", variance[2], "%)")

p4 <- ggplot(pcoa_after_df, aes(x = Axis1, y = Axis2, fill = Batch)) +
  geom_point(shape = 21, size = 3, color = "black", stroke = 0.3) +
  stat_ellipse(color = NA, fill = NA, size = 0) + 
  labs(
    x = xlab_txt,
    y = ylab_txt,
    fill = "Batch"
  ) +
  scale_fill_manual(values = c(
    "Vilà-Quintana_2024" = "#7DABD1",
    "Kartal-1_2022"      = "#C7D3E5", 
    "Matsukawa_2021"     = "#FAADD5",
    "Nagata_2022"        = "#BC8CC3"
  )) +
  theme_classic() +
  theme(
    legend.position      = c(1, 1),
    legend.justification = c("right", "top"),
    legend.background    = element_blank(),
    legend.key           = element_blank(),
    legend.title         = element_blank(),
    legend.text          = element_text(size = 14),
    legend.key.size      = unit(0.9, "lines"),
    plot.title           = element_text(hjust = 0.5),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold")
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = "After~correction*\",\"~R^2==\"1.8\"*'%'",
    hjust = -0.05, vjust = 1.5,
    parse = TRUE,
    size = 7.5,
    fontface = "bold"
  )
print(p4)

p3 <- p3 + stat_ellipse(level = 0.95, show.legend = F,linetype = 2)
p4 <- p4 + stat_ellipse(level = 0.95, show.legend = F,linetype = 2)
print(p3)
print(p4)

p5 <- plot_grid(
  p1, p2, p3, p4, 
  # labels = c("C", "D", "E", "F"),     
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

p34 <- plot_grid(
  p3, p4, 
  # labels = c("C", "D", "E", "F"),     
  # label_size = 20,          
  ncol = 2,
  align = "hv"
)
print(p34)

ggsave("fig/batch_stool.pdf", plot = p12, width = 20, height = 10)
ggsave("fig/batch_saliva.pdf", plot = p34, width = 20, height = 10)

ggsave("fig/batch.png", plot = p5, width = 15, height = 10, dpi = 500)
ggsave("fig/batch_all.pdf", plot = p5, width = 25, height = 20)
