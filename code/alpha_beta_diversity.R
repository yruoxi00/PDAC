library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(zCompositions)
library(MMUPHin)
library(ggplot2)
library(ggpubr)
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

## batch effect correction 

meta_data_stool <- stool_data %>%
  column_to_rownames("ID") %>%
  dplyr::select(-Group) %>%
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

species_adj_stool <- adjust_batch(
  feature_abd = t(species),
  batch       = "Data",
  covariates  = c("Group_new"),
  data        = meta_data_stool,
  control     = list(verbose = TRUE, diagnostics = TRUE)
)
species_adj_matrix_stool <- species_adj_stool$feature_abd_adj
species_adj_matrix_stool <- sweep(species_adj_matrix_stool, 2, colSums(species_adj_matrix_stool), "/")
summary(colSums(species_adj_matrix_stool))  # should all be 1

## validation (PERMANOVA and PCoA)

dist_before <- vegdist(species, method = "bray")
dist_after <- vegdist(t(species_adj_matrix_stool), method = "bray")

adonis_before <- adonis2(dist_before ~ Data, data = meta_data_stool)
adonis_after <- adonis2(dist_after ~ Data, data = meta_data_stool)
print(adonis_before)
print(adonis_after)

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
  column_to_rownames("ID") %>%
  dplyr::select(-Group) %>%
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

species_adj_saliva <- adjust_batch(
  feature_abd = t(species),
  batch       = "Data",
  covariates  = c("Group_new"),
  data        = meta_data_saliva,
  control     = list(verbose = TRUE, diagnostics = TRUE)
)
species_adj_matrix_saliva <- species_adj_saliva$feature_abd_adj
species_adj_matrix_saliva <- sweep(species_adj_matrix_saliva, 2, colSums(species_adj_matrix_saliva), "/")
summary(colSums(species_adj_matrix_saliva))  # should all be 1

## validation (PERMANOVA and PCoA)

dist_before <- vegdist(species, method = "bray")
dist_after <- vegdist(t(species_adj_matrix_saliva), method = "bray")

adonis_before <- adonis2(dist_before ~ Data, data = meta_data_saliva)
adonis_after <- adonis2(dist_after ~ Data, data = meta_data_saliva)
print(adonis_before)
print(adonis_after)

# Alpha diversity (Shannon + Observed richness)
alpha_diversity <- function(species_adj_matrix, metadata) {
  X <- t(species_adj_matrix)
  
  common_ids <- intersect(rownames(X), rownames(metadata))
  if (length(common_ids) == 0) stop("No overlapping sample IDs between abundance and metadata.")
  X <- X[common_ids, , drop = FALSE]
  md <- metadata[common_ids, , drop = FALSE]
  
  shannon  <- diversity(X, index = "shannon")
  observed <- specnumber(X)
  
  out <- data.frame(
    Sample   = common_ids,
    Shannon  = as.numeric(shannon),
    Observed = as.numeric(observed),
    Data     = md$Data,
    Group    = md$Group_new,
    row.names = NULL
  )
  out
}

# Beta Diversity (CLR + Euclidean)
beta_clr_pseudocount <- function(abundance_matrix, metadata, pseudocount = 1e-6) {
  
  common_ids <- intersect(rownames(metadata), colnames(abundance_matrix))
  if (length(common_ids) < 2) stop("Fewer than 2 overlapping samples.")
  X  <- abundance_matrix[, common_ids, drop = FALSE] + pseudocount
  md <- metadata[common_ids, , drop = FALSE]
  
  clr_mat <- apply(X, 2, function(x) log(x) - mean(log(x)))  
  
  D <- dist(t(clr_mat), method = "euclidean")  
  
  list(dist = D, metadata = md, clr = clr_mat)
}

grp_cols <- c(
  "Control"              = "#4C72B0",
  "Pancreatic Cancer"    = "#ED6A6B"
)

# Alpha diversity (stool)
alpha_stool <- alpha_diversity(species_adj_matrix_stool, meta_data_stool)

alpha_stool$Group <- factor(
  alpha_stool$Group, 
  levels = c("Control", "Pancreatic Cancer")
)

sample_sizes <- alpha_stool %>%
  group_by(Group) %>%
  summarise(n = n())

p1 <- ggplot(alpha_stool, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(width = 0.6, color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 4,
    hide.ns = FALSE,
    comparisons = list(
      c("Control", "Pancreatic Cancer")
    ), 
    step.increase = 0.12
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#4C72B0",
      "Pancreatic Cancer" = "#ED6A6B"
    )
  ) +
  labs(
    title = "Alpha Diversity: Shannon Index (Stool)",
    x = NULL,
    y = "Shannon Diversity Index"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  geom_text(
    data = sample_sizes,
    aes(
      x = Group,
      y = min(alpha_stool$Shannon) - 0.1 * diff(range(alpha_stool$Shannon)),
      label = paste0("n = ", n)
    ),
    size = 4, vjust = 1
  )
print(p1)

# Beta diversity (stool)
meta_data_stool$Group_new <- factor(
  meta_data_stool$Group_new,
  levels = c("Control", "Pancreatic Cancer")
)

stool_res <- beta_clr_pseudocount(species_adj_matrix_stool, meta_data_stool, pseudocount = 1e-6)

adonis_aitch_stool <- adonis2(stool_res$dist ~ Data + Group_new, 
                              data = stool_res$metadata, 
                              by = "terms")
print(adonis_aitch_stool)

perm_p_stool  <- adonis_aitch_stool["Group_new", "Pr(>F)"]
perm_r2_stool <- adonis_aitch_stool["Group_new", "R2"]

label_stool <- paste0("PERMANOVA: p = ", signif(perm_p_stool, 3))

pcoa_stool <- cmdscale(stool_res$dist, k = 2, eig = TRUE)
var_stool <- round(pcoa_stool$eig / sum(pcoa_stool$eig) * 100, 2)
pcoa_stool_df <- data.frame(
  Axis1 = pcoa_stool$points[,1],
  Axis2 = pcoa_stool$points[,2],
  Group = stool_res$metadata$Group_new,
  Batch = stool_res$metadata$Data,
  row.names = rownames(stool_res$metadata)
)

p2 <- ggplot(pcoa_stool_df, aes(Axis1, Axis2, fill = Group)) +
  geom_point(shape = 21, size = 2, color = "black", stroke = 0.3) +
  scale_fill_manual(values = grp_cols, drop = FALSE) +
  labs(
    title = "PCoA (Robust Aitchison, Stool)",
    x = paste0("PCoA1 (", var_stool[1], "%)"),
    y = paste0("PCoA2 (", var_stool[2], "%)")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11)
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.5,
    label = label_stool,
    size = 4
  )
print(p2)

# Alpha diversity (saliva)
alpha_saliva <- alpha_diversity(species_adj_matrix_saliva, meta_data_saliva)

alpha_saliva$Group <- factor(
  alpha_saliva$Group, 
  levels = c("Control", "Pancreatic Cancer")
)

sample_sizes <- alpha_saliva %>%
  group_by(Group) %>%
  summarise(n = n())

p3 <- ggplot(alpha_saliva, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(width = 0.6, color = "black", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 4,
    hide.ns = FALSE,
    comparisons = list(
      c("Control", "Pancreatic Cancer")
    ), 
    step.increase = 0.12
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#4C72B0",
      "Pancreatic Cancer" = "#ED6A6B"
    )
  ) +
  labs(
    title = "Alpha Diversity: Shannon Index (Saliva)",
    x = NULL,
    y = "Shannon Diversity Index"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  geom_text(
    data = sample_sizes,
    aes(
      x = Group,
      y = min(alpha_stool$Shannon) - 0.1 * diff(range(alpha_stool$Shannon)),
      label = paste0("n = ", n)
    ),
    size = 4, vjust = 1
  )
print(p3)

# Beta diversity (saliva)
meta_data_saliva$Group_new <- factor(
  meta_data_saliva$Group_new,
  levels = c("Control", "Pancreatic Cancer")
)

saliva_res <- beta_clr_pseudocount(species_adj_matrix_saliva, meta_data_saliva, pseudocount = 1e-6)

adonis_aitch_saliva <- adonis2(saliva_res$dist ~ Data + Group_new, 
                               data = saliva_res$metadata, 
                               by = "term")
print(adonis_aitch_saliva)

perm_p_saliva  <- adonis_aitch_saliva["Group_new", "Pr(>F)"]
perm_r2_saliva <- adonis_aitch_saliva["Group_new", "R2"]

label_saliva <- paste0("PERMANOVA: p = ", signif(perm_p_saliva, 3))

pcoa_saliva <- cmdscale(saliva_res$dist, k = 2, eig = TRUE)
var_saliva <- round(pcoa_saliva$eig / sum(pcoa_saliva$eig) * 100, 2)
pcoa_saliva_df <- data.frame(
  Axis1 = pcoa_saliva$points[,1],
  Axis2 = pcoa_saliva$points[,2],
  Group = saliva_res$metadata$Group_new,
  Batch = saliva_res$metadata$Data,
  row.names = rownames(saliva_res$metadata)
)

p4 <- ggplot(pcoa_saliva_df, aes(Axis1, Axis2, fill = Group)) +
  geom_point(shape = 21, size = 2, color = "black", stroke = 0.3) +
  scale_fill_manual(values = grp_cols, drop = FALSE) +
  labs(
    title = "PCoA (Robust Aitchison, Saliva)",
    x = paste0("PCoA1 (", var_stool[1], "%)"),
    y = paste0("PCoA2 (", var_stool[2], "%)")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11)
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.05, vjust = 1.5,
    label = label_saliva,
    size = 4
  )
print(p4)

p5 <- plot_grid(
  p1, p3, p2, p4, 
  labels = c("A", "B", "C", "D"),     
  label_size = 20,          
  ncol = 2,
  align = "hv"
)
print(p5)

ggsave("fig/alpha_beta.png", plot = p5, width = 15, height = 10, dpi = 500)
ggsave("fig/alpha_beta.pdf", plot = p5, width = 15, height = 10)

### Meta Data 
meta_data <- read_csv("data/meta_data.csv")
species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

# species data 

species_matrix <- species %>%
  column_to_rownames(var = "ID")
species <- as.data.frame(t(species_matrix))
species <- species[meta_data$ID, ]
species <- species / 100

abundance_threshold <- 2e-4   # 0.0002
prevalence_threshold <- 0.10  # 10%

taxa <- lapply(split(meta_data$ID, meta_data$Data), function(sample_ids) {
  species_sub <- species[sample_ids, , drop = FALSE]  
  apply(species_sub, 2, function(x) {
    mean(x > abundance_threshold) >= prevalence_threshold
  })
})

taxa <- Reduce("|", taxa)
cat("Original number of taxa:", ncol(species), "\n")
species <- species[, taxa]
cat("Filtered number of taxa:", ncol(species), "\n")

meta_data <- meta_data %>%
  mutate(
    ID   = as.character(ID),
    Geo  = factor(Geo),
    Sex  = factor(Sex),
    Type = factor(Type),
    Age  = as.numeric(Age)
  ) %>%
  column_to_rownames("ID")

clr_transform <- function(x) log(x) - mean(log(x))

species_clr <- apply(t(species + 1e-6), 1, clr_transform)  
dist_matrix <- dist(species_clr, method = "euclidean")      
samples <- attr(dist_matrix, "Labels")
md <- meta_data[samples, , drop = FALSE]

md <- md %>%
  mutate(
    Geo  = factor(Geo),
    Sex  = factor(Sex),
    Type = factor(Type),
    Age  = suppressWarnings(as.numeric(Age))
  )

host_vars <- c("Geo","Type","Age","Sex")
keep <- complete.cases(md[, host_vars])
cat("Total samples:", nrow(md), "\n")
cat("Dropping due to NA in host vars:", sum(!keep), "\n")

md_cc <- md[keep, , drop = FALSE]

Dmat <- as.matrix(dist_matrix)
D_cc <- as.dist(Dmat[keep, keep])

### PERMANOVA
adonis_res <- adonis2(
  D_cc ~ Geo + Age + Sex + Type,
  data = md_cc,
  permutations = 999,
  by = "margin"
)
print(adonis_res)
