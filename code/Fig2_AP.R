library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(MMUPHin)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(ggtext)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/ap.csv")

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

species <- species_raw[, species_union, drop = FALSE]

# ===================== Batch Effect ===================== #

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

# validation (PERMANOVA and PCoA)

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
    "Liu_2024" = "#7DABD1",
    "Yu_2021" = "#F6D963",
    "Song_2025" = "#FAADD5",
    "Kartal-1_2022" = "#BC8CC3"
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
    label = "Before~correction*\",\"~R^2==\"4.6\"*'%'",
    hjust = -0.05, vjust = 1.5,
    parse = TRUE,
    size = 7.5, 
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
    "Liu_2024" = "#7DABD1",
    "Yu_2021" = "#F6D963",
    "Song_2025" = "#FAADD5",
    "Kartal-1_2022" = "#BC8CC3"
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
    label = "After~correction*\",\"~R^2==\"2.2\"*'%'",
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

# ===================== Species Distribution ===================== #

plain <- function(x, ...) format(x, ..., scientific = FALSE, drop0trailing = TRUE)

studycol <- c(
  "Liu_2024" = "#7DABD1",
  "Yu_2021" = "#F6D963",
  "Song_2025" = "#FAADD5",
  "Kartal-1_2022" = "#BC8CC3"
)

map_study_names <- function(df, col = "study") {
  df %>%
    mutate(
      !!col := case_when(
        .data[[col]] == "_relabDat1" ~ "Liu_2024",
        .data[[col]] == "_relabDat3" ~ "Yu_2021",
        .data[[col]] == "_relabDat4" ~ "Song_2025",
        .data[[col]] == "_relabDat5" ~ "Kartal-1_2022",
        TRUE ~ .data[[col]]
      )
    ) %>%
    mutate(
      !!col := factor(
        .data[[col]],
        levels = c("Liu_2024","Yu_2021","Song_2025","Kartal-1_2022")
      )
    )
}

get_top_species_labels <- function(stool_data_sub, species_all, top_n = 25) {
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
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n)
  
  overlap <- spe_all %>%
    filter(present != n_studies & present != 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n)
  
  solo <- spe_all %>%
    filter(present == 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n)
  
  spe_plot <- bind_rows(universal, overlap, solo)
  
  labels <- spe_plot$species %>%
    gsub("^s__", "", .) %>%
    gsub("_", " ", .)
  
  unique(labels)
}

# ---- subsets ----
stool_PC <- stool_data %>%
  filter(Group_new == "AP")

stool_con <- stool_data %>%
  filter(Group_new == "Control")

# ---- overlapping species (using new filtered species matrix) ----
pc_labels  <- get_top_species_labels(stool_PC, species, top_n = 25)
con_labels <- get_top_species_labels(stool_con, species, top_n = 25)
overlap_labels <- intersect(pc_labels, con_labels)

# ---- main plotting function ----
make_species_dist_plot <- function(stool_data_sub, species_all, title_text,
                                   top_n = 25, overlap_species = NULL) {
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
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    mutate(group = "Universal")
  
  overlap <- spe_all %>%
    filter(present != n_studies & present != 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    mutate(group = "Overlapping")
  
  solo <- spe_all %>%
    filter(present == 1) %>%
    mutate(mean_abund = rowMeans(across(contains("_relab")), na.rm = TRUE)) %>%
    arrange(desc(mean_abund)) %>%
    slice_head(n = top_n) %>%
    mutate(group = "Singular")
  
  spe_plot <- bind_rows(universal, overlap, solo)
  
  spe_plot2 <- spe_plot %>%
    select(-mean_abund) %>%
    pivot_longer(cols = contains("_relab"),
                 names_to = "study", values_to = "relab") %>%
    filter(relab > 0) %>%
    select(species, group, study, relab) %>%
    mutate(
      species_clean = gsub("^s__", "", species),
      species_clean = gsub("_", " ", species_clean)
    )
  
  spe_order <- spe_plot2 %>%
    group_by(group, species_clean) %>%
    summarise(mean = mean(relab), .groups = "drop_last") %>%
    arrange(group, desc(mean))
  
  spe_plot2$species_clean <- factor(
    spe_plot2$species_clean,
    levels = spe_order$species_clean
  )
  spe_plot2$group <- factor(
    spe_plot2$group,
    levels = c("Universal","Overlapping","Singular")
  )
  
  spe_plot2 <- map_study_names(spe_plot2, col = "study")
  
  label_fun <- function(x) {
    if (is.null(overlap_species)) {
      paste0("<i>", x, "</i>")
    } else {
      ifelse(
        x %in% overlap_species,
        paste0("<b><i>", x, "</i></b>"),
        paste0("<i>", x, "</i>")
      )
    }
  }
  
  ggplot(spe_plot2, aes(x = species_clean, y = relab)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
    geom_point(aes(fill = study), shape = 21, size = 2.6, alpha = 0.85,
               position = position_jitter(width = 0.08, height = 0)) +
    scale_fill_manual(values = studycol, drop = FALSE) +
    ylab("Relative abundance (%)") + xlab("") +
    scale_y_log10(breaks = c(10, 1, 0.1, 0.01), labels = plain) +
    scale_x_discrete(labels = label_fun) +
    facet_wrap(~ group, scales = "free_x") +
    theme(
      strip.text.x = element_text(hjust = 0.5, color = "black", size = 14),
      axis.text.x  = ggtext::element_markdown(
        angle = 90, vjust = 0.5, hjust = 1, size = 11, color = "black"
      ),
      axis.text.y  = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = c(.9, .85),
      legend.title = element_blank(),
      legend.text = element_text(size = 10, color = "black"),
      legend.background = element_blank(),
      legend.key.height = unit(.3, 'cm'),
      legend.key.width  = unit(.3, 'cm'),
      legend.key = element_rect(fill = "white")
    ) +
    guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
    ggtitle(title_text)
}

p3 <- make_species_dist_plot(
  stool_data_sub  = stool_PC,
  species_all     = species,
  title_text      = "Stool Species Distribution (AP)",
  top_n           = 25,
  overlap_species = overlap_labels
)

p4 <- make_species_dist_plot(
  stool_data_sub  = stool_con,
  species_all     = species,
  title_text      = "Stool Species Distribution (Control)",
  top_n           = 25,
  overlap_species = overlap_labels
)

print(p3)
print(p4)

# ===================== Top 30 Species ===================== #

get_top30_labels_bar <- function(group_name, stool_data, n_top = 30) {
  if (group_name == "AP") {
    species_mat <- species_pdac
  } else if (group_name == "Control") {
    species_mat <- species_control
  } else {
    stop("Unknown group_name")
  }
  
  idx_ids <- stool_data$ID[stool_data$Group_new == group_name]
  species_sub <- species_mat[idx_ids, , drop = FALSE]
  
  species_means <- colMeans(species_sub, na.rm = TRUE)
  
  df <- data.frame(
    species        = names(species_means),
    mean_abundance = species_means
  )
  
  remove_pattern <- "(SGB|GGB|unknown|uncultured|unclassified|metagenome)"
  df <- df %>%
    dplyr::filter(!grepl(remove_pattern, species, ignore.case = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    slice(1:n_top)
  
  labels <- df$species %>%
    gsub("^s__", "", .) %>%
    gsub("_", " ", .)
  
  unique(labels)
}

pc_labels_bar  <- get_top30_labels_bar("AP", stool_data, n_top = 30)
con_labels_bar <- get_top30_labels_bar("Control", stool_data, n_top = 30)
overlap_bar    <- intersect(pc_labels_bar, con_labels_bar)

make_top30_barplot <- function(group_name, stool_data,
                               n_top = 30, plot_title = NULL,
                               overlap_species = NULL) {
  if (group_name == "AP") {
    species_mat <- species_pdac
  } else if (group_name == "Control") {
    species_mat <- species_control
  } else {
    stop("Unknown group_name")
  }
  
  idx_ids <- stool_data$ID[stool_data$Group_new == group_name]
  species_sub <- species_mat[idx_ids, , drop = FALSE]
  
  species_means <- colMeans(species_sub, na.rm = TRUE)
  
  df <- data.frame(
    species        = names(species_means),
    mean_abundance = species_means
  )
  
  remove_pattern <- "(SGB|GGB|unknown|uncultured|unclassified|metagenome)"
  df <- df %>%
    dplyr::filter(!grepl(remove_pattern, species, ignore.case = TRUE)) %>%
    arrange(desc(mean_abundance)) %>%
    slice(1:n_top)
  
  df$species <- df$species %>%
    gsub("^s__", "", .) %>%
    gsub("_", " ", .)
  
  df$species <- factor(df$species, levels = rev(df$species))
  df$mean_abundance <- df$mean_abundance * 100
  
  label_fun <- function(x) {
    if (is.null(overlap_species)) {
      paste0("<i>", x, "</i>")
    } else {
      ifelse(
        x %in% overlap_species,
        paste0("<b><i>", x, "</i></b>"),
        paste0("<i>", x, "</i>")
      )
    }
  }
  
  ggplot(df, aes(x = mean_abundance, y = species)) +
    geom_bar(stat = "identity", aes(fill = mean_abundance)) +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    scale_y_discrete(labels = label_fun) +
    xlab("Mean relative abundance (%)") +
    ylab("") +
    ggtitle(plot_title) +
    theme_classic() +
    theme(
      axis.text.y  = ggtext::element_markdown(size = 10),
      axis.text.x  = element_text(size = 12),
      axis.title.x = element_text(size = 13, face = "bold"),
      plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
}

p5 <- make_top30_barplot(
  group_name      = "AP",
  stool_data      = stool_data,
  n_top           = 30,
  plot_title      = "Top 30 Stool Species (AP)",
  overlap_species = overlap_bar
)

p6 <- make_top30_barplot(
  group_name      = "Control",
  stool_data      = stool_data,
  n_top           = 30,
  plot_title      = "Top 30 Stool Species (Control)",
  overlap_species = overlap_bar
)

print(p5)
print(p6)

# ===================== Genera Diversity ===================== #

control_df <- stool_data %>% filter(Group_new == "Control")
pdac_df    <- stool_data %>% filter(Group_new == "AP")

compute_genus_per_sample <- function(group_df, species_grp) {
  species_sub <- species_grp[group_df$ID, , drop = FALSE]
  nonzero_species <- colSums(species_sub > 0, na.rm = TRUE) > 0
  species_sub <- species_sub[, nonzero_species, drop = FALSE]
  species_names <- colnames(species_sub)
  genus_vec <- species_names %>%
    str_replace("^s__", "") %>%
    str_replace("^s_", "") %>%
    str_replace("^s ", "") %>%
    str_split_fixed(pattern = "_", n = 2) %>%
    .[, 1]
  keep_idx <- genus_vec != "" &
    !is.na(genus_vec) &
    !str_detect(genus_vec, regex("SGB|GGB", ignore_case = TRUE)) &
    !str_detect(genus_vec, regex("unknown|unclassified|uncultured", ignore_case = TRUE))
  species_sub <- species_sub[, keep_idx, drop = FALSE]
  genus_vec   <- genus_vec[keep_idx]
  map_df <- tibble(
    species = colnames(species_sub),
    genus   = genus_vec
  )
  long_df <- as.data.frame(species_sub) %>%
    rownames_to_column("ID") %>%
    pivot_longer(-ID, names_to = "species", values_to = "abundance") %>%
    left_join(map_df, by = "species")
  per_sample_genus <- long_df %>%
    filter(abundance > 0) %>%
    group_by(ID, genus) %>%
    summarise(n_species = n_distinct(species), .groups = "drop")
  per_sample_genus
}

control_genus_sample <- compute_genus_per_sample(control_df, species_control)
pdac_genus_sample    <- compute_genus_per_sample(pdac_df, species_pdac)

control_genus_summary <- control_genus_sample %>%
  group_by(genus) %>%
  summarise(mean_n_species = mean(n_species), .groups = "drop") %>%
  arrange(desc(mean_n_species))

pdac_genus_summary <- pdac_genus_sample %>%
  group_by(genus) %>%
  summarise(mean_n_species = mean(n_species), .groups = "drop") %>%
  arrange(desc(mean_n_species))

overlap_genera <- intersect(control_genus_summary$genus,
                            pdac_genus_summary$genus)

combined_overlap <- control_genus_summary %>%
  filter(genus %in% overlap_genera) %>%
  rename(mean_control = mean_n_species) %>%
  left_join(
    pdac_genus_summary %>%
      filter(genus %in% overlap_genera) %>%
      rename(mean_pdac = mean_n_species),
    by = "genus"
  ) %>%
  mutate(mean_overall = (mean_control + mean_pdac) / 2) %>%
  arrange(desc(mean_overall))

top30_overlap_genus <- combined_overlap %>%
  slice_head(n = 30) %>%
  pull(genus)

pval_df <- lapply(top30_overlap_genus, function(g) {
  x <- control_genus_sample %>% filter(genus == g) %>% pull(n_species)
  y <- pdac_genus_sample %>% filter(genus == g) %>% pull(n_species)
  p <- wilcox.test(x, y)$p.value
  tibble(genus = g, p_value = p)
}) %>% bind_rows()

pval_df <- pval_df %>%
  left_join(combined_overlap, by = "genus") %>%
  mutate(
    is_sig = p_value < 0.05,
    
    direction_control = case_when(
      is_sig & mean_control > mean_pdac ~ "+",
      is_sig & mean_control < mean_pdac ~ "-",
      TRUE ~ ""
    ),
    sign_control = ifelse(
      is_sig & direction_control != "",
      paste0(" (", direction_control, ")"),
      ""
    ),
    
    direction_pdac = case_when(
      direction_control == "+" ~ "-",
      direction_control == "-" ~ "+",
      TRUE ~ ""
    ),
    sign_pdac = ifelse(
      is_sig & direction_pdac != "",
      paste0(" (", direction_pdac, ")"),
      ""
    ),
    
    label_control = ifelse(
      is_sig,
      paste0("\\* ", "**", genus, "**", sign_control),
      genus
    ),
    
    label_pdac = ifelse(
      is_sig,
      paste0("\\* ", "**", genus, "**", sign_pdac),
      genus
    )
  )

level_vec <- rev(top30_overlap_genus)

control_plot_df <- control_genus_sample %>%
  filter(genus %in% top30_overlap_genus) %>%
  mutate(genus = factor(genus, levels = level_vec))

pdac_plot_df <- pdac_genus_sample %>%
  filter(genus %in% top30_overlap_genus) %>%
  mutate(genus = factor(genus, levels = level_vec))

control_means <- control_genus_summary %>%
  filter(genus %in% top30_overlap_genus) %>%
  mutate(genus = factor(genus, levels = level_vec))

pdac_means <- pdac_genus_summary %>%
  filter(genus %in% top30_overlap_genus) %>%
  mutate(genus = factor(genus, levels = level_vec))

label_map_control <- setNames(pval_df$label_control, pval_df$genus)
label_map_pdac    <- setNames(pval_df$label_pdac,    pval_df$genus)

p7 <- ggplot() +
  geom_col(
    data = control_means,
    aes(x = mean_n_species, y = genus),
    width = 0.6, fill = "#1f78b4"
  ) +
  geom_point(
    data = control_plot_df,
    aes(x = n_species, y = genus),
    alpha = 0.4, size = 0.7,
    position = position_jitter(height = 0)
  ) +
  scale_y_discrete(labels = label_map_control) +
  labs(
    title = "Top 30 Overlapping Stool Genera (Control)",
    x     = "Average number of species",
    y     = ""
  ) +
  theme_classic() +
  theme(
    axis.text.y  = ggtext::element_markdown(size = 10),
    axis.text.x  = element_text(size = 12),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

p8 <- ggplot() +
  geom_col(
    data = pdac_means,
    aes(x = mean_n_species, y = genus),
    width = 0.6, fill = "#1f78b4"
  ) +
  geom_point(
    data = pdac_plot_df,
    aes(x = n_species, y = genus),
    alpha = 0.4, size = 0.7,
    position = position_jitter(height = 0)
  ) +
  scale_y_discrete(labels = label_map_pdac) +
  labs(
    title = "Top 30 Overlapping Stool Genera (AP)",
    x     = "Average number of species",
    y     = ""
  ) +
  theme_classic() +
  theme(
    axis.text.y  = ggtext::element_markdown(size = 10),
    axis.text.x  = element_text(size = 12),
    axis.title.x = element_text(size = 13, face = "bold"),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

print(p7)
print(p8)

# ===================== Final Plot ===================== #

panel <- ggplot() + theme_void()

row1 <- plot_grid(
  panel, panel, p1, p2,
  labels = c("a", "b", "c", "d"),
  label_size = 20,
  ncol = 4,
  align = "hv"
)

row2 <- plot_grid(
  p4, p3,
  labels = c("e", "f"),
  label_size = 20,
  ncol = 2,
  align = "hv"
)

row3 <- plot_grid(
  p7, p8, p6, p5, 
  labels = c("g", "h", "i", "j"),
  label_size = 20,
  ncol = 4,
  align = "hv"
)

p <- plot_grid(
  row1, row2, row3,
  ncol = 1
)
print(p)

ggsave("fig/Fig2_AP.pdf", p, width = 30, height = 20, limitsize = FALSE)
