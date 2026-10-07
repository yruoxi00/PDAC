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

setwd("/Users/yrx/Downloads/pancreatitis")

## ===================== STOOL SECTION ===================== ##

## --------------------- Data prep --------------------- ##

pdac    <- read_csv("data/PDAC.csv")
species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

set.seed(123)

stool_data <- pdac %>%
  filter(Type == "Stool")

cohort_var <- "Study"

stool_data <- stool_data %>%
  mutate(
    cohort    = .data[[cohort_var]],
    Group_new = factor(Group_new, levels = c("Control", "Pancreatic Cancer"))
  )

## species matrix
species_matrix <- species %>%
  column_to_rownames(var = "ID")
species <- as.data.frame(t(species_matrix))
species <- species[stool_data$ID, ]
species <- species / 100

abundance_threshold  <- 2e-4   # 0.0002
prevalence_threshold <- 0.10   # 10%

taxa <- lapply(split(stool_data$ID, stool_data$Data), function(sample_ids) {
  species_sub <- species[sample_ids, , drop = FALSE]
  apply(species_sub, 2, function(x) {
    mean(x > abundance_threshold) >= prevalence_threshold
  })
})

taxa <- Reduce("|", taxa)
species_stool <- species[, taxa]

stopifnot(all(rownames(species_stool) == stool_data$ID))

meta_stool <- stool_data %>%
  column_to_rownames("ID")

## ---------- Prevalence (for prevalence-integrated importance) ---------- ##
# Prevalence across all stool samples for each species
prev_stool <- colMeans(species_stool > 0)  # names = taxa

## --------------------- RF helper (prevalence-integrated) --------------------- ##

run_rf_one_prev <- function(X, y, cohort_name, prev_vec, top_n = 50) {
  y <- factor(y, levels = c("Control", "Pancreatic Cancer"))
  
  rf_model <- randomForest(
    x          = X,
    y          = y,
    ntree      = 1000,
    importance = TRUE
  )
  
  imp <- as.data.frame(importance(rf_model))
  imp$feature <- rownames(imp)
  imp$cohort  <- cohort_name
  
  # attach prevalence for each feature
  imp$prevalence     <- prev_vec[imp$feature]
  imp$Gini_prev2.5   <- imp$MeanDecreaseGini * (imp$prevalence ^ 2.5)
  
  imp %>%
    arrange(desc(Gini_prev2.5)) %>%
    slice(1:top_n)
}

## --------------------- RF (All + by cohort) --------------------- ##

X_all <- species_stool
y_all <- stool_data$Group_new

rf_all <- run_rf_one_prev(X_all, y_all, "All_stool", prev_vec = prev_stool, top_n = 50)

cohort_levels <- sort(unique(stool_data$cohort))

rf_by_cohort <- map_dfr(cohort_levels, function(cc) {
  idx <- stool_data$cohort == cc
  X   <- X_all[idx, , drop = FALSE]
  y   <- y_all[idx]
  
  run_rf_one_prev(X, y, cohort_name = cc, prev_vec = prev_stool, top_n = 50)
})

rf_top50 <- bind_rows(rf_all, rf_by_cohort)

rf_top50 <- rf_top50 %>%
  filter(
    (!grepl("unknown", feature, ignore.case = TRUE) |
       grepl("SGB", feature, ignore.case = TRUE)) &
      !grepl("GGB", feature, ignore.case = TRUE)
  )

feature_overlap <- rf_top50 %>%
  group_by(feature) %>%
  summarise(n_cohort = n_distinct(cohort), .groups = "drop")

rf_top50 <- rf_top50 %>%
  left_join(feature_overlap, by = "feature") %>%
  mutate(
    overlap_group = ifelse(
      n_cohort >= 3,
      "≥3 cohorts",
      "<3 cohorts"
    )
  )

## --------------------- Cohort labels --------------------- ##

cohort_counts <- stool_data %>%
  group_by(cohort, Group_new) %>%
  tally(name = "n") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    label = sprintf("%s\nC=%s, PC=%s",
                    cohort,
                    Control,
                    `Pancreatic Cancer`)
  )

all_counts <- stool_data %>%
  count(Group_new) %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    cohort = "All_stool",
    label  = sprintf("All_stool\nC=%s, PC=%s",
                     Control,
                     `Pancreatic Cancer`)
  )

cohort_labels_tbl <- bind_rows(
  all_counts %>% select(cohort, label),
  cohort_counts %>% select(cohort, label)
)

stool_x_labels <- setNames(cohort_labels_tbl$label, cohort_labels_tbl$cohort)

cohort_levels_rf <- sort(unique(rf_top50$cohort[rf_top50$cohort != "All_stool"]))

if ("Matsukawa_2021" %in% cohort_levels_rf) {
  cohort_order <- c(
    "All_stool",
    setdiff(cohort_levels_rf, "Matsukawa_2021"),
    "Matsukawa_2021"
  )
} else {
  cohort_order <- c("All_stool", cohort_levels_rf)
}

## --------------------- Mean differences & direction --------------------- ##

stool_long <- stool_data %>%
  select(cohort, Group_new) %>%
  bind_cols(as.data.frame(species_stool)) %>%
  pivot_longer(
    cols      = -c(cohort, Group_new),
    names_to  = "feature",
    values_to = "abundance"
  )

mean_diff_tbl <- stool_long %>%
  group_by(cohort, feature, Group_new) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = mean_abundance
  ) %>%
  mutate(
    mean_diff = `Pancreatic Cancer` - Control,
    direction = case_when(
      mean_diff >  0  ~ "Enriched",
      mean_diff <  0  ~ "Depleted",
      TRUE            ~ "No difference"
    )
  ) %>%
  select(cohort, feature, mean_diff, direction)

mean_diff_all <- stool_data %>%
  select(Group_new) %>%
  bind_cols(as.data.frame(species_stool)) %>%
  pivot_longer(
    cols      = -Group_new,
    names_to  = "feature",
    values_to = "abundance"
  ) %>%
  group_by(feature, Group_new) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = mean_abundance
  ) %>%
  mutate(
    cohort    = "All_stool",
    mean_diff = `Pancreatic Cancer` - Control,
    direction = case_when(
      mean_diff >  0  ~ "Enriched",
      mean_diff <  0  ~ "Depleted",
      TRUE            ~ "No difference"
    )
  ) %>%
  select(cohort, feature, mean_diff, direction)

mean_diff_tbl2 <- bind_rows(mean_diff_tbl, mean_diff_all)

rf_top50 <- rf_top50 %>%
  left_join(mean_diff_tbl2, by = c("cohort", "feature"))

dir_consistency <- rf_top50 %>%
  group_by(feature) %>%
  summarise(
    n_dir_cohort = n_distinct(cohort[direction != "No difference"]),
    n_pos        = sum(direction == "Enriched", na.rm = TRUE),
    n_neg        = sum(direction == "Depleted", na.rm = TRUE),
    n_dir_types  = (n_pos > 0) + (n_neg > 0),
    consistent_ge4 = (n_dir_cohort >= 3 & n_dir_types == 1),
    .groups = "drop"
  )

## --------------------- Order of features --------------------- ##

top50_pooled_order <- rf_all %>%
  arrange(desc(Gini_prev2.5)) %>%   # use prevalence-integrated importance
  pull(feature) %>%
  unique()

df_plot_rf <- rf_top50 %>%
  filter(feature %in% top50_pooled_order) %>%
  mutate(
    cohort  = factor(cohort,  levels = cohort_order),
    feature = factor(feature, levels = rev(top50_pooled_order)),
    overlap_group = factor(overlap_group,
                           levels = c("<3 cohorts", "≥3 cohorts")),
    direction = factor(
      direction,
      levels = c("Enriched", "Depleted", "No difference")
    )
  ) %>%
  filter(
    !is.na(Gini_prev2.5),
    !is.na(feature),
    !is.na(cohort)
  )

feature_overlap_plot <- df_plot_rf %>%
  distinct(cohort, feature) %>% 
  group_by(feature) %>%
  summarise(n_cohort = n_distinct(cohort), .groups = "drop")

## --------------------- Y-axis labels (clean + colored) --------------------- ##

dir_all <- mean_diff_tbl2 %>%
  filter(cohort == "All_stool") %>%
  select(feature, direction)

stool_rf_label_tbl <- feature_overlap_plot %>%
  left_join(dir_all, by = "feature") %>%
  left_join(dir_consistency %>% select(feature, consistent_ge4),
            by = "feature") %>%
  mutate(
    feature_clean = feature %>%
      gsub("^s__", "", .) %>%
      gsub("^s_",  "", .) %>%
      gsub("^s ",  "", .) %>%
      gsub("_", " ", .),
    
    highlight = (n_cohort >= 3 & !is.na(consistent_ge4) & consistent_ge4),
    
    label = case_when(
      highlight & direction == "Enriched" ~
        sprintf("<span style='color:#E41A1C;'><i>%s</i></span>", feature_clean),
      highlight & direction == "Depleted" ~
        sprintf("<span style='color:#377EB8;'><i>%s</i></span>", feature_clean),
      highlight & direction == "No difference" ~
        sprintf("<span style='color:grey40;'><i>%s</i></span>", feature_clean),
      TRUE ~ sprintf("<i>%s</i>", feature_clean)
    )
  ) %>%
  select(feature, label) %>%
  { setNames(.$label, .$feature) }

## --------------------- RF Plot (stool) --------------------- ##

p1 <- ggplot(
  df_plot_rf,
  aes(x = cohort, y = feature)
) +
  geom_point(
    aes(
      size  = Gini_prev2.5,
      color = direction
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    name  = "Prev-integrated RF importance",
    range = c(2, 6)
  ) +
  scale_color_manual(
    values = c(
      "Enriched" = "#E41A1C",
      "Depleted" = "#377EB8",
      "No difference" = "grey70"
    ),
    name = "Direction"
  ) +
  scale_x_discrete(labels = function(x) stool_x_labels[x]) +
  scale_y_discrete(labels = function(x) stool_rf_label_tbl[x]) +
  theme_classic() +
  labs(
    x = "Cohort (C = controls, PC = pancreatic cancer)",
    y = "Top 50 important species",
    title = "Prevalence-integrated RF – top 50 species across PDAC stool cohorts"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

## --------------------- ROC for original RF (stool) --------------------- ##
# (unchanged; prevalence affects only feature importance, not the classifier)

rf_roc_one_orig <- function(X, y, cohort_name) {
  y <- factor(y, levels = c("Control", "Pancreatic Cancer"))
  
  rf_model <- randomForest(
    x          = X,
    y          = y,
    ntree      = 1000,
    importance = FALSE
  )
  
  prob_pc <- predict(rf_model, type = "prob")[, "Pancreatic Cancer"]
  
  roc_obj <- pROC::roc(
    response  = y,
    predictor = prob_pc,
    levels    = c("Control", "Pancreatic Cancer"),
    direction = "<"
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  coords_df <- pROC::coords(
    roc_obj,
    x   = "all",
    ret = c("specificity", "sensitivity"),
    transpose = FALSE
  )
  
  df <- as.data.frame(coords_df)
  df$fpr    <- 1 - df$specificity
  df$tpr    <- df$sensitivity
  df$cohort <- cohort_name
  df$AUC    <- auc_val
  
  df[, c("fpr", "tpr", "cohort", "AUC")]
}

X_all_stool <- species_stool
y_all_stool <- stool_data$Group_new

roc_stool_all <- rf_roc_one_orig(X_all_stool, y_all_stool, "All_stool")

stool_cohort_levels <- sort(unique(stool_data$cohort))

roc_stool_by_cohort <- purrr::map_dfr(stool_cohort_levels, function(cc) {
  idx <- stool_data$cohort == cc
  X   <- X_all_stool[idx, , drop = FALSE]
  y   <- y_all_stool[idx]
  
  rf_roc_one_orig(X, y, cohort_name = cc)
})

roc_stool_df <- dplyr::bind_rows(roc_stool_all, roc_stool_by_cohort)

stool_auc_labels <- roc_stool_df %>%
  dplyr::distinct(cohort, AUC) %>%
  dplyr::mutate(
    cohort = factor(cohort, levels = cohort_order)
  ) %>%
  dplyr::arrange(cohort) %>%
  dplyr::mutate(
    cohort_label = sprintf("%s (AUC = %.2f)", as.character(cohort), AUC),
    cohort_label = factor(cohort_label, levels = cohort_label)
  )

roc_stool_df <- roc_stool_df %>%
  dplyr::left_join(stool_auc_labels, by = c("cohort", "AUC"))

p2 <- ggplot(roc_stool_df,
             aes(x = fpr, y = tpr, color = cohort_label)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Random Forest ROC – PDAC vs Control (stool)",
    x = "Specificity",
    y = "Sensitivity",
    color = "Cohort"
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

## ===================== SALIVA SECTION ===================== ##

pdac    <- read_csv("data/PDAC.csv")
species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

set.seed(123)

saliva_data <- pdac %>%
  filter(Type == "Saliva")

cohort_var <- "Study"

saliva_data <- saliva_data %>%
  mutate(
    cohort    = .data[[cohort_var]],
    Group_new = factor(Group_new, levels = c("Control", "Pancreatic Cancer"))
  )

## species matrix
species_matrix <- species %>%
  column_to_rownames(var = "ID")
species <- as.data.frame(t(species_matrix))
species <- species[saliva_data$ID, ]
species <- species / 100

abundance_threshold  <- 2e-4   # 0.0002
prevalence_threshold <- 0.10   # 10%

taxa <- lapply(split(saliva_data$ID, saliva_data$Data), function(sample_ids) {
  species_sub <- species[sample_ids, , drop = FALSE]
  apply(species_sub, 2, function(x) {
    mean(x > abundance_threshold) >= prevalence_threshold
  })
})

taxa <- Reduce("|", taxa)
species_saliva <- species[, taxa]

stopifnot(all(rownames(species_saliva) == saliva_data$ID))

meta_saliva <- saliva_data %>%
  column_to_rownames("ID")

## ---------- Prevalence (saliva) ---------- ##
prev_saliva <- colMeans(species_saliva > 0)

## --------------------- RF (prevalence-integrated) --------------------- ##

X_all <- species_saliva
y_all <- saliva_data$Group_new

rf_all <- run_rf_one_prev(X_all, y_all, "All_saliva", prev_vec = prev_saliva, top_n = 50)

cohort_levels <- sort(unique(saliva_data$cohort))

valid_cohorts <- saliva_data %>%
  group_by(cohort) %>%
  filter(n_distinct(Group_new) >= 2) %>%
  pull(cohort) %>%
  unique() %>%
  sort()

rf_by_cohort <- map_dfr(valid_cohorts, function(cc) {
  idx <- saliva_data$cohort == cc
  X   <- X_all[idx, , drop = FALSE]
  y   <- y_all[idx]
  
  run_rf_one_prev(X, y, cohort_name = cc, prev_vec = prev_saliva, top_n = 50)
})

rf_top50 <- bind_rows(rf_all, rf_by_cohort)

rf_top50 <- rf_top50 %>%
  filter(
    (!grepl("unknown", feature, ignore.case = TRUE) |
       grepl("SGB", feature, ignore.case = TRUE)) &
      !grepl("GGB", feature, ignore.case = TRUE)
  )

feature_overlap <- rf_top50 %>%
  group_by(feature) %>%
  summarise(n_cohort = n_distinct(cohort), .groups = "drop")

rf_top50 <- rf_top50 %>%
  left_join(feature_overlap, by = "feature") %>%
  mutate(
    overlap_group = ifelse(
      n_cohort >= 3,
      "≥3 cohorts",
      "<3 cohorts"
    )
  )

## --------------------- Cohort labels --------------------- ##

cohort_counts <- saliva_data %>%
  group_by(cohort, Group_new) %>%
  tally(name = "n") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    label = sprintf("%s\nC=%s, PC=%s",
                    cohort,
                    Control,
                    `Pancreatic Cancer`)
  )

all_counts <- saliva_data %>%
  count(Group_new) %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    cohort = "All_saliva",
    label  = sprintf("All_saliva\nC=%s, PC=%s",
                     Control,
                     `Pancreatic Cancer`)
  )

cohort_labels_tbl <- bind_rows(
  all_counts %>% select(cohort, label),
  cohort_counts %>% select(cohort, label)
)

saliva_x_labels <- setNames(cohort_labels_tbl$label, cohort_labels_tbl$cohort)

cohort_levels_rf <- sort(unique(rf_top50$cohort[rf_top50$cohort != "All_saliva"]))

cohort_order <- c("All_saliva", cohort_levels_rf)

## --------------------- Mean differences & direction --------------------- ##

saliva_long <- saliva_data %>%
  select(cohort, Group_new) %>%
  bind_cols(as.data.frame(species_saliva)) %>%
  pivot_longer(
    cols      = -c(cohort, Group_new),
    names_to  = "feature",
    values_to = "abundance"
  )

mean_diff_tbl <- saliva_long %>%
  group_by(cohort, feature, Group_new) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = mean_abundance
  ) %>%
  mutate(
    mean_diff = `Pancreatic Cancer` - Control,
    direction = case_when(
      mean_diff >  0  ~ "Enriched",
      mean_diff <  0  ~ "Depleted",
      TRUE            ~ "No difference"
    )
  ) %>%
  select(cohort, feature, mean_diff, direction)

mean_diff_all <- saliva_data %>%
  select(Group_new) %>%
  bind_cols(as.data.frame(species_saliva)) %>%
  pivot_longer(
    cols      = -Group_new,
    names_to  = "feature",
    values_to = "abundance"
  ) %>%
  group_by(feature, Group_new) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = Group_new,
    values_from = mean_abundance
  ) %>%
  mutate(
    cohort    = "All_saliva",
    mean_diff = `Pancreatic Cancer` - Control,
    direction = case_when(
      mean_diff >  0  ~ "Enriched",
      mean_diff <  0  ~ "Depleted",
      TRUE            ~ "No difference"
    )
  ) %>%
  select(cohort, feature, mean_diff, direction)

mean_diff_tbl2 <- bind_rows(mean_diff_tbl, mean_diff_all)

rf_top50 <- rf_top50 %>%
  left_join(mean_diff_tbl2, by = c("cohort", "feature"))

dir_consistency <- rf_top50 %>%
  group_by(feature) %>%
  summarise(
    n_dir_cohort = n_distinct(cohort[direction != "No difference"]),
    n_pos        = sum(direction == "Enriched", na.rm = TRUE),
    n_neg        = sum(direction == "Depleted", na.rm = TRUE),
    n_dir_types  = (n_pos > 0) + (n_neg > 0),
    consistent_ge4 = (n_dir_cohort >= 3 & n_dir_types == 1),
    .groups = "drop"
  )

## --------------------- Order of features --------------------- ##

top50_pooled_order <- rf_all %>%
  arrange(desc(Gini_prev2.5)) %>%
  pull(feature) %>%
  unique()

df_plot_rf <- rf_top50 %>%
  filter(feature %in% top50_pooled_order) %>%
  mutate(
    cohort  = factor(cohort,  levels = cohort_order),
    feature = factor(feature, levels = rev(top50_pooled_order)),
    overlap_group = factor(overlap_group,
                           levels = c("<3 cohorts", "≥3 cohorts")),
    direction = factor(
      direction,
      levels = c("Enriched", "Depleted", "No difference")
    )
  ) %>%
  filter(
    !is.na(Gini_prev2.5),
    !is.na(feature),
    !is.na(cohort)
  )

feature_overlap_plot <- df_plot_rf %>%
  distinct(cohort, feature) %>% 
  group_by(feature) %>%
  summarise(n_cohort = n_distinct(cohort), .groups = "drop")

## --------------------- Y-axis labels (clean + colored) --------------------- ##

dir_all <- mean_diff_tbl2 %>%
  filter(cohort == "All_saliva") %>%
  select(feature, direction)

saliva_rf_label_tbl <- feature_overlap_plot %>%
  left_join(dir_all, by = "feature") %>%
  left_join(dir_consistency %>% select(feature, consistent_ge4),
            by = "feature") %>%
  mutate(
    feature_clean = feature %>%
      gsub("^s__", "", .) %>%
      gsub("^s_",  "", .) %>%
      gsub("^s ",  "", .) %>%
      gsub("_", " ", .),
    
    highlight = (n_cohort >= 3 & !is.na(consistent_ge4) & consistent_ge4),
    
    label = case_when(
      highlight & direction == "Enriched" ~
        sprintf("<span style='color:#E41A1C;'><i>%s</i></span>", feature_clean),
      highlight & direction == "Depleted" ~
        sprintf("<span style='color:#377EB8;'><i>%s</i></span>", feature_clean),
      highlight & direction == "No difference" ~
        sprintf("<span style='color:grey40;'><i>%s</i></span>", feature_clean),
      TRUE ~ sprintf("<i>%s</i>", feature_clean)
    )
  ) %>%
  select(feature, label) %>%
  { setNames(.$label, .$feature) }

## --------------------- RF Plot (saliva) --------------------- ##

p3 <- ggplot(
  df_plot_rf,
  aes(x = cohort, y = feature)
) +
  geom_point(
    aes(
      size  = Gini_prev2.5,
      color = direction
    ),
    alpha = 0.9
  ) +
  scale_size_continuous(
    name  = "Prev-integrated RF importance",
    range = c(2, 6)
  ) +
  scale_color_manual(
    values = c(
      "Enriched" = "#E41A1C",
      "Depleted" = "#377EB8",
      "No difference" = "grey70"
    ),
    name = "Direction"
  ) +
  scale_x_discrete(labels = function(x) saliva_x_labels[x]) +
  scale_y_discrete(labels = function(x) saliva_rf_label_tbl[x]) +
  theme_classic() +
  labs(
    x = "Cohort (C = controls, PC = pancreatic cancer)",
    y = "Top 50 important species",
    title = "Prevalence-integrated RF – top 50 species across PDAC saliva cohorts"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

## --------------------- ROC for original RF (saliva) --------------------- ##

X_all_saliva <- species_saliva
y_all_saliva <- saliva_data$Group_new

roc_saliva_all <- rf_roc_one_orig(X_all_saliva, y_all_saliva, "All_saliva")

saliva_cohort_levels <- sort(unique(saliva_data$cohort))

roc_saliva_by_cohort <- purrr::map_dfr(saliva_cohort_levels, function(cc) {
  idx <- saliva_data$cohort == cc
  X   <- X_all_saliva[idx, , drop = FALSE]
  y   <- y_all_saliva[idx]
  
  if (dplyr::n_distinct(y) < 2L) {
    message("Skipping cohort ", cc,
            " (only one class: ", paste(unique(y), collapse = ", "), ")")
    return(NULL)
  }
  
  rf_roc_one_orig(X, y, cohort_name = cc)
})

roc_saliva_df <- dplyr::bind_rows(roc_saliva_all, roc_saliva_by_cohort)

saliva_auc_labels <- roc_saliva_df %>%
  dplyr::distinct(cohort, AUC) %>%
  dplyr::mutate(
    cohort = factor(cohort, levels = cohort_order)
  ) %>%
  dplyr::arrange(cohort) %>%
  dplyr::mutate(
    cohort_label = sprintf("%s (AUC = %.2f)", as.character(cohort), AUC),
    cohort_label = factor(cohort_label, levels = cohort_label)
  )

roc_saliva_df <- roc_saliva_df %>%
  dplyr::left_join(saliva_auc_labels, by = c("cohort", "AUC"))

p4 <- ggplot(roc_saliva_df,
             aes(x = fpr, y = tpr, color = cohort_label)) +
  geom_line(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_equal() +
  theme_classic() +
  labs(
    title = "Random Forest ROC – PDAC vs Control (saliva)",
    x = "Specificity",
    y = "Sensitivity",
    color = "Cohort"
  ) +
  theme(
    plot.title      = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

## --------------------- Combined plot --------------------- ##

p <- plot_grid(
  p2, p1, p4, p3, 
  labels = c("A", "B", "C", "D"),
  label_size = 20,
  ncol = 2,
  align = "hv"
)
print(p)

ggsave("fig/rf_prevalence.pdf", p, width = 20, height = 10)
