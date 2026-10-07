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

pdac <- read_csv("data/PDAC.csv")

### ===================== STOOL DATA ===================== ###

stool_data <- pdac %>%
  filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

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
species_stool <- species[, taxa]

stopifnot(all(rownames(species_stool) == stool_data$ID))

meta_stool <- stool_data %>%
  column_to_rownames("ID")

### ===================== SALIVA DATA ===================== ###

saliva_data <- pdac %>%
  filter(Type == "Saliva")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

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
species_saliva <- species[, taxa]

stopifnot(all(rownames(species_saliva) == saliva_data$ID))

meta_saliva <- saliva_data %>%
  column_to_rownames("ID")

### ===================== MaAsLin3 RESULTS (READ ONCE) ===================== ###

res_stool_all  <- read_tsv("result/Maaslin3_stool/all_results.tsv")
res_saliva_all <- read_tsv("result/Maaslin3_saliva/all_results.tsv")

res_stool_dat2  <- read_tsv("result/Maaslin3_stool_dat2/all_results.tsv")
res_stool_dat5  <- read_tsv("result/Maaslin3_stool_dat5/all_results.tsv")
res_stool_dat6  <- read_tsv("result/Maaslin3_stool_dat6/all_results.tsv")
res_stool_dat8  <- read_tsv("result/Maaslin3_stool_dat8/all_results.tsv")
res_stool_dat9  <- read_tsv("result/Maaslin3_stool_dat9/all_results.tsv")

res_saliva_dat2 <- read_tsv("result/Maaslin3_saliva_dat2/all_results.tsv")
res_saliva_dat5 <- read_tsv("result/Maaslin3_saliva_dat5/all_results.tsv")
res_saliva_dat9 <- read_tsv("result/Maaslin3_saliva_dat9/all_results.tsv")

## split by model once
res_stool_all_abund   <- res_stool_all   %>% filter(model == "abundance")
res_stool_all_prev    <- res_stool_all   %>% filter(model == "prevalence")
res_saliva_all_abund  <- res_saliva_all  %>% filter(model == "abundance")
res_saliva_all_prev   <- res_saliva_all  %>% filter(model == "prevalence")

res_stool_dat2_abund  <- res_stool_dat2  %>% filter(model == "abundance")
res_stool_dat2_prev   <- res_stool_dat2  %>% filter(model == "prevalence")
res_stool_dat5_abund  <- res_stool_dat5  %>% filter(model == "abundance")
res_stool_dat5_prev   <- res_stool_dat5  %>% filter(model == "prevalence")
res_stool_dat6_abund  <- res_stool_dat6  %>% filter(model == "abundance")
res_stool_dat6_prev   <- res_stool_dat6  %>% filter(model == "prevalence")
res_stool_dat8_abund  <- res_stool_dat8  %>% filter(model == "abundance")
res_stool_dat8_prev   <- res_stool_dat8  %>% filter(model == "prevalence")
res_stool_dat9_abund  <- res_stool_dat9  %>% filter(model == "abundance")
res_stool_dat9_prev   <- res_stool_dat9  %>% filter(model == "prevalence")

res_saliva_dat2_abund <- res_saliva_dat2 %>% filter(model == "abundance")
res_saliva_dat2_prev  <- res_saliva_dat2 %>% filter(model == "prevalence")
res_saliva_dat5_abund <- res_saliva_dat5 %>% filter(model == "abundance")
res_saliva_dat5_prev  <- res_saliva_dat5 %>% filter(model == "prevalence")
res_saliva_dat9_abund <- res_saliva_dat9 %>% filter(model == "abundance")
res_saliva_dat9_prev  <- res_saliva_dat9 %>% filter(model == "prevalence")

### helper for both abundance & prevalence
q_cut <- 0.25

add_cohort <- function(res, cohort_label) {
  res %>%
    filter(value == "Pancreatic Cancer") %>%
    mutate(cohort = cohort_label) %>%
    filter(qval_individual < q_cut) %>%      # <<< use qval_individual
    select(cohort, feature, coef, qval_individual, model)
}

clean_species <- function(x) {
  x <- gsub("^s__", "", x)
  x <- gsub("^s_",  "", x)
  x <- gsub("^s ",  "", x)
  x <- gsub("_", " ", x)
  x
}

### ===================== ABUNDANCE: SIG TABLES ===================== ###

sig_stool_abund <- bind_rows(
  add_cohort(res_stool_all_abund,   "All"),
  add_cohort(res_stool_dat2_abund,  "Vilà-Quintana_2024"),
  add_cohort(res_stool_dat5_abund,  "Kartal-1_2022"),
  add_cohort(res_stool_dat6_abund,  "Kartal-2_2022"),
  add_cohort(res_stool_dat8_abund,  "Matsukawa_2021"),
  add_cohort(res_stool_dat9_abund,  "Nagata_2022")
) %>%
  mutate(direction = if_else(coef > 0, "Enriched", "Depleted"))

sig_saliva_abund <- bind_rows(
  add_cohort(res_saliva_all_abund,  "All"),
  add_cohort(res_saliva_dat2_abund, "Vilà-Quintana_2024"),
  add_cohort(res_saliva_dat5_abund, "Kartal-1_2022"),
  add_cohort(res_saliva_dat9_abund, "Nagata_2022")
) %>%
  mutate(direction = if_else(coef > 0, "Enriched", "Depleted"))

overlap_stool_abund <- sig_stool_abund %>%
  distinct(cohort, feature) %>%
  count(feature) %>%
  filter(n >= 3)

overlap_saliva_abund <- sig_saliva_abund %>%
  distinct(cohort, feature) %>%
  count(feature) %>%
  filter(n >= 3)

df_plot_stool_abund <- sig_stool_abund %>%
  inner_join(overlap_stool_abund, by = "feature") %>%
  mutate(
    feature = factor(feature, levels = sort(unique(feature))),
    cohort  = factor(
      cohort,
      levels = c("All", "Vilà-Quintana_2024", "Kartal-1_2022",
                 "Kartal-2_2022", "Nagata_2022", "Matsukawa_2021")
    ),
    direction = factor(direction, levels = c("Enriched", "Depleted"))
  ) %>%
  filter(!str_detect(feature, "SGB"))

df_plot_saliva_abund <- sig_saliva_abund %>%
  inner_join(overlap_saliva_abund, by = "feature") %>%
  mutate(
    feature = factor(feature, levels = sort(unique(feature))),
    cohort  = factor(
      cohort,
      levels = c("All", "Vilà-Quintana_2024", "Kartal-1_2022", "Nagata_2022")
    ),
    direction = factor(direction, levels = c("Enriched", "Depleted"))
  ) %>%
  filter(!str_detect(feature, "SGB"))

df_plot_stool_abund <- df_plot_stool_abund %>%
  mutate(feature = factor(feature, levels = sort(unique(feature))))

df_plot_saliva_abund <- df_plot_saliva_abund %>%
  mutate(feature = factor(feature, levels = sort(unique(feature))))

### LABELS (ABUNDANCE) ###

stool_label_tbl_abund <- df_plot_stool_abund %>%
  group_by(feature) %>%
  summarise(
    n_cohort      = n_distinct(cohort),
    dir_consistent = n_distinct(direction) == 1,
    dir_first      = first(direction),
    .groups = "drop"
  ) %>%
  mutate(
    feature_clean = clean_species(feature),
    label_colored = case_when(
      dir_consistent & dir_first == "Enriched" & n_cohort >= 3 ~
        sprintf("<span style='color:#E41A1C;'><i>%s</i></span>", feature_clean),
      dir_consistent & dir_first == "Depleted" & n_cohort >= 3 ~
        sprintf("<span style='color:#377EB8;'><i>%s</i></span>", feature_clean),
      TRUE ~ sprintf("<i>%s</i>", feature_clean)
    )
  ) %>%
  { setNames(.$label_colored, .$feature) }

saliva_label_tbl_abund <- df_plot_saliva_abund %>%
  group_by(feature) %>%
  summarise(
    n_cohort      = n_distinct(cohort),
    dir_consistent = n_distinct(direction) == 1,
    dir_first      = first(direction),
    .groups = "drop"
  ) %>%
  mutate(
    label_colored = case_when(
      dir_consistent & dir_first == "Enriched" & n_cohort >= 3 ~
        sprintf("<span style='color:#E41A1C;'>%s</span>", feature),
      dir_consistent & dir_first == "Depleted" & n_cohort >= 3 ~
        sprintf("<span style='color:#377EB8;'>%s</span>", feature),
      TRUE ~ as.character(feature)
    )
  ) %>%
  { setNames(.$label_colored, .$feature) }

### COHORT SIZES (COMMON) ###

cohort_sizes_stool <- tibble(
  cohort = c("All", "Vilà-Quintana_2024", "Kartal-1_2022",
             "Kartal-2_2022", "Nagata_2022", "Matsukawa_2021"),
  n_case    = c(189, 11, 67, 44, 43, 24),
  n_control = c(351, 18, 53, 32, 235, 13)
)

stool_x_labels <- cohort_sizes_stool %>%
  mutate(
    label = sprintf("%s\nCases=%s, Controls=%s",
                    cohort, n_case, n_control)
  ) %>%
  { setNames(.$label, .$cohort) }

cohort_sizes_saliva <- tibble(
  cohort = c("All", "Vilà-Quintana_2024", "Kartal-1_2022", "Nagata_2022"),
  n_case    = c(116, 11, 43, 47),
  n_control = c(304, 21, 48, 235)
)

saliva_x_labels <- cohort_sizes_saliva %>%
  mutate(
    label = sprintf("%s\nCases=%s, Controls=%s",
                    cohort, n_case, n_control)
  ) %>%
  { setNames(.$label, .$cohort) }

### RANKING + TOP 50 (STOOL, ABUNDANCE) ###

stool_rank_abund <- df_plot_stool_abund %>%
  group_by(feature) %>%
  summarise(
    mean_abs_coef = mean(abs(coef), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_coef))

top_n <- 50

top_stool_features_abund <- stool_rank_abund %>%
  slice_head(n = top_n) %>%
  pull(feature)

df_plot_stool_abund <- df_plot_stool_abund %>%
  filter(feature %in% top_stool_features_abund) %>%
  mutate(
    feature = factor(feature, levels = rev(top_stool_features_abund))
  )

### PLOTS (ABUNDANCE) ###

p1 <- ggplot(df_plot_stool_abund,
             aes(x = cohort, y = feature, color = direction)) +
  geom_point(aes(size = abs(coef)), alpha = 0.9) +
  scale_color_manual(values = c("Enriched" = "#E41A1C",
                                "Depleted" = "#377EB8")) +
  scale_size_continuous(name = "Effect Size", range = c(2, 6)) +
  scale_x_discrete(labels = function(x) stool_x_labels[x]) +
  scale_y_discrete(labels = function(x) stool_label_tbl_abund[x]) +
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Top 50 overlapped species",
    color = "Direction",
    title = "Overlapped Species Across PDAC Stool Cohorts (Abundance)"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

p2 <- ggplot(df_plot_saliva_abund,
             aes(x = cohort, y = feature, color = direction)) +
  geom_point(aes(size = abs(coef)), alpha = 0.9) +
  scale_color_manual(values = c("Enriched" = "#E41A1C",
                                "Depleted" = "#377EB8")) +
  scale_size_continuous(name = "Effect Size", range = c(2, 6)) +
  scale_x_discrete(labels = function(x) saliva_x_labels[x]) +
  scale_y_discrete(labels = function(x) saliva_label_tbl_abund[x]) +
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Overlapped species",
    color = "Direction",
    title = "Overlapped Species Across PDAC Saliva Cohorts (Abundance)"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

### ===================== PREVALENCE: SIG TABLES ===================== ###

sig_stool_prev <- bind_rows(
  add_cohort(res_stool_all_prev,   "All"),
  add_cohort(res_stool_dat2_prev,  "Vilà-Quintana_2024"),
  add_cohort(res_stool_dat5_prev,  "Kartal-1_2022"),
  add_cohort(res_stool_dat6_prev,  "Kartal-2_2022"),
  add_cohort(res_stool_dat8_prev,  "Matsukawa_2021"),
  add_cohort(res_stool_dat9_prev,  "Nagata_2022")
) %>%
  mutate(direction = if_else(coef > 0, "Enriched", "Depleted"))

sig_saliva_prev <- bind_rows(
  add_cohort(res_saliva_all_prev,  "All"),
  add_cohort(res_saliva_dat2_prev, "Vilà-Quintana_2024"),
  add_cohort(res_saliva_dat5_prev, "Kartal-1_2022"),
  add_cohort(res_saliva_dat9_prev, "Nagata_2022")
) %>%
  mutate(direction = if_else(coef > 0, "Enriched", "Depleted"))

overlap_stool_prev <- sig_stool_prev %>%
  distinct(cohort, feature) %>%
  count(feature) %>%
  filter(n >= 3)

overlap_saliva_prev <- sig_saliva_prev %>%
  distinct(cohort, feature) %>%
  count(feature) %>%
  filter(n >= 3)

df_plot_stool_prev <- sig_stool_prev %>%
  inner_join(overlap_stool_prev, by = "feature") %>%
  mutate(
    feature = factor(feature, levels = sort(unique(feature))),
    cohort  = factor(
      cohort,
      levels = c("All", "Vilà-Quintana_2024", "Kartal-1_2022",
                 "Kartal-2_2022", "Nagata_2022", "Matsukawa_2021")
    ),
    direction = factor(direction, levels = c("Enriched", "Depleted"))
  ) %>%
  filter(!str_detect(feature, "SGB"))

df_plot_saliva_prev <- sig_saliva_prev %>%
  inner_join(overlap_saliva_prev, by = "feature") %>%
  mutate(
    feature = factor(feature, levels = sort(unique(feature))),
    cohort  = factor(
      cohort,
      levels = c("All", "Vilà-Quintana_2024", "Kartal-1_2022", "Nagata_2022")
    ),
    direction = factor(direction, levels = c("Enriched", "Depleted"))
  ) %>%
  filter(!str_detect(feature, "SGB"))

df_plot_stool_prev <- df_plot_stool_prev %>%
  mutate(feature = factor(feature, levels = sort(unique(feature))))

df_plot_saliva_prev <- df_plot_saliva_prev %>%
  mutate(feature = factor(feature, levels = sort(unique(feature))))

### LABELS (PREVALENCE) ###

stool_label_tbl_prev <- df_plot_stool_prev %>%
  group_by(feature) %>%
  summarise(
    n_cohort      = n_distinct(cohort),
    dir_consistent = n_distinct(direction) == 1,
    dir_first      = first(direction),
    .groups = "drop"
  ) %>%
  mutate(
    feature_clean = clean_species(feature),
    label_colored = case_when(
      dir_consistent & dir_first == "Enriched" & n_cohort >= 3 ~
        sprintf("<span style='color:#E41A1C;'><i>%s</i></span>", feature_clean),
      dir_consistent & dir_first == "Depleted" & n_cohort >= 3 ~
        sprintf("<span style='color:#377EB8;'><i>%s</i></span>", feature_clean),
      TRUE ~ sprintf("<i>%s</i>", feature_clean)
    )
  ) %>%
  { setNames(.$label_colored, .$feature) }

saliva_label_tbl_prev <- df_plot_saliva_prev %>%
  group_by(feature) %>%
  summarise(
    n_cohort      = n_distinct(cohort),
    dir_consistent = n_distinct(direction) == 1,
    dir_first      = first(direction),
    .groups = "drop"
  ) %>%
  mutate(
    label_colored = case_when(
      dir_consistent & dir_first == "Enriched" & n_cohort >= 3 ~
        sprintf("<span style='color:#E41A1C;'>%s</span>", feature),
      dir_consistent & dir_first == "Depleted" & n_cohort >= 3 ~
        sprintf("<span style='color:#377EB8;'>%s</span>", feature),
      TRUE ~ as.character(feature)
    )
  ) %>%
  { setNames(.$label_colored, .$feature) }

### RANKING + TOP 50 (STOOL, PREVALENCE) ###

stool_rank_prev <- df_plot_stool_prev %>%
  group_by(feature) %>%
  summarise(
    mean_abs_coef = mean(abs(coef), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_coef))

top_stool_features_prev <- stool_rank_prev %>%
  slice_head(n = top_n) %>%
  pull(feature)

df_plot_stool_prev <- df_plot_stool_prev %>%
  filter(feature %in% top_stool_features_prev) %>%
  mutate(
    feature = factor(feature, levels = rev(top_stool_features_prev))
  )

### PLOTS (PREVALENCE) ###

p3 <- ggplot(df_plot_stool_prev,
             aes(x = cohort, y = feature, color = direction)) +
  geom_point(aes(size = abs(coef)), alpha = 0.9) +
  scale_color_manual(values = c("Enriched" = "#E41A1C",
                                "Depleted" = "#377EB8")) +
  scale_size_continuous(name = "Effect Size", range = c(2, 6)) +
  scale_x_discrete(labels = function(x) stool_x_labels[x]) +
  scale_y_discrete(labels = function(x) stool_label_tbl_prev[x]) +
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Top 50 overlapped species",
    color = "Direction",
    title = "Overlapped Species Across PDAC Stool Cohorts (Prevalence)"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

p4 <- ggplot(df_plot_saliva_prev,
             aes(x = cohort, y = feature, color = direction)) +
  geom_point(aes(size = abs(coef)), alpha = 0.9) +
  scale_color_manual(values = c("Enriched" = "#E41A1C",
                                "Depleted" = "#377EB8")) +
  scale_size_continuous(name = "Effect Size", range = c(2, 6)) +
  scale_x_discrete(labels = function(x) saliva_x_labels[x]) +
  scale_y_discrete(labels = function(x) saliva_label_tbl_prev[x]) +
  theme_classic() +
  labs(
    x = "Cohort",
    y = "Overlapped species",
    color = "Direction",
    title = "Overlapped Species Across PDAC Saliva Cohorts (Prevalence)"
  ) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_markdown(size = 7),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

### COMBINED FIGURE ###

p <- plot_grid(
  p1, p2, p3, p4,
  labels = c("A", "B", "C", "D"),
  label_size = 20,
  ncol = 2,
  align = "hv"
)

print(p)

ggsave("fig/overlapped_species.pdf", p, width = 20, height = 10)
