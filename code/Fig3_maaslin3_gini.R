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
library(pheatmap)
library(RColorBrewer)
library(gridExtra)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/PDAC.csv")

## ===================== MAASLIN3 HEATMAP ===================== ##

### stool  

stool_data <- pdac %>% filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")
taxonomy <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "metaphlan4")

clean_species <- function(x) {
  x %>%
    str_replace("^s__", "") %>%
    str_replace("^s_",  "") %>%
    str_replace("^s ",  "") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

species_matrix <- species %>% column_to_rownames(var = "ID")
species_raw <- as.data.frame(t(species_matrix))
species_raw <- species_raw[stool_data$ID, , drop = FALSE]
species_raw <- species_raw / 100

abundance_threshold  <- 2e-4
prevalence_threshold <- 0.10

filter_group_species <- function(group_name, stool_df, species_mat,
                                 abundance_threshold = 2e-4,
                                 prevalence_threshold = 0.10) {
  group_df <- stool_df %>% filter(Group_new == group_name)
  species_group <- species_mat[group_df$ID, , drop = FALSE]
  taxa_list <- lapply(split(group_df$ID, group_df$Data), function(sample_ids) {
    species_sub <- species_group[sample_ids, , drop = FALSE]
    apply(species_sub, 2, function(x) mean(x > abundance_threshold) >= prevalence_threshold)
  })
  taxa_keep <- Reduce("|", taxa_list)
  species_group[, taxa_keep, drop = FALSE]
}

species_control <- filter_group_species("Control", stool_data, species_raw, abundance_threshold, prevalence_threshold)
species_pdac <- filter_group_species("PDAC", stool_data, species_raw, abundance_threshold, prevalence_threshold)

species_union <- union(colnames(species_control), colnames(species_pdac))
species_use <- species_raw[, species_union, drop = FALSE]

coerce_maaslin_types <- function(df) {
  df %>%
    mutate(
      coef       = readr::parse_number(as.character(coef)),
      qval_joint = readr::parse_number(as.character(qval_joint))
    )
}

res_all <- read_tsv("result/Maaslin3_stool/significant_results.tsv") %>% coerce_maaslin_types()
res2    <- read_tsv("result/Maaslin3_stool_dat2/significant_results.tsv") %>% coerce_maaslin_types()
res5    <- read_tsv("result/Maaslin3_stool_dat5/significant_results.tsv") %>% coerce_maaslin_types()
res6    <- read_tsv("result/Maaslin3_stool_dat6/significant_results.tsv") %>% coerce_maaslin_types()
res8    <- read_tsv("result/Maaslin3_stool_dat8/significant_results.tsv") %>% coerce_maaslin_types()
res9    <- read_tsv("result/Maaslin3_stool_dat9/significant_results.tsv") %>% coerce_maaslin_types()

filter_maaslin <- function(df) {
  df %>%
    filter(model %in% c("abundance", "prevalence"), metadata == "Group_new", value == "PDAC") %>%
    group_by(feature) %>%
    summarise(
      n_models  = n(),
      n_pos     = sum(coef > 0, na.rm = TRUE),
      n_neg     = sum(coef < 0, na.rm = TRUE),
      consistent = (n_pos == n_models | n_neg == n_models),
      mean_coef  = ifelse(consistent, mean(coef, na.rm = TRUE), NA_real_),
      qval_joint = min(qval_joint, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    filter(consistent, !is.na(mean_coef)) %>%
    filter(!str_detect(feature, regex("SGB|unknown|uncultured|unclassified", ignore_case = TRUE)))
}

res_all_f <- filter_maaslin(res_all)
res2_f    <- filter_maaslin(res2)
res5_f    <- filter_maaslin(res5)
res6_f    <- filter_maaslin(res6)
res8_f    <- filter_maaslin(res8)
res9_f    <- filter_maaslin(res9)

res <- list(
  all  = res_all_f,
  dat2 = res2_f,
  dat5 = res5_f,
  dat6 = res6_f,
  dat8 = res8_f,
  dat9 = res9_f
) %>% bind_rows(.id = "cohort")

study_map <- c(
  "all"  = "All Cohorts",
  "dat2" = "Vilà-Quintana_2024",
  "dat5" = "Kartal-1_2022",
  "dat6" = "Kartal-2_2022",
  "dat8" = "Matsukawa_2021",
  "dat9" = "Nagata_2022"
)

res$study <- study_map[res$cohort]

heat_df <- res %>%
  select(feature, study, mean_coef) %>%
  group_by(feature, study) %>%
  summarise(mean_coef = mean(mean_coef, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = study, values_from = mean_coef) %>%
  mutate(species_clean = clean_species(feature))

cohort_order_pretty <- c(
  "Nagata_2022",
  "Matsukawa_2021",
  "Kartal-2_2022",
  "Kartal-1_2022",
  "Vilà-Quintana_2024",
  "All Cohorts"
)

cohort_cols_maaslin <- cohort_order_pretty[cohort_order_pretty %in% colnames(heat_df)]

top50_maaslin <- heat_df %>%
  mutate(avg_abs = rowMeans(abs(across(all_of(cohort_cols_maaslin))), na.rm = TRUE)) %>%
  arrange(desc(avg_abs)) %>%
  slice_head(n = 50)

mat1 <- as.matrix(top50_maaslin[, cohort_cols_maaslin, drop = FALSE])
rownames(mat1) <- top50_maaslin$species_clean
mat1[is.na(mat1)] <- 0

## ===================== RANDOM FOREST HEATMAP ===================== ##

set.seed(123)

stool_data <- stool_data %>% mutate(cohort = Study)

X_all <- species_use
idx <- match(stool_data$ID, rownames(X_all))
keep <- !is.na(idx)

stool_data_rf <- stool_data[keep, , drop = FALSE]
X_all <- X_all[idx[keep], , drop = FALSE]
y_all <- factor(stool_data_rf$Group_new, levels = c("Control", "PDAC"))

stopifnot(all(stool_data_rf$ID == rownames(X_all)))

run_rf_gini_tbl <- function(X, y, cohort_name) {
  rf_model <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
  imp <- as.data.frame(importance(rf_model))
  imp$feature <- rownames(imp)
  imp$cohort  <- cohort_name
  imp %>%
    transmute(
      cohort,
      feature,
      MeanDecreaseGini = as.numeric(MeanDecreaseGini)
    ) %>%
    filter(!grepl("unknown|uncultured|unclassified|SGB|GGB", feature, ignore.case = TRUE))
}

cohort_levels <- sort(unique(stool_data_rf$cohort))

rf_by_cohort <- purrr::map_dfr(cohort_levels, function(cc) {
  idxc <- stool_data_rf$cohort == cc
  run_rf_gini_tbl(X_all[idxc, , drop = FALSE], y_all[idxc], cohort_name = cc)
})

rf_all <- run_rf_gini_tbl(X_all, y_all, cohort_name = "All Cohorts")

rf_imp <- bind_rows(rf_by_cohort, rf_all) %>%
  mutate(species_clean = clean_species(feature))

rf_imp2 <- rf_imp %>%
  group_by(species_clean, cohort) %>%
  summarise(MeanDecreaseGini = mean(MeanDecreaseGini, na.rm = TRUE), .groups = "drop")

cohort_order_rf <- c(
  "Nagata_2022",
  "Matsukawa_2021",
  "Kartal-2_2022",
  "Kartal-1_2022",
  "Vilà-Quintana_2024",
  "All Cohorts"
)

cohort_cols <- cohort_order_rf[cohort_order_rf %in% unique(rf_imp2$cohort)]

gini_abs_rank <- rf_imp2 %>%
  filter(cohort %in% cohort_cols) %>%
  group_by(species_clean) %>%
  summarise(avg_abs_gini = mean(abs(MeanDecreaseGini), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_abs_gini))

rf_top50_names <- gini_abs_rank %>%
  slice_head(n = 50) %>%
  pull(species_clean)

rf_wide_raw <- rf_imp2 %>%
  filter(cohort %in% cohort_cols, species_clean %in% rf_top50_names) %>%
  select(species_clean, cohort, MeanDecreaseGini) %>%
  pivot_wider(names_from = cohort, values_from = MeanDecreaseGini, values_fill = 0) %>%
  mutate(species_clean = factor(species_clean, levels = rf_top50_names)) %>%
  arrange(species_clean)

mat2 <- as.matrix(rf_wide_raw[, cohort_cols, drop = FALSE])
rownames(mat2) <- as.character(rf_wide_raw$species_clean)
colnames(mat2) <- cohort_cols

stool_long <- stool_data_rf %>%
  select(cohort, Group_new) %>%
  bind_cols(as.data.frame(X_all)) %>%
  pivot_longer(cols = -c(cohort, Group_new), names_to = "feature", values_to = "abundance") %>%
  mutate(species_clean = clean_species(feature))

diff_tbl <- stool_long %>%
  filter(species_clean %in% rf_top50_names) %>%
  group_by(cohort, species_clean, Group_new) %>%
  summarise(m = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group_new, values_from = m, values_fill = 0) %>%
  mutate(diff = PDAC - Control) %>%
  select(cohort, species_clean, diff)

diff_all <- stool_data_rf %>%
  select(Group_new) %>%
  bind_cols(as.data.frame(X_all)) %>%
  pivot_longer(cols = -Group_new, names_to = "feature", values_to = "abundance") %>%
  mutate(species_clean = clean_species(feature)) %>%
  filter(species_clean %in% rf_top50_names) %>%
  group_by(species_clean, Group_new) %>%
  summarise(m = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group_new, values_from = m, values_fill = 0) %>%
  mutate(cohort = "All Cohorts", diff = PDAC - Control) %>%
  select(cohort, species_clean, diff)

diff_wide <- bind_rows(diff_tbl, diff_all) %>%
  pivot_wider(names_from = cohort, values_from = diff, values_fill = 0) %>%
  mutate(species_clean = factor(species_clean, levels = rf_top50_names)) %>%
  arrange(species_clean)

mat2_diff <- as.matrix(diff_wide[, cohort_cols, drop = FALSE])
rownames(mat2_diff) <- as.character(diff_wide$species_clean)

mat2_sign <- ifelse(mat2_diff > 0, "+", ifelse(mat2_diff < 0, "-", ""))

## ===================== BOLD OVERLAPPING SPECIES ===================== ## 

maaslin_species <- clean_species(rownames(mat1)) 
rf_species <- clean_species(rownames(mat2)) 
common_species <- intersect(maaslin_species, rf_species) 
label_expr_maaslin <- ifelse( clean_species(rownames(mat1)) %in% 
                                common_species, sprintf("bolditalic('%s')", 
                                                        rownames(mat1)), 
                              sprintf("italic('%s')", rownames(mat1)) ) 
labels_row_maaslin <- parse(text = label_expr_maaslin) 
label_expr_rf <- ifelse( clean_species(rownames(mat2)) %in% 
                           common_species, sprintf("bolditalic('%s')", 
                                                   rownames(mat2)), 
                         sprintf("italic('%s')", rownames(mat2)) ) 
labels_row_rf <- parse(text = label_expr_rf) 

lim_maaslin  <- 3
brks_maaslin <- seq(-lim_maaslin, lim_maaslin, length.out = 201)
mat1c <- pmax(pmin(mat1, lim_maaslin), -lim_maaslin)

lim_rf  <- 3
brks_rf <- seq(0, lim_rf, length.out = 201)
mat2c   <- pmax(pmin(mat2, lim_rf), -lim_rf)

hmcols <- colorRampPalette(c("#313695", "white", "#A50026"))(200)
hmcols_rf <- colorRampPalette(c("white", "#A50026"))(200)

CELL <- 10

p1 <- pheatmap(
  mat1c,
  cellwidth = CELL, cellheight = CELL,
  color = hmcols,
  breaks = brks_maaslin,
  na_col = "white",
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = "black",
  labels_row = labels_row_maaslin,
  legend_breaks = seq(-lim_maaslin, lim_maaslin, by = 1),
  silent = TRUE
)

p2 <- pheatmap(
  mat2c,
  cellwidth = CELL, cellheight = CELL,
  color = hmcols_rf,
  breaks = brks_rf,
  na_col = "white",
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = "black",
  labels_row = labels_row_rf,
  display_numbers = mat2_sign,
  number_color = "black",
  fontsize_number = 10,
  legend_breaks = pretty(c(0, lim_rf), n = 5),
  silent = TRUE
)

g1 <- p1$gtable
g2 <- p2$gtable

gg1 <- ggplotify::as.ggplot(g1)
gg2 <- ggplotify::as.ggplot(g2)

p <- plot_grid(
  gg1, gg2,
  labels = c("a", "b"),
  label_size = 20,
  ncol = 2,
  align = "hv",
  axis = "tb"
)
print(p)

ggsave("PDAC_Fig3_maaslin3_gini.pdf", p, width = 10, height = 9, limitsize = FALSE)
