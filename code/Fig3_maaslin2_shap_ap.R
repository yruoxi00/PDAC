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
library(ranger)
library(fastshap)
library(stringr)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/AP.csv")

## ===================== MAASLIN2 HEATMAP ===================== ##

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
species_pdac <- filter_group_species("AP", stool_data, species_raw, abundance_threshold, prevalence_threshold)

species_union <- union(colnames(species_control), colnames(species_pdac))
species_use <- species_raw[, species_union, drop = FALSE]

coerce_maaslin_types <- function(df) {
  df %>%
    mutate(
      coef = readr::parse_number(as.character(coef)),
      qval = readr::parse_number(as.character(qval))
    )
}

res_all <- read_tsv("result/Maaslin2_stool_ap/significant_results.tsv") %>% coerce_maaslin_types()
res1    <- read_tsv("result/Maaslin2_stool_dat1_ap/significant_results.tsv") %>% coerce_maaslin_types()
res3    <- read_tsv("result/Maaslin2_stool_dat3_ap/significant_results.tsv") %>% coerce_maaslin_types()
res4    <- read_tsv("result/Maaslin2_stool_dat4_ap/significant_results.tsv") %>% coerce_maaslin_types()
res5    <- read_tsv("result/Maaslin2_stool_dat5_ap/significant_results.tsv") %>% coerce_maaslin_types()

filter_maaslin <- function(df) {
  df %>%
    filter(
      metadata == "Group_new",
      value == "AP"
    ) %>%
    transmute(
      feature,
      coef,
      qval
    ) %>%
    filter(
      !str_detect(
        feature,
        regex("SGB|unknown|uncultured|unclassified", ignore_case = TRUE)
      )
    )
}

res_all_f <- filter_maaslin(res_all)
res1_f    <- filter_maaslin(res1)
res3_f    <- filter_maaslin(res3)
res4_f    <- filter_maaslin(res4)
res5_f    <- filter_maaslin(res5)

res <- list(
  all  = res_all_f,
  dat1 = res1_f,
  dat3 = res3_f,
  dat4 = res4_f,
  dat5 = res5_f
) %>% bind_rows(.id = "cohort")

study_map <- c(
  "all"  = "All Cohorts",
  "dat1" = "Liu_2024",
  "dat3" = "Yu_2021",
  "dat4" = "Song_2025",
  "dat5" = "Kartal-1_2022"
)

res$study <- study_map[res$cohort]

heat_df <- res %>%
  select(feature, study, coef) %>%
  group_by(feature, study) %>%
  summarise(coef = mean(coef, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = study, values_from = coef) %>%
  mutate(species_clean = clean_species(feature))

cohort_order_pretty <- c(
  "Kartal-1_2022",
  "Song_2025",
  "Yu_2021",
  "Liu_2024",
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
idx  <- match(stool_data$ID, rownames(X_all))
keep <- !is.na(idx)

stool_data_rf <- stool_data[keep, , drop = FALSE]
X_all <- X_all[idx[keep], , drop = FALSE]
y_all <- factor(stool_data_rf$Group_new, levels = c("Control", "AP"))

stopifnot(all(stool_data_rf$ID == rownames(X_all)))

prev <- colMeans(X_all > 0, na.rm = TRUE)
top_features <- names(sort(prev, decreasing = TRUE))[1:100]
X_all <- X_all[, top_features, drop = FALSE]

stopifnot(all(stool_data_rf$ID == rownames(X_all)))

run_rf_shap_tbl <- function(X, y, cohort_name, ntrees = 100, nsim = 5, seed = 123) {
  set.seed(seed)
  y <- factor(y, levels = c("Control", "AP"))
  dat <- data.frame(y = y, as.data.frame(X))
  rf <- ranger(
    y ~ .,
    data = dat,
    num.trees = ntrees,
    probability = TRUE,
    num.threads = parallel::detectCores()
  )
  pred_fun <- function(object, newdata) {
    p <- predict(object, data = newdata)$predictions
    as.numeric(p[, "AP"])
  }
  shap_mat <- fastshap::explain(
    object = rf,
    X = dat[, -1, drop = FALSE],
    pred_wrapper = pred_fun,
    nsim = nsim,
    adjust = TRUE
  )
  as.data.frame(shap_mat) %>%
    pivot_longer(cols = everything(), names_to = "feature", values_to = "shap") %>%
    group_by(feature) %>%
    summarise(
      mean_shap = mean(shap, na.rm = TRUE),
      mean_abs  = mean(abs(shap), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(cohort = cohort_name) %>%
    filter(!grepl("unknown|uncultured|unclassified|SGB|GGB", feature, ignore.case = TRUE))
}

cohort_levels <- sort(unique(stool_data_rf$cohort))

shap_by_cohort <- purrr::map_dfr(cohort_levels, function(cc) {
  idxc <- stool_data_rf$cohort == cc
  run_rf_shap_tbl(X_all[idxc, , drop = FALSE], y_all[idxc], cohort_name = cc)
})

shap_all <- run_rf_shap_tbl(X_all, y_all, cohort_name = "All Cohorts")

shap_imp <- bind_rows(shap_by_cohort, shap_all) %>%
  mutate(species_clean = clean_species(feature))

shap_imp2 <- shap_imp %>%
  group_by(species_clean, cohort) %>%
  summarise(
    mean_shap = mean(mean_shap, na.rm = TRUE),
    mean_abs  = mean(mean_abs,  na.rm = TRUE),
    .groups = "drop"
  )

cohort_order_rf <- c(
  "Kartal-1_2022",
  "Song_2025",
  "Yu_2021",
  "Liu_2024",
  "All Cohorts"
)

cohort_cols <- cohort_order_rf[cohort_order_rf %in% unique(shap_imp2$cohort)]

top50_names <- shap_imp2 %>%
  filter(cohort %in% cohort_cols) %>%
  group_by(species_clean) %>%
  summarise(avg_abs = mean(mean_abs, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(avg_abs)) %>%
  slice_head(n = 50) %>%
  pull(species_clean)

shap_wide_signed <- shap_imp2 %>%
  filter(cohort %in% cohort_cols, species_clean %in% top50_names) %>%
  select(species_clean, cohort, mean_shap) %>%
  pivot_wider(names_from = cohort, values_from = mean_shap, values_fill = 0) %>%
  mutate(species_clean = factor(species_clean, levels = top50_names)) %>%
  arrange(species_clean)

mat2 <- as.matrix(shap_wide_signed[, cohort_cols, drop = FALSE])
rownames(mat2) <- as.character(shap_wide_signed$species_clean)
colnames(mat2) <- cohort_cols

stool_long <- stool_data_rf %>%
  select(cohort, Group_new) %>%
  bind_cols(as.data.frame(X_all)) %>%
  pivot_longer(
    cols = -c(cohort, Group_new),
    names_to = "feature",
    values_to = "abundance"
  ) %>%
  mutate(species_clean = clean_species(feature))

diff_tbl <- stool_long %>%
  filter(species_clean %in% top50_names) %>%
  group_by(cohort, species_clean, Group_new) %>%
  summarise(m = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group_new, values_from = m, values_fill = 0) %>%
  mutate(diff = AP - Control) %>%
  select(cohort, species_clean, diff)

diff_all <- stool_data_rf %>%
  select(Group_new) %>%
  bind_cols(as.data.frame(X_all)) %>%
  pivot_longer(cols = -Group_new, names_to = "feature", values_to = "abundance") %>%
  mutate(species_clean = clean_species(feature)) %>%
  filter(species_clean %in% top50_names) %>%
  group_by(species_clean, Group_new) %>%
  summarise(m = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group_new, values_from = m, values_fill = 0) %>%
  mutate(cohort = "All Cohorts", diff = AP - Control) %>%
  select(cohort, species_clean, diff)

diff_wide <- bind_rows(diff_tbl, diff_all) %>%
  pivot_wider(names_from = cohort, values_from = diff, values_fill = 0) %>%
  mutate(species_clean = factor(species_clean, levels = top50_names)) %>%
  arrange(species_clean)

mat2_diff <- as.matrix(diff_wide[, cohort_cols, drop = FALSE])
rownames(mat2_diff) <- as.character(diff_wide$species_clean)

mat2_sign <- ifelse(mat2_diff > 0, "+", ifelse(mat2_diff < 0, "-", ""))

## ===================== BOLD OVERLAPPING SPECIES ===================== ## 

maaslin_species <- clean_species(rownames(mat1)) 
rf_species      <- clean_species(rownames(mat2)) 
common_species  <- intersect(maaslin_species, rf_species) 

label_expr_maaslin <- ifelse(
  clean_species(rownames(mat1)) %in% common_species,
  sprintf("bolditalic('%s')", rownames(mat1)),
  sprintf("italic('%s')", rownames(mat1))
)
labels_row_maaslin <- parse(text = label_expr_maaslin)

label_expr_rf <- ifelse(
  clean_species(rownames(mat2)) %in% common_species,
  sprintf("bolditalic('%s')", rownames(mat2)),
  sprintf("italic('%s')", rownames(mat2))
)
labels_row_rf <- parse(text = label_expr_rf)

lim_maaslin  <- 3
brks_maaslin <- seq(-lim_maaslin, lim_maaslin, length.out = 201)
mat1c <- pmax(pmin(mat1, lim_maaslin), -lim_maaslin)

lim_rf  <- max(mat2, na.rm = TRUE)
brks_rf <- seq(-lim_rf, lim_rf, length.out = 201)
mat2c <- pmax(pmin(mat2, lim_rf), -lim_rf)

hmcols <- colorRampPalette(c("#313695", "white", "#A50026"))(200)

CELL <- 10

p1 <- pheatmap(
  mat1c,
  cellwidth = CELL, cellheight = CELL,
  color = hmcols, breaks = brks_maaslin,
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
  color = hmcols, breaks = brks_rf,
  na_col = "white",
  cluster_rows = FALSE, cluster_cols = FALSE,
  border_color = "black",
  labels_row = labels_row_rf,
  display_numbers = mat2_sign,
  number_color = "black",
  fontsize_number = 10,
  legend_breaks = pretty(c(-lim_rf, lim_rf), n = 5),
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

ggsave("AP_Fig3_maaslin2_shap.pdf", p, width = 10, height = 9, limitsize = FALSE)
