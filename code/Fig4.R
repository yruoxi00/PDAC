library(readxl)
library(dplyr)
library(tidyverse)
library(tibble)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gghalves)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/AP.csv")
stool_data  <- pdac %>% filter(Type == "Stool")
saliva_data <- pdac %>% filter(Type == "Saliva")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

species_matrix <- species %>% column_to_rownames(var = "ID")
species_raw <- as.data.frame(t(species_matrix)) / 100

clean_species_name <- function(x) {
  x %>%
    str_replace("^s__", "") %>%
    str_replace("^s_",  "") %>%
    str_replace("^s\\s+", "") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

is_unknown_sgb <- function(x) {
  str_detect(x, regex("SGB|GGB|unknown|unclassified|uncultured", ignore_case = TRUE)) |
    x == "" | is.na(x)
}

new_names <- clean_species_name(colnames(species_raw))

if (any(duplicated(new_names))) {
  colnames(species_raw) <- new_names
  species_raw <- species_raw %>%
    as_tibble(rownames = "ID") %>%
    pivot_longer(-ID, names_to = "species", values_to = "abund") %>%
    group_by(ID, species) %>%
    summarise(abund = sum(abund, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = species, values_from = abund, values_fill = 0) %>%
    column_to_rownames("ID") %>%
    as.data.frame()
} else {
  colnames(species_raw) <- new_names
}

species_raw <- species_raw[, !is_unknown_sgb(colnames(species_raw)), drop = FALSE]

abundance_threshold  <- 2e-4
prevalence_threshold <- 0.10

get_group_species_by_study <- function(df, species_mat, group_value,
                                       abundance_threshold = 2e-4,
                                       prevalence_threshold = 0.10,
                                       study_col = "Study",
                                       id_col = "ID",
                                       group_col = "Group_new") {
  df2 <- df %>% filter(.data[[group_col]] == group_value)
  ids <- df2[[id_col]]
  ids <- ids[ids %in% rownames(species_mat)]
  df2 <- df2 %>% filter(.data[[id_col]] %in% ids)
  
  taxa_list <- lapply(split(df2[[id_col]], df2[[study_col]]), function(sample_ids) {
    Xsub <- species_mat[sample_ids, , drop = FALSE]
    colMeans(Xsub > abundance_threshold, na.rm = TRUE) >= prevalence_threshold
  })
  
  keep <- Reduce("|", taxa_list)
  names(keep)[keep]
}

case_label    <- "AP"
control_label <- "Control"

stool_case   <- get_group_species_by_study(stool_data,  species_raw, case_label,    abundance_threshold, prevalence_threshold)
saliva_case  <- get_group_species_by_study(saliva_data, species_raw, case_label,    abundance_threshold, prevalence_threshold)
stool_ctrl   <- get_group_species_by_study(stool_data,  species_raw, control_label, abundance_threshold, prevalence_threshold)
saliva_ctrl  <- get_group_species_by_study(saliva_data, species_raw, control_label, abundance_threshold, prevalence_threshold)

overlap_case  <- intersect(stool_case, saliva_case)
overlap_ctrl  <- intersect(stool_ctrl, saliva_ctrl)
final_overlap <- intersect(overlap_case, overlap_ctrl)

writeLines(overlap_case,  "result/overlap1_case_species_ap.txt")
writeLines(overlap_ctrl,  "result/overlap2_ctrl_species_ap.txt")
writeLines(final_overlap, "result/overlap_final_species_ap.txt")

# ===================== Read FINAL overlap list =====================

overlap_final <- readLines("result/overlap_final_species_ap.txt") %>%
  clean_species_name() %>%
  unique()

overlap_final <- intersect(overlap_final, colnames(species_raw))

stool_ids  <- intersect(stool_data$ID,  rownames(species_raw))
saliva_ids <- intersect(saliva_data$ID, rownames(species_raw))

stopifnot(length(stool_ids) > 0, length(saliva_ids) > 0)

mean_stool <- colMeans(species_raw[stool_ids,  , drop = FALSE], na.rm = TRUE)
mean_oral  <- colMeans(species_raw[saliva_ids, , drop = FALSE], na.rm = TRUE)

prev_stool <- colMeans(species_raw[stool_ids,  , drop = FALSE] > 0, na.rm = TRUE)
prev_oral  <- colMeans(species_raw[saliva_ids, , drop = FALSE] > 0, na.rm = TRUE)

df_scatter <- tibble(
  species    = names(mean_stool),
  mean_stool = as.numeric(mean_stool),
  mean_oral  = as.numeric(mean_oral),
  prev_stool = as.numeric(prev_stool),
  prev_oral  = as.numeric(prev_oral)
) %>%
  mutate(group = if_else(species %in% overlap_final, "Overlap", "Other"))

eps_abund <- 1e-10
eps_prev  <- 1e-4

df_scatter <- df_scatter %>%
  mutate(
    mean_stool_plot = pmax(mean_stool, eps_abund),
    mean_oral_plot  = pmax(mean_oral,  eps_abund),
    prev_stool_plot = pmax(prev_stool, eps_prev),
    prev_oral_plot  = pmax(prev_oral,  eps_prev)
  )

p1 <- ggplot(
  df_scatter,
  aes(
    x = mean_oral_plot,
    y = mean_stool_plot,
    color = group,
    shape = group
  )
) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("Other" = "grey70", "Overlap" = "red3")) +
  scale_shape_manual(values = c("Other" = 16, "Overlap" = 17)) +  
  labs(
    x = "Mean relative abundance (across all saliva samples)",
    y = "Mean relative abundance (across all stool samples)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

print(p1)

p2 <- ggplot(
  df_scatter,
  aes(
    x = prev_oral_plot,
    y = prev_stool_plot,
    color = group,
    shape = group
  )
) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_x_log10(limits = c(eps_prev, 1)) +
  scale_y_log10(limits = c(eps_prev, 1)) +
  scale_color_manual(values = c("Other" = "grey70", "Overlap" = "red3")) +
  scale_shape_manual(values = c("Other" = 16, "Overlap" = 17)) +
  labs(
    x = "Prevalence (across all saliva samples)",
    y = "Prevalence (across all stool samples)"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")

print(p2)

# ----------------- ORAL → GUT TRANSLOCATION BY COHORT (STOOL ONLY) -----------------

meta_use <- stool_data

overlap_use <- intersect(overlap_final, colnames(species_raw))
stopifnot(length(overlap_use) > 0)

keep_ids <- intersect(meta_use$ID, rownames(species_raw))

df_frac_cohort <- meta_use %>%
  filter(ID %in% keep_ids, Group_new %in% c("Control", "AP")) %>%
  transmute(
    ID,
    cohort = Study,
    Group_new = factor(Group_new, levels = c("Control", "AP"))
  ) %>%
  filter(!is.na(cohort), cohort != "") %>%
  rowwise() %>%
  mutate(
    fraction_oral_gut = mean(species_raw[ID, overlap_use] > 0, na.rm = TRUE)
  ) %>%
  ungroup()

df_frac_overall <- df_frac_cohort %>% mutate(cohort = "Overall")

df_frac_plot <- bind_rows(df_frac_overall, df_frac_cohort) %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c("Overall", sort(unique(as.character(df_frac_cohort$cohort))))
    )
  )

# p3 <- ggplot(df_frac_plot, aes(x = Group_new, y = fraction_oral_gut, fill = Group_new)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.7) +
#   geom_jitter(width = 0.18, size = 1.8, alpha = 0.55) +
#   stat_compare_means(method = "wilcox.test", label = "p.format", label.y = 0.98) +
#   scale_y_continuous(limits = c(0, 1)) +
#   scale_fill_manual(values = c("Control" = "grey70", "PDAC" = "red3")) +
#   facet_wrap(~ cohort, ncol = length(levels(df_frac_plot$cohort))) +
#   labs(x = NULL, y = "Fraction of oral-origin species detected in stool") +
#   theme_classic(base_size = 14) +
#   theme(
#     legend.position = "none",
#     strip.background = element_blank(),
#     strip.text = element_text(face = "bold")
#   )
# 
# print(p3)

df_frac_plot2 <- df_frac_plot %>%
  mutate(
    x_center = as.numeric(Group_new),
    x_cloud  = x_center + 0.22
  )

p3 <- ggplot(df_frac_plot2, aes(y = fraction_oral_gut)) +
  
geom_half_violin(
  aes(x = x_cloud, fill = Group_new),
  side = "r",
  trim = FALSE,
  alpha = 0.45,
  width = 0.9, 
  color = "black"
) +

geom_point(
  aes(x = x_cloud),
  position = position_jitter(width = 0.06, height = 0.015),
  size = 1.6,
  alpha = 0.6,
  color = "black"
  ) + 
  
geom_boxplot(
  aes(x = x_center, fill = Group_new),
  width = 0.18,
  outlier.shape = NA,
  alpha = 0.9, 
  color = "black"
) +
  
stat_compare_means(
  aes(x = x_center, group = Group_new),
  method = "wilcox.test",
  label = "p.format",
  label.y = 0.98
) +
  
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Control", "PDAC")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  
scale_fill_manual(
  values = c(
    "Control" = "#9ecae1",
    "AP"    = "red3" 
  )
) +
  scale_color_manual(
    values = c(
      "Control" = "#6baed6",
      "AP"    = "red3"
    )
  ) +
  
  facet_wrap(~ cohort, ncol = length(levels(df_frac_plot$cohort))) +
  labs(
    x = NULL,
    y = "Fraction of oral-origin species detected in stool"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p3)

# -------------------- HEATMAP (STOOL): 25 overlap species --------------------

meta_hm <- stool_data %>%
  filter(Group_new %in% c("Control", "AP")) %>%
  mutate(
    cohort = Study,
    Group_new = factor(Group_new, levels = c("Control", "AP"))
  )

keep_ids <- intersect(meta_hm$ID, rownames(species_raw))
meta_hm  <- meta_hm %>% filter(ID %in% keep_ids)

overlap_use <- intersect(overlap_final, colnames(species_raw))
stopifnot(length(overlap_use) == 25)

df_long <- species_raw[meta_hm$ID, overlap_use, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "species", values_to = "abund") %>%
  left_join(meta_hm %>% select(ID, cohort, Group_new), by = "ID")

df_heat <- df_long %>%
  group_by(Group_new, cohort, species) %>%
  summarise(
    mean_abund = mean(abund, na.rm = TRUE),
    prev       = mean(abund > 0, na.rm = TRUE),
    .groups = "drop"
  )

eps <- 1e-10
df_heat <- df_heat %>%
  mutate(value = log10(mean_abund + eps))

species_order <- df_heat %>%
  group_by(species) %>%
  summarise(overall = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall) %>%
  pull(species)

df_heat <- df_heat %>%
  mutate(
    species = factor(species, levels = species_order),
    cohort  = factor(cohort, levels = sort(unique(cohort)))
  )

p4 <- ggplot(df_heat, aes(x = cohort, y = species, fill = value)) +
  geom_tile() +
  facet_grid(. ~ Group_new) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "white",
    high = "red3",
    midpoint = median(df_heat$value, na.rm = TRUE),
    name = "log10(mean RA)"
  ) +
  labs(
    x = "Cohort",
    y = "Overlap species (n = 25)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
print(p4)

# -------------------- HEATMAP (SALIVA): 25 overlap species --------------------

meta_hm_s <- saliva_data %>%
  filter(Group_new %in% c("Control", "AP")) %>%
  mutate(
    cohort = Study,
    Group_new = factor(Group_new, levels = c("Control", "AP"))
  )

keep_ids_s <- intersect(meta_hm_s$ID, rownames(species_raw))
meta_hm_s  <- meta_hm_s %>% filter(ID %in% keep_ids_s)

overlap_use_s <- intersect(overlap_final, colnames(species_raw))
stopifnot(length(overlap_use_s) == 25)

df_long_s <- species_raw[meta_hm_s$ID, overlap_use_s, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "species", values_to = "abund") %>%
  left_join(meta_hm_s %>% select(ID, cohort, Group_new), by = "ID")

df_heat_s <- df_long_s %>%
  group_by(Group_new, cohort, species) %>%
  summarise(
    mean_abund = mean(abund, na.rm = TRUE),
    prev       = mean(abund > 0, na.rm = TRUE),
    .groups = "drop"
  )

eps <- 1e-10
df_heat_s <- df_heat_s %>%
  mutate(value = log10(mean_abund + eps))

species_order_s <- df_heat_s %>%
  group_by(species) %>%
  summarise(overall = mean(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(overall) %>%
  pull(species)

df_heat_s <- df_heat_s %>%
  mutate(
    species = factor(species, levels = species_order_s),
    cohort  = factor(cohort, levels = sort(unique(cohort)))
  )

p5 <- ggplot(df_heat_s, aes(x = cohort, y = species, fill = value)) +
  geom_tile() +
  facet_grid(. ~ Group_new) +
  scale_fill_gradient2(
    low = "blue4",
    mid = "white",
    high = "red3",
    midpoint = median(df_heat_s$value, na.rm = TRUE),
    name = "log10(mean RA)"
  ) +
  labs(
    x = "Cohort",
    y = "Overlap species (n = 25)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p5)

# -------------------- HEATMAP (ALL): Stool + Saliva, 25 overlap species --------------------

species_levels <- rev(sort(overlap_final))

eps <- 1e-10

meta_hm_all <- pdac %>%
  filter(Type %in% c("Stool", "Saliva")) %>%
  filter(Group_new %in% c("Control", "AP")) %>%
  mutate(
    cohort = Study,
    Group_new = factor(Group_new, levels = c("Control", "AP")),
    Type = factor(Type, levels = c("Stool", "Saliva"))
  )

keep_ids_all <- intersect(meta_hm_all$ID, rownames(species_raw))
meta_hm_all  <- meta_hm_all %>% filter(ID %in% keep_ids_all)

overlap_use_all <- intersect(overlap_final, colnames(species_raw))
stopifnot(length(overlap_use_all) == 25)

df_long_all <- species_raw[meta_hm_all$ID, overlap_use_all, drop = FALSE] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID") %>%
  pivot_longer(-ID, names_to = "species", values_to = "abund") %>%
  left_join(meta_hm_all %>% select(ID, cohort, Group_new, Type), by = "ID")

df_heat_all <- df_long_all %>%
  group_by(Type, Group_new, cohort, species) %>%
  summarise(
    mean_abund = mean(abund, na.rm = TRUE),
    prev       = mean(abund > 0, na.rm = TRUE),
    .groups = "drop"
  )

species_levels <- rev(sort(overlap_final))

df_heat_all_plot <- df_heat_all %>%
  mutate(
    value = log10(mean_abund + eps),
    species = factor(species, levels = species_levels),
    Type = factor(Type, levels = c("Stool", "Saliva")),
    Group_new = factor(Group_new, levels = c("Control", "AP")),
    panel = factor(
      paste(Type, Group_new, sep = " | "),
      levels = c("Stool | Control", "Stool | AP", "Saliva | Control", "Saliva | AP")
    )
  )

df_heat_all_plot <- df_heat_all %>%
  mutate(
    value = log10(mean_abund + eps),
    species = factor(species, levels = species_levels),
    Type = factor(Type, levels = c("Stool", "Saliva")),
    Group_new = factor(Group_new, levels = c("Control", "AP")),
    panel = factor(
      paste(Type, Group_new, sep = " | "),
      levels = c("Stool | Control", "Stool | AP", "Saliva | Control", "Saliva | AP")
    )
  )

df_heat_all_plot <- df_heat_all_plot %>%
  complete(panel, species, cohort, fill = list(value = NA_real_)) %>%
  group_by(panel, cohort) %>%
  mutate(keep_col = any(!is.na(value))) %>%
  ungroup() %>%
  filter(keep_col) %>%
  select(-keep_col)

mid_all <- median(df_heat_all_plot$value, na.rm = TRUE)

p45 <- ggplot(df_heat_all_plot, aes(x = cohort, y = species, fill = value)) +
  geom_tile() +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_x_discrete(drop = TRUE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient2(
    low = "blue4", mid = "white", high = "red3",
    midpoint = mid_all,
    na.value = "white",
    name = "log10(mean RA)"
  ) +
  labs(x = "Cohort", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10, face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(0.5, "lines")
  )

print(p45)

# ===================== GLOBAL THEME SETTINGS =====================

AXIS_TITLE_SIZE <- 15
AXIS_TEXT_SIZE  <- 12

theme_all <- theme_classic(base_size = AXIS_TEXT_SIZE) +
  theme(
    axis.title = element_text(size = AXIS_TITLE_SIZE),
    axis.text  = element_text(size = AXIS_TEXT_SIZE)
  )

p1 <- p1 + theme_all
p2 <- p2 + theme_all
p3 <- p3 + theme_all

row1 <- plot_grid(p1, p2, labels = c("a", "b"), label_size = 20, ncol = 2, align = "hv")
row2 <- plot_grid(p3, labels = c("c"), label_size = 20, ncol = 1, align = "hv")
row3 <- plot_grid(p45, labels = c("d"), label_size = 20, ncol = 1, align = "hv")

p <- plot_grid(row1, row2, row3, ncol = 1, rel_heights = c(0.9, 0.75, 0.9))
print(p)

ggsave("Fig4_AP.pdf", p, width = 15, height = 18)
