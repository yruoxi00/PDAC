library(tidyverse)
library(dplyr)
library(anpan)
library(ape)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(posterior)
library(ggnewscale)
library(ggsci)
library(tidytree)
library(purrr)
library(phangorn)
library(cowplot)

setwd("/Users/yrx/Downloads/pancreatitis")

### Fusicatenibacter saccharivorans

meta <- readr::read_tsv("data/PDAC.tsv", show_col_types = FALSE)
tr <- read.tree("out4874/trees/RAxML_bestTree.t__SGB4874.StrainPhlAn4.tre")

fit <- anpan_pglmm(
  meta_file = "data/PDAC.tsv",
  tree_file = "out4874/trees/RAxML_bestTree.t__SGB4874.StrainPhlAn4.tre",
  outcome   = "PDAC",
  family    = "binomial"
)

draws_summary <- as.data.frame(fit$pglmm_fit$draws("phylo_effect") %>% posterior::summarize_draws())
s_ids <- fit$model_input$sample_id

phylo_effect_df_clean <- data.frame(
  sample_id = as.character(s_ids), 
  effect_mean = draws_summary$mean,
  stringsAsFactors = FALSE
) %>% filter(!is.na(effect_mean))

plot_meta_clean <- meta %>% 
  filter(sample_id %in% tr$tip.label) %>%
  left_join(phylo_effect_df_clean, by = "sample_id") %>%
  filter(!is.na(Study), !is.na(PDAC), !is.na(effect_mean))

keep_tips <- plot_meta_clean$sample_id
tr_trimmed <- keep.tip(tr, keep_tips)

tree_obj_clean <- as_tibble(tr_trimmed) %>%
  left_join(plot_meta_clean, by = c("label" = "sample_id")) %>%
  as.treedata()

study_col <- c(
  "Vilà-Quintana_2024" = "#7DABD1",
  "Kartal-1_2022"      = "#C7D3E5", 
  "Kartal-2_2022"      = "#F6D963",
  "Matsukawa_2021"     = "#FAADD5",
  "Nagata_2022"        = "#BC8CC3"
)
pdac_col_val <- c("0" = "#336fac", "1" = "#b82929")
col_sec <- pal_d3('category20')(13)

target_clades <- data.frame(
  node = c(329, 554, 395, 473, 452), 
  color = col_sec[c(1, 12, 3, 5, 7)],
  stringsAsFactors = FALSE
)

p_f <- ggtree(tree_obj_clean, layout = "fan", size = 0.15, open.angle = 10) +
  xlim(-0.1, 0.6) 

for(i in 1:nrow(target_clades)) {
  p_f <- p_f + 
    geom_hilight(
      node = target_clades$node[i], 
      fill = target_clades$color[i], 
      alpha = 0.2, 
      extendto = 0.45
    )
}

p_f <- p_f +
  geom_tippoint(aes(color = as.factor(PDAC)), size = 0.6) +
  scale_color_manual(values = pdac_col_val, name = "PDAC Status") +
  
  geom_fruit(
    geom = geom_tile, 
    mapping = aes(fill = Study), 
    width = 0.02, 
    offset = 0.08, 
    inherit.aes = FALSE
  ) + 
  scale_fill_manual(values = study_col, name = "Study Group") +
  
  new_scale_fill() +
  
  geom_fruit(
    geom = geom_tile, 
    mapping = aes(fill = effect_mean), 
    width = 0.02, 
    offset = 0.14, 
    inherit.aes = FALSE
  ) +  
  scale_fill_viridis_c(option = "G", direction = -1, name = "Phylo Effect Mean")

p_f <- p_f + 
  ggtitle("Fusicatenibacter saccharivorans") +
  theme(
    legend.position = "right", 
    legend.text = element_text(size = 7), 
    plot.title = element_text(hjust = 0.5, face = "italic", size = 12)
  ) + 
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2),
    fill_new = guide_colorbar(
      title.position = "top", 
      barwidth = 0.8, 
      barheight = 6,
      order = 3
    )
  )
print(p_f)

### Eubacterium rectale

meta <- readr::read_tsv("data/PDAC.tsv", show_col_types = FALSE)
tr <- read.tree("out4933/trees/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre")

fit <- anpan_pglmm(
  meta_file = "data/PDAC.tsv",
  tree_file = "out4933/trees/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre",
  outcome   = "PDAC",
  family    = "binomial"
)

draws_summary <- as.data.frame(fit$pglmm_fit$draws("phylo_effect") %>% posterior::summarize_draws())
s_ids <- fit$model_input$sample_id

phylo_effect_df_clean <- data.frame(
  sample_id = as.character(s_ids), 
  effect_mean = draws_summary$mean,
  stringsAsFactors = FALSE
) %>% filter(!is.na(effect_mean))

plot_meta_clean <- meta %>% 
  filter(sample_id %in% tr$tip.label) %>%
  left_join(phylo_effect_df_clean, by = "sample_id") %>%
  filter(!is.na(Study), !is.na(PDAC), !is.na(effect_mean))

keep_tips <- plot_meta_clean$sample_id
tr_trimmed <- keep.tip(tr, keep_tips)

tree_obj_clean <- as_tibble(tr_trimmed) %>%
  left_join(plot_meta_clean, by = c("label" = "sample_id")) %>%
  as.treedata()

study_col <- c(
  "Vilà-Quintana_2024" = "#7DABD1",
  "Kartal-1_2022"      = "#C7D3E5", 
  "Kartal-2_2022"      = "#F6D963",
  "Matsukawa_2021"     = "#FAADD5",
  "Nagata_2022"        = "#BC8CC3"
)
pdac_col_val <- c("0" = "#336fac", "1" = "#b82929")
col_sec <- pal_d3('category20')(13)

target_clades <- data.frame(
  node = c(303, 352, 485, 531), 
  color = col_sec[c(1, 12, 3, 7)],
  stringsAsFactors = FALSE
)

p_e <- ggtree(tree_obj_clean, layout = "fan", size = 0.15, open.angle = 10) +
  xlim(-0.1, 0.6) 

for(i in 1:nrow(target_clades)) {
  p_e <- p_e + 
    geom_hilight(
      node = target_clades$node[i], 
      fill = target_clades$color[i], 
      alpha = 0.2, 
      extendto = 0.45
    )
}

p_e <- p_e +
  geom_tippoint(aes(color = as.factor(PDAC)), size = 0.6) +
  scale_color_manual(values = pdac_col_val, name = "PDAC Status") +
  
  geom_fruit(
    geom = geom_tile, 
    mapping = aes(fill = Study), 
    width = 0.02, 
    offset = 0.08, 
    inherit.aes = FALSE
  ) + 
  scale_fill_manual(values = study_col, name = "Study Group") +
  
  new_scale_fill() +
  
  geom_fruit(
    geom = geom_tile, 
    mapping = aes(fill = effect_mean), 
    width = 0.02, 
    offset = 0.14, 
    inherit.aes = FALSE
  ) +  
  scale_fill_viridis_c(option = "G", direction = -1, name = "Phylo Effect Mean")

p_e <- p_e + 
  ggtitle("Eubacterium rectale") +
  theme(
    legend.position = "right", 
    legend.text = element_text(size = 7), 
    plot.title = element_text(hjust = 0.5, face = "italic", size = 12)
  ) + 
  guides(
    color = guide_legend(order = 1),
    fill = guide_legend(order = 2),
    fill_new = guide_colorbar(
      title.position = "top", 
      barwidth = 0.8, 
      barheight = 6,
      order = 3
    )
  )
print(p_e)

p <- plot_grid(p_f, p_e, 
               labels = c("b", ""), 
               ncol = 2, 
               align = 'hv')
print(p)

ggsave("Fig5b.pdf", p, width = 10, height = 11)

### subclade 

plot_meta_clean$PDAC_num <- as.numeric(as.character(plot_meta_clean$PDAC))
total_pdac <- sum(plot_meta_clean$PDAC_num == 1, na.rm = TRUE)
total_samples <- nrow(plot_meta_clean)

clades_list <- prop.part(tr_trimmed)
n_tips <- ape::Ntip(tr_trimmed)

analyze_node_native <- function(node_idx, clades, tr, meta) {
  curr_node <- n_tips + node_idx
  tip_indices <- clades[[node_idx]]
  descendant_tips <- tr$tip.label[tip_indices]
  n_samples <- length(descendant_tips)
  
  if (n_samples < 10 || n_samples > (n_tips * 0.5)) return(NULL)
  
  sub_meta <- meta %>% filter(sample_id %in% descendant_tips)
  pdac_in_clade <- sum(sub_meta$PDAC_num == 1, na.rm = TRUE)
  pdac_ratio <- pdac_in_clade / n_samples
  
  mat <- matrix(c(
    pdac_in_clade, total_pdac - pdac_in_clade,
    n_samples - pdac_in_clade, (total_samples - n_samples) - (total_pdac - pdac_in_clade)
  ), nrow = 2)
  p_val <- fisher.test(mat)$p.value
  
  return(data.frame(node = curr_node, n_samples = n_samples, 
                    pdac_ratio = pdac_ratio, p_value = p_val))
}

node_results <- map_dfr(1:length(clades_list), ~analyze_node_native(.x, clades_list, tr_trimmed, plot_meta_clean))

candidates <- node_results %>%
  filter(p_value < 0.05) %>%
  arrange(desc(pdac_ratio)) 

final_nodes <- c()
for (i in 1:nrow(candidates)) {
  curr_node <- candidates$node[i]
  
  desc_nodes <- tidytree::offspring(as_tibble(tr_trimmed), curr_node)$node
  anc_nodes <- phangorn::Ancestors(tr_trimmed, curr_node, type = "all")
  
  if (!any(c(desc_nodes, anc_nodes) %in% final_nodes)) {
    final_nodes <- c(final_nodes, curr_node)
  }
  
  if (length(final_nodes) >= 5) break
}

target_clades <- candidates %>%
  filter(node %in% final_nodes) %>%
  mutate(
    tag = LETTERS[1:n()],
    color = col_sec[1:n()]
  ) %>%
  select(node, tag, color)

print(target_clades)
