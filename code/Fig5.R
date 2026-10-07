library(anpan)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(ape)
library(ggplot2)
library(egg)
library(data.table)
library(ggpubr)
library(cowplot)
library(grid)
library(tidyverse)
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

setwd("/Users/yrx/Downloads/pancreatitis")

meta_file <- "data/PDAC.tsv"
meta <- readr::read_tsv(meta_file, show_col_types = FALSE)

tree_map <- tibble::tribble(
  ~species,   ~tree_file,
  "Veillonella_atypica", "out6936/trees/RAxML_bestTree.t__SGB6936.StrainPhlAn4.tre",
  "Veillonella_parvula", "out6939/trees/RAxML_bestTree.t__SGB6939.StrainPhlAn4.tre",
  "Clostridium_symbiosum", "out4699/trees/RAxML_bestTree.t__SGB4699.StrainPhlAn4.tre", 
  "Akkermansia_muciniphila", "out9226/trees/RAxML_bestTree.t__SGB9226.StrainPhlAn4.tre", 
  "Hungatella_hathewayi", "out4742/trees/RAxML_bestTree.t__SGB4742.StrainPhlAn4.tre", 
  
  "Eubacterium_ventriosum", "out5045/trees/RAxML_bestTree.t__SGB5045.StrainPhlAn4.tre", 
  "Eubacterium_rectale", "out4933/trees/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre",
  "Agathobaculum_butyriciproducens", "out14993/trees/RAxML_bestTree.t__SGB14993.StrainPhlAn4.tre", 
  "Fusicatenibacter_saccharivorans", "out4874/trees/RAxML_bestTree.t__SGB4874.StrainPhlAn4.tre", 
  "Anaerobutyricum_hallii", "out4532/trees/RAxML_bestTree.t__SGB4532.StrainPhlAn4.tre", 
  "Ruminococcus_bromii", "out4285/trees/RAxML_bestTree.t__SGB4285.StrainPhlAn4.tre", 
  "Clostridiaceae_bacterium", "out4721/trees/RAxML_bestTree.t__SGB4721.StrainPhlAn4.tre", 
  "Lachnospira_pectinoschiza", "out5089/trees/RAxML_bestTree.t__SGB5089.StrainPhlAn4.tre"
)

run_one <- function(species, tree_file){
  tr <- ape::read.tree(tree_file)
  
  fit <- anpan::anpan_pglmm(
    meta_file    = meta_file,
    tree_file    = tree_file,
    trim_pattern = NULL,
    covariates   = NULL,
    outcome      = "PDAC",
    family       = "binomial",
    verbose      = TRUE,
    cores        = 4
  )
  
  comp <- as.data.frame(fit$loo$comparison)
  base_row <- comp["base_fit", , drop = FALSE]
  
  ids <- tr$tip.label
  n_tbl <- meta %>%
    filter(sample_id %in% ids) %>%
    count(PDAC) %>%
    mutate(Group = ifelse(PDAC == 1, "PDAC", "Control")) %>%
    select(Group, n)
  
  tibble(
    species   = species,
    tree_file = tree_file,
    elpd_diff = base_row$elpd_diff,
    se_diff   = base_row$se_diff,
    n_total   = sum(n_tbl$n),
    n_pdac    = n_tbl$n[n_tbl$Group=="PDAC"]    %>% {ifelse(length(.)==0,0,.)},
    n_ctrl    = n_tbl$n[n_tbl$Group=="Control"] %>% {ifelse(length(.)==0,0,.)}
  )
}

re_pglmm_t <- purrr::pmap_dfr(tree_map, run_one)
re_pglmm_t

# PDAC status species
pdac_depleted_spe <- c("s__Eubacterium_ventriosum", #SGB5045
                       "s__Eubacterium_rectale", #SGB4933
                       "s__Agathobaculum_butyriciproducens", #SGB14993
                       "s__Fusicatenibacter_saccharivorans", #SGB4874
                       "s__Anaerobutyricum_hallii", #SGB4532
                       "s__Ruminococcus_bromii", #SGB4285
                       "s__Clostridiaceae_bacterium", #SGB4721
                       "s__Lachnospira_pectinoschiza") #SGB5089
length(pdac_depleted_spe) #8
pdac_enriched_spe <- c("s__Veillonella_atypica", #SGB6936
                       "s__Veillonella_parvula", #SGB6939
                       "s__Clostridium_symbiosum", #SGB4699
                       "s__Akkermansia_muciniphila", #SGB9226
                       "s__Hungatella_hathewayi") #SGB4742
length(pdac_enriched_spe) #5
spe_pdac_status <- data.table(
  species=c(pdac_depleted_spe,pdac_enriched_spe),
  pdac_status=c(rep("PDAC depleted",8),
                rep("PDAC enriched",5))
)
spe_pdac_status
spe_pdac_status$species <- gsub("s__","",spe_pdac_status$species)

re_plot <- re_pglmm_t %>%
  mutate(elpd_2se = abs(elpd_diff) - 2*se_diff,
         sig = ifelse(elpd_2se > 0, "*", "")) %>%
  left_join(spe_pdac_status, by="species") %>%
  mutate(pdac_status = ifelse(is.na(pdac_status), "Non-significant", pdac_status),
         species_update = gsub("_"," ", species)) %>%
  arrange(desc(elpd_diff))

re_plot$species_update <- factor(re_plot$species_update, levels = re_plot$species_update)

re_samplesize <- re_plot %>%
  select(species, species_update, n_pdac, n_ctrl) %>%
  tidyr::pivot_longer(cols = c(n_ctrl, n_pdac),
                      names_to = "Group", values_to = "n") %>%
  mutate(Group = recode(Group, n_ctrl="Control", n_pdac="PDAC"),
         species_update = factor(species_update, levels = levels(re_plot$species_update)))

p_pglmm <- ggplot(re_plot, aes(x = -elpd_diff, y = species_update)) +
  geom_errorbar(aes(xmin = -(elpd_diff - se_diff), xmax = -(elpd_diff + se_diff)),
                width = 0.2, color = "#1F4690") +
  geom_point(color = "#1F4690") +
  theme_bw() +
  labs(x = expression(Delta * "ELPD"))+
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank())

p_status <- ggplot(re_plot, aes(x=1, y=species_update, fill=pdac_status)) +
  geom_col(width=1) +
  scale_fill_manual(
    breaks = c("Non-significant","PDAC depleted","PDAC enriched"),
    values = c("grey","#5050FFFF","#CE3D32FF")
  ) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, face="italic", color="black"),
        axis.ticks.x = element_blank(),
        legend.position="none")

p_n <- ggplot(re_samplesize, aes(y=species_update, x=n, fill=Group)) +
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=c("Control"="#72b0da", "PDAC"="#C75127")) +
  labs(x="Sample size") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        legend.position="none")

blank <- ggplot() + theme_void()

p <- egg::ggarrange(
  p_status, blank, p_pglmm, p_n,
  nrow=1, widths=c(0.15, -0.15, 3, 1.5)
)

p

p_status_lgd <- p_status +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "Species status"))

p_n_lgd <- p_n +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title = "PDAC status"))

lgd_status <- ggpubr::get_legend(p_status_lgd)
lgd_group  <- ggpubr::get_legend(p_n_lgd)

lgd_all <- cowplot::plot_grid(
  lgd_status, lgd_group,
  ncol = 1)

p_with_lgd <- cowplot::plot_grid(
  p, lgd_all,
  ncol = 2,
  rel_widths = c(1, 0.2)
)

p_with_lgd

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

p1 <- plot_grid(p_f, p_e, 
               labels = c("b", ""), 
               ncol = 2, 
               align = 'hv')
p <- plot_grid(p_with_lgd, p1, 
               labels = c("a", ""), 
               ncol = 1, 
               align = 'hv')
p

ggsave("Fig5.pdf", p, width = 10, height = 12)
