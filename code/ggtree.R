rm(list = ls())

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(tidyr)
  library(Biostrings)
})

tree_path <- "/Users/yrx/Downloads/pancreatitis/out6939/trees/RAxML_bestTree.t__SGB6939.StrainPhlAn4.tre"
aln_path  <- "/Users/yrx/Downloads/pancreatitis/out6939/trees/t__SGB6939.StrainPhlAn4_concatenated.aln"
meta_path <- "/Users/yrx/Downloads/pancreatitis/data/PDAC.txt"

out_dir <- "/Users/yrx/Downloads/strainphlan_vis_SGB6939"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

slice_start <- 1     
slice_len   <- 100   
dx          <- 0.02  
tile_h      <- 0.6

tr <- read.tree(tree_path)

meta <- read.delim(meta_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

meta2 <- meta %>%
  transmute(
    label     = trimws(SampleID),     
    Group_new = trimws(Group_new),
    Study     = trimws(Study),
    Type      = trimws(Type)
  )

meta2 <- meta2 %>%
  filter(!is.na(Group_new) & Group_new %in% c("Control", "PDAC"))

keep <- intersect(tr$tip.label, meta2$label)
tr2  <- keep.tip(tr, keep)

aln <- readDNAStringSet(aln_path)
aln_ids <- names(aln)

tips <- tr2$tip.label
common_ids <- Reduce(intersect, list(tips, meta2$label, aln_ids))

if (length(common_ids) < 3) {
  stop("Too few common samples among tree tips, metadata, and alignment. Check IDs / trimming.")
}

tr2 <- keep.tip(tr2, common_ids)

meta2 <- meta2 %>%
  filter(label %in% common_ids)

aln2 <- aln[common_ids]
aln2 <- aln2[tr2$tip.label]  

stopifnot(all(tr2$tip.label == names(aln2)))

aln_mat <- as.matrix(aln2)
n_pos <- ncol(aln_mat)

slice_end <- min(n_pos, slice_start + slice_len - 1)
pos_idx <- slice_start:slice_end

snp_slice <- aln_mat[, pos_idx, drop = FALSE]

snp_df <- as.data.frame(snp_slice, stringsAsFactors = FALSE)
snp_df$label <- rownames(snp_df)

snp_long <- snp_df %>%
  pivot_longer(
    cols = -label,
    names_to = "pos",
    values_to = "nt"
  ) %>%
  mutate(
    pos_num = suppressWarnings(as.integer(gsub("^V", "", pos))),
    label   = trimws(label),
    nt      = as.character(nt)
  )

if (all(is.na(snp_long$pos_num))) {
  snp_long <- snp_long %>%
    group_by(pos) %>%
    mutate(pos_num = cur_group_id()) %>%
    ungroup()
}

p_base <- ggtree(tr2) %<+% meta2 +
  geom_tippoint(aes(color = Group_new), size = 0.6) +
  scale_color_manual(values = c(Control = "#336fac", PDAC = "#b82929"),
                     na.translate = FALSE, name = "Group") +
  theme_tree2() +
  ggtitle(sprintf("StrainPhlAn tree + SNP slice (%d-%d)", slice_start, slice_end))

tip_y <- p_base$data %>%
  dplyr::filter(isTip) %>%
  dplyr::select(label, y)

snp_long2 <- snp_long %>%
  left_join(tip_y, by = "label") %>%
  filter(!is.na(y))   

max_sites   <- 20   
min_var_n   <- 10   
use_maf     <- FALSE
maf_min     <- 0.05 
maf_max     <- 0.50 

valid_nt <- c("A","C","G","T","-")

snp_f <- snp_long2 %>%
  mutate(nt2 = ifelse(nt %in% valid_nt, nt, NA_character_)) %>%
  filter(!is.na(nt2))

site_stat <- snp_f %>%
  group_by(pos_num) %>%
  summarise(
    n_nonmiss = sum(!is.na(nt2)),
    major = names(sort(table(nt2), decreasing = TRUE))[1],
    major_freq = max(table(nt2)) / n_nonmiss,
    var_n = n_nonmiss - max(table(nt2)),
    maf = 1 - major_freq,
    entropy = {
      p <- table(nt2) / n_nonmiss
      -sum(p * log(p))
    },
    .groups = "drop"
  )

if (!use_maf) {
  site_keep <- site_stat %>%
    filter(var_n >= min_var_n) %>%
    arrange(desc(var_n), desc(entropy)) %>%
    slice_head(n = max_sites) %>%
    pull(pos_num)
} else {
  site_keep <- site_stat %>%
    filter(maf >= maf_min, maf <= maf_max) %>%
    arrange(desc(maf), desc(entropy)) %>%
    slice_head(n = max_sites) %>%
    pull(pos_num)
}

cat("Keeping SNP sites:", length(site_keep), "\n")

snp_long_keep <- snp_f %>%
  filter(pos_num %in% site_keep) %>%
  mutate(
    pos_rank = as.integer(factor(pos_num, levels = sort(unique(site_keep))))
  )

x0 <- max(p_base$data$x, na.rm = TRUE)

dx2 <- 0.05

p_clean <- p_base +
  xlim(0, x0 + dx2 * (max(snp_long_keep$pos_rank) + 5)) +
  geom_tile(
    data = snp_long_keep,
    aes(
      x = x0 + dx2 * pos_rank,
      y = y,
      fill = nt2
    ),
    inherit.aes = FALSE,
    height = tile_h
  ) +
  scale_fill_manual(
    values = c(A="#4daf4a", C="#377eb8", G="#ff7f00", T="#e41a1c", "-"="grey80"),
    na.value = "grey90",
    name = "SNP"
  ) +
  ggtitle(sprintf("StrainPhlAn tree + high-info SNPs (kept %d sites)", length(site_keep)))
p_clean

ggsave("/Users/yrx/Downloads/pancreatitis/ggtree.pdf", p_clean, width = 12, height = 10, limitsize = FALSE)
