library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

tree_file <- "/Users/yrx/Downloads/pancreatitis/out4874/trees/RAxML_bestTree.t__SGB4874.StrainPhlAn4.tre"
meta_file <- "/Users/yrx/Downloads/pancreatitis/data/PDAC.txt"

id_col_name    <- NULL        
group_col_name <- "Group_new"
type_col_name  <- "Type" 
study_col_name <- "Study"

plot_title <- "Fusicatenibacter saccharivorans"

tr0   <- ape::read.tree(tree_file)
meta0 <- readr::read_tsv(meta_file, show_col_types = FALSE)

if (is.null(id_col_name)) id_col_name <- names(meta0)[1]

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s+", "", x)
  x <- sub("\\|.*$", "", x)
  x <- sub("(_ST.*$)$", "", x)
  x
}

tr <- tr0
tr$tip.label <- clean_id(tr$tip.label)

meta <- meta0 %>%
  mutate(
    label = clean_id(.data[[id_col_name]]),
    Group = trimws(as.character(.data[[group_col_name]])),
    Type  = tolower(trimws(as.character(.data[[type_col_name]]))),
    Study = trimws(as.character(.data[[study_col_name]]))
  ) %>%
  select(label, Group, Type, Study)

meta <- meta %>%
  mutate(
    Group = case_when(
      tolower(Group) %in% c("control", "con", "ctrl", "0", "false") ~ "Control",
      tolower(Group) %in% c("pdac", "case", "1", "true")            ~ "PDAC",
      TRUE ~ Group
    ),
    Type = case_when(
      Type %in% c("stool", "feces", "faeces") ~ "stool",
      Type %in% c("saliva", "oral")          ~ "saliva",
      TRUE ~ Type
    ),
    Study = ifelse(is.na(Study) | Study == "", NA_character_, Study)
  )

ann0 <- meta %>% filter(label %in% tr$tip.label)

keep <- ann0 %>%
  filter(!is.na(Group) & !is.na(Type) & !is.na(Study)) %>%
  pull(label) %>%
  unique()

tr2  <- ape::keep.tip(tr, keep)
ann2 <- ann0 %>% filter(label %in% keep)

group_col <- c(Control = "#336fac", PDAC = "#b82929")
type_col  <- c(saliva = "#5DB1DD", stool = "#C75127")

study_levels <- sort(unique(ann2$Study))
study_pal <- setNames(
  grDevices::hcl.colors(length(study_levels), palette = "Dark 3"),
  study_levels
)

tree_line <- 0.5
dot_size  <- 0.8

ring_w  <- 0.05
offset1 <- 0.01
offset2 <- offset1 + ring_w - 0.001

get_node_tips <- function(tr, node) {
  ape::extract.clade(tr, node)$tip.label
}

internal_nodes <- function(tr) {
  (Ntip(tr) + 1):(Ntip(tr) + tr$Nnode)
}

find_enriched_subtrees <- function(tr, ann, class_col = "Group",
                                   purity = 0.8, min_tips = 10,
                                   target_levels = NULL) {
  stopifnot("label" %in% names(ann))
  stopifnot(class_col %in% names(ann))
  
  cls <- setNames(as.character(ann[[class_col]]), ann$label)
  
  cand <- lapply(internal_nodes(tr), function(nd){
    tips <- get_node_tips(tr, nd)
    tips <- tips[tips %in% names(cls)]
    if (length(tips) < min_tips) return(NULL)
    
    tab <- table(cls[tips])
    if (length(tab) == 0) return(NULL)
    
    top_level <- names(which.max(tab))
    top_prop  <- max(tab) / sum(tab)
    
    if (!is.null(target_levels) && !(top_level %in% target_levels)) return(NULL)
    if (top_prop < purity) return(NULL)
    
    data.frame(
      node = nd,
      n_tips = length(tips),
      top_level = top_level,
      top_prop = top_prop,
      stringsAsFactors = FALSE
    )
  })
  cand <- bind_rows(cand)
  if (nrow(cand) == 0) return(list(nodes = integer(0), info = cand))
  
  cand <- cand[order(-cand$n_tips, -cand$top_prop), ]
  
  chosen <- integer(0)
  chosen_tips <- list()
  
  for (i in seq_len(nrow(cand))) {
    nd <- cand$node[i]
    tips <- get_node_tips(tr, nd)
    
    overlap <- FALSE
    for (j in seq_along(chosen)) {
      tj <- chosen_tips[[j]]
      if (length(intersect(tips, tj)) > 0) {
        overlap <- TRUE
        break
      }
    }
    if (overlap) next
    
    chosen <- c(chosen, nd)
    chosen_tips[[length(chosen_tips)+1]] <- tips
  }
  
  list(nodes = chosen, info = cand)
}

enr <- find_enriched_subtrees(
  tr2, ann2,
  class_col = "Group",
  purity = 0.8,
  min_tips = 10,
  target_levels = c("PDAC","Control")
)

ann2$Subclade <- "Other"
if (length(enr$nodes) > 0) {
  for (k in seq_along(enr$nodes)) {
    nd <- enr$nodes[k]
    tips_k <- get_node_tips(tr2, nd)
    top_group <- enr$info$top_level[enr$info$node == enr$nodes[k]]
    
    ann2$Subclade[ann2$label %in% tips_k] <-
      paste0(top_group, "_Clade_", k)
  }
}

subclade_levels <- c("Other", sort(unique(ann2$Subclade[ann2$Subclade!="Other"])))
subclade_pal <- setNames(
  c("grey85", grDevices::hcl.colors(length(subclade_levels)-1, "Set 3")),
  subclade_levels
)

p <- ggtree(tr2, layout = "circular") %<+% ann2 +
  geom_tree(size = tree_line, color = "grey20") +
  geom_tippoint(aes(color = Group), size = dot_size, alpha = 0.85) +
  scale_color_manual(values = group_col, name = "Group")

p1 <- p +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Type),
    width = ring_w, offset = offset1,
    colour = NA
  ) +
  scale_fill_manual(values = type_col, name = "Type")

p2 <- p1 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Study),
    width = ring_w, offset = offset2,
    colour = NA
  ) +
  scale_fill_manual(values = study_pal, name = "Study") +
  ggtitle(plot_title)

offset3 <- offset2 + ring_w - 0.002

p3 <- p2 +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = Subclade),
    width = ring_w, offset = offset3,
    colour = NA
  ) +
  scale_fill_manual(values = subclade_pal, name = "Subclade") +
  ggtitle("Fusicatenibacter saccharivorans") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "italic",
      size = 14
    )
  )

p3

ggsave(file = "/Users/yrx/Downloads/pancreatitis/Fusicatenibacter_saccharivorans.pdf", p3, width = 5, height = 8)
