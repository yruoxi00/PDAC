library(anpan)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(ape)
library(ggplot2)
library(egg)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/PDAC.csv")
pdac <- pdac %>%
  rename(sample_id = ID)
pdac <- pdac %>%
  mutate(
    PDAC = ifelse(Group_new == "PDAC", 1, 0)
  )
readr::write_tsv(pdac, "data/PDAC.tsv")

# Veillonella atypica
meta_file = "data/PDAC.tsv"
tree_file = "out6936/trees/RAxML_bestTree.t__SGB6936.StrainPhlAn4.tre"

anpan_pglmm(meta_file,
            tree_file,
            trim_pattern = NULL,
            covariates = NULL,
            outcome = "PDAC",
            verbose = TRUE,
            cores = 4)

# Veillonella parvula
meta_file = "data/PDAC.tsv"
tree_file = "out6939/trees/RAxML_bestTree.t__SGB6939.StrainPhlAn4.tre"

anpan_pglmm(meta_file,
            tree_file,
            trim_pattern = NULL,
            covariates = NULL,
            outcome = "PDAC",
            verbose = TRUE,
            cores = 4)

# Eubacterium rectale
meta_file = "data/PDAC.tsv"
tree_file = "out4933/trees/RAxML_bestTree.t__SGB4933.StrainPhlAn4.tre"

anpan_pglmm(meta_file,
            tree_file,
            trim_pattern = NULL,
            covariates = NULL,
            outcome = "PDAC",
            verbose = TRUE,
            cores = 4)

# Agathobaculum_butyriciproducens
meta_file = "data/PDAC.tsv"
tree_file = "out14993/trees/RAxML_bestTree.t__SGB14993.StrainPhlAn4.tre"

anpan_pglmm(meta_file,
            tree_file,
            trim_pattern = NULL,
            covariates = NULL,
            outcome = "PDAC",
            verbose = TRUE,
            cores = 4)
