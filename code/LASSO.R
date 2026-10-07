library(readxl)
library(dplyr)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(glmnet)
library(pROC)
library(stringr)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/PDAC.csv")

### stool 

stool_data <- pdac %>%
  filter(Type == "Stool")

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

# species data 

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
cat("Original number of taxa:", ncol(species), "\n")
species <- species[, taxa]
cat("Filtered number of taxa:", ncol(species), "\n")

X <- as.matrix(species)
X <- log10(X + 1e-6)

y <- factor(stool_data$Group_new, levels = c("Control", "PDAC"))

set.seed(123)

cvfit <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",   
  alpha = 1,             
  nfolds = 5,            
  type.measure = "deviance"
)

coef_lasso <- coef(cvfit, s = "lambda.1se")

coef_lasso <- as.matrix(coef_lasso)

selected_taxa <- rownames(coef_lasso)[coef_lasso[,1] != 0]
selected_taxa <- setdiff(selected_taxa, "(Intercept)")

length(selected_taxa)
selected_taxa

lasso_df <- data.frame(
  Species = selected_taxa,
  Coef    = coef_lasso[selected_taxa, 1]
) %>%
  arrange(desc(abs(Coef)))

pred_prob <- predict(cvfit, newx = X, s = "lambda.1se", type = "response")

roc_obj <- roc(as.numeric(y), as.numeric(pred_prob))
auc(roc_obj)

lasso_plot_df <- lasso_df %>%
  mutate(
    Species_clean = Species %>%
      str_replace("^s__", "") %>%
      str_replace_all("_", " ")
  ) %>%
  
  filter(
    !str_detect(Species_clean, "GGB|SGB| sp\\.| sp |unclassified|unknown")
  ) %>%
  
  mutate(
    Direction = ifelse(Coef > 0, "Enriched in PDAC", "Enriched in Control"),
    absCoef   = abs(Coef)
  ) %>%
  
  arrange(desc(absCoef)) %>%
  
  mutate(
    Species_clean = factor(Species_clean, levels = rev(Species_clean))
  )

p <- ggplot(lasso_plot_df,
            aes(x = Coef, y = Species_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(aes(color = Direction), size = 2.6) +
  labs(
    x = "LASSO coefficient (log10 abundance)",
    y = NULL,
    title = "LASSO-selected taxa"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(face = "italic"),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )
print(p)

ggsave(filename = "LASSO.pdf", p, width = 7.5, height = 7)
