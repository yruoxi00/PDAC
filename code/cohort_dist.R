library(ggplot2)
library(webr)
library(dplyr)
library(readr)
library(cowplot)
library(grid)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/AP.csv")

stool_data  <- pdac %>% filter(Type == "Stool")
saliva_data <- pdac %>% filter(Type == "Saliva")

stool_counts  <- stool_data  %>% count(Study, Group_new, name = "n")
saliva_counts <- saliva_data %>% count(Study, Group_new, name = "n")

total_stool  <- sum(stool_counts$n)
total_saliva <- sum(saliva_counts$n)

## Stool
PieDonut(
  stool_counts,
  aes(Study, Group_new, count = n),
  title       = "Stool: AP vs Control by Cohort",
  showPieName = FALSE
)
grid.text(
  paste0("N = ", total_stool),
  x  = 0.5, y = 0.5,
  gp = gpar(fontface = "bold", cex = 1.2)
)

## Saliva
PieDonut(
  saliva_counts,
  aes(Study, Group_new, count = n),
  title       = "Saliva: AP vs Control by Cohort",
  showPieName = FALSE
)
grid.text(
  paste0("N = ", total_saliva),
  x  = 0.5, y = 0.5,
  gp = gpar(fontface = "bold", cex = 1.2)
)

png("fig/stool_piedonut_ap.pdf", res = 300)

PieDonut(
  stool_counts,
  aes(Study, Group_new, count = n),
  title       = "Stool: AP vs Control by Cohort",
  showPieName = FALSE
)
grid.text(
  paste0("N = ", total_stool),
  x  = 0.5, y = 0.5,
  gp = gpar(fontface = "bold", cex = 1.2)
)

dev.off()

png("fig/saliva_piedonut_ap.pdf", res = 300)

PieDonut(
  saliva_counts,
  aes(Study, Group_new, count = n),
  title       = "Saliva: AP vs Control by Cohort",
  showPieName = FALSE
)
grid.text(
  paste0("N = ", total_saliva),
  x  = 0.5, y = 0.5,
  gp = gpar(fontface = "bold", cex = 1.2)
)

dev.off()

pdac %>%
  count(Study, name = "N_total") %>%
  arrange(desc(N_total))
