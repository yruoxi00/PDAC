library(readxl)
library(ggplot2)
library(webr)
library(dplyr)
library(readr)
library(cowplot)
library(grid)

setwd("/Users/yrx/Downloads/pancreatitis")

pdac <- read_csv("data/PDAC.csv")

stool_data  <- pdac %>% filter(Type == "Stool")
saliva_data <- pdac %>% filter(Type == "Saliva")

stool_counts  <- stool_data  %>% count(Study, Group_new, name = "n")
saliva_counts <- saliva_data %>% count(Study, Group_new, name = "n")

total_stool  <- sum(stool_counts$n)
total_saliva <- sum(saliva_counts$n)

species <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Species")

paired_ids <- intersect(stool_data$ID, saliva_data$ID)

num_paired <- length(paired_ids)

print(paste("Number of individuals with paired oral-gut samples:", num_paired))
