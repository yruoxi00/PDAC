library(readxl)
library(dplyr)

setwd("/Users/yrx/Downloads/pancreatitis")

pancreatitis_raw <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "Pancreatitis")

pdac_raw <- read_excel("data/Meta_PA_microbiome_new.xlsx", sheet = "PDAC")

data <- bind_rows(pancreatitis_raw, pdac_raw) %>%
  mutate(
    Age = case_when(
      Data == "Dat9" & Age == ">70" ~ 75,
      Data == "Dat9" & Age == "<70" ~ 65,
      TRUE ~ suppressWarnings(as.numeric(Age))
    )
  )

pdac <- data %>%
  filter(Data %in% c("Dat2", "Dat5", "Dat6", "Dat8", "Dat9")) %>%
  mutate(
    Group_new = case_when(
      Group %in% c("PC", "PDAC") ~ "PDAC",
      Group %in% c("Con", "Control") ~ "Control",
      TRUE ~ NA_character_
    ),
    Type = case_when(
      Type %in% c("Stool", "stool") ~ "Stool",
      Type %in% c("Saliva", "saliva") ~ "Saliva",
      TRUE ~ Type
    )
  ) %>%
  filter(Group_new %in% c("PDAC", "Control"))

write.csv(pdac, "data/PDAC.csv", row.names = FALSE)
