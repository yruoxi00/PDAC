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

ap <- data %>%
  filter(Data %in% c("Dat1", "Dat3", "Dat4", "Dat5")) %>%
  mutate(
    Group_new = case_when(
      Group %in% c("AP", "Pancreatitis") ~ "AP",
      Group %in% c("Con", "Control") ~ "Control",
      TRUE ~ NA_character_
    ),
    Type = case_when(
      Type %in% c("Stool", "stool") ~ "Stool",
      Type %in% c("Saliva", "saliva") ~ "Saliva",
      TRUE ~ Type
    )
  ) %>%
  filter(Group_new %in% c("AP", "Control"))

write.csv(ap, "data/AP.csv", row.names = FALSE)
