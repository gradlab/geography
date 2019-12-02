#!/usr/bin/env Rscript --vanilla

library(tidyverse)

raw <- read_csv('lfsi_emp_a_1_Data.csv', col_types = cols(Value = 'n'))

data <- raw %>%
  filter(
    INDIC_EM == 'Active population',
    !str_detect(GEO, 'Euro area'),
    !str_detect(GEO, 'European Union'),
    GEO != 'France'
  ) %>%
  mutate(
    GEO = recode(GEO,
      `Germany (until 1990 former territory of the FRG)` = 'Germany',
      # exclude "France", which includes Guadaloupe, etc. and keep only mainland
      `France (metropolitan)` = 'France',
    ),
    Value = Value * 1e3
  ) %>%
  select_if(~ length(unique(.)) > 1) %>%
  group_by(GEO) %>%
  summarize_at('Value', sum) %>%
  select(country = GEO, n_workers = Value)

write_tsv(data, 'workers.tsv')
