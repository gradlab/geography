#!/usr/bin/env Rscript --vanilla

library(tidyverse)

raw <- read_csv('ilc_di15_1_Data.csv', col_types = cols(Value = 'n'))

data <- raw %>%
  # drop uninformative columns
  select_if(~ length(unique(.)) > 1) %>%
  filter(
    TIME == 2015,
    !str_detect(GEO, 'European Union'),
    !str_detect(GEO, 'Euro area')
  ) %>%
  mutate(GEO = recode(GEO, `Germany (until 1990 former territory of the FRG)` = 'Germany')) %>%
  arrange(GEO) %>%
  select(country = GEO, income = Value)

# check all the income values are good
stopifnot(all(!is.na(data$income)))

write_tsv(data, '../income.tsv')
