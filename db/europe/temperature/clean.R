#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(countrycode)

raw <- read_tsv('cckp_historical_data_0.txt')

data <- raw %>%
  # left_join(codes) %>%
  mutate(country = countrycode(ISO_3DIGIT, 'iso3c', 'country.name', custom_match = c('KSV' = 'Kosovo'))) %>%
  select(country, temperature = Annual_temp) %>%
  bind_rows(tibble(country = 'Malta', temperature = 19.5))

write_tsv(data, '../temperature.tsv')
