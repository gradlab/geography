#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(countrycode)

raw <- read_tsv('estat_tps00001_filtered.tsv')

data <- raw %>%
  separate(`freq,indic_de,geo\\TIME_PERIOD`, c('freq', 'indic', 'geo'), sep = ',') %>%
  select(geo, value = `2015`) %>%
  filter(value != ':') %>%
  mutate(value = as.numeric(str_replace(value, ' (e|b)$', '')))

write_tsv('population.tsv')
