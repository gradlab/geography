#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(countrycode)

raw <- read_tsv('tps00003.tsv')

data <- raw %>%
  select(`metric,code` = `unit,geo\\time`, density = `2015`) %>%
  separate(`metric,code`, c('metric', 'country_code'), sep = ',') %>%
  filter(
    # ":" encodes a missing value, weirdly
    density != ':',
    !(country_code %in% c('EU27', 'EU28'))
  ) %>%
  mutate(
    country = countrycode(country_code, 'iso2c', 'country.name', custom_match = c('UK' = 'United Kingdom', 'EL' = 'Greece'))
  )

stopifnot(all(data$metric == 'PER_KM2'))

data %>%
  select(country, density) %>%
  write_tsv('../density.tsv')
