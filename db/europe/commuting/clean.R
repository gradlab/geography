#!/usr/bin/env Rscript --vanilla

library(tidyverse)
library(countrycode)

raw <- read_tsv('table-64.tsv')

data <- raw %>%
  gather('country_of_work', 'n_commuters', -country_of_residence) %>%
  filter(country_of_residence != country_of_work) %>%
  mutate(
    n_commuters = as.numeric(case_when(
      # "-1" for HR/NL is, I think, just 1
      .$n_commuters == '-1' ~ '1',
      # "." will be approximated as zero
      .$n_commuters == '.' ~ '0',
      # strip the parens
      str_detect(.$n_commuters, '^\\(') ~ str_match(.$n_commuters, '^\\((\\d+)\\)$')[, 2],
      # NA is just zero
      is.na(.$n_commuters) ~ '0',
      TRUE ~ .$n_commuters
    ))
  ) %>%
  mutate_at(vars(matches('country')), ~ countrycode(., 'iso2c', 'country.name', custom_match = c('UK' = 'United Kingdom'))) %>%
  mutate_at('n_commuters', ~ . * 1e3)

# I compared the results from these row sums against the "Total" column in the
# original table. The results are so-so, but as noted in the table footer,
# they aren't expected to be very good.

# data %>%
#   group_by(country_of_residence) %>%
#   summarize_at('n_commuters', sum) %>%
#   print(n = Inf)

# Read in the population data, and get commuting fraction, as a fraction of total population
workers <- read_tsv('employment/workers.tsv')

commuting <- data %>%
  mutate(
    country1 = pmin(country_of_residence, country_of_work),
    country2 = pmax(country_of_residence, country_of_work)
  ) %>%
  group_by(country1, country2) %>%
  summarize_at('n_commuters', sum) %>%
  ungroup() %>%
  left_join(rename(workers, within1 = n_workers), by = c('country1' = 'country')) %>%
  left_join(rename(workers, within2 = n_workers), by = c('country2' = 'country')) %>%
  mutate(f_commuting = n_commuters / (within1 + within2)) %>%
  select(country1, country2, f_commuting) %>%
  # symmetrize
  bind_rows(rename(., country1 = country2, country2 = country1))

write_tsv(commuting, '../commuting.tsv')
