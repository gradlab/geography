#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(readr)

# Antibiotic use data ---------------------------------------------------------

use_data = read_tsv('esac.tsv') %>%
  # get only the years, setting, and antibiotics of interest
  mutate(antibiotic = recode(abx_atc,
    J01C = 'beta_lactam',
    J01D = 'beta_lactam',
    J01M = 'quinolone',
    J01FA = 'macrolide'
  )) %>%
  filter(
    between(year, 2011, 2015),
    toc == 'AC',
    antibiotic %in% c('beta_lactam', 'quinolone', 'macrolide')
  ) %>%
  # add together the two types of beta-lactam use
  group_by(year, country, antibiotic) %>%
  summarize(did = sum(did)) %>%
  ungroup() %>%
  # average over the years
  group_by(country, antibiotic) %>%
  summarize(did = mean(did)) %>%
  ungroup()

# Antibiotic resistance data --------------------------------------------------

resistance_data = read_csv('ECDC_surveillance_data_Antimicrobial_resistance.csv') %>%
  separate(Population, c('pathogen', 'antibiotic'), '\\|') %>%
  select(
    pathogen,
    antibiotic,
    metric = Indicator,
    year = Time,
    country = RegionName,
    value = NumValue
  ) %>%
  # keep only the years, pathogen/antibiotic combinations, and metrics of interest
  mutate(
    pathogen = recode(pathogen,
      `Escherichia coli` = 'E. coli',
      `Streptococcus pneumoniae` = 'S. pneumoniae'
    ),
    antibiotic = recode(antibiotic,
      Fluoroquinolones = 'quinolone',
      Macrolides = 'macrolide',
      Penicillins = 'beta_lactam'
    ),
    metric = recode(metric,
      `Non-susceptible (I and R) isolates` = 'n_nonsusceptible',
      `Total tested isolates` = 'n_isolates'
    )
  ) %>%
  filter(
    between(year, 2011, 2015),
    pathogen %in% c('E. coli', 'S. pneumoniae'),
    antibiotic %in% c('quinolone', 'macrolide', 'beta_lactam'),
    metric %in% c('n_nonsusceptible', 'n_isolates')
  ) %>%
  mutate(value = as.integer(value)) %>%
  spread(metric, value) %>%
  filter(n_isolates > 0) %>%
  # add up numerator and denominator isolates over years
  group_by(country, pathogen, antibiotic) %>%
  summarize(
    n_resistant = sum(n_nonsusceptible),
    n_isolates = sum(n_isolates)
  ) %>%
  ungroup()

# Combined data ---------------------------------------------------------------

# Combine the two data sets
combined_data = inner_join(use_data, resistance_data, by = c('country', 'antibiotic')) %>%
  rename(unit = country)

# Save the combined data
write_tsv(combined_data, '../ecdc_data.tsv')
