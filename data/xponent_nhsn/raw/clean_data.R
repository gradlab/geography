#!/usr/bin/env Rscript --vanilla

library(tidyverse)

# Antibiotic use data ---------------------------------------------------------

use_data <- read_csv('xponent/Data_OAU.csv') %>%
  filter(
    Antibiotic_Class == 'Fluoroquinolones',
    Location %in% state.abb, # excludes Puerto Rico and DC
    between(Year, 2011, 2014)
  ) %>%
  # average over years
  group_by(Location) %>%
  summarize(Rate = mean(as.numeric(Rate))) %>%
  # "Rate" is prescriptions per 1,000 people per year, we want prescriptions per person per year
  mutate(rx_person_year = Rate / 1000) %>%
  select(state = Location, rx_person_year)

# Antibiotic resistance data --------------------------------------------------

resistance_data <- read_csv('nhsn/PSA_States.csv') %>%
  filter(
    State %in% state.abb,
    Phenotype == 'E.coli R to fluoroquinolones',
    EventType == 'CAUTI',
    EventYear == 'All Years',
    AgeCategory == 'All Ages'
  ) %>%
  select(state = State, n_isolates = NumberTested, n_resistant = NumberResistant) %>%
  mutate_at(c('n_isolates', 'n_resistant'), as.integer)

# Combined data ---------------------------------------------------------------

# Check that the two data sets have the same number of rows (i.e., one per state)
stopifnot(nrow(use_data) == nrow(resistance_data))

# Combine the two data sets
combined_data <- inner_join(use_data, resistance_data, by = 'state') %>%
  mutate(state = state.name[match(state, state.abb)])

# Check that the result has the same number of rows
stopifnot(nrow(use_data) == nrow(combined_data))

# Save the combined data
write_tsv(combined_data, '../xponent_nhsn_data.tsv')
