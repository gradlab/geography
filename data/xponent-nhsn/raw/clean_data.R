#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

# Antibiotic use data ---------------------------------------------------------

use_data = read_csv('xponent/Data_OAU.csv') %>%
  filter(Antibiotic_Class == 'Fluoroquinolones', Location != 'National', between(Year, 2011, 2014)) %>%
  # average over years
  group_by(Location) %>%
  summarize(Rate = mean(as.numeric(Rate))) %>%
  select(
    state = Location,
    rx_1k_year = Rate # prescriptions per 1,000 people per year
  ) %>%
  filter(state %in% state.abb) # excludes Puerto Rico and DC

# Antibiotic resistance data --------------------------------------------------

resistance_data = read_csv('nhsn/PSA_States.csv') %>%
  filter(
    State %in% state.abb,
    Phenotype == 'E.coli R to fluoroquinolones',
    EventType == 'CAUTI',
    EventYear == 'All Years',
    AgeCategory == 'All Ages'
  ) %>%
  mutate_at(vars(starts_with('Number')), as.integer) %>%
  mutate_at(vars(PercentResistant, matches('95CI')), as.numeric) %>%
  select(state = State, n_isolates = NumberTested, n_resistant = NumberResistant)

# Combined data ---------------------------------------------------------------

# Check that the two data sets have the same number of rows (i.e., one per state)
stopifnot(nrow(use_data) == nrow(resistance_data))

# Combine the two data sets
combined_data = inner_join(use_data, resistance_data, by = 'state')

# Check that the result has the same number of rows
stopifnot(nrow(use_data) == nrow(combined_data))

# Save the combined data
write_tsv(combined_data, '../use_resistance_data.tsv')
