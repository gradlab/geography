#!/usr/bin/env Rscript

read_csv('ACS_15_5YR_S1903_with_ann.csv', skip = 1) %>%
  select(state = `Geography`, income = `Median income (dollars); Estimate; Households`) %>%
  filter(state %in% state.name) %>%
  arrange(state) %>%
  write_tsv('state_income.tsv')
