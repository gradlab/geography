#!/usr/bin/env Rscript

units <- read_tsv('../../europe_units.tsv') %$%
  unit

codes <- read_csv('COW country codes.csv') %>%
  select(country = StateNme, code = CCode) %>%
  distinct()

pairs <- read_csv('DirectContiguity320/contdird.csv') %>%
  select(code1 = state1no, code2 = state2no, year, contiguity_type = conttype)

adjacency <- pairs %>%
  left_join(codes, by = c('code1' = 'code')) %>% rename(unit1 = country) %>%
  left_join(codes, by = c('code2' = 'code')) %>% rename(unit2 = country) %>%
  # require land/river contiguity (not large water bodies)
  filter(contiguity_type == 1) %>%
  filter(between(year, 2011, 2015)) %>%
  select(year, unit1, unit2)

# make sure we got all the units, except for the islands
my_units <- adjacency %$%
  c(unit1, unit2) %>%
  unique()

missing_units <- setdiff(units, my_units)
stopifnot(setequal(missing_units, c('Iceland', 'Malta')))

adjacency %>%
  filter(unit1 %in% units, unit2 %in% units) %>%
  select(-year) %>%
  distinct() %>%
  write_tsv('../../europe_adjacency.tsv')
