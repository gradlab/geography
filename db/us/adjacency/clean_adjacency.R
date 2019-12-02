#!/usr/bin/env Rscript

raw_adjacency = read_tsv(
  'county_adjacency.txt',
  col_names = c('county1', 'state_county_fips1', 'county2', 'state_county_fips2')
)

raw_fips = read_csv(
  'national_county.txt',
  col_names = c('state_abb', 'state_fips', 'county_fips', 'county', 'fips_class_code')
)

fips = raw_fips %>%
  select(state_abb, state_fips) %>%
  distinct()

recode_many <- function(x, from_values, to_values) {
  to_values[match(x, from_values)]
}

fips_to_state = function(x) {
  x %>%
    recode_many(fips$state_fips, fips$state_abb) %>%
    recode_many(state.abb, state.name)
}

adjacency = raw_adjacency %>%
  mutate(
    fips1 = zoo::na.locf(str_sub(state_county_fips1, 1, 2)),
    fips2 = str_sub(state_county_fips2, 1, 2)
  ) %>%
  select(fips1, fips2) %>%
  distinct() %>%
  mutate(
    state1 = fips_to_state(fips1),
    state2 = fips_to_state(fips2),
  ) %>%
  select(state1, state2) %>%
  filter(!is.na(state1), !is.na(state2))

# there should be 50 states in the "from" column
stopifnot(length(unique(adjacency$state1)) == 50)

# the result should be symmetrical
stopifnot(all_equal(adjacency, rename(adjacency, state1 = state2, state2 = state1)))

write_tsv(adjacency, '../adjacency.tsv')
