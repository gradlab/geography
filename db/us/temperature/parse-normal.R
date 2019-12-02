#!/usr/bin/env Rscript

# Temperatures are annual averages for 1981-2010

station_col_positions = fwf_positions(
  start = c(1, 13, 22, 32, 39, 42, 73, 77, 81, 87),
  end = c(11, 20, 30, 37, 40, 71, 75, 79, 85, 99),
  col_names = c('station_id', 'lat', 'long', 'elevation', 'state_abb', 'name', 'gsn_flag', 'hcn_flag', 'wmo_id', 'method')
)

stations = read_fwf('temp-inventory.txt', col_positions = station_col_positions)

temp_col_positions = fwf_positions(
  start = c(1, 19, 24),
  end = c(11, 23, 24),
  col_names = c('station_id', 'value', 'completeness_flag')
)

raw = read_fwf('ann-tavg-normal.txt', col_positions = temp_col_positions)

temperature = raw %>%
  # value is 10ths of degrees F
  mutate(temperature = value / 10) %>%
  # include station info
  left_join(stations, by = 'station_id') %>%
  # get state names
  left_join(data_frame(state_abb = state.abb, state = state.name), by = 'state_abb') %>%
  # remove, e.g., AS = American Samoa
  filter(!is.na(state)) %>%
  group_by(state) %>%
  summarize(temperature = mean(temperature))

stopifnot(nrow(temperature) == 50)

write_tsv(temperature, 'state_temperature.tsv')
