#!/usr/bin/env Rscript --vanilla

library(tidyverse)

# Mapping from DID to TPY
#  DID = DDD per 1,000 Inhabitants per Day
#  DDD = Defined Daily Dose
#  TPY = Preatment per Person per Year
did_tpy_map = data_frame(
  antibiotic = c('beta_lactam', 'quinolone', 'macrolide'),
  ddd_per_tx = c(10, 10, 7),
  tpy_per_did = 365 / (1e3 * ddd_per_tx)
)

# Load use/resistance data ---------------------------------------------

# Load Xponent/NHSN use/resistance data
xponent_nhsn_data <- read_tsv('data/xponent_nhsn/xponent_nhsn_data.tsv') %>%
  rename(use = rx_person_year) %>%
  mutate(
    pathogen = 'E. coli',
    antibiotic = 'quinolone',
    n_susceptible = n_isolates - n_resistant,
    f_resistant = n_resistant / n_isolates # proportion resistant
  )

# Load ECDC use/resistance data
ecdc_data <- read_tsv('data/ecdc/ecdc_data.tsv') %>%
  mutate(
    n_susceptible = n_isolates - n_resistant,
    f_resistant = n_resistant / n_isolates # proportion resistant
  ) %>%
  left_join(did_tpy_map, by = 'antibiotic') %>%
  mutate(use = did * tpy_per_did) # treatments per person per day

eu_units <- ecdc_data %>% pull(country) %>% unique()

# Load state/country metadata ------------------------------------------

# Load data about US states' temperature, income, density
us_temp <- read_tsv('db/us/temperature.tsv')
us_income <- read_tsv('db/us/income.tsv')
us_density <- read_tsv('db/us/density.tsv')

# Load data about European countries' temperature, income density
eu_temp <- read_tsv('db/europe/temperature.tsv')
eu_income <- read_tsv('db/europe/income.tsv')
eu_density <- read_tsv('db/europe/density.tsv')

# Load data about the adjacency of states/countries
us_adjacency <- read_tsv('db/us/adjacency.tsv') %>%
  mutate(adjacent = TRUE) %>%
  complete(state1, state2) %>%
  replace_na(list(adjacent = FALSE)) %>%
  rename(unit1 = state1, unit2 = state2)

eu_adjacency <- read_tsv('db/europe/adjacency.tsv') %>%
  mutate(adjacent = TRUE) %>%
  right_join(
    crossing(country1 = eu_units, country2 = eu_units),
    by = c('country1', 'country2')
  ) %>%
  replace_na(list(adjacent = FALSE)) %>%
  rename(unit1 = country1, unit2 = country2)

adjacency_db <- bind_rows(
  'US' = us_adjacency,
  'Europe' = eu_adjacency,
  .id = 'dataset'
)

# Load data about commuting
us_commuting <- read_tsv('db/us/commuting.tsv') %>%
  rename_all(~ str_replace(., '^state', 'unit'))

eu_commuting <- read_tsv('db/europe/commuting.tsv') %>%
  rename_all(~ str_replace(., '^country', 'unit')) %>%
  filter(unit1 %in% eu_units, unit2 %in% eu_units)

commuting_db <- bind_rows(
  'US' = us_commuting,
  'Europe' = eu_commuting,
  .id = 'dataset'
)


# Combine use/resistance and other state characteristics --------------

us_data <- xponent_nhsn_data %>%
  left_join(us_temp, by = 'state') %>%
  left_join(us_income, by = 'state') %>%
  left_join(us_density, by = 'state') %>%
  rename(unit = state)

eu_data <- ecdc_data %>%
  left_join(eu_temp, by = 'country') %>%
  left_join(eu_income, by = 'country') %>%
  left_join(eu_density, by = 'country') %>%
  rename(unit = country)

data <- bind_rows(
  'Xponent/NHSN' = us_data,
  'ECDC' = eu_data,
  .id = 'dataset'
) %>%
  nest(-dataset, -pathogen, -antibiotic, .key = 'unit_data') %>%
  # Make dataset/pathogen/antibiotic labels like "ECDC Ec/q"
  mutate(
    bug_drug = case_when(
      .$pathogen == 'E. coli' & .$antibiotic == 'quinolone' ~ 'Ec/q',
      .$pathogen == 'S. pneumoniae' & .$antibiotic == 'macrolide' ~ 'Sp/m',
      .$pathogen == 'S. pneumoniae' & .$antibiotic == 'beta_lactam' ~ 'Sp/bl'
    ),
    label = str_c(dataset, ' ', bug_drug)
  )

# Regional analysis ------------------------------------------------------------

# Aggregate the units (US state or European country) into regions. Aggregate use
# data by population-weighted means. Aggregate resistance data by adding up the
# numbers of susceptible and resistant isolates.

aggregate_at = function(unit_data, group_chr) {
  stopifnot(class(group_chr) == 'character')
  group = as.symbol(group_chr)

  group_by_at(unit_data, group_chr, .add = TRUE) %>%
    summarize(
      use = weighted.mean(use, w = population),
      n_resistant = sum(n_resistant),
      n_isolates = sum(n_isolates)
    ) %>%
    ungroup() %>%
    rename(group := !!group) %>%
    mutate(
      group_type = group_chr,
      n_susceptible = n_isolates - n_resistant
    )
}

# Function to run logistic regressions predicting resistance from use rates
model_f = function(df) glm(
  cbind(n_resistant, n_susceptible) ~ use,
  family = 'binomial', data = df
)

# Bootstrap the data, making new aggregate data and models with each resampling
aggregate_boot_f = function(unit_data, group_chr) {
  # For each bootstrapped dataset, run the logistic regression, and use the
  # regression to make predictions over the range of data. The 2.5% and 97.5%
  # percentiles of these predictions will be the 95% confidence intervals.
  pred_use = seq(min(unit_data$use), max(unit_data$use), length.out = pred_n)
  newdata = data_frame(use = pred_use)

  # Run the model on the base (unbootstrapped) data
  model = unit_data %>%
    aggregate_at(group_chr) %>%
    model_f()

  # Get the predictions from that model, which are the solid line in the plot
  model_pred = data_frame(
    use = pred_use,
    f_resistant = predict(model, newdata = newdata, type = 'response')
  )

  results = data_frame(i = 1:n_trials) %>%
    mutate(
      # Get bootstrapped data sets
      boot_unit_data = map(i, ~ sample_n(unit_data, size = nrow(unit_data), replace = TRUE)),
      boot_aggregate_data = map(boot_unit_data, ~ aggregate_at(., group_chr)),
      # Run the model and predict the outcomes
      model = map(boot_aggregate_data, model_f),
      slope = map_dbl(model, ~ coef(.)['use']),
      pred_res = map(model, ~ predict(., newdata = newdata, type = 'response')), # predicted resistance
      pred_df = map(pred_res, ~ data_frame(use = pred_use, f_resistant = .)) # data frame of those predictions
    )

  # At each point of the pred_n points pred_use, get the 95% CIs for the predictions
  pred_data = results %>%
    select(i, pred_df) %>% unnest() %>%
    group_by(use) %>%
    summarize(
      ymin = quantile(f_resistant, alpha / 2),
      ymax = quantile(f_resistant, 1 - alpha / 2)
    ) %>%
    left_join(model_pred, by = 'use')

  list(
    slope_cil = quantile(results$slope, alpha / 2),
    slope_ciu = quantile(results$slope, 1 - alpha / 2),
    pred_data = pred_data
  )
}

data_region_combinations = data_frame(
  dataset = c(rep('Xponent/NHSN', 3), rep('ECDC', 2)),
  level = c('unit', 'us_division', 'us_region', 'unit', 'eu_region')
)

regional_results = left_join(data_region_combinations, data, by = 'dataset') %>%
  mutate(
    # Make a nice display for level across US and Europe
    display_level = factor(
      recode(level, unit = 'Unit', us_division = 'Division', us_region = 'Region', eu_region = 'Region'),
      levels = c('Unit', 'Division', 'Region')
    ),
    aggregated_data = map2(unit_data, level, aggregate_at),
    # Fit models to each aggregate type
    model = map(aggregated_data, model_f),
    # Extract the use/resistance regression coefficient ("slope") from the base model
    slope = map_dbl(model, ~ coef(.)['use']),
    # Bootstrap the 95% CIs for that regression coefficient
    boot_results = map2(unit_data, level, aggregate_boot_f),
    slope_cil = map_dbl(boot_results, ~ .$slope_cil),
    slope_ciu = map_dbl(boot_results, ~ .$slope_ciu),
    pred_data = map(boot_results, ~ .$pred_data)
  )


# Show regional plots ----------------------------------------------------------

my_palette = c('black', '#1b9e77', '#d95f02')

# Pull out the data and the predictions separately
aggregate_plot_data = regional_results %>%
  select(label, display_level, aggregated_data) %>%
  unnest() %>%
  mutate(f_resistant = n_resistant / n_isolates)

aggregate_pred_data = regional_results %>%
  select(label, display_level, pred_data) %>%
  unnest()

supp_figure_3 = ggplot(
    data = NULL,
    aes(use, f_resistant, color = display_level, fill = display_level)
  ) +
  facet_wrap(~ label, scales = 'free') +
  geom_ribbon(
    data = aggregate_pred_data,
    aes(ymin = ymin, ymax = ymax, color = NULL),
    alpha = 0.25
  ) +
  geom_line(data = aggregate_pred_data, size = 1) +
  geom_point(data = aggregate_plot_data) +
  xlab('antibiotic use (annual treatments per capita)') +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  scale_y_continuous(
    'antibiotic resistance (% of isolates)',
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

ggsave('fig/supp_figure_3.pdf', plot = supp_figure_3)

supp_figure_4 = regional_results %>%
  ggplot(aes(display_level, slope, ymin = slope_cil, ymax = slope_ciu, color = display_level)) +
  facet_wrap(~ label, scales = 'free') +
  geom_point() +
  geom_errorbar() +
  scale_color_manual(values = my_palette) +
  guides(color = 'none') +
  xlab('') +
  ylab('Regression coefficient')

ggsave('fig/supp_figure_4.pdf', plot = supp_figure_4)


# Distance analysis -------------------------------------------------------------

odds = function(p) p / (1 - p)
log_odds_ratio = function(p, q) log(odds(p)) - log(odds(q))

# Function to get information about pairs of states, to be used multiple times
# when bootstrapping lists of states.
cross_units = function(unit_data, adjacency_data) {
  unit_data %$%
    crossing(unit1 = unit, unit2 = unit) %>%
    left_join(rename_all(unit_data, ~ str_c(., '1')), by = 'unit1') %>%
    left_join(rename_all(unit_data, ~ str_c(., '2')), by = 'unit2') %>%
    filter(use1 > use2) %>%
    mutate(
      d_resistant = f_resistant1 - f_resistant2,
      d_use = log_odds_ratio(use1, use2),
      d_income = income1 - income2,
      d_density = density1 - density2,
      d_temperature = temperature1 - temperature2
    ) %>%
    select(unit1, unit2, starts_with('d_')) %>%
    left_join(adjacency_data, by = c('unit1', 'unit2')) %>%
    # pairs of states not in the adjacency list are not adjacent
    replace_na(list(adjacent = 0))
}

# Function to do Tukey's bisquare regression
bisquare = partial(MASS::rlm, psi = MASS::psi.bisquare)

# Function to predict difference in resistant from difference in use (and other
# factors) as per equation in Methods, "Use-resistance relationships by
# adjacency"
distance_model_f = function(df) bisquare(d_resistant ~ d_use + d_use:adjacent + d_income + d_density + d_temperature, data = df)

distance_boot_f = function(unit_data, adjacency_data) {
  pair_data = cross_units(unit_data, adjacency_data)
  base_distance_model = distance_model_f(pair_data)

  # Bootstrap lists of states and run the distance model on those bootstrapped data sets
  base_d_use = pair_data$d_use
  pred_d_use = seq(min(base_d_use), max(base_d_use), length.out = pred_n)
  newdata = crossing(
    d_use = pred_d_use,
    adjacent = c(0, 1),
    d_income = mean(pair_data$d_income),
    d_density = mean(pair_data$d_density),
    d_temperature = mean(pair_data$d_temperature)
  )

  results = data_frame(i = 1:n_trials) %>%
    mutate(
      # Get bootstrapped data sets
      boot_unit_data = map(i, ~ sample_n(unit_data, size = nrow(unit_data), replace = TRUE)),
      boot_pair_data = map(boot_unit_data, ~ cross_units(., adjacency_data)),
      # Run the model and predict the outcomes
      model = map(boot_pair_data, distance_model_f),
      inter = map_dbl(model, ~ coef(.)['d_use:adjacent']), # interaction term
      pred_y = map(model, ~ predict(., newdata = newdata, type = 'response')),
      pred_df = map(pred_y, ~ mutate(newdata, d_resistant = .)) # data frame of predictions
    )

  # Get the base (non-bootstrapped) model predictions
  base_pred_df = newdata %>%
    mutate(d_resistant = predict(base_distance_model, newdata = newdata, type = 'response')) %>%
    select(d_use, adjacent, d_resistant)

  # At each point pred_use, get the CIs for the predictions
  pred_data = results %>%
    select(i, pred_df) %>% unnest() %>%
    group_by(d_use, adjacent) %>%
    summarize(
      ymin = quantile(d_resistant, alpha / 2),
      ymax = quantile(d_resistant, 1 - alpha / 2)
    ) %>%
    ungroup() %>%
    # Merge in the base model predictions to use as the solid line
    left_join(base_pred_df, by = c('d_use', 'adjacent'))

  list(
    inter_cil = quantile(results$inter, alpha / 2),
    inter_ciu = quantile(results$inter, 1 - alpha / 2),
    pred_data = pred_data
  )
}

cat('Note: Some robust regressions may fail during the bootstrapped runs. This is not abnormal.\n')

distance_results = data %>%
  mutate(
    pair_data = map2(unit_data, adjacency_data, cross_units),
    model = map(pair_data, distance_model_f),
    # interaction term
    inter = map_dbl(model, ~ coef(.)['d_use:adjacent']),
    # bootstrap
    boot_results = map2(unit_data, adjacency_data, distance_boot_f),
    inter_cil = map_dbl(boot_results, ~ .$inter_cil),
    inter_ciu = map_dbl(boot_results, ~ .$inter_ciu),
    pred_data = map(boot_results, ~.$pred_data)
  )


## Plot distance results -------------------------------------------------------

distance_plot_data = distance_results %>%
  select(label, pair_data) %>%
  unnest()

distance_pred_data = distance_results %>%
  select(label, pred_data) %>%
  unnest()

supp_figure_5 = ggplot(
    data = NULL,
    aes(d_use, d_resistant, color = as.logical(adjacent), fill = as.logical(adjacent))
  ) +
  facet_wrap(~ label, scales = 'free') +
  geom_ribbon(
    data = distance_pred_data,
    aes(ymin = ymin, ymax = ymax, color = NULL),
    alpha = 0.25
  ) +
  geom_point(data = distance_plot_data, shape = 1) +
  geom_line(data = distance_pred_data) +
  scale_color_manual(values = c('black', 'red')) +
  scale_fill_manual(values = c('black', 'red')) +
  xlab('Difference in antibiotic use (annual treatments per capita)') +
  ylab('Difference in resistance (log odds ratio)')

ggsave('fig/supp_figure_5.pdf', plot = supp_figure_5)


distance_results %>%
  select(dataset, pathogen, antibiotic, inter, inter_ciu, inter_cil) %>%
  mutate_if(is.numeric, ~ round(., 2)) %>%
  write_tsv('fig/supp_table_2.tsv')
