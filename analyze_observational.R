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

adjacency_db <- bind_rows(us_adjacency, eu_adjacency)

# Load data about commuting
us_commuting <- read_tsv('db/us/commuting.tsv') %>%
  rename_all(~ str_replace(., '^state', 'unit'))

eu_commuting <- read_tsv('db/europe/commuting.tsv') %>%
  rename_all(~ str_replace(., '^country', 'unit')) %>%
  filter(unit1 %in% eu_units, unit2 %in% eu_units)

commuting_db <- bind_rows(us_commuting, eu_commuting)


# Combine use/resistance and other state characteristics --------------

us_data <- xponent_nhsn_data %>%
  left_join(us_temp, by = 'state') %>%
  left_join(us_income, by = 'state') %>%
  left_join(us_density, by = 'state') %>%
  select(
    unit = state, pathogen, antibiotic, use, f_resistant,
    density, income, temperature
  )

eu_data <- ecdc_data %>%
  left_join(eu_temp, by = 'country') %>%
  left_join(eu_income, by = 'country') %>%
  left_join(eu_density, by = 'country') %>%
  select(
    unit = country, pathogen, antibiotic, use, f_resistant,
    density, income, temperature
  )

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


# Plot the use/resistance data ----------------------------------------

round_up <- function(x, digits) ceiling(x * 10 ** digits) / 10 ** digits
round_down <- function(x, digits) floor(x * 10 ** digits) / 10 ** digits

breaker <- function(digits) {
  function(x) {
    lower <- round_up(x[1], digits)
    upper <- round_down(x[2], digits)
    round(seq(lower, upper, length.out = 3), digits)
  }
}

figure2 <- data %>%
  unnest() %>%
  ggplot(aes(use, f_resistant * 100)) +
  facet_wrap(~ dataset + bug_drug, scales = 'free') +
  geom_smooth(method = 'lm', color = 'gray50') +
  geom_point() +
  scale_x_continuous(
    'Treatments per person per year',
    breaks = breaker(2)
  ) +
  scale_y_continuous(
    name = 'Resistance (% isolates nonsusceptible)'
  )

ggsave('fig/figure_2.pdf', plot = figure2)

# Set up "crossed" data -----------------------------------------------

odds <- function(p) p / (1 - p)
log_odds_ratio <- function(p, q) log(odds(p)) - log(odds(q))

cross_units <- function(df) {
  df %>%
    with(crossing(unit1 = unit, unit2 = unit)) %>%
    filter(unit1 < unit2) %>%
    left_join(rename_all(df, ~ str_c(., '1')), by = 'unit1') %>%
    left_join(rename_all(df, ~ str_c(., '2')), by = 'unit2') %>%
    mutate(
      d_resistant = f_resistant1 - f_resistant2,
      d_use = use1 - use2,
      d_income = income1 - income2,
      d_density = density1 - density2,
      d_temperature = temperature1 - temperature2,
      dr_du = d_resistant / d_use,
      lor_resistant = log_odds_ratio(f_resistant1, f_resistant2),
      lorr_du = lor_resistant / d_use
    ) %>%
    select(unit1, unit2, starts_with('d_'), dr_du, starts_with('lor')) %>%
    left_join(adjacency_db, by = c('unit1', 'unit2')) %>%
    left_join(commuting_db, by = c('unit1', 'unit2'))
}

leave_one_out_from_cross <- function(df) {
  units <- unique(c(df$unit1, df$unit2))
  map(units, ~ filter(df, unit1 != ., unit2 != .))
}

cross_data <- data %>%
  mutate(
    cross_data = map(unit_data, cross_units),
    l1o_cross_data = map(cross_data, leave_one_out_from_cross)
  )

# Show the commuting histogram ----------------------------------------

commuting_histogram <- cross_data %>%
  select(dataset, cross_data) %>%
  unnest() %>%
  select(dataset, unit1, unit2, adjacent, f_commuting) %>%
  filter(
    unit1 < unit2,
    !is.na(f_commuting)
  ) %>%
  ggplot(aes(log10(f_commuting), fill = adjacent)) +
  facet_wrap(~ dataset, scales = 'free') +
  geom_histogram(color = 'black') +
  scale_fill_manual(
    limits = c('TRUE', 'FALSE'),
    labels = c('adjacent', 'not adj.'),
    values = c('black', 'white')
  ) +
  xlab('Commuting fraction (log 10)') +
  ylab('Number of states or countries') +
  theme(
    legend.position = c(0.85, 0.75),
    legend.title = element_blank()
  )

ggsave('fig/supplemental_figure_1.pdf', plot = commuting_histogram)


# Do the pairwise analyses --------------------------------------------

jackknife.sd <- function(x) {
  n <- length(x)
  sqrt((n - 1) / n * sum((x - mean(x)) ** 2))
}

analysis_f <- function(model_f, coef_f, ratio_f) {
  base_results <- cross_data %>%
    mutate(
      model = map(cross_data, model_f),
      coef = map_dbl(model, coef_f),
      ratio = map_dbl(model, ratio_f)
    ) %>%
    select(dataset, coef, ratio)

  l1o_results <- cross_data %>%
    select(dataset, l1o_cross_data) %>%
    unnest() %>%
    mutate(
      model = map(l1o_cross_data, model_f),
      coef = map_dbl(model, coef_f),
      ratio = map_dbl(model, ratio_f)
    ) %>%
    group_by(dataset) %>%
    summarize(
      coef_se = jackknife.sd(coef),
      ratio_se = jackknife.sd(ratio)
    )

  base_results %>%
    left_join(l1o_results, by = 'dataset')
}

rlm <- MASS::rlm

# First, for simple adjacency
adjacency_results <- analysis_f(
  function(df) with(df, {
    list(
      adj_med = median(df$dr_du[df$adjacent]),
      nonadj_med = median(df$dr_du[!df$adjacent])
    )
  }),
  function(model) with(model, { adj_med - nonadj_med }),
  function(model) with(model, { (adj_med - nonadj_med) / nonadj_med })
)

# Adjacency, replacing dr/du with LOR(r)/du
adjacency_lorru_results <- analysis_f(
  function(df) with(df, {
    list(
      adj_med = median(df$lorr_du[df$adjacent]),
      nonadj_med = median(df$lorr_du[!df$adjacent])
    )
  }),
  function(model) with(model, { adj_med - nonadj_med }),
  function(model) with(model, { (adj_med - nonadj_med) / nonadj_med })
)

# Adjancency, using a robust regression
adjacency_rlm_results <- analysis_f(
  function(df) rlm(dr_du ~ adjacent, data = df),
  function(model) coef(model)['adjacentTRUE'],
  function(model) coef(model) %>% { .['adjacentTRUE'] / .['(Intercept)'] }
)

# Using a robust regression with covariates
adjacency_covariates_results <- analysis_f(
  function(df) rlm(dr_du ~ adjacent + d_income + d_temperature + d_density, data = df),
  function(model) coef(model)['adjacentTRUE'],
  function(model) coef(model) %>% { .['adjacentTRUE'] / .['(Intercept)'] }
)

# Analysis of commuting
commuting_results <- analysis_f(
  function(df) rlm(dr_du ~ f_commuting, data = df),
  function(model) coef(model)['f_commuting'] * 1e-4,
  function(model) coef(model) %>% { .['f_commuting'] / .['(Intercept)'] * 1e-4 }
)

# Print out the tables ------------------------------------------------

sigfig <- function(x, n = 2) {
  formatC(signif(x, digits = n), digits = n, format = "fg", flag = "#")
}

save_results <- function(df, fn) {
  df %>%
    mutate(
      coef_hci = 1.96 * coef_se * 0.5,
      coef_cil = coef - coef_hci,
      coef_ciu = coef + coef_hci,
      coef_star = coef_cil > 0 | coef_ciu < 0,
      ratio_hci = 1.96 * ratio_se * 0.5,
      ratio_cil = ratio - ratio_hci,
      ratio_ciu = ratio + ratio_hci,
      ratio_star = ratio_cil > 0 | ratio_ciu < 0
    ) %>%
    mutate_if(is.numeric, sigfig) %>%
    mutate_at(vars(ends_with('star')), ~ recode(as.numeric(.), `1` = '*', `0` = '')) %>%
    mutate(
      coef_display = str_glue('{coef} ({coef_cil} to {coef_ciu}){coef_star}'),
      ratio_display = str_glue('{ratio} ({ratio_cil} to {ratio_ciu}){ratio_star}')
    ) %>%
    select(dataset, coef_display, ratio_display) %>%
    write_tsv(fn)
}

save_results(
  adjacency_results,
  'fig/supplemental_table_2.tsv'
)

save_results(
  adjacency_lorru_results,
  'fig/supplemental_table_3.tsv'
)

save_results(
  adjacency_rlm_results,
  'fig/supplemental_table_4.tsv'
)

show_results(
  commuting_results,
  'fig/supplemental_table_5.tsv'
)

# Plots ---------------------------------------------------------------

boxplot_data_f <- function(df, ymin, ymax) {
  df %>%
    nest(-adjacent) %>%
    mutate(
      y = map(data, ~ .$dr_du),
      box = map(y, ~ boxplot.stats(.)$stats),
      boxplot_data = map(box, ~ tibble(
        ymin = max(.[1], ymin),
        lower = .[2],
        middle = .[3],
        upper = .[4],
        ymax = min(.[5], ymax)
      ))
    ) %>%
    select(adjacent, boxplot_data) %>%
    unnest()
}

boxplot_f <- function(cross_data, f_to_keep) {
  half_drop <- (1 - f_to_keep) / 2

  plot_data <- cross_data %>%
    mutate(
      y = map(cross_data, ~ .$dr_du),
      ymin = map_dbl(y, ~ quantile(., half_drop)),
      ymax = map_dbl(y, ~ quantile(., 1 - half_drop)),
      boxplot_data = pmap(list(cross_data, ymin, ymax), boxplot_data_f),
      point_data = pmap(list(cross_data, ymin, ymax), ~ filter(..1, between(dr_du, ..2, ..3)))
    )

  point_data <- plot_data %>%
    select(dataset, point_data) %>%
    unnest()

  boxplot_data <- plot_data %>%
    select(dataset, boxplot_data) %>%
    unnest()

  plot <- ggplot(data = NULL, aes(x = factor(adjacent))) +
    facet_wrap(~ dataset, scales = 'free_y') +
    geom_boxplot(
      data = boxplot_data,
      aes(lower = lower, upper = upper, middle = middle, ymin = ymin, ymax = ymax),
      stat = 'identity'
    ) +
    geom_jitter(data = point_data, aes(y = dr_du), size = 0.1, width = 0.2) +
    scale_x_discrete(
      '',
      labels = c(`TRUE` = 'Adjacent', `FALSE` = 'Not adj.')
    ) +
    ylab(expression(paste('Use-resistance association ', (Delta * rho / Delta * tau)))) +
    theme_cowplot() +
    theme(strip.background = element_blank())

  plot
}

adjacency_plot <- boxplot_f(cross_data, 0.90)
ggsave('fig/figure_3.pdf', plot = adjacency_plot)


commute_plot <- cross_data %>%
  select(dataset, cross_data) %>%
  unnest() %>%
  group_by(dataset) %>%
  mutate(x = case_when(
    f_commuting == 0 ~ -6,
    TRUE ~ log10(f_commuting)
  )) %>%
  filter(between(ecdf(dr_du)(dr_du), 0.025, 0.975)) %>%
  ungroup() %>%
  ggplot(aes(x)) +
  facet_wrap(~ dataset, scales = 'free') +
  geom_point(aes(y = dr_du), shape = 1, size = 0.5) +
  xlab('Commuting fraction (log10)') +
  ylab(expression(paste('Use-resistance association ', (Delta * rho / Delta * tau)))) +
  theme_cowplot() +
  theme(strip.background = element_blank())

ggsave('fig/figure_4.pdf', plot = commute_plot)
