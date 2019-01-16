library(dplyr) # for manipulating data frames
library(purrr) # for map and related functions
library(ggplot2) # for plots
library(tidyr) # for "crossing"
library(stringr) # for "str_c"

# Within-host neutral (WHN) model ----------------------------------------------

# Parameters that are shared by all the WHN models
whn_base_parms = list(
  beta = 4.0,
  u = 1.0,
  cost = 0.1683543
)


## WHN single-population model -------------------------------------------------

# The ODE simulations use a vector encode the state of the simulation (i.e.,
# the size of each compartment). However, it's easier to express the ODEs in
# terms of a data frame, so I make two helper functions that let me pack the
# state data frame into a vector, feed that to the ODE function, and then
# unpack it again afterward.

# helper function that unpacks state vector to data frame
whn_unpack = function(x) {
  matrix(x, ncol = 4) %>%
    as_tibble() %>%
    setNames(c('X', 'S', 'R', 'D'))
}

# inverse helper function: pack data frame to vector
whn_pack = function(df) {
  stopifnot(names(df) == c('X', 'S', 'R', 'D'))
  unlist(df)
}

# ODE function that advances the simulation state for all the WHN simulations
whn_ode_func = function(t, state_vector, parms) {
  state = whn_unpack(state_vector)

  with(c(state, parms), {
    n_pop = nrow(state) # number of populations
    stopifnot(dim(betaij) == rep(n_pop, 2)) # betaij is defined in the next function

    N = X + S + R + D # total number of individuals

    # Compare to the equations in the Supplement under "WHN model with multiple populations"
    dS = (betaij %*% ((S + D) / N)) * X - (u + taui) * S - (1 - cost) * (betaij %*% (R / N)) * S
    dR = (1 - cost) * (betaij %*% (R / N)) * X - u * R + taui * D
    dD = (1 - cost) * (betaij %*% (R / N)) * S - (u + taui) * D
    dX = -(dS + dR + dD)

    list(c(dX, dS, dR, dD))
  })
}

# Run the simulation
# taui: a vector of the antibiotic use rates in each population
whn_sim = function(taui, epsilon) {
  n_pop = length(taui)

  parms = with(whn_base_parms, {
    # Compare to the equations in the Supplement under "WHN model with multiple populations"
    betaij = matrix(beta * epsilon / n_pop, ncol = n_pop, nrow = n_pop)
    diag(betaij) <- beta * (1 - epsilon * (1 - 1 / n_pop))

    stopifnot(all(rowSums(betaij) == beta))
    stopifnot(all(colSums(betaij) == beta))

    c(whn_base_parms, list(
      taui = taui, epsilon = epsilon,
      betaij = betaij
    ))
  })

  # Initial state of the simulation
  state = data_frame(
    X = rep(0.990 / n_pop, n_pop),
    S = rep(0.005 / n_pop, n_pop),
    R = S,
    D = 0
  )

  state_vector = whn_pack(state)

  # Run the simulation
  result = rootSolve::runsteady(state_vector, func = whn_ode_func, parms = parms)

  # Unpack the output state and clean it up
  result$y %>%
    whn_unpack() %>%
    mutate(
      tau = taui,
      pop = seq_along(taui),
      rho = R / (S + R + D)
    )
}


# To show that the simulation is working, this code reproduces Davies et al.
# Figure 2g with k = 1.

whn_repro = whn_sim(taui = seq(0, 3 / 12, length.out = 16), epsilon = 0.0)

whn_repro_plot = whn_repro %>%
  ggplot(aes(tau, rho)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(
    name = 'antibiotic tx per month',
    # top axis shows treatments per year
    sec.axis = sec_axis(~ . * 12)
  ) +
  ylab('proportion resistant')

ggsave('fig/davies_figure_2g_k1.pdf', plot = whn_repro_plot)


## WHN two-population model ----------------------------------------------------

# Compare to the parameter values in Methods, "Dynamical model of antibiotic resistance"
twopop_taui_df = data_frame(
  size = c('small', 'large'), # intervention magnitude
  taui = list(c(0.10, 0.15), c(0.05, 0.20))
)

whn2 = twopop_taui_df %>%
  crossing(epsilon = c(0, 0.01, 0.1, 0.25, 0.5)) %>%
  mutate(results = map2(taui, epsilon, whn_sim))

whn2_plot = whn2 %>%
  unnest() %>%
  mutate_at('pop', factor) %>%
  ggplot(aes(x = factor(epsilon), y = rho, fill = pop)) +
  facet_wrap(~ size, nrow = 2) +
  geom_col(position = 'dodge', color = 'black') +
  scale_fill_manual(
    name = '',
    labels = c('intervention', 'control'),
    values = c('white', 'black')
  ) +
  scale_x_discrete(name = expression('interaction strength ' (epsilon))) +
  scale_y_continuous(
    name = 'resistance (%)',
    labels = scales::percent_format(accuracy = 1, suffix = ''),
    limits = c(0, 1.0),
    expand = c(0, 0)
  )

# Plot the two-population results
ggsave('fig/figure_1b.pdf', plot = whn2_plot)


## WHN nested populations ------------------------------------------------------

whn_nest_sim = function(taui, eps_super, eps_sub) {
  n_pop = length(taui) # total number of populations
  pop_per = sqrt(n_pop) # subpopulations per super-population

  # The code is currently limited only to allow N subpopulations per N
  # superpopulations
  if (pop_per != as.integer(pop_per)) {
    stop('length of taui must be a perfect square')
  }

  # Assign a superpopulation to each subpopulation
  superpop = rep(1:pop_per, each = pop_per)
  stopifnot(length(superpop) == length(taui))

  # The "contact" matrix is beta_ij without the beta
  # Compare Supplemental Methods, "Nested population structure"
  contact = outer(1:n_pop, 1:n_pop, Vectorize(function(i, j) {
    if (i == j) {
      (1 - eps_super * (1 - 1 / pop_per)) * (1 - eps_sub * (1 - 1 / pop_per))
    } else if (superpop[i] == superpop[j]) {
      (1 - eps_super * (1 - 1 / pop_per)) * eps_sub / pop_per
    } else {
      eps_super / pop_per * 1 / pop_per
    }
  }))

  stopifnot(all.equal(rowSums(contact), rep(1.0, n_pop)))
  stopifnot(all.equal(colSums(contact), rep(1.0, n_pop)))

  parms = with(whn_base_parms, {
    betaij = beta * contact

    c(whn_base_parms, list(
      taui = taui, eps_super = eps_super, eps_sub = eps_sub,
      betaij = betaij
    ))
  })

  state = data_frame(
    X = rep(0.990 / n_pop, n_pop),
    S = rep(0.005 / n_pop, n_pop),
    R = S,
    D = 0
  )

  state_vector = whn_pack(state)

  result = rootSolve::runsteady(state_vector, func = whn_ode_func, parms = parms)

  result$y %>%
    whn_unpack() %>%
    mutate(pop = seq_along(taui),
           superpop = superpop[pop],
           tau = taui,
           rho = R / (S + R + D))
}

# Compare the taui values with those in Supplemental Table 1
whn_nest_taui = rep(seq(0.07, 0.18, length.out = 4), each = 4) + 0.02 * rep(seq(-1, 1, length.out = 4), 4)

# Run the nested simulations
whn_nest = data_frame(
  eps_super = c(0.0, 0.0, 0.5),
  eps_sub = c(0.0, 0.5, 0.5),
  panel = c(
    'list(epsilon[super] == 0.0, epsilon[sub] == 0.0)',
    'list(epsilon[super] == 0.0, epsilon[sub] == 0.5)',
    'list(epsilon[super] == 0.5, epsilon[sub] == 0.5)'
  )
) %>%
  mutate(data = map2(eps_super, eps_sub, ~ whn_nest_sim(whn_nest_taui, .x, .y))) %>%
  unnest()

# Average the data over superpopulations as in Figure 2b
whn_nest_region = whn_nest %>%
  group_by(panel, superpop) %>%
  summarize_at(c('tau', 'rho'), mean) %>%
  ungroup()

# Get the differences between populations, as shown in Figure 2c
cross_nest = function(df) {
  pops = unique(df$pop)

  crossing(pop1 = pops, pop2 = pops) %>%
    left_join(rename_all(df, ~ str_c(., '1')), by = 'pop1') %>%
    left_join(rename_all(df, ~ str_c(., '2')), by = 'pop2') %>%
    filter(tau1 > tau2)
}

whn_nest_contrast = whn_nest %>%
  nest(-panel) %>%
  mutate(data = map(data, cross_nest)) %>%
  unnest() %>%
  mutate(
    d_rho = rho1 - rho2,
    d_tau = tau1 - tau2,
    adjacent = superpop1 == superpop2
  )

figure_2a = whn_nest %>%
  ggplot(aes(tau, rho, color = factor(superpop))) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_smooth(method = 'lm', se = FALSE) +
  geom_point(shape = 1) +
  guides(color = 'none') +
  scale_x_continuous(
    'antibiotic use (monthly treatments per capita)',
    limits = c(0.05, 0.2)
  ) +
  scale_y_continuous(
    'resistance (%)',
    limits = c(-0.02, 1.02),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

figure_2b = whn_nest_region %>%
  ggplot(aes(tau, rho)) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_smooth(method = 'lm', se = FALSE, color = 'black') +
  geom_point(aes(color = factor(superpop)), size = 2) +
  guides(color = 'none') +
  scale_x_continuous(
    'antibiotic use (monthly treatments per capita)',
    limits = c(0.05, 0.2)
  ) +
  scale_y_continuous(
    'resistance (%)',
    limits = c(-0.0, 1.0),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

figure_2c = whn_nest_contrast %>%
  ggplot(aes(d_tau, d_rho, color = adjacent)) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(
    name = '',
    values = c('black', 'red'),
    labels = c('not adjacent', 'adjacent')
  ) +
  xlab('difference in antibiotic use (monthly treatments per capita)') +
  scale_y_continuous(
    'resistance difference (%)    ',
    labels = scales::percent_format(accuracy = 1, suffix = '')
  ) +
  guides(color = 'none')

ggsave('fig/figure_2a.pdf', plot = figure_2a)
ggsave('fig/figure_2b.pdf', plot = figure_2b)
ggsave('fig/figure_2c.pdf', plot = figure_2c)


# D-types model ----------------------------------------------------------------

# Parameters common to all D-types simulations
dtypes_base_parms = list(
  n_d = 16,
  min_mu = 0.5,
  max_mu = 2.0,
  beta = 2.0,
  c_beta = 1.0,
  c_mu = 1.1,
  k = 15.0
)

dtypes_base_parms$mu_d = with(dtypes_base_parms, seq(min_mu, max_mu, length = n_d))


## Single-population model -----------------------------------------------------

# pack the list into a vector
dtypes_pack = function(lst) {
  stopifnot(setequal(names(lst), c('X_i', 'S_id', 'R_id')))
  stopifnot(all(dim(lst$S_id) == dim(lst$R_id)))
  stopifnot(dim(lst$S_id)[1] == length(lst$X_i))
  c(lst$X_i, as.vector(lst$S_id), as.vector(lst$R_id))
}

# unpack the vector into a list
dtypes_unpack = function(x, n_pop, n_d) {
  stopifnot(length(x) == n_pop + 2 * (n_d * n_pop))

  x_end = n_pop
  s_start = x_end + 1
  s_end = x_end + n_d * n_pop
  r_start = s_end + 1
  r_end = s_end + n_d * n_pop
  stopifnot(r_end == length(x))

  list(X_i = x[1:x_end],
       S_id = matrix(x[s_start:s_end], nrow = n_pop, ncol = n_d),
       R_id = matrix(x[r_start:r_end], nrow = n_pop, ncol = n_d))
}

# for initially packing the data frame into the list
dtypes_df_to_list = function(df) {
  # check we have the right column names
  stopifnot(setequal(names(df), c('population', 'phenotype', 'dtype', 'value')))
  # check we have the right phenotypes
  stopifnot(setequal(df$phenotype, c('S', 'R', 'X')))

  # extract the no. populations and D-types
  n_pop = length(unique(df$population))
  n_d = length(unique(df$dtype))

  # there should be S and R for each pop/D-type combo
  stopifnot(nrow(filter(df, phenotype == 'S')) == n_pop * n_d)
  stopifnot(nrow(filter(df, phenotype == 'R')) == n_pop * n_d)
  # but X for only each pop
  stopifnot(nrow(filter(df, phenotype == 'X')) == n_pop)

  list(X_i = df %>% filter(phenotype == 'X') %$% values,
       S_id = df %>% filter(phenotype == 'S') %$% value %>% matrix(nrow = n_pop, ncol = n_d),
       R_id = df %>% filter(phenotype == 'R') %$% value %>% matrix(nrow = n_pop, ncol = n_d))
}

# for finally unpacking the list into a data frame
dtypes_list_to_df = function(lst) {
  # check for names
  stopifnot(setequal(names(lst), c('X_i', 'S_id', 'R_id')))

  # get dimensions
  n_pop = length(lst$X_i)
  n_d = dim(lst$S_id)[2]

  # check shapes
  stopifnot(all(dim(lst$S_id) == c(n_pop, n_d)))
  stopifnot(all(dim(lst$R_id) == c(n_pop, n_d)))

  X_rows = data_frame(pop = 1:n_pop, phenotype = 'X', dtype = NA, value = lst$X_i)
  S_rows = crossing(dtype = 1:n_d, pop = 1:n_pop) %>%
    mutate(phenotype = 'S', value = as.vector(lst$S_id))
  R_rows = crossing(dtype = 1:n_d, pop = 1:n_pop) %>%
    mutate(phenotype = 'R', value = as.vector(lst$R_id))

  bind_rows(X_rows, S_rows, R_rows) %>%
    arrange(pop, phenotype, dtype)
}

dtypes_ode_func = function(t, state_vector, parms) {
  state = dtypes_unpack(state_vector, n_pop = parms$n_pop, n_d = parms$n_d)

  with(c(state, parms), {
    stopifnot(length(X_i) == n_pop)
    stopifnot(dim(S_id) == c(n_pop, n_d))
    stopifnot(dim(R_id) == c(n_pop, n_d))

    # rowSums are over D-types (second index)
    v_id = (1.0 - ((S_id + R_id) / rowSums(S_id + R_id) - 1.0 / n_d)) ** k
    stopifnot(dim(v_id) == c(n_pop, n_d))

    N_i = X_i + rowSums(S_id + R_id)
    stopifnot(length(N_i) == n_pop)

    # "tau_i * S_id" does sum_i { tau_i S_id }, which is length P vector
    # to multiply by rows, need to do some fancy footwork: "mat %*% diag(row)"
    dS_id = v_id * (beta_ij %*% (S_id / N_i)) * X_i - tau_i * S_id - S_id %*% diag(mu_d)
    dR_id = v_id * (beta_ij %*% (R_id / N_i)) * X_i - c_mu * R_id %*% diag(mu_d)
    dX_i = -rowSums(dS_id + dR_id)

    list(dtypes_pack(list(X_i = dX_i, S_id = dS_id, R_id = dR_id)))
  })
}

dtypes_sim = function(tau_i, epsilon = 0) {
  n_pop = length(tau_i)

  parms = with(dtypes_base_parms, {
    beta_ij = matrix(beta * epsilon / n_pop, ncol = n_pop, nrow = n_pop)
    diag(beta_ij) <- beta * (1 - epsilon * (1 - 1 / n_pop))

    stopifnot(all(rowSums(beta_ij) == beta))
    stopifnot(all(colSums(beta_ij) == beta))

    c(dtypes_base_parms, list(
      n_pop = n_pop, tau_i = tau_i, epsilon = epsilon, beta_ij = beta_ij
    ))
  })

  n_d = parms$n_d

  state = list(X_i = rep(0.9 / n_pop, n_pop),
               S_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d),
               R_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d))

  state_vector = dtypes_pack(state)
  stopifnot(all.equal(sum(state_vector), 1))

  result = rootSolve::runsteady(state_vector, func = dtypes_ode_func, parms = parms,
                                stol = 1e-8 / n_pop, rtol = 1e-6 / n_pop, atol = 1e-6 / n_pop)

  result$y %>%
    dtypes_unpack(n_pop = n_pop, n_d = n_d) %>%
    dtypes_list_to_df() %>%
    left_join(data_frame(pop = 1:n_pop, tau = tau_i), by = 'pop')
}

# This reproduces Figure 3 from Lehtinen et al.
dtypes_repro = dtypes_sim(0.075, 0)

dtypes_repro_plot = dtypes_repro %>%
  mutate(mu_d = dtypes_base_parms$mu_d[dtype]) %>%
  ggplot(aes(mu_d, value, fill = phenotype)) +
  geom_col() +
  xlab('Clearance rate') +
  ylab('Strain prevalence')

ggsave('fig/lehtinen_figure_3.pdf', plot = dtypes_repro_plot)


## D-types two-population ------------------------------------------------------

dtypes2 = twopop_taui_df %>%
  crossing(epsilon = c(0.01, 0.10, 0.5)) %>%
  mutate(results = map2(taui, epsilon, dtypes_sim))

dtypes_simplify_results = function(df) {
  df %>%
    group_by(pop, phenotype) %>%
    summarize(
      value = sum(value),
      tau = unique(tau)
    ) %>%
    ungroup() %>%
    spread(phenotype, value) %>%
    mutate(rho = R / (R + S))
}

dtypes2_plot = dtypes2 %>%
  mutate(results = map(results, dtypes_simplify_results)) %>%
  unnest() %>%
  group_by(epsilon, size) %>%
  mutate(delta_tau = max(tau) - min(tau)) %>%
  ungroup() %>%
  mutate_at('pop', factor) %>%
  ggplot(aes(x = factor(epsilon), y = rho, fill = pop)) +
  facet_grid(
    rows = vars(delta_tau),
    labeller = label_bquote(Delta * tau == .(delta_tau))
  ) +
  geom_col(position = 'dodge', color = 'black') +
  geom_hline(yintercept = 0) +
  scale_fill_manual(
    name = '',
    labels = c('intervention', 'control'),
    values = c('white', 'black')
  ) +
  scale_x_discrete(name = expression('interaction strength ' (epsilon))) +
  scale_y_continuous(
    name = 'resistance (% of hosts)',
    labels = scales::percent_format(accuracy = 1, suffix = ''),
    limits = c(0, 1.0),
    expand = c(0, 0)
  )

ggsave('fig/supp_figure_1.pdf', plot = dtypes2_plot)

## D-types nested --------------------------------------------------------------

dtypes_nest_sim = function(tau_i, eps_super, eps_sub) {
  n_pop = length(tau_i)
  pop_per = sqrt(n_pop)

  # a lot of this code is reproduced from the WHN functions and should be refactored
  if (pop_per != as.integer(pop_per)) {
    stop('length of taui must be a perfect square')
  }

  superpop = rep(1:pop_per, each = pop_per)
  stopifnot(length(superpop) == length(tau_i))

  contact = outer(1:n_pop, 1:n_pop, Vectorize(function(i, j) {
    if (i == j) {
      (1 - eps_super * (1 - 1 / pop_per)) * (1 - eps_sub * (1 - 1 / pop_per))
    } else if (superpop[i] == superpop[j]) {
      (1 - eps_super * (1 - 1 / pop_per)) * eps_sub / pop_per
    } else {
      eps_super / pop_per * 1 / pop_per
    }
  }))

  stopifnot(all.equal(rowSums(contact), rep(1.0, n_pop)))
  stopifnot(all.equal(colSums(contact), rep(1.0, n_pop)))

  parms = with(dtypes_base_parms, {
    c(dtypes_base_parms, list(
      beta_ij = beta * contact,
      n_pop = n_pop, tau_i = tau_i,
      eps_super = eps_super, eps_sub = eps_sub
    ))
  })

  n_d = parms$n_d

  state = list(X_i = rep(0.9 / n_pop, n_pop),
               S_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d),
               R_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d))

  state_vector = dtypes_pack(state)
  stopifnot(all.equal(sum(state_vector), 1))

  result = rootSolve::runsteady(state_vector, func = dtypes_ode_func, parms = parms,
                                stol = 1e-8 / n_pop, rtol = 1e-6 / n_pop, atol = 1e-6 / n_pop)

  result$y %>%
    dtypes_unpack(n_pop = n_pop, n_d = n_d) %>%
    dtypes_list_to_df() %>%
    left_join(data_frame(pop = 1:n_pop, superpop = superpop[pop], tau = tau_i), by = 'pop')
}

dtypes_nest_taui = rep(seq(0.07, 0.18, length.out = 4), each = 4) + 0.02 * rep(seq(-1, 1, length.out = 4), 4)

dtypes_nest = data_frame(
  eps_super = c(0.0, 0.0, 0.5),
  eps_sub = c(0.0, 0.5, 0.5),
  panel = c(
    'list(epsilon[super] == 0.0, epsilon[sub] == 0.0)',
    'list(epsilon[super] == 0.0, epsilon[sub] == 0.5)',
    'list(epsilon[super] == 0.5, epsilon[sub] == 0.5)'
  )
) %>%
  mutate(data = map2(eps_super, eps_sub, ~ dtypes_nest_sim(dtypes_nest_taui, .x, .y))) %>%
  unnest() %>%
  group_by(panel, eps_super, eps_sub, superpop, pop, tau, phenotype) %>%
  summarize(value = sum(value)) %>%
  ungroup() %>%
  spread(phenotype, value) %>%
  mutate(rho = R / (X + S + R))

# Data for Supplemental Figure 2b
dtypes_nest_region = dtypes_nest %>%
  group_by(panel, superpop) %>%
  summarize_at(c('tau', 'rho'), mean) %>%
  ungroup()

# Data for Supplemental Figure 2c
dtypes_nest_contrast = dtypes_nest %>%
  nest(-panel) %>%
  mutate(data = map(data, cross_nest)) %>%
  unnest() %>%
  mutate(
    d_rho = rho1 - rho2,
    d_tau = tau1 - tau2,
    adjacent = superpop1 == superpop2
  )

supp_figure_2a = dtypes_nest %>%
  ggplot(aes(tau, rho, color = factor(superpop))) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_smooth(method = 'lm', se = FALSE) +
  geom_point(shape = 1) +
  guides(color = 'none') +
  scale_x_continuous(
    'antibiotic use (monthly treatments per capita)',
    limits = c(0.05, 0.2)
  ) +
  scale_y_continuous(
    'resistance (%)',
    limits = c(0, 0.5),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

supp_figure_2b = dtypes_nest_region %>%
  ggplot(aes(tau, rho)) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_smooth(method = 'lm', se = FALSE, color = 'black') +
  geom_point(aes(color = factor(superpop)), size = 2) +
  guides(color = 'none') +
  scale_x_continuous(
    'antibiotic use (monthly treatments per capita)',
    limits = c(0.05, 0.2)
  ) +
  scale_y_continuous(
    'resistance (%)',
    limits = c(0, 0.5),
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

supp_figure_2c = dtypes_nest_contrast %>%
  ggplot(aes(d_tau, d_rho, color = adjacent)) +
  facet_grid(cols = vars(panel), labeller = label_parsed) +
  geom_point(shape = 1, size = 2) +
  geom_smooth(method = 'lm', se = FALSE) +
  scale_color_manual(values = c('black', 'red')) +
  guides(color = 'none') +
  xlab('difference in antibiotic use (monthly treatments per capita)') +
  scale_y_continuous(
    'resistance difference (%)    ',
    labels = scales::percent_format(accuracy = 1, suffix = '')
  )

ggsave('fig/supp_figure_2a.pdf', supp_figure_2a)
ggsave('fig/supp_figure_2b.pdf', supp_figure_2b)
ggsave('fig/supp_figure_2c.pdf', supp_figure_2c)
