#!/usr/bin/env Rscript --vanilla

library(tidyverse)

# Within-host neutral (WHN) model --------------------------------------

# Parameters that are shared by all the WHN models
whn_base_parms <- list(
  beta = 4.0,
  u = 1.0,
  cost = 0.1683543
)


## WHN single-population model -----------------------------------------

# The ODE simulations use a vector encode the state of the simulation (i.e.,
# the size of each compartment). However, it's easier to express the ODEs in
# terms of a data frame, so I make two helper functions that let me pack the
# state data frame into a vector, feed that to the ODE function, and then
# unpack it again afterward.

# helper function that unpacks state vector to data frame
whn_unpack <- function(x) {
  matrix(x, ncol = 4) %>%
    as_tibble() %>%
    setNames(c('X', 'S', 'R', 'D'))
}

# inverse helper function: pack data frame to vector
whn_pack <- function(df) {
  stopifnot(names(df) == c('X', 'S', 'R', 'D'))
  unlist(df)
}

# ODE function that advances the simulation state for all the WHN simulations
whn_ode_func <- function(t, state_vector, parms) {
  state <- whn_unpack(state_vector)

  with(c(state, parms), {
    n_pop <- nrow(state) # number of populations
    stopifnot(dim(betaij) == rep(n_pop, 2)) # betaij is defined in the next function

    N <- X + S + R + D # total number of individuals

    # Compare to the equations in the Supplement under "WHN model with multiple populations"
    dS <- (betaij %*% ((S + D) / N)) * X - (u + taui) * S - (1 - cost) * (betaij %*% (R / N)) * S
    dR <- (1 - cost) * (betaij %*% (R / N)) * X - u * R + taui * D
    dD <- (1 - cost) * (betaij %*% (R / N)) * S - (u + taui) * D
    dX <- -(dS + dR + dD)

    list(c(dX, dS, dR, dD))
  })
}

# Helper function to check that a transmission matrix:
# - Has nonnegative entries
# - All rows add up to the overall constant
# - All columns add up to the overall constant
check_beta <- function(mat, const) {
  stopifnot(all(mat >= 0))
  stopifnot(all(rowSums(mat) == const))
  stopifnot(all(colSums(mat) == const))
}

# Run the simulation
# taui: a vector of the antibiotic use rates in each population
whn_sim <- function(taui, epsilon) {
  n_pop <- length(taui)

  # Check that epsilon is in the right range
  stopifnot(between(epsilon, 0, 1 - 1 / n_pop))

  parms <- with(whn_base_parms, {
    # Compare to the equations in the Supplement under "WHN model with 2 populations"
    betaij = matrix(beta * epsilon, ncol = n_pop, nrow = n_pop)
    diag(betaij) <- beta * (1 - epsilon)

    check_beta(betaij, beta)

    # Add transmission matrix beta_ij to the parameter list
    c(whn_base_parms, list(
      taui = taui, epsilon = epsilon,
      betaij = betaij
    ))
  })

  # Initial state of the simulation
  state <- tibble(
    X = rep(0.990 / n_pop, n_pop),
    S = rep(0.005 / n_pop, n_pop),
    R = S,
    D = 0
  )

  state_vector <- whn_pack(state)

  # Run the simulation
  result <- rootSolve::runsteady(state_vector, func = whn_ode_func, parms = parms)

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

whn_repro <- whn_sim(taui = seq(0, 3 / 12, length.out = 16), epsilon = 0.0)

whn_repro_plot <- whn_repro %>%
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


## WHN two-population model --------------------------------------------

whn_epsilon_values <- c(0, 1e-4, 0.01, 0.025, 0.050, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5)

whn2 <- tibble(
  delta_tau = round(seq(0.0, 0.15, length.out = 7), 3),
  taui = map(delta_tau, ~ 0.125 + c(-0.5, 0.5) * .)
) %>%
  crossing(epsilon = whn_epsilon_values) %>%
  mutate(
    results = map2(taui, epsilon, whn_sim),
    delta_rho = map_dbl(results, ~ max(.$rho) - min(.$rho)),
    drho_dtau = delta_rho / delta_tau
  )

# Figure 1bc

figure_1bc <- whn2 %>%
  mutate(pop = map(1:n(), ~ c('intervention', 'control'))) %>%
  unnest() %>%
  filter(
    epsilon %in% c(0, 0.01, 0.1),
    delta_tau %in% c(0.05, 0.1)
  ) %>%
  ggplot(aes(x = factor(epsilon), y = rho, fill = pop)) +
  facet_wrap(~ delta_tau) +
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

ggsave('fig/figure_1bc.pdf', plot = figure_1bc)

# Figure 1d

figure_1d <- whn2 %>%
  filter(epsilon %in% c(0, 0.01, 0.1)) %>%
  ggplot(aes(delta_tau, delta_rho, group = factor(epsilon))) +
  geom_point() +
  geom_line()

ggsave('fig/figure_1d.pdf', plot = figure_1d)

# Figure 1e

figure_1e <- whn2 %>%
  filter(delta_tau %in% c(0.05, 0.1)) %>%
  ggplot(aes(epsilon, drho_dtau, group = factor(delta_tau))) +
  geom_point() +
  geom_line()

ggsave('fig/figure_1e.pdf', plot = figure_1e)




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

dtypes_base_parms$mu_d <- with(
  dtypes_base_parms,
  seq(min_mu, max_mu, length = n_d)
)


## Single-population model -----------------------------------------------------

# pack the list into a vector
dtypes_pack <- function(lst) {
  stopifnot(setequal(names(lst), c('X_i', 'S_id', 'R_id')))
  stopifnot(all(dim(lst$S_id) == dim(lst$R_id)))
  stopifnot(dim(lst$S_id)[1] == length(lst$X_i))
  c(lst$X_i, as.vector(lst$S_id), as.vector(lst$R_id))
}

# unpack the vector into a list
dtypes_unpack <- function(x, n_pop, n_d) {
  stopifnot(length(x) == n_pop + 2 * (n_d * n_pop))

  x_end <- n_pop
  s_start <- x_end + 1
  s_end <- x_end + n_d * n_pop
  r_start <- s_end + 1
  r_end <- s_end + n_d * n_pop
  stopifnot(r_end == length(x))

  list(
    X_i = x[1:x_end],
    S_id = matrix(x[s_start:s_end], nrow = n_pop, ncol = n_d),
    R_id = matrix(x[r_start:r_end], nrow = n_pop, ncol = n_d)
  )
}

# for initially packing the data frame into the list
dtypes_df_to_list <- function(df) {
  # check we have the right column names
  stopifnot(setequal(names(df), c('population', 'phenotype', 'dtype', 'value')))
  # check we have the right phenotypes
  stopifnot(setequal(df$phenotype, c('S', 'R', 'X')))

  # extract the no. populations and D-types
  n_pop <- length(unique(df$population))
  n_d <- length(unique(df$dtype))

  # there should be S and R for each pop/D-type combo
  stopifnot(nrow(filter(df, phenotype == 'S')) == n_pop * n_d)
  stopifnot(nrow(filter(df, phenotype == 'R')) == n_pop * n_d)
  # but X for only each pop
  stopifnot(nrow(filter(df, phenotype == 'X')) == n_pop)

  list(
    X_i = df %>% filter(phenotype == 'X') %>% pull(value),
    S_id = df %>% filter(phenotype == 'S') %>% pull(value) %>% matrix(nrow = n_pop, ncol = n_d),
    R_id = df %>% filter(phenotype == 'R') %>% pull(value) %>% matrix(nrow = n_pop, ncol = n_d))
}

# for finally unpacking the list into a data frame
dtypes_list_to_df <- function(lst) {
  # check for names
  stopifnot(setequal(names(lst), c('X_i', 'S_id', 'R_id')))

  # get dimensions
  n_pop <- length(lst$X_i)
  n_d <- dim(lst$S_id)[2]

  # check shapes
  stopifnot(all(dim(lst$S_id) == c(n_pop, n_d)))
  stopifnot(all(dim(lst$R_id) == c(n_pop, n_d)))

  X_rows <- tibble(pop = 1:n_pop, phenotype = 'X', dtype = NA, value = lst$X_i)
  S_rows <- crossing(dtype = 1:n_d, pop = 1:n_pop) %>%
    mutate(phenotype = 'S', value = as.vector(lst$S_id))
  R_rows <- crossing(dtype = 1:n_d, pop = 1:n_pop) %>%
    mutate(phenotype = 'R', value = as.vector(lst$R_id))

  bind_rows(X_rows, S_rows, R_rows) %>%
    arrange(pop, phenotype, dtype)
}

dtypes_ode_func <- function(t, state_vector, parms) {
  state <- dtypes_unpack(state_vector, n_pop = parms$n_pop, n_d = parms$n_d)

  with(c(state, parms), {
    stopifnot(length(X_i) == n_pop)
    stopifnot(dim(S_id) == c(n_pop, n_d))
    stopifnot(dim(R_id) == c(n_pop, n_d))

    # rowSums are over D-types (second index)
    v_id <- (1.0 - ((S_id + R_id) / rowSums(S_id + R_id) - 1.0 / n_d)) ** k
    stopifnot(dim(v_id) == c(n_pop, n_d))

    N_i <- X_i + rowSums(S_id + R_id)
    stopifnot(length(N_i) == n_pop)

    # "tau_i * S_id" does sum_i { tau_i S_id }, which is length P vector
    # to multiply by rows, need to do some fancy footwork: "mat %*% diag(row)"
    dS_id <- v_id * (beta_ij %*% (S_id / N_i)) * X_i - tau_i * S_id - S_id %*% diag(mu_d)
    dR_id <- v_id * (beta_ij %*% (R_id / N_i)) * X_i - c_mu * R_id %*% diag(mu_d)
    dX_i <- -rowSums(dS_id + dR_id)

    list(dtypes_pack(list(X_i = dX_i, S_id = dS_id, R_id = dR_id)))
  })
}

dtypes_sim <- function(tau_i, epsilon = 0) {
  n_pop <- length(tau_i)
  stopifnot(between(epsilon, 0, 1 - 1 / n_pop))

  parms <- with(dtypes_base_parms, {
    beta_ij = matrix(beta * epsilon, ncol = n_pop, nrow = n_pop)
    diag(beta_ij) <- beta * (1 - epsilon)

    check_beta(beta_ij, beta)

    c(dtypes_base_parms, list(
      n_pop = n_pop, tau_i = tau_i, epsilon = epsilon, beta_ij = beta_ij
    ))
  })

  n_d <- parms$n_d

  state <- list(
    X_i = rep(0.9 / n_pop, n_pop),
    S_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d),
    R_id = matrix(0.05 / (n_pop * n_d), nrow = n_pop, ncol = n_d)
  )

  state_vector <- dtypes_pack(state)
  stopifnot(all.equal(sum(state_vector), 1))

  result <- rootSolve::runsteady(
    state_vector,
    func = dtypes_ode_func,
    parms = parms,
    stol = 1e-8 / n_pop,
    rtol = 1e-6 / n_pop,
    atol = 1e-6 / n_pop
  )

  result$y %>%
    dtypes_unpack(n_pop = n_pop, n_d = n_d) %>%
    dtypes_list_to_df() %>%
    left_join(tibble(pop = 1:n_pop, tau = tau_i), by = 'pop')
}

# This reproduces Figure 3 from Lehtinen et al.
dtypes_repro <- dtypes_sim(0.075, 0)

dtypes_repro_plot <- dtypes_repro %>%
  mutate(mu_d = dtypes_base_parms$mu_d[dtype]) %>%
  ggplot(aes(mu_d, value, fill = phenotype)) +
  geom_col() +
  xlab('Clearance rate') +
  ylab('Strain prevalence')

ggsave('fig/lehtinen_figure_3.pdf', plot = dtypes_repro_plot)


## D-types two-population ------------------------------------------------------

dtypes_epsilon_values <- c(0, 1e-4, 0.001, 0.005, 0.01, 0.0175, 0.025, 0.05, 0.075, 0.1)
dtypes2 <- tibble(
  delta_tau = c(0.05, 0.10),
  taui = map(delta_tau, ~ 0.125 + c(-0.5, 0.5) * .)
) %>%
  crossing(epsilon = dtypes_epsilon_values) %>%
  mutate(results = map2(taui, epsilon, dtypes_sim))


dtypes_simplify_results <- function(df) {
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

dtypes2_plot <- dtypes2 %>%
  mutate(
    results = map(results, dtypes_simplify_results),
    delta_rho = map_dbl(results, ~ max(.$rho) - min(.$rho)),
    dr_du = delta_rho / delta_tau
  ) %>%
  ggplot(aes(epsilon, dr_du, group = factor(delta_tau))) +
  geom_point(aes(shape = factor(delta_tau))) +
  geom_line() +
  geom_blank(data = tibble(epsilon = 0, dr_du = 0, delta_tau = 0.01)) +
  scale_shape_manual(
    values = c(1, 16),
    labels = c(
      expression(Delta * tau == 0.05),
      expression(Delta * tau == 0.10)
    )
  ) +
  guides(shape = guide_legend(title = "", label.hjust = 0)) +
  scale_x_continuous(
    expression(epsilon),
    limits = c(0, 0.1),
    expand = c(0.02, 0, 0.05, 0)
  ) +
  scale_y_continuous(
    expression(Delta * rho / Delta * tau),
    limits = c(0, 6.5),
    expand = c(0, 0.1)
  ) +
  theme(legend.position = c(0.5, 0.75))

ggsave('fig/supp_figure_2.pdf', plot = dtypes2_plot)
