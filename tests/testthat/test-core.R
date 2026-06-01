test_that("quadrature weights are standard-normal weights", {
  quad <- make_quadrature(7)

  expect_equal(sum(quad$weight), 1, tolerance = 1e-10)
  expect_equal(sum(quad$theta * quad$weight), 0, tolerance = 1e-10)
})

test_that("posterior weights and expected counts have expected dimensions", {
  pars <- data.frame(item = paste0("I", 1:3), a = c(0.8, 1.1, 1.3),
                     d = c(-0.2, 0.1, 0.4))
  resp <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
  colnames(resp) <- pars$item

  weights <- posterior_weights_2pl(resp, pars, n_quad = 5)
  counts <- summarize_expected_counts(resp, weights)

  expect_equal(dim(weights), c(2, 5))
  expect_equal(rowSums(weights), c(1, 1), tolerance = 1e-10)
  expect_equal(dim(counts$N), c(3, 5))
  expect_equal(dim(counts$R), c(3, 5))
})

test_that("mixed-subjects fitting works with supplied starting parameters", {
  set.seed(11)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(60), pars)
  predicted <- observed
  generated <- simulate_2pl(rnorm(120), pars)

  fit <- fit_mixed_subjects(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = 0.5,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 50)
  )

  expect_s3_class(fit, "mixedsubjects_fit")
  expect_equal(nrow(fit$item_pars), 4)
  expect_true(all(is.finite(fit$item_pars$a)))
  expect_true(all(is.finite(fit$item_pars$d)))
})

test_that("paired missingness is enforced by default", {
  pars <- data.frame(item = paste0("I", 1:3), a = c(0.8, 1.1, 1.3),
                     d = c(-0.2, 0.1, 0.4))
  observed <- matrix(c(1, NA, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
  predicted <- observed
  predicted[1, 2] <- 1
  generated <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 2, byrow = TRUE)
  colnames(observed) <- colnames(predicted) <- colnames(generated) <- pars$item

  expect_error(
    fit_mixed_subjects(
      observed = observed,
      predicted = predicted,
      generated = generated,
      lambda = 0.5,
      initial_pars = pars,
      n_quad = 5,
      control = list(maxit = 20)
    ),
    "same missingness pattern"
  )

  expect_s3_class(
    fit_mixed_subjects(
      observed = observed,
      predicted = predicted,
      generated = generated,
      lambda = 0.5,
      initial_pars = pars,
      n_quad = 5,
      paired_missing = "allow",
      control = list(maxit = 20)
    ),
    "mixedsubjects_fit"
  )
})

test_that("split fitting pools held-out counts", {
  set.seed(12)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(60), pars)
  predicted <- observed
  generated <- simulate_2pl(rnorm(120), pars)
  split_id <- rep(1:2, each = 30)

  fit <- fit_mixed_subjects_split(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = 0.5,
    split_id = split_id,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 50)
  )

  expect_s3_class(fit, "mixedsubjects_fit")
  expect_equal(fit$split$n_splits, 2)
  expect_equal(fit$q_observed$counts$n, nrow(observed))
  expect_equal(fit$q_generated$counts$n, nrow(generated))
})

test_that("from-quadrature interface matches raw workflow", {
  set.seed(13)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(60), pars)
  predicted <- observed
  generated <- simulate_2pl(rnorm(120), pars)
  quad <- make_quadrature(7)

  q_observed <- mixed_subjects_quadrature(observed, item_pars = pars,
                                          quadrature = quad)
  q_predicted <- mixed_subjects_quadrature(
    predicted,
    item_pars = pars,
    quadrature = quad,
    weights = NULL
  )
  q_predicted$weights <- q_observed$weights
  q_predicted$counts <- summarize_expected_counts(predicted, q_observed$weights)
  q_predicted$quad <- mixedsubjectsirt:::counts_to_quad(q_predicted$counts)
  q_generated <- mixed_subjects_quadrature(generated, item_pars = pars,
                                           quadrature = quad)

  raw_fit <- fit_mixed_subjects(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = 0.5,
    initial_pars = pars,
    quadrature = quad,
    control = list(maxit = 50)
  )
  count_fit <- fit_mixed_subjects_from_quadrature(
    q_observed = q_observed,
    q_predicted = q_predicted,
    q_generated = q_generated,
    lambda = 0.5,
    initial_pars = pars,
    control = list(maxit = 50)
  )

  expect_equal(raw_fit$par, count_fit$par, tolerance = 1e-8)
})

test_that("analytic expected-count gradient matches finite differences", {
  set.seed(14)
  n_items <- 3
  n_nodes <- 5
  theta <- seq(-2, 2, length.out = n_nodes)
  item_names <- paste0("I", seq_len(n_items))
  N <- matrix(runif(n_items * n_nodes, 5, 30), nrow = n_items)
  p <- matrix(runif(n_items * n_nodes, 0.15, 0.85), nrow = n_items)
  counts <- list(
    N = N,
    R = N * p,
    n = 100,
    n_items = n_items,
    n_nodes = n_nodes,
    item_names = item_names,
    theta = theta,
    weight = rep(1 / n_nodes, n_nodes),
    node_count = colSums(N) / n_items
  )
  class(counts) <- "mixedsubjects_counts"

  pars <- data.frame(item = item_names, a = c(0.9, 1.1, 1.3),
                     d = c(-0.4, 0.2, 0.5))
  par <- c(pars$a, pars$d)
  analytic <- mixedsubjectsirt:::gradient_expected_counts(counts, pars)

  eps <- 1e-6
  numeric <- numeric(length(par))
  for (i in seq_along(par)) {
    plus <- minus <- par
    plus[i] <- plus[i] + eps
    minus[i] <- minus[i] - eps
    numeric[i] <- (
      mixedsubjectsirt:::loss_expected_counts(
        counts,
        mixedsubjectsirt:::item_pars_from_vector(plus, item_names)
      ) -
        mixedsubjectsirt:::loss_expected_counts(
          counts,
          mixedsubjectsirt:::item_pars_from_vector(minus, item_names)
        )
    ) / (2 * eps)
  }

  expect_equal(analytic, numeric, tolerance = 1e-5)
})

test_that("lambda zero removes the correction term", {
  set.seed(15)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(50), pars)
  predicted <- 1 - observed
  generated <- simulate_2pl(rnorm(200), pars)

  fit0 <- fit_mixed_subjects(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = 0,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )
  fit0_changed <- fit_mixed_subjects(
    observed = observed,
    predicted = observed,
    generated = 1 - generated,
    lambda = 0,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )

  expect_equal(fit0$par, fit0_changed$par, tolerance = 1e-8)
})

test_that("informative generated responses can improve recovery", {
  set.seed(16)
  pars <- data.frame(item = paste0("I", 1:5), a = c(0.8, 1.0, 1.2, 1.4, 1.6),
                     d = c(-0.8, -0.3, 0.1, 0.4, 0.9))
  observed <- simulate_2pl(rnorm(80), pars)
  predicted <- observed
  generated <- simulate_2pl(rnorm(2000), pars)

  fit0 <- fit_mixed_subjects(observed, predicted, generated,
                             lambda = 0, initial_pars = pars, n_quad = 9,
                             control = list(maxit = 100))
  fit_high <- fit_mixed_subjects(observed, predicted, generated,
                                 lambda = 0.9, initial_pars = pars, n_quad = 9,
                                 control = list(maxit = 100))

  rmse0 <- sqrt(mean((fit0$item_pars$d - pars$d)^2))
  rmse_high <- sqrt(mean((fit_high$item_pars$d - pars$d)^2))

  expect_lt(rmse_high, rmse0)
})

test_that("diagnostic lambda grid flags uninformative correction", {
  set.seed(17)
  pars <- data.frame(item = paste0("I", 1:5), a = c(0.8, 1.0, 1.2, 1.4, 1.6),
                     d = c(-0.8, -0.3, 0.1, 0.4, 0.9))
  observed <- simulate_2pl(rnorm(120), pars)
  predicted <- matrix(rbinom(length(observed), 1, 0.5), nrow = nrow(observed))
  generated <- matrix(rbinom(500 * ncol(observed), 1, 0.5), nrow = 500)
  colnames(predicted) <- colnames(generated) <- colnames(observed)

  diagnostic <- diagnose_lambda_grid(
    lambda_grid = c(0, 0.5, 1),
    observed = observed,
    predicted = predicted,
    generated = generated,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )

  expect_equal(diagnostic$lowest_observed_loss_lambda, 0)
  expect_warning(
    tune_lambda_grid(
      lambda_grid = c(0, 0.5),
      observed = observed,
      predicted = predicted,
      generated = generated,
      initial_pars = pars,
      n_quad = 5,
      control = list(maxit = 30)
    ),
    "diagnostic"
  )
})

test_that("split-specific zero lambda ignores that fold's predictions", {
  set.seed(18)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(80), pars)
  predicted <- observed
  generated <- simulate_2pl(rnorm(400), pars)
  split_id <- rep(1:2, each = 40)

  changed <- predicted
  changed[split_id == 2, ] <- 1 - changed[split_id == 2, ]

  fit_original <- fit_mixed_subjects_split(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = c(0.8, 0),
    split_id = split_id,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )
  fit_changed <- fit_mixed_subjects_split(
    observed = observed,
    predicted = changed,
    generated = generated,
    lambda = c(0.8, 0),
    split_id = split_id,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )

  expect_equal(fit_original$par, fit_changed$par, tolerance = 1e-8)
})

test_that("sandwich covariance and ability risk have expected shapes", {
  set.seed(19)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(80), pars)
  generated <- simulate_2pl(rnorm(200), pars)

  fit <- fit_mixed_subjects(
    observed = observed,
    predicted = observed,
    generated = generated,
    lambda = 0.5,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )
  Sigma <- vcov_mixed_subjects(fit)
  risk <- ability_risk(observed, fit, vcov = Sigma)

  expect_equal(dim(Sigma), c(8, 8))
  expect_equal(Sigma, t(Sigma), tolerance = 1e-10)
  expect_equal(nrow(risk$details), nrow(observed))
  expect_true(is.finite(risk$summary$mean_total_risk))
})

test_that("ability gradient matches finite differences of score_theta", {
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  resp <- matrix(c(1, 0, 1, 0), nrow = 1)
  colnames(resp) <- pars$item

  theta <- score_theta(resp, pars)
  analytic <- ability_gradient(resp, pars, theta = theta)[1, ]
  par <- c(pars$a, pars$d)
  eps <- 1e-5
  numeric <- numeric(length(par))

  for (i in seq_along(par)) {
    plus <- minus <- par
    plus[i] <- plus[i] + eps
    minus[i] <- minus[i] - eps
    plus_pars <- mixedsubjectsirt:::item_pars_from_vector(plus, pars$item)
    minus_pars <- mixedsubjectsirt:::item_pars_from_vector(minus, pars$item)
    numeric[i] <- (
      score_theta(resp, plus_pars) -
        score_theta(resp, minus_pars)
    ) / (2 * eps)
  }

  expect_equal(unname(analytic), numeric, tolerance = 1e-3)
})

test_that("ability-risk tuning chooses a candidate lambda and crossfits", {
  set.seed(20)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed <- simulate_2pl(rnorm(80), pars)
  generated <- simulate_2pl(rnorm(240), pars)

  tuned <- tune_lambda_ability(
    lambda_grid = c(0, 0.5),
    observed = observed,
    predicted = observed,
    generated = generated,
    initial_pars = pars,
    n_quad = 7,
    control = list(maxit = 80)
  )

  expect_true(tuned$best_lambda %in% c(0, 0.5))
  expect_equal(nrow(tuned$summary), 2)
  expect_true(all(is.finite(tuned$summary$mean_total_risk)))

  split_tuned <- tune_lambda_ability_crossfit(
    lambda_grid = c(0, 0.5),
    observed = observed,
    predicted = observed,
    generated = generated,
    split_id = rep(1:2, each = 40),
    initial_pars = pars,
    n_quad = 5,
    control = list(maxit = 60)
  )

  expect_length(split_tuned$lambda_by_split, 2)
  expect_s3_class(split_tuned$final_fit, "mixedsubjects_fit")
})
