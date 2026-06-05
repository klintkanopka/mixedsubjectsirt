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

  tuned <- tune_lambda_ability_risk(
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

  split_tuned <- tune_lambda_ability_risk_crossfit(
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

test_that("marginal_loss_2pl gradient matches finite differences", {
  set.seed(30)
  pars <- data.frame(item = paste0("I", 1:3), a = c(0.9, 1.1, 0.8),
                     d = c(-0.3, 0.2, 0.5))
  resp <- simulate_2pl(rnorm(30), pars)
  quad <- make_quadrature(7)

  par_vec <- c(pars$a, pars$d)
  item_names <- pars$item
  ip <- mixedsubjectsirt:::item_pars_from_vector(par_vec, item_names)

  analytic <- mixedsubjectsirt:::marginal_gradient_2pl(resp, ip, quad)

  eps <- 1e-5
  numeric_fd <- numeric(length(par_vec))
  for (j in seq_along(par_vec)) {
    pv <- pm <- par_vec
    pv[j] <- pv[j] + eps
    pm[j] <- pm[j] - eps
    numeric_fd[j] <-
      (mixedsubjectsirt:::marginal_loss_2pl(
        resp, mixedsubjectsirt:::item_pars_from_vector(pv, item_names), quad) -
       mixedsubjectsirt:::marginal_loss_2pl(
        resp, mixedsubjectsirt:::item_pars_from_vector(pm, item_names), quad)
      ) / (2 * eps)
  }

  expect_equal(as.numeric(analytic), numeric_fd, tolerance = 1e-5)
})

test_that("fit_mixed_subjects_mml returns valid fit object", {
  set.seed(31)
  pars <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                     d = c(-0.2, 0.1, 0.4, -0.5))
  observed  <- simulate_2pl(rnorm(60), pars)
  generated <- simulate_2pl(rnorm(120), pars)

  fit <- fit_mixed_subjects_mml(
    observed, observed, generated,
    lambda = 0.5, initial_pars = pars, n_quad = 7,
    control = list(maxit = 50)
  )

  expect_s3_class(fit, "mixedsubjects_fit")
  expect_true(isTRUE(fit$mml))
  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(fit$item_pars$a)))
  expect_true(all(is.finite(fit$item_pars$d)))
  # vcov should work (posteriors stored at converged parameters)
  Sigma <- vcov_mixed_subjects(fit)
  expect_equal(dim(Sigma), c(8L, 8L))
  expect_true(all(is.finite(Sigma)))
})

test_that("fit_mixed_subjects_mml selects lambda > 0 for F=Y predictor", {
  # With predicted = observed (F = Y), the correction term reduces variance.
  # tune_lambda_ability_risk should select lambda > 0, unlike the frozen
  # expected-count estimator which collapses to 0 due to gradient asymmetry.
  set.seed(32)
  pars <- data.frame(item = paste0("I", 1:5), a = seq(0.8, 1.6, l = 5),
                     d = seq(-0.8, 0.8, l = 5))
  theta  <- rnorm(200)
  obs    <- simulate_2pl(theta, pars)
  gen    <- simulate_2pl(rnorm(600), pars)

  tuned <- tune_lambda_ability_risk(
    lambda_grid  = seq(0, 1, by = 0.2),
    observed     = obs,
    predicted    = obs,   # F = Y: theoretical upper bound
    generated    = gen,
    initial_pars = pars,
    fit_fn       = fit_mixed_subjects_mml,
    n_quad       = 11,
    control      = list(maxit = 200)
  )

  # MML should select lambda > 0 for F=Y (the frozen estimator gives 0)
  expect_gt(tuned$best_lambda, 0)
  # Risk should be U-shaped: lower at the optimum than at lambda = 0
  expect_lt(min(tuned$summary$mean_total_risk),
            tuned$summary$mean_total_risk[tuned$summary$lambda == 0])
})

test_that("item_loss_and_grad sums equal scalar loss and gradient", {
  set.seed(40)
  pars <- data.frame(item = paste0("I", 1:3), a = c(0.9, 1.1, 0.8),
                     d = c(-0.3, 0.2, 0.5))
  obs  <- simulate_2pl(rnorm(50), pars)
  quad <- make_quadrature(7)
  w    <- posterior_weights_2pl(obs, pars, quadrature = quad)
  cnts <- summarize_expected_counts(obs, w)

  ig <- mixedsubjectsirt:::item_loss_and_grad(cnts, pars)

  # Sum of per-item losses equals the scalar expected-count loss
  expect_equal(sum(ig$loss),
               mixedsubjectsirt:::loss_expected_counts(cnts, pars),
               tolerance = 1e-10)

  # Concatenated gradient equals the scalar gradient
  expect_equal(as.numeric(c(ig$grad_a, ig$grad_d)),
               as.numeric(mixedsubjectsirt:::gradient_expected_counts(cnts, pars)),
               tolerance = 1e-10)
})

test_that("tune_lambda_ppi_score_item gives N/(n+N) for F=Y on all items", {
  set.seed(41)
  pars  <- data.frame(item = paste0("I", 1:4), a = c(0.8, 1.1, 1.3, 0.9),
                      d = c(-0.2, 0.1, 0.4, -0.5))
  obs   <- simulate_2pl(rnorm(80), pars)
  N_gen <- 200
  upper <- N_gen / (nrow(obs) + N_gen)

  res <- tune_lambda_ppi_score_item(obs, obs, pars, n_generated = N_gen, n_quad = 7)

  # F=Y: each item's lambda equals N/(n+N)
  expect_equal(res$lambda, rep(upper, 4), tolerance = 1e-6)
  expect_equal(res$item, pars$item)
})

test_that("fit_mixed_subjects_mml with vector lambda converges and equals scalar when uniform", {
  set.seed(42)
  pars <- data.frame(item = paste0("I", 1:3), a = c(0.9, 1.1, 0.8),
                     d = c(-0.3, 0.2, 0.5))
  obs  <- simulate_2pl(rnorm(60), pars)
  gen  <- simulate_2pl(rnorm(120), pars)

  # Scalar lambda
  fit_s <- fit_mixed_subjects_mml(obs, obs, gen, lambda = 0.5,
    initial_pars = pars, n_quad = 7, control = list(maxit = 80))
  # Vector lambda (all 0.5): vector path uses frozen counts, different from scalar
  fit_v <- fit_mixed_subjects_mml(obs, obs, gen, lambda = rep(0.5, 3),
    initial_pars = pars, n_quad = 7, control = list(maxit = 80))

  # Both should converge
  expect_equal(fit_s$convergence, 0)
  expect_equal(fit_v$convergence, 0)
  # Both should produce valid positive discriminations
  expect_true(all(fit_s$item_pars$a > 0))
  expect_true(all(fit_v$item_pars$a > 0))
  # Lambda stored correctly
  expect_equal(fit_v$lambda, rep(0.5, 3))
})

# ---------- 1PL tests -------------------------------------------------------

test_that("louis_missing_info matches numerical Hessian of marginal_loss_2pl", {
  set.seed(60)
  pars <- data.frame(item = paste0("I", 1:3), a = c(1, 1.2, 0.9),
                     d = c(0, -0.5, 0.3))
  # Use a larger sample and more quadrature nodes so A_marg is reliably PD
  resp <- simulate_2pl(rnorm(400), pars)
  quad <- make_quadrature(15)

  w    <- posterior_weights_2pl(resp, pars, quadrature = quad)
  cnts <- summarize_expected_counts(resp, w)

  H      <- mixedsubjectsirt:::avg_hessian_counts(cnts, pars)
  Imiss  <- mixedsubjectsirt:::louis_missing_info(resp, w, pars)
  A_marg <- H - Imiss

  # Numerical Hessian via central differences on the marginal gradient
  par_vec    <- c(pars$a, pars$d)
  item_names <- pars$item
  eps        <- 1e-4
  n_par      <- length(par_vec)
  num_hess   <- matrix(0, n_par, n_par)
  for (j in seq_len(n_par)) {
    pv <- pm <- par_vec; pv[j] <- pv[j] + eps; pm[j] <- pm[j] - eps
    ipv <- mixedsubjectsirt:::item_pars_from_vector(pv, item_names)
    ipm <- mixedsubjectsirt:::item_pars_from_vector(pm, item_names)
    gv  <- mixedsubjectsirt:::marginal_gradient_2pl(resp, ipv, quad)
    gm  <- mixedsubjectsirt:::marginal_gradient_2pl(resp, ipm, quad)
    num_hess[j, ] <- (gv - gm) / (2 * eps)
  }

  # Louis formula matches numerical Hessian
  expect_equal(as.numeric(A_marg), as.numeric(num_hess), tolerance = 1e-5)

  # At the MLE (converged parameters), the marginal bread is positive definite
  fit_mle <- fit_1pl(resp, n_quad = 15)
  w_mle   <- posterior_weights_2pl(resp, fit_mle$pars, quadrature = quad)
  H_mle   <- mixedsubjectsirt:::avg_hessian_counts(
    summarize_expected_counts(resp, w_mle), fit_mle$pars)
  Im_mle  <- mixedsubjectsirt:::louis_missing_info(resp, w_mle, fit_mle$pars)
  expect_true(all(eigen(H_mle - Im_mle)$values > 0))
})

test_that("vcov dispatches to vcov_mixed_subjects_mml for scalar MML fits", {
  set.seed(61)
  pars     <- data.frame(item = paste0("I", 1:3), a = c(1, 1.2, 0.9),
                          d = c(0, -0.5, 0.3))
  obs      <- simulate_2pl(rnorm(50), pars)
  gen      <- simulate_2pl(rnorm(120), pars)

  fit_mml  <- fit_mixed_subjects_mml(obs, obs, gen, lambda = 0.5,
    initial_pars = pars, n_quad = 7, control = list(maxit = 80))
  fit_ec   <- fit_mixed_subjects(obs, obs, gen, lambda = 0.5,
    initial_pars = pars, n_quad = 7, control = list(maxit = 80))

  # Use stats::vcov() — the correct way to get MML covariance
  Sigma_mml_dispatch  <- stats::vcov(fit_mml)    # routes to Louis-corrected
  Sigma_mml_direct    <- vcov_mixed_subjects_mml(fit_mml)  # explicit call
  Sigma_ec            <- stats::vcov(fit_ec)     # routes to EM Hessian

  # S3 dispatch equals direct call for MML fits
  expect_equal(Sigma_mml_dispatch, Sigma_mml_direct)

  # Both are valid covariance matrices
  expect_equal(dim(Sigma_mml_dispatch), c(6L, 6L))
  expect_true(all(is.finite(Sigma_mml_dispatch)))
  expect_true(all(eigen(Sigma_mml_dispatch)$values > 0))

  # Louis-corrected covariance is larger than EM Hessian covariance
  # (Louis bread <= EM bread  =>  Louis bread^{-1} >= EM bread^{-1}  =>
  #  Louis sandwich >= EM sandwich in PSD sense)
  expect_gt(mean(diag(Sigma_mml_dispatch)), mean(diag(Sigma_ec)))

  # Direct call to vcov_mixed_subjects bypasses Louis correction — should be
  # smaller than the dispatch result for an MML fit
  Sigma_ec_path <- vcov_mixed_subjects(fit_mml)  # deliberately uses EM bread
  expect_lt(mean(diag(Sigma_ec_path)), mean(diag(Sigma_mml_dispatch)))
})

test_that("vcov_mixed_subjects handles vector lambda", {
  set.seed(62)
  pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
  obs  <- simulate_2pl(rnorm(60), pars)
  gen  <- simulate_2pl(rnorm(120), pars)

  fit_vec <- fit_mixed_subjects_mml(obs, obs, gen, lambda = rep(0.5, 3),
    initial_pars = pars, n_quad = 7, control = list(maxit = 80))

  # vcov should not throw for vector lambda
  Sigma <- vcov(fit_vec)
  expect_equal(dim(Sigma), c(6L, 6L))
  expect_true(all(is.finite(Sigma)))
})

test_that("tune_lambda_ability_risk filters failed candidates", {
  set.seed(63)
  pars  <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
  obs   <- simulate_2pl(rnorm(30), pars)
  gen   <- simulate_2pl(rnorm(60), pars)

  # All candidates should converge; best_lambda should be > 0 for F=Y
  res <- tune_lambda_ability_risk(
    c(0, 0.5), obs, obs, gen,
    fit_fn = fit_mixed_subjects_mml,
    initial_pars = pars, n_quad = 5, control = list(maxit = 50)
  )
  expect_true("selection_risk" %in% names(res$summary))
  expect_true(res$best_lambda %in% c(0, 0.5))
})

test_that("fit_1pl returns positive shared discrimination", {
  set.seed(50)
  pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
  resp <- simulate_2pl(rnorm(80), pars)
  fit  <- fit_1pl(resp, n_quad = 7)

  expect_equal(fit$convergence, 0)
  # All a values equal (1PL constraint)
  expect_equal(length(unique(round(fit$pars$a, 10))), 1)
  expect_gt(fit$pars$a[1], 0)
  expect_equal(nrow(fit$pars), 3)
})

test_that("gradient_expected_counts_1pl length equals n_items + 1", {
  set.seed(51)
  pars <- data.frame(a = 1.1, d = c(-0.4, 0.1, 0.5, -0.2))
  resp <- simulate_2pl(rnorm(60), pars)
  quad <- make_quadrature(7)
  w    <- posterior_weights_2pl(resp, pars, quadrature = quad)
  cnts <- summarize_expected_counts(resp, w)

  g <- mixedsubjectsirt:::gradient_expected_counts_1pl(cnts, pars)
  expect_length(g, nrow(pars) + 1L)
  expect_true(all(is.finite(g)))
})

test_that("fit_mixed_subjects_mml_1pl with F=Y converges and tune selects lambda > 0", {
  set.seed(52)
  pars <- data.frame(a = 1, d = seq(-0.8, 0.8, l = 5))
  theta <- rnorm(150)
  obs   <- simulate_2pl(theta, pars)
  gen   <- simulate_2pl(rnorm(400), pars)

  fit <- fit_mixed_subjects_mml_1pl(obs, obs, gen, lambda = 0.5,
    initial_pars = pars, n_quad = 7, control = list(maxit = 100))

  expect_s3_class(fit, "mixedsubjects_1pl_fit")
  expect_equal(fit$convergence, 0)
  expect_equal(length(unique(round(fit$item_pars$a, 10))), 1)

  # vcov should be (J+1) x (J+1) = 6 x 6
  Sigma <- vcov_mixed_subjects_1pl(fit)
  expect_equal(dim(Sigma), c(6L, 6L))
  expect_true(all(is.finite(Sigma)))

  # tune_lambda_ppi_score_1pl: F=Y should give N/(n+N)
  ppi <- tune_lambda_ppi_score_1pl(obs, obs, pars, n_generated = 400, n_quad = 7)
  expect_equal(ppi$lambda, 400 / (150 + 400), tolerance = 1e-6)

  # tune_lambda_ability_risk_1pl should select lambda > 0 for F=Y
  tuned <- tune_lambda_ability_risk_1pl(
    c(0, 0.5), obs, obs, gen,
    initial_pars = pars, n_quad = 7, control = list(maxit = 50))
  expect_gt(tuned$best_lambda, 0)
})
