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
