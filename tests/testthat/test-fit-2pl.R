test_that("fit_2pl recovers item parameters by marginal maximum likelihood", {
  # fit_2pl() is the human-only baseline every mixed-subjects fit starts from,
  # so a scale bias here (e.g. a mis-scaled quadrature) would propagate into
  # every downstream estimate.
  set.seed(707)
  pars <- data.frame(
    item = paste0("Item", 1:8),
    a = c(1, 1.2, 0.9, 1.1, 0.8, 1.5, 0.7, 1.3),
    d = c(0, 0.5, -0.5, 0.2, -0.3, 0.8, -0.9, 0.1)
  )
  resp <- simulate_2pl(rnorm(2000), pars)

  fit <- fit_2pl(resp)

  expect_equal(fit$convergence, 0)
  expect_equal(fit$pars$item, pars$item)
  expect_equal(fit$pars$a, pars$a, tolerance = 0.2)
  expect_equal(fit$pars$d, pars$d, tolerance = 0.2)
  # No systematic scale bias in the discriminations
  expect_gt(mean(fit$pars$a / pars$a), 0.9)
  expect_lt(mean(fit$pars$a / pars$a), 1.1)
  expect_equal(fit$pars$b, -fit$pars$d / fit$pars$a)
  expect_null(fit$model)
})

test_that("fit_2pl handles missing responses and honors initial_pars/quadrature", {
  set.seed(708)
  pars <- data.frame(item = paste0("I", 1:5), a = seq(0.8, 1.4, l = 5),
                     d = seq(-0.8, 0.8, l = 5))
  resp <- simulate_2pl(rnorm(800), pars)
  resp[cbind(sample(nrow(resp), 100), sample(ncol(resp), 100, replace = TRUE))] <- NA

  fit <- fit_2pl(resp)
  expect_equal(fit$convergence, 0)
  expect_equal(fit$pars$item, pars$item)
  expect_true(all(is.finite(fit$pars$a)))

  # Starting from the truth and starting from the default should agree
  from_truth <- fit_2pl(resp, initial_pars = pars)
  expect_equal(from_truth$pars$a, fit$pars$a, tolerance = 1e-3)
  expect_equal(from_truth$pars$d, fit$pars$d, tolerance = 1e-3)

  # An explicitly supplied quadrature is used
  quad <- make_quadrature(21)
  expect_equal(fit_2pl(resp, quadrature = quad)$pars$a,
               fit_2pl(resp, n_quad = 21)$pars$a)

  # Discriminations respect the lower bound
  expect_true(all(fit_2pl(resp, slope_lower = 0.5)$pars$a >= 0.5 - 1e-8))
})

test_that("fit_2pl agrees with the mirt 2PL MLE", {
  # mirt is a suggested dependency only; skip where it is unavailable or its
  # own dependencies fail to load.
  skip_if_not_installed("mirt")
  set.seed(709)
  pars <- data.frame(item = paste0("Item", 1:8),
                     a = c(1, 1.2, 0.9, 1.1, 0.8, 1.5, 0.7, 1.3),
                     d = c(0, 0.5, -0.5, 0.2, -0.3, 0.8, -0.9, 0.1))
  resp <- simulate_2pl(rnorm(2000), pars)

  native <- fit_2pl(resp)$pars
  ref <- as.data.frame(mirt::coef(
    mirt::mirt(as.data.frame(resp), 1, "2PL", verbose = FALSE),
    simplify = TRUE
  )$items)

  expect_equal(native$a, ref$a1, tolerance = 0.05)
  expect_equal(native$d, ref$d, tolerance = 0.05)
})
