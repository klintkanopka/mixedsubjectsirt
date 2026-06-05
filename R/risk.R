avg_hessian_counts <- function(counts, item_pars) {
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = counts$n_items,
    item_names = counts$item_names
  )

  p <- 2 * counts$n_items
  H <- matrix(0, nrow = p, ncol = p)

  for (j in seq_len(counts$n_items)) {
    eta <- item_pars$d[j] + item_pars$a[j] * counts$theta
    pr <- stats::plogis(eta)
    w <- counts$N[j, ] * pr * (1 - pr) / counts$n

    idx_a <- j
    idx_d <- counts$n_items + j
    H[idx_a, idx_a] <- sum(w * counts$theta^2)
    H[idx_a, idx_d] <- sum(w * counts$theta)
    H[idx_d, idx_a] <- H[idx_a, idx_d]
    H[idx_d, idx_d] <- sum(w)
  }

  H
}

person_scores_2pl <- function(resp, weights, item_pars) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  weights <- as.matrix(weights)

  if (nrow(resp) != nrow(weights)) {
    stop("resp and weights must have the same number of rows.", call. = FALSE)
  }

  theta <- attr(weights, "theta")
  if (is.null(theta)) {
    stop("weights must have a theta attribute.", call. = FALSE)
  }

  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )

  n <- nrow(resp)
  n_items <- ncol(resp)
  scores <- matrix(0, nrow = n, ncol = 2 * n_items)

  eta <- outer(theta, item_pars$a, `*`) +
    matrix(item_pars$d, nrow = length(theta), ncol = n_items, byrow = TRUE)
  prob <- stats::plogis(eta)

  for (j in seq_len(n_items)) {
    observed <- !is.na(resp[, j])
    if (any(observed)) {
      W <- weights[observed, , drop = FALSE]
      y <- resp[observed, j]
      ep <- as.vector(W %*% prob[, j])
      ep_theta <- as.vector(W %*% (prob[, j] * theta))
      etheta <- as.vector(W %*% theta)

      scores[observed, j] <- ep_theta - y * etheta
      scores[observed, n_items + j] <- ep - y
    }
  }

  colnames(scores) <- c(
    paste0("a_", item_pars$item),
    paste0("d_", item_pars$item)
  )
  scores
}

stable_inverse <- function(x, ridge = 1e-8, max_tries = 8) {
  if (nrow(x) != ncol(x)) {
    stop("x must be square.", call. = FALSE)
  }

  for (i in 0:max_tries) {
    ridge_i <- ridge * 10^i
    out <- tryCatch(
      solve(x + diag(ridge_i, nrow(x))),
      error = function(e) NULL
    )
    if (!is.null(out) && all(is.finite(out))) {
      return(out)
    }
  }

  stop("Could not invert matrix, even after ridge regularization.",
       call. = FALSE)
}

safe_cov <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) <= 1) {
    return(matrix(0, ncol(x), ncol(x)))
  }
  stats::cov(x)
}

extract_q_for_vcov <- function(fit) {
  if (!inherits(fit, "mixedsubjects_fit")) {
    stop("fit must be a mixedsubjects_fit object.", call. = FALSE)
  }

  needed <- c("q_observed", "q_predicted", "q_generated")
  if (!all(needed %in% names(fit))) {
    stop("fit does not contain quadrature summaries.", call. = FALSE)
  }

  summaries <- list(fit$q_observed, fit$q_predicted, fit$q_generated)
  has_data <- vapply(
    summaries,
    function(q) is.list(q) && !is.null(q$resp) && !is.null(q$weights),
    logical(1)
  )

  if (!all(has_data)) {
    stop(
      "fit must contain raw response matrices and posterior weights. ",
      "Refit with fit_mixed_subjects() or supply complete quadrature summaries.",
      call. = FALSE
    )
  }

  summaries
}

#' Sandwich covariance for a mixed-subjects fit
#'
#' Estimates the full sandwich covariance matrix for item parameters from the
#' fixed-posterior expected-count estimating equations. The parameter order is
#' all discriminations followed by all intercepts, matching `fit$par`.
#'
#' @param object A fitted object returned by [fit_mixed_subjects()] or
#'   [fit_mixed_subjects_from_quadrature()] with response matrices and posterior
#'   weights available in its quadrature summaries.
#' @param ridge Small ridge value used when inverting the Hessian.
#' @param ... Unused; included for method compatibility.
#'
#' @return A covariance matrix with attributes `bread` and `meat`.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' fit <- fit_mixed_subjects(
#'   observed, observed, simulate_2pl(rnorm(80), pars),
#'   lambda = 0.5, initial_pars = pars, n_quad = 7
#' )
#' dim(vcov_mixed_subjects(fit))
vcov_mixed_subjects <- function(object, ridge = 1e-8, ...) {
  summaries <- extract_q_for_vcov(object)
  q_observed <- summaries[[1]]
  q_predicted <- summaries[[2]]
  q_generated <- summaries[[3]]

  item_pars <- object$item_pars
  n_items <- nrow(item_pars)
  lambda <- validate_lambda(object$lambda)

  H_observed <- avg_hessian_counts(q_observed$counts, item_pars)
  H_predicted <- avg_hessian_counts(q_predicted$counts, item_pars)
  H_generated <- avg_hessian_counts(q_generated$counts, item_pars)
  bread <- H_observed + lambda * (H_generated - H_predicted)

  S_observed <- person_scores_2pl(q_observed$resp, q_observed$weights, item_pars)
  S_predicted <- person_scores_2pl(q_predicted$resp, q_predicted$weights, item_pars)
  S_generated <- person_scores_2pl(q_generated$resp, q_generated$weights, item_pars)

  if (nrow(S_observed) != nrow(S_predicted)) {
    stop("Observed and paired predicted score matrices must have the same rows.",
         call. = FALSE)
  }

  S_labeled <- S_observed - lambda * S_predicted
  meat <- safe_cov(S_labeled) / nrow(S_labeled) +
    lambda^2 * safe_cov(S_generated) / nrow(S_generated)

  bread_inv <- stable_inverse(bread, ridge = ridge)
  Sigma <- bread_inv %*% meat %*% t(bread_inv)
  Sigma <- (Sigma + t(Sigma)) / 2

  names <- c(paste0("a_", item_pars$item), paste0("d_", item_pars$item))
  dimnames(Sigma) <- list(names, names)
  attr(Sigma, "bread") <- bread
  attr(Sigma, "meat") <- meat
  Sigma
}

#' @export
vcov.mixedsubjects_fit <- function(object, ...) {
  vcov_mixed_subjects(object, ...)
}

#' Gradient of ML ability scores with respect to item parameters
#'
#' Computes the implicit derivative of bounded maximum-likelihood ability scores
#' with respect to 2PL item parameters. The column order is all discriminations
#' followed by all intercepts.
#'
#' @param resp Response matrix with rows for subjects and columns for items.
#' @param item_pars Item parameters in slope-intercept form, or a
#'   `"mixedsubjects_fit"` object.
#' @param theta Optional precomputed ability estimates. If omitted,
#'   [score_theta()] is used.
#' @param bounds Bounds passed to [score_theta()] when `theta` is omitted.
#' @param eps Tolerance used to mark near-zero test information as undefined.
#'
#' @return A matrix with one row per response pattern and one column per item
#'   parameter.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' ability_gradient(resp, pars)
ability_gradient <- function(resp, item_pars, theta = NULL,
                             bounds = c(-6, 6), eps = 1e-10) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  if (inherits(item_pars, "mixedsubjects_fit")) {
    item_pars <- item_pars$item_pars
  }
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )

  if (is.null(theta)) {
    theta <- score_theta(resp, item_pars, bounds = bounds)
  }
  if (length(theta) != nrow(resp)) {
    stop("theta must have length nrow(resp).", call. = FALSE)
  }

  n_items <- ncol(resp)
  grad <- matrix(NA_real_, nrow = nrow(resp), ncol = 2 * n_items)
  colnames(grad) <- c(paste0("a_", item_pars$item), paste0("d_", item_pars$item))

  for (i in seq_len(nrow(resp))) {
    observed <- which(!is.na(resp[i, ]))
    if (length(observed) == 0 || !is.finite(theta[i])) {
      next
    }

    eta <- item_pars$d[observed] + item_pars$a[observed] * theta[i]
    pr <- stats::plogis(eta)
    pq <- pr * (1 - pr)
    a <- item_pars$a[observed]
    d_score_d_theta <- -sum(a^2 * pq)

    if (abs(d_score_d_theta) < eps) {
      next
    }

    for (r in seq_along(observed)) {
      j <- observed[r]
      d_score_d_a <- resp[i, j] - pr[r] - item_pars$a[j] * pq[r] * theta[i]
      d_score_d_d <- -item_pars$a[j] * pq[r]

      grad[i, j] <- -d_score_d_a / d_score_d_theta
      grad[i, n_items + j] <- -d_score_d_d / d_score_d_theta
    }
    missing_cols <- setdiff(seq_len(n_items), observed)
    if (length(missing_cols) > 0) {
      grad[i, missing_cols] <- 0
      grad[i, n_items + missing_cols] <- 0
    }
  }

  grad
}

#' Propagated ability risk from item-parameter uncertainty
#'
#' Computes `g_i' Sigma g_i` for each response pattern, where `g_i` is the
#' gradient of the ability estimate with respect to item parameters. If
#' `theta_true` is supplied, the returned total risk also includes squared
#' ability estimation error.
#'
#' @param resp Target response matrix.
#' @param fit_or_pars A `"mixedsubjects_fit"` object or item-parameter data
#'   frame.
#' @param vcov Optional covariance matrix. Required when `fit_or_pars` is item
#'   parameters rather than a fitted mixed-subjects object.
#' @param theta_true Optional true theta values for simulation studies.
#' @param bounds Bounds passed to [score_theta()].
#'
#' @return A list with `summary` and per-pattern `details`.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- simulate_2pl(rnorm(30), pars)
#' Sigma <- diag(0.01, 4)
#' ability_risk(resp, pars, vcov = Sigma)$summary
ability_risk <- function(resp, fit_or_pars, vcov = NULL, theta_true = NULL,
                         bounds = c(-6, 6)) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  if (inherits(fit_or_pars, "mixedsubjects_fit")) {
    item_pars <- fit_or_pars$item_pars
    if (is.null(vcov)) {
      vcov <- vcov_mixed_subjects(fit_or_pars)
    }
  } else {
    item_pars <- fit_or_pars
    if (is.null(vcov)) {
      stop("vcov is required when fit_or_pars is not a fitted object.",
           call. = FALSE)
    }
  }

  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )
  theta_hat <- score_theta(resp, item_pars, bounds = bounds)
  grad <- ability_gradient(resp, item_pars, theta = theta_hat, bounds = bounds)

  param_var <- rep(NA_real_, nrow(resp))
  for (i in seq_len(nrow(resp))) {
    g <- grad[i, ]
    if (all(is.finite(g))) {
      param_var[i] <- max(as.numeric(t(g) %*% vcov %*% g), 0)
    }
  }

  squared_error <- rep(NA_real_, nrow(resp))
  if (!is.null(theta_true)) {
    if (length(theta_true) != nrow(resp)) {
      stop("theta_true must have length nrow(resp).", call. = FALSE)
    }
    squared_error <- (theta_hat - theta_true)^2
  }

  total_risk <- param_var
  if (!is.null(theta_true)) {
    total_risk <- total_risk + squared_error
  }

  details <- data.frame(
    theta_hat = theta_hat,
    param_var = param_var,
    squared_error = squared_error,
    total_risk = total_risk
  )

  list(
    summary = data.frame(
      mean_param_var = mean(param_var, na.rm = TRUE),
      mean_squared_error = if (is.null(theta_true)) NA_real_ else
        mean(squared_error, na.rm = TRUE),
      mean_total_risk = mean(total_risk, na.rm = TRUE)
    ),
    details = details,
    gradient = grad
  )
}

#' Plug-in PPI++ optimal tuning parameter
#'
#' Implements the closed-form estimator from Proposition 2 of Angelopoulos,
#' Duchi and Zrnic (2023) for the lambda that minimises the trace of the
#' asymptotic item-parameter covariance matrix `Tr(Sigma_gamma)`.
#'
#' **This is the item-parameter variance objective, not the psychometric
#' scoring objective.** For IRT applications where accurate ability scoring
#' is the goal, use [tune_lambda_ability()] or [tune_lambda_ability_crossfit()]
#' instead. Those functions directly minimise the propagated ability-score
#' risk `E[g' Sigma_gamma g]` — the quantity that matters for test scoring —
#' rather than item-parameter estimation efficiency. `tune_lambda_ppi_score()`
#' is provided as a theoretical diagnostic and to facilitate method validation.
#'
#' The formula uses the **same** human posterior weights for both the human and
#' paired-LLM score vectors. This symmetry is required for the PPI++
#' unbiasedness condition `E[grad_gen] = E[grad_pred]` at the true parameters.
#'
#' @param observed Human response matrix.
#' @param predicted Paired LLM responses for the same rows as `observed`.
#' @param item_pars Item parameters in slope-intercept form at which to
#'   evaluate the score vectors. Typically the human 2PL MLE from [fit_2pl()].
#' @param n_generated Number of generated (unpaired) LLM subjects, used to
#'   compute `r = n / n_generated`.
#' @param quadrature Optional quadrature grid. If omitted, a standard-normal
#'   grid with `n_quad` nodes is created.
#' @param n_quad Number of quadrature nodes when `quadrature` is omitted.
#'
#' @return A list with elements `lambda` (the plug-in estimate, clipped to
#'   \[0, 1\]), `n`, `n_generated`, `r`, and the intermediate matrices `C_hf`
#'   (cross-covariance of human and paired-LLM score vectors) and `V_f`
#'   (variance of paired-LLM score vectors).
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' predicted <- observed
#' tune_lambda_ppi_score(observed, predicted, pars, n_generated = 100, n_quad = 7)$lambda
tune_lambda_ppi_score <- function(observed, predicted, item_pars, n_generated,
                                   quadrature = NULL, n_quad = 31) {
  observed  <- validate_response_matrix(observed,  "observed",
                                        allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted",
                                        allow_fractional = TRUE)
  check_same_items(observed, predicted, "observed", "predicted")
  if (nrow(observed) != nrow(predicted)) {
    stop("observed and predicted must have the same number of rows.",
         call. = FALSE)
  }
  if (!is.numeric(n_generated) || length(n_generated) != 1 ||
      !is.finite(n_generated) || n_generated <= 0) {
    stop("n_generated must be a single positive finite number.", call. = FALSE)
  }

  item_pars <- standardize_item_pars(item_pars, n_items = ncol(observed),
                                      item_names = colnames(observed))
  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  # Compute human posterior weights from observed responses. The SAME weights
  # are used for both score vectors so that the cross-covariance C_hf correctly
  # represents the PPI++ correction variance, satisfying the unbiasedness
  # condition E[grad_gen] = E[grad_pred] at the true item parameters.
  weights <- posterior_weights_2pl(observed, item_pars, quadrature = quadrature)

  S_h <- person_scores_2pl(observed,  weights, item_pars)
  S_f <- person_scores_2pl(predicted, weights, item_pars)

  n <- nrow(observed)
  r <- n / n_generated

  counts <- summarize_expected_counts(observed, weights)
  H      <- avg_hessian_counts(counts, item_pars)
  H_inv  <- stable_inverse(H)

  mu_h <- colMeans(S_h)
  mu_f <- colMeans(S_f)
  C_hf <- t(S_h - mu_h) %*% (S_f - mu_f) / (n - 1)
  V_f  <- t(S_f - mu_f) %*% (S_f - mu_f) / (n - 1)

  # Proposition 2, Angelopoulos, Duchi & Zrnic (2023):
  # lambda* = Tr(H^{-1}(C_hf + C_hf')H^{-1}) / (2(1+r) Tr(H^{-1} V_f H^{-1}))
  sym_cov     <- C_hf + t(C_hf)
  numerator   <- sum(diag(H_inv %*% sym_cov   %*% H_inv))
  denominator <- 2 * (1 + r) * sum(diag(H_inv %*% V_f %*% H_inv))

  lambda <- if (denominator <= 0) 0 else max(0, min(1, numerator / denominator))

  list(
    lambda      = lambda,
    n           = n,
    n_generated = n_generated,
    r           = r,
    C_hf        = C_hf,
    V_f         = V_f
  )
}

#' Tune lambda by downstream ability-score risk
#'
#' Fits candidate mixed-subjects calibrations, estimates the item-parameter
#' sandwich covariance for each, and chooses the lambda that minimises average
#' propagated ability-score risk on a target response matrix.
#'
#' This function minimises `E[g' Sigma_gamma g]` — the propagated ability-score
#' risk — which is the appropriate objective for IRT applications where accurate
#' test scoring is the goal. This is **distinct** from [tune_lambda_ppi_score()],
#' which minimises the trace of the item-parameter covariance matrix
#' `Tr(Sigma_gamma)` (the PPI++ theoretical objective). The two criteria
#' generally yield different lambda values:
#'
#' - `tune_lambda_ability_risk()` asks: which lambda produces the most accurate
#'   ability scores for the target population? Use this for operational scoring.
#' - [tune_lambda_ppi_score()] asks: which lambda minimises item-parameter
#'   estimation variance? Use this for method validation and diagnostics.
#'
#' Diagnostic note: if `tune_lambda_ability_risk()` selects `lambda = 0` for a
#' misaligned LLM (one whose item parameters differ from the human calibration),
#' this is the correct mathematical outcome under the current fixed-posterior
#' expected-count implementation. The frozen posteriors create a gradient
#' asymmetry that inflates item parameters at any `lambda > 0`, increasing
#' ability risk. This is not a bug in the risk function; it is a property of the
#' estimating equations. See [fit_mixed_subjects_mml()] for a marginal-likelihood
#' implementation that removes this asymmetry.
#'
#' @param lambda_grid Numeric vector of candidate lambda values in `[0, 1]`.
#' @param observed,predicted,generated Response matrices passed to
#'   [fit_mixed_subjects()].
#' @param target_resp Response matrix defining the target scoring population. If
#'   omitted, `observed` is used.
#' @param theta_true Optional true theta values for `target_resp`, used in
#'   simulation studies to add squared scoring error to the risk. When omitted,
#'   `mean_squared_error` in the summary is `NA`; only `mean_param_var` is
#'   computed.
#' @param n_quad Number of quadrature nodes.
#' @param initial_pars Optional starting item parameters.
#' @param fit_fn Fitting function to use. Defaults to [fit_mixed_subjects()]
#'   (frozen expected-count). Pass [fit_mixed_subjects_mml()] to use the
#'   marginal-likelihood PPI++ estimator, which eliminates gradient asymmetry
#'   and should select `lambda > 0` for genuinely informative predictors.
#' @param bounds Bounds passed to [score_theta()].
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to `fit_fn`.
#'
#' @return A list with `summary`, `best_lambda`, `best_fit`, and all fitted
#'   candidate objects.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' generated <- simulate_2pl(rnorm(100), pars)
#' tuned <- tune_lambda_ability_risk(
#'   c(0, 0.5), observed, observed, generated,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$best_lambda
#' @seealso [tune_lambda_ppi_score()] for the PPI++ theoretical lambda that
#'   minimises the trace of the item-parameter covariance matrix;
#'   [fit_mixed_subjects_mml()] for the marginal-likelihood estimator.
tune_lambda_ability_risk <- function(lambda_grid, observed, predicted, generated,
                                     target_resp = NULL, theta_true = NULL,
                                     n_quad = 31, initial_pars = NULL,
                                     fit_fn = fit_mixed_subjects,
                                     bounds = c(-6, 6),
                                     control = list(maxit = 500), ...) {
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0) {
    stop("lambda_grid must be a non-empty numeric vector.", call. = FALSE)
  }
  lambda_grid <- sort(unique(lambda_grid))
  vapply(lambda_grid, validate_lambda, numeric(1))

  observed <- validate_response_matrix(observed, "observed",
                                       allow_fractional = FALSE)
  if (is.null(target_resp)) {
    target_resp <- observed
  } else {
    target_resp <- validate_response_matrix(target_resp, "target_resp",
                                            allow_fractional = TRUE)
    check_same_items(observed, target_resp, "observed", "target_resp")
  }

  rows <- vector("list", length(lambda_grid))
  fits <- vector("list", length(lambda_grid))
  risks <- vector("list", length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    fits[[i]] <- fit_fn(
      observed = observed,
      predicted = predicted,
      generated = generated,
      lambda = lambda,
      n_quad = n_quad,
      initial_pars = initial_pars,
      control = control,
      ...
    )

    Sigma <- vcov_mixed_subjects(fits[[i]])
    risks[[i]] <- ability_risk(
      target_resp,
      fits[[i]],
      vcov = Sigma,
      theta_true = theta_true,
      bounds = bounds
    )

    rows[[i]] <- data.frame(
      lambda = lambda,
      mean_param_var = risks[[i]]$summary$mean_param_var,
      mean_squared_error = risks[[i]]$summary$mean_squared_error,
      mean_total_risk = risks[[i]]$summary$mean_total_risk,
      convergence = fits[[i]]$convergence
    )
  }

  summary <- do.call(rbind, rows)
  best_idx <- which.min(summary$mean_total_risk)

  list(
    summary = summary,
    best_lambda = summary$lambda[best_idx],
    best_fit = fits[[best_idx]],
    fits = fits,
    risks = risks
  )
}

#' Cross-fit ability-score-risk lambda tuning
#'
#' Estimates lambda separately for each held-out split using only the remaining
#' labeled rows, then fits the final split-sample mixed-subjects estimator with
#' those fold-specific lambda values.
#'
#' @inheritParams tune_lambda_ability_risk
#' @param n_splits Number of sample splits.
#' @param split_id Optional integer split assignment for labeled rows.
#' @param seed Optional seed used when `split_id` is omitted.
#'
#' @return A list with fold-specific lambda values, fold tuning objects, and the
#'   final split-sample fit.
#' @export
tune_lambda_ability_risk_crossfit <- function(lambda_grid, observed, predicted,
                                              generated, target_resp = NULL,
                                              theta_true = NULL, n_splits = 2,
                                              split_id = NULL, seed = NULL,
                                              n_quad = 31, initial_pars = NULL,
                                              bounds = c(-6, 6),
                                              control = list(maxit = 500), ...) {
  observed <- validate_response_matrix(observed, "observed",
                                       allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted",
                                        allow_fractional = TRUE)
  generated <- validate_response_matrix(generated, "generated",
                                        allow_fractional = TRUE)
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")

  if (nrow(observed) != nrow(predicted)) {
    stop("observed and predicted must have the same number of rows.",
         call. = FALSE)
  }

  if (is.null(split_id)) {
    split_id <- make_split_id(nrow(observed), n_splits, seed = seed)
  } else {
    split_id <- as.integer(split_id)
    if (length(split_id) != nrow(observed)) {
      stop("split_id must have length nrow(observed).", call. = FALSE)
    }
  }

  if (is.null(target_resp)) {
    target_resp <- observed
  } else {
    target_resp <- validate_response_matrix(target_resp, "target_resp",
                                            allow_fractional = TRUE)
    check_same_items(observed, target_resp, "observed", "target_resp")
  }

  split_values <- sort(unique(split_id))
  lambda_by_split <- numeric(length(split_values))
  fold_tuning <- vector("list", length(split_values))

  for (s in seq_along(split_values)) {
    fold <- split_values[s]
    train <- split_id != fold

    theta_train <- if (is.null(theta_true)) {
      NULL
    } else {
      theta_true[train]
    }

    fold_tuning[[s]] <- tune_lambda_ability_risk(
      lambda_grid = lambda_grid,
      observed = observed[train, , drop = FALSE],
      predicted = predicted[train, , drop = FALSE],
      generated = generated,
      target_resp = target_resp[train, , drop = FALSE],
      theta_true = theta_train,
      n_quad = n_quad,
      initial_pars = initial_pars,
      bounds = bounds,
      control = control,
      ...
    )
    lambda_by_split[s] <- fold_tuning[[s]]$best_lambda
  }

  final_fit <- fit_mixed_subjects_split(
    observed = observed,
    predicted = predicted,
    generated = generated,
    lambda = lambda_by_split,
    split_id = split_id,
    n_quad = n_quad,
    initial_pars = initial_pars,
    control = control,
    ...
  )

  list(
    lambda_by_split = lambda_by_split,
    split_id = split_id,
    fold_tuning = fold_tuning,
    final_fit = final_fit
  )
}

#' Per-item PPI++ optimal tuning parameters
#'
#' Applies the PPI++ Proposition 2 plug-in formula independently for each item,
#' producing a vector of item-specific lambda values `λ_j ∈ [0, 1]`.
#'
#' The global [tune_lambda_ppi_score()] uses the full parameter covariance matrix
#' `Tr(Σ_γ)` as the objective. This function instead applies the same formula
#' using only the 2×2 diagonal block of the inverse Hessian for item `j`, and
#' the 2D sub-vectors of the human and paired-LLM score vectors. The result is
#' the λ that minimises the marginal variance of `(a_j, d_j)` independently for
#' each item.
#'
#' **Use case.** When a single global λ is forced to zero because a few items
#' have poor LLM predictions, per-item λ_j allows well-predicted items to still
#' benefit from the LLM data. Pass the returned vector to
#' [fit_mixed_subjects_mml()] as the `lambda` argument.
#'
#' This is a **theoretical diagnostic**: it minimises item-parameter variance,
#' not ability-score risk. For operational scoring use
#' [tune_lambda_ability_risk_item()] instead.
#'
#' @param observed Human response matrix.
#' @param predicted Paired LLM responses for the same rows as `observed`.
#' @param item_pars Item parameters at which to evaluate the score vectors.
#' @param n_generated Number of generated (unpaired) LLM subjects.
#' @param quadrature Optional quadrature grid.
#' @param n_quad Number of quadrature nodes when `quadrature` is omitted.
#'
#' @return A list with `lambda` (numeric vector of length `n_items`), `item`
#'   (item names), `n`, `n_generated`, and `r = n / n_generated`.
#' @export
#'
#' @seealso [tune_lambda_ppi_score()] for the global version;
#'   [fit_mixed_subjects_mml()] to fit with a per-item lambda vector.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' tune_lambda_ppi_score_item(observed, observed, pars, n_generated = 100, n_quad = 7)$lambda
tune_lambda_ppi_score_item <- function(observed, predicted, item_pars, n_generated,
                                        quadrature = NULL, n_quad = 31) {
  observed  <- validate_response_matrix(observed,  "observed",
                                        allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted",
                                        allow_fractional = TRUE)
  check_same_items(observed, predicted, "observed", "predicted")
  if (nrow(observed) != nrow(predicted)) {
    stop("observed and predicted must have the same number of rows.", call. = FALSE)
  }
  if (!is.numeric(n_generated) || length(n_generated) != 1 ||
      !is.finite(n_generated) || n_generated <= 0) {
    stop("n_generated must be a single positive finite number.", call. = FALSE)
  }

  item_pars  <- standardize_item_pars(item_pars, n_items = ncol(observed),
                                      item_names = colnames(observed))
  quadrature <- check_quadrature(quadrature, n_quad = n_quad)
  n_items    <- nrow(item_pars)

  # Human posterior weights — same for both score vectors (symmetry condition)
  weights <- posterior_weights_2pl(observed, item_pars, quadrature = quadrature)
  S_h     <- person_scores_2pl(observed,  weights, item_pars)
  S_f     <- person_scores_2pl(predicted, weights, item_pars)

  n   <- nrow(observed)
  r   <- n / n_generated

  # Full Hessian; use H^{-1} block for each item
  counts <- summarize_expected_counts(observed, weights)
  H      <- avg_hessian_counts(counts, item_pars)
  H_inv  <- stable_inverse(H)

  lambda_j <- numeric(n_items)

  for (j in seq_len(n_items)) {
    idx <- c(j, n_items + j)      # (a_j, d_j) columns in score matrices

    s_h_j <- S_h[, idx, drop = FALSE]   # n × 2
    s_f_j <- S_f[, idx, drop = FALSE]

    H_inv_j <- H_inv[idx, idx, drop = FALSE]   # 2 × 2

    mu_h_j  <- colMeans(s_h_j)
    mu_f_j  <- colMeans(s_f_j)
    C_hf_j  <- t(s_h_j - mu_h_j) %*% (s_f_j - mu_f_j) / (n - 1)
    V_f_j   <- t(s_f_j - mu_f_j) %*% (s_f_j - mu_f_j) / (n - 1)

    sym_cov_j   <- C_hf_j + t(C_hf_j)
    numerator_j   <- sum(diag(H_inv_j %*% sym_cov_j   %*% H_inv_j))
    denominator_j <- 2 * (1 + r) * sum(diag(H_inv_j %*% V_f_j %*% H_inv_j))

    lambda_j[j] <- if (denominator_j <= 0) 0 else
      max(0, min(1, numerator_j / denominator_j))
  }

  list(lambda = lambda_j, item = item_pars$item, n = n,
       n_generated = n_generated, r = r)
}

#' Per-item ability-risk lambda tuning via coordinate descent
#'
#' Finds a per-item vector of lambda values `λ_j ∈ [0, 1]` that minimises
#' propagated ability-score risk `E[g' Σ_γ g]` using coordinate descent on the
#' items. Each coordinate step selects the `λ_j` in `lambda_grid` that gives
#' the smallest mean ability risk while holding all other `λ_{j'}` fixed.
#'
#' Uses [fit_mixed_subjects_mml()] at each candidate, so posteriors are
#' recomputed from the current parameters (no frozen-posterior bias). The
#' resulting lambda vector can then be used directly with
#' [fit_mixed_subjects_mml()].
#'
#' **Computational cost.** Each pass evaluates `n_items × length(lambda_grid)`
#' fits. For `n_items = 8` and a 5-point grid this is 40 fits per pass. Use
#' `n_pass = 1` (the default) for a single greedy sweep, which is usually
#' sufficient.
#'
#' @param lambda_grid Numeric vector of candidate λ values in `[0, 1]` to try
#'   for each item independently.
#' @param observed,predicted,generated Response matrices passed to
#'   [fit_mixed_subjects_mml()].
#' @param target_resp Target scoring population. If omitted, `observed` is used.
#' @param theta_true Optional true theta values, used to add squared scoring
#'   error to the risk.
#' @param n_quad Number of quadrature nodes.
#' @param initial_pars Optional starting item parameters.
#' @param n_pass Number of coordinate-descent passes (default 1).
#' @param init_lambda Starting lambda vector for coordinate descent. Supply the
#'   global scalar optimum from [tune_lambda_ability_risk()] (e.g.
#'   `init_lambda = 0.5`) to start the search around a useful operating point
#'   rather than all-zeros. A scalar is broadcast to all items; a vector of
#'   length `n_items` sets per-item starting values.
#' @param bounds Bounds passed to [score_theta()].
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to [fit_mixed_subjects_mml()].
#'
#' @return A list with `lambda` (per-item vector), `item` (item names),
#'   `n_pass`, and `final_fit` (the [fit_mixed_subjects_mml()] fit at the
#'   selected lambda).
#' @export
#'
#' @seealso [tune_lambda_ppi_score_item()] for the faster PPI++-score version;
#'   [tune_lambda_ability_risk()] for the global scalar version.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed  <- simulate_2pl(rnorm(40), pars)
#' generated <- simulate_2pl(rnorm(100), pars)
#' tuned <- tune_lambda_ability_risk_item(
#'   c(0, 0.5), observed, observed, generated,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$lambda
tune_lambda_ability_risk_item <- function(lambda_grid, observed, predicted, generated,
                                          target_resp = NULL, theta_true = NULL,
                                          n_quad = 31, initial_pars = NULL,
                                          n_pass = 1,
                                          init_lambda = 0,
                                          bounds = c(-6, 6),
                                          control = list(maxit = 300), ...) {
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0) {
    stop("lambda_grid must be a non-empty numeric vector.", call. = FALSE)
  }
  lambda_grid <- sort(unique(lambda_grid))
  vapply(lambda_grid, validate_lambda, numeric(1))

  observed <- validate_response_matrix(observed, "observed", allow_fractional = FALSE)
  n_items  <- ncol(observed)
  if (is.null(target_resp)) target_resp <- observed

  # Initialise lambda vector. init_lambda can be a scalar (same for all items)
  # or a length-n_items vector. Passing the global scalar optimum from
  # tune_lambda_ability_risk() is recommended: the coordinate descent then
  # refines which items benefit from that level of correction.
  init_lambda <- validate_lambda_vector(init_lambda, n = n_items)
  lambda_vec  <- if (length(init_lambda) == 1) rep(init_lambda, n_items) else init_lambda

  for (pass in seq_len(n_pass)) {
    lambda_prev <- lambda_vec

    for (j in seq_len(n_items)) {
      best_risk_j <- Inf
      best_lam_j  <- lambda_vec[j]

      for (lam in lambda_grid) {
        lambda_try    <- lambda_vec
        lambda_try[j] <- lam

        fit_j <- tryCatch(
          fit_mixed_subjects_mml(
            observed     = observed,
            predicted    = predicted,
            generated    = generated,
            lambda       = lambda_try,
            n_quad       = n_quad,
            initial_pars = initial_pars,
            control      = control,
            ...
          ),
          error = function(e) NULL
        )
        if (is.null(fit_j)) next

        risk_j <- tryCatch({
          Sigma_j <- vcov_mixed_subjects(fit_j)
          ability_risk(target_resp, fit_j, vcov = Sigma_j,
                       theta_true = theta_true, bounds = bounds)$summary$mean_total_risk
        }, error = function(e) Inf)

        if (is.finite(risk_j) && risk_j < best_risk_j) {
          best_risk_j <- risk_j
          best_lam_j  <- lam
        }
      }
      lambda_vec[j] <- best_lam_j
    }

    if (max(abs(lambda_vec - lambda_prev)) < 1e-4) break
  }

  final_fit <- tryCatch(
    fit_mixed_subjects_mml(
      observed     = observed,
      predicted    = predicted,
      generated    = generated,
      lambda       = lambda_vec,
      n_quad       = n_quad,
      initial_pars = initial_pars,
      control      = control,
      ...
    ),
    error = function(e) NULL
  )

  list(
    lambda    = lambda_vec,
    item      = colnames(observed),
    n_pass    = n_pass,
    final_fit = final_fit
  )
}

#' @describeIn tune_lambda_ability_risk Deprecated name; use
#'   [tune_lambda_ability_risk()] instead.
#' @export
tune_lambda_ability <- function(lambda_grid, observed, predicted, generated,
                                target_resp = NULL, theta_true = NULL,
                                n_quad = 31, initial_pars = NULL,
                                bounds = c(-6, 6),
                                control = list(maxit = 500), ...) {
  .Deprecated(
    new   = "tune_lambda_ability_risk",
    msg   = paste(
      "'tune_lambda_ability' has been renamed 'tune_lambda_ability_risk'",
      "to make its objective (ability-score risk, not PPI++ score) explicit.",
      "Please update your code."
    )
  )
  tune_lambda_ability_risk(
    lambda_grid = lambda_grid,
    observed = observed, predicted = predicted, generated = generated,
    target_resp = target_resp, theta_true = theta_true,
    n_quad = n_quad, initial_pars = initial_pars,
    bounds = bounds, control = control, ...
  )
}

#' @describeIn tune_lambda_ability_risk_crossfit Deprecated name; use
#'   [tune_lambda_ability_risk_crossfit()] instead.
#' @export
tune_lambda_ability_crossfit <- function(lambda_grid, observed, predicted,
                                         generated, target_resp = NULL,
                                         theta_true = NULL, n_splits = 2,
                                         split_id = NULL, seed = NULL,
                                         n_quad = 31, initial_pars = NULL,
                                         bounds = c(-6, 6),
                                         control = list(maxit = 500), ...) {
  .Deprecated(
    new   = "tune_lambda_ability_risk_crossfit",
    msg   = paste(
      "'tune_lambda_ability_crossfit' has been renamed",
      "'tune_lambda_ability_risk_crossfit'.",
      "Please update your code."
    )
  )
  tune_lambda_ability_risk_crossfit(
    lambda_grid = lambda_grid,
    observed = observed, predicted = predicted, generated = generated,
    target_resp = target_resp, theta_true = theta_true,
    n_splits = n_splits, split_id = split_id, seed = seed,
    n_quad = n_quad, initial_pars = initial_pars,
    bounds = bounds, control = control, ...
  )
}
