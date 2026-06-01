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

#' Tune lambda by downstream ability risk
#'
#' Fits candidate mixed-subjects calibrations, estimates the item-parameter
#' sandwich covariance for each, and chooses the lambda minimizing average
#' propagated ability risk on a target response matrix.
#'
#' @param lambda_grid Numeric vector of candidate lambda values in `[0, 1]`.
#' @param observed,predicted,generated Response matrices passed to
#'   [fit_mixed_subjects()].
#' @param target_resp Response matrix defining the target scoring population. If
#'   omitted, `observed` is used.
#' @param theta_true Optional true theta values for `target_resp`, used in
#'   simulation studies to add squared scoring error to the risk.
#' @param n_quad Number of quadrature nodes.
#' @param initial_pars Optional starting item parameters.
#' @param bounds Bounds passed to [score_theta()].
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to [fit_mixed_subjects()].
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
#' tuned <- tune_lambda_ability(
#'   c(0, 0.5), observed, observed, generated,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$best_lambda
tune_lambda_ability <- function(lambda_grid, observed, predicted, generated,
                                target_resp = NULL, theta_true = NULL,
                                n_quad = 31, initial_pars = NULL,
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
    fits[[i]] <- fit_mixed_subjects(
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

#' Cross-fit ability-risk lambda tuning
#'
#' Estimates lambda separately for each held-out split using only the remaining
#' labeled rows, then fits the final split-sample mixed-subjects estimator with
#' those fold-specific lambda values.
#'
#' @inheritParams tune_lambda_ability
#' @param n_splits Number of sample splits.
#' @param split_id Optional integer split assignment for labeled rows.
#' @param seed Optional seed used when `split_id` is omitted.
#'
#' @return A list with fold-specific lambda values, fold tuning objects, and the
#'   final split-sample fit.
#' @export
tune_lambda_ability_crossfit <- function(lambda_grid, observed, predicted,
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

    fold_tuning[[s]] <- tune_lambda_ability(
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
