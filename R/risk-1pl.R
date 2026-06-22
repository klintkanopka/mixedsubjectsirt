# ---------------------------------------------------------------------------- #
# 1PL risk and lambda-tuning functions
#
# The 1PL parameter vector is c(a_shared, d_1, ..., d_J) — length J+1.
# All sandwich-covariance and ability-risk functions mirror their 2PL
# counterparts; the only structural differences are:
#   - Score vectors are (J+1)-dimensional rather than 2J-dimensional.
#   - The Hessian is (J+1) × (J+1).
#   - The ability gradient is (J+1)-dimensional.
# ---------------------------------------------------------------------------- #

# ---------- Internal helpers ------------------------------------------------

person_scores_1pl <- function(resp, weights, item_pars) {
  # Per-person (J+1)-dimensional score vectors for the 1PL model.
  # Component 1  : shared discrimination score = sum over all items of a_j scores
  # Components 2..J+1: per-item difficulty scores (same as 2PL d-scores)
  scores_2pl <- person_scores_2pl(resp, weights, item_pars)
  n_items <- nrow(item_pars)

  score_a <- rowSums(scores_2pl[, seq_len(n_items), drop = FALSE])
  score_d <- scores_2pl[, n_items + seq_len(n_items), drop = FALSE]

  result <- cbind(a_shared = score_a, score_d)
  colnames(result) <- c("a_shared", paste0("d_", item_pars$item))
  result
}

avg_hessian_counts_1pl <- function(counts, item_pars) {
  # (J+1) × (J+1) expected information matrix for the 1PL model.
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = counts$n_items,
    item_names = counts$item_names
  )

  p <- counts$n_items + 1L
  H <- matrix(0, nrow = p, ncol = p)

  a <- item_pars$a[1] # shared discrimination

  for (j in seq_len(counts$n_items)) {
    eta <- item_pars$d[j] + a * counts$theta
    pr <- stats::plogis(eta)
    w <- counts$N[j, ] * pr * (1 - pr) / counts$n

    idx_d <- j + 1L

    H[1L, 1L] <- H[1L, 1L] + sum(w * counts$theta^2) # H[a, a]
    H[1L, idx_d] <- sum(w * counts$theta) # H[a, d_j]
    H[idx_d, 1L] <- H[1L, idx_d]
    H[idx_d, idx_d] <- sum(w) # H[d_j, d_j]
  }

  H
}


# ---------- Exported risk and covariance functions --------------------------

#' Sandwich covariance for a 1PL mixed-subjects fit
#'
#' Estimates the `(J+1) × (J+1)` sandwich covariance matrix for the shared
#' discrimination and per-item intercepts of a 1PL mixed-subjects calibration.
#'
#' @note **Bread approximation.** The bread uses `avg_hessian_counts_1pl()`,
#'   the EM complete-data Hessian for the 1PL model, rather than the Louis
#'   (1982) marginal observed-information correction implemented for 2PL in
#'   [vcov_mixed_subjects_mml()]. The EM bread over-states efficiency by
#'   ignoring missing information about theta. A Louis-corrected 1PL bread is
#'   planned for a future release.
#'
#' @param object A `"mixedsubjects_1pl_fit"` object from [fit_mixed_subjects_1pl()]
#'   or [fit_mixed_subjects_mml_1pl()].
#' @param ridge Ridge regularization for Hessian inversion.
#' @param ... Unused.
#'
#' @return A `(J+1) × (J+1)` covariance matrix.  Row/column names are
#'   `"a_shared"` and `"d_Item1"`, `"d_Item2"`, etc.
#' @export
vcov_mixed_subjects_1pl <- function(object, ridge = 1e-8, ...) {
  summaries <- extract_q_for_vcov(object)
  q_obs <- summaries[[1]]
  q_pred <- summaries[[2]]
  q_gen <- summaries[[3]]

  item_pars <- object$item_pars
  lambda <- validate_lambda(object$lambda)

  H_obs <- avg_hessian_counts_1pl(q_obs$counts, item_pars)
  H_pred <- avg_hessian_counts_1pl(q_pred$counts, item_pars)
  H_gen <- avg_hessian_counts_1pl(q_gen$counts, item_pars)
  bread <- H_obs + lambda * (H_gen - H_pred)

  S_obs <- person_scores_1pl(q_obs$resp, q_obs$weights, item_pars)
  S_pred <- person_scores_1pl(q_pred$resp, q_pred$weights, item_pars)
  S_gen <- person_scores_1pl(q_gen$resp, q_gen$weights, item_pars)

  if (nrow(S_obs) != nrow(S_pred)) {
    stop(
      "Observed and predicted score matrices must have the same rows.",
      call. = FALSE
    )
  }

  S_labeled <- S_obs - lambda * S_pred
  meat <- safe_cov(S_labeled) /
    nrow(S_labeled) +
    lambda^2 * safe_cov(S_gen) / nrow(S_gen)

  bread_inv <- stable_inverse(bread, ridge = ridge)
  Sigma <- bread_inv %*% meat %*% t(bread_inv)
  Sigma <- (Sigma + t(Sigma)) / 2

  nms <- c("a_shared", paste0("d_", item_pars$item))
  dimnames(Sigma) <- list(nms, nms)
  attr(Sigma, "bread") <- bread
  attr(Sigma, "meat") <- meat
  Sigma
}

#' @export
vcov.mixedsubjects_1pl_fit <- function(object, ...) {
  vcov_mixed_subjects_1pl(object, ...)
}

#' Gradient of ML ability scores w.r.t. 1PL item parameters
#'
#' Computes the implicit derivative of bounded maximum-likelihood ability
#' scores with respect to the 1PL parameters `(a_shared, d_1, ..., d_J)`.
#'
#' The gradient for the shared discrimination is the sum of the per-item
#' discrimination gradients:
#' `da_shared = sum_j da_j` (chain rule via the constraint `a_j = a_shared`).
#'
#' @param resp Response matrix.
#' @param item_pars Item parameters with all `a` equal (1PL), or a
#'   `"mixedsubjects_1pl_fit"` object.
#' @param theta Optional precomputed ability estimates.
#' @param bounds Bounds passed to [score_theta()].
#' @param eps Tolerance for near-zero test information.
#'
#' @return A matrix with one row per response pattern and `J + 1` columns
#'   (`a_shared`, then one column per item's `d_j`).
#' @export
ability_gradient_1pl <- function(
  resp,
  item_pars,
  theta = NULL,
  bounds = c(-6, 6),
  eps = 1e-10
) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  if (inherits(item_pars, "mixedsubjects_1pl_fit")) {
    item_pars <- item_pars$item_pars
  }
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )

  # Compute the standard 2PL gradient (w.r.t. all 2J parameters)
  grad_2pl <- ability_gradient(
    resp,
    item_pars,
    theta = theta,
    bounds = bounds,
    eps = eps
  )
  n_items <- nrow(item_pars)

  # Map to 1PL parameterization: da_shared = sum_j da_j
  grad_a_shared <- rowSums(grad_2pl[, seq_len(n_items), drop = FALSE])
  grad_d <- grad_2pl[, n_items + seq_len(n_items), drop = FALSE]

  result <- cbind(a_shared = grad_a_shared, grad_d)
  colnames(result) <- c("a_shared", paste0("d_", item_pars$item))
  result
}

#' Propagated ability risk for a 1PL fit
#'
#' Computes `g_i' Sigma_1pl g_i` for each response pattern, where `g_i` is
#' the `(J+1)`-dimensional gradient of the ability estimate with respect to
#' `(a_shared, d_1, ..., d_J)` and `Sigma_1pl` is the sandwich covariance
#' from [vcov_mixed_subjects_1pl()].
#'
#' @param resp Target response matrix.
#' @param fit_or_pars A `"mixedsubjects_1pl_fit"` object or item-parameter
#'   data frame.
#' @param vcov Optional `(J+1) × (J+1)` covariance matrix.  Required when
#'   `fit_or_pars` is not a fitted object.
#' @param theta_true Optional true theta values for simulation studies.
#' @param bounds Bounds passed to [score_theta()].
#'
#' @return A list with `summary` and per-pattern `details`, the same structure
#'   as [ability_risk()].
#' @export
ability_risk_1pl <- function(
  resp,
  fit_or_pars,
  vcov = NULL,
  theta_true = NULL,
  bounds = c(-6, 6)
) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)

  if (inherits(fit_or_pars, "mixedsubjects_1pl_fit")) {
    item_pars <- fit_or_pars$item_pars
    if (is.null(vcov)) {
      vcov <- vcov_mixed_subjects_1pl(fit_or_pars)
    }
  } else {
    item_pars <- fit_or_pars
    if (is.null(vcov)) {
      stop(
        "vcov is required when fit_or_pars is not a fitted object.",
        call. = FALSE
      )
    }
  }

  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )
  theta_hat <- score_theta(resp, item_pars, bounds = bounds)
  grad <- ability_gradient_1pl(
    resp,
    item_pars,
    theta = theta_hat,
    bounds = bounds
  )

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

  list(
    summary = data.frame(
      mean_param_var = mean(param_var, na.rm = TRUE),
      mean_squared_error = if (is.null(theta_true)) {
        NA_real_
      } else {
        mean(squared_error, na.rm = TRUE)
      },
      mean_total_risk = mean(total_risk, na.rm = TRUE)
    ),
    details = data.frame(
      theta_hat = theta_hat,
      param_var = param_var,
      squared_error = squared_error,
      total_risk = total_risk
    ),
    gradient = grad
  )
}


# ---------- Lambda tuning ---------------------------------------------------

#' Plug-in PPI++ optimal tuning parameter for a 1PL model
#'
#' Applies the PPI++ Proposition 2 formula using `(J+1)`-dimensional score
#' vectors for the 1PL parameterization `(a_shared, d_1, ..., d_J)`.
#'
#' This is the **item-parameter variance** objective — it minimizes
#' `Tr(Sigma_1pl)`. For practical scoring applications use
#' [tune_lambda_ability_risk_1pl()] instead.
#'
#' @inheritParams tune_lambda_ppi_score
#'
#' @return A list with `lambda`, `n`, `n_generated`, `r`, `C_hf`, `V_f`.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
#' obs  <- simulate_2pl(rnorm(40), pars)
#' tune_lambda_ppi_score_1pl(obs, obs, pars, n_generated = 100, n_quad = 7)$lambda
tune_lambda_ppi_score_1pl <- function(
  observed,
  predicted,
  item_pars,
  n_generated,
  quadrature = NULL,
  n_quad = 31
) {
  observed <- validate_response_matrix(
    observed,
    "observed",
    allow_fractional = FALSE
  )
  predicted <- validate_response_matrix(
    predicted,
    "predicted",
    allow_fractional = FALSE
  )
  check_same_items(observed, predicted, "observed", "predicted")
  if (nrow(observed) != nrow(predicted)) {
    stop(
      "observed and predicted must have the same number of rows.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(n_generated) ||
      length(n_generated) != 1 ||
      !is.finite(n_generated) ||
      n_generated <= 0
  ) {
    stop("n_generated must be a single positive finite number.", call. = FALSE)
  }

  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(observed),
    item_names = colnames(observed)
  )
  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  weights <- posterior_weights_2pl(observed, item_pars, quadrature = quadrature)
  S_h <- person_scores_1pl(observed, weights, item_pars)
  S_f <- person_scores_1pl(predicted, weights, item_pars)

  n <- nrow(observed)
  r <- n / n_generated

  counts <- summarize_expected_counts(observed, weights)
  H <- avg_hessian_counts_1pl(counts, item_pars)
  H_inv <- stable_inverse(H)

  mu_h <- colMeans(S_h)
  mu_f <- colMeans(S_f)
  C_hf <- t(S_h - mu_h) %*% (S_f - mu_f) / (n - 1)
  V_f <- t(S_f - mu_f) %*% (S_f - mu_f) / (n - 1)

  sym_cov <- C_hf + t(C_hf)
  numerator <- sum(diag(H_inv %*% sym_cov %*% H_inv))
  denominator <- 2 * (1 + r) * sum(diag(H_inv %*% V_f %*% H_inv))

  lambda <- if (denominator <= 0) 0 else max(0, min(1, numerator / denominator))

  list(
    lambda = lambda,
    n = n,
    n_generated = n_generated,
    r = r,
    C_hf = C_hf,
    V_f = V_f
  )
}

#' Tune lambda by downstream ability-score risk for a 1PL model
#'
#' Selects the lambda minimizing `E[g' Sigma_1pl g]` — the propagated
#' ability-score risk in the 1PL parameterization — using
#' [fit_mixed_subjects_mml_1pl()] by default. As in the 2PL
#' [tune_lambda_ability_risk()], lambda is chosen by direct 1-D optimization
#' (`method = "optimize"`, the default) or over `lambda_grid`
#' (`method = "grid"`).
#'
#' Passes `fit_fn` to allow switching between the frozen expected-count
#' estimator ([fit_mixed_subjects_1pl()]) and the marginal-MML estimator
#' ([fit_mixed_subjects_mml_1pl()]).
#'
#' @inheritParams tune_lambda_ability_risk
#' @param fit_fn Fitting function. Defaults to [fit_mixed_subjects_mml_1pl()].
#'
#' @return A list with `summary`, `best_lambda`, `best_fit`, `fits`, `risks`.
#' @export
#'
#' @seealso [tune_lambda_ability_risk()] for the 2PL version;
#'   [tune_lambda_ppi_score_1pl()] for the PPI++ score diagnostic.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
#' obs  <- simulate_2pl(rnorm(40), pars)
#' gen  <- simulate_2pl(rnorm(100), pars)
#' tuned <- tune_lambda_ability_risk_1pl(
#'   c(0, 0.5), obs, obs, gen,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$best_lambda
tune_lambda_ability_risk_1pl <- function(
  lambda_grid = seq(0, 1, by = 0.1),
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
  n_quad = 31,
  initial_pars = NULL,
  fit_fn = fit_mixed_subjects_mml_1pl,
  method = c("optimize", "grid"),
  bounds = c(-6, 6),
  max_discrimination = 10,
  control = list(maxit = 500),
  ...
) {
  method <- match.arg(method)
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0) {
    stop("lambda_grid must be a non-empty numeric vector.", call. = FALSE)
  }
  lambda_grid <- sort(unique(lambda_grid))
  vapply(lambda_grid, validate_lambda, numeric(1))

  observed <- validate_response_matrix(
    observed,
    "observed",
    allow_fractional = FALSE
  )
  if (is.null(target_resp)) {
    target_resp <- observed
  } else {
    target_resp <- validate_response_matrix(
      target_resp,
      "target_resp",
      allow_fractional = TRUE
    )
    check_same_items(observed, target_resp, "observed", "target_resp")
  }

  inf_risk <- list(summary = data.frame(mean_param_var = Inf,
                                        mean_squared_error = Inf,
                                        mean_total_risk = Inf))

  # Memoized per-lambda evaluation, robust to failed 1PL fits (tryCatch).
  cache <- new.env(parent = emptyenv())
  evaluate <- function(lambda) {
    lambda <- validate_lambda(lambda)
    key <- sprintf("%.12f", lambda)
    hit <- get0(key, envir = cache, inherits = FALSE, ifnotfound = NULL)
    if (!is.null(hit)) return(hit)

    fit <- tryCatch(
      fit_fn(observed = observed, predicted = predicted, generated = generated,
             lambda = lambda, n_quad = n_quad, initial_pars = initial_pars,
             control = control, ...),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      risk <- inf_risk; conv <- 99L; max_disc <- Inf
    } else {
      Sigma <- tryCatch(vcov_mixed_subjects_1pl(fit), error = function(e) NULL)
      risk <- if (is.null(Sigma)) inf_risk else tryCatch(
        ability_risk_1pl(target_resp, fit, vcov = Sigma,
                         theta_true = theta_true, bounds = bounds),
        error = function(e) inf_risk)
      conv <- fit$convergence
      max_disc <- max(abs(fit$item_pars$a))
    }
    sel <- risk$summary$mean_total_risk
    if (!is.finite(sel) || conv != 0 || max_disc > max_discrimination) sel <- Inf
    out <- list(lambda = lambda, fit = fit, risk = risk, max_disc = max_disc,
                convergence = conv, selection_risk = sel)
    assign(key, out, envir = cache)
    out
  }

  lo <- min(lambda_grid)
  hi <- max(lambda_grid)

  if (method == "grid") {
    invisible(lapply(lambda_grid, evaluate))
  } else {
    penalty <- 1e12
    obj <- function(lambda) {
      s <- evaluate(lambda)$selection_risk
      if (is.finite(s)) s else penalty
    }
    evaluate(lo)
    evaluate(hi)
    if (hi > lo) stats::optimize(obj, interval = c(lo, hi), tol = 1e-3)
  }

  evals <- mget(ls(cache), envir = cache)
  rows <- lapply(evals, function(e) data.frame(
    lambda = e$lambda,
    mean_param_var = e$risk$summary$mean_param_var,
    mean_squared_error = e$risk$summary$mean_squared_error,
    mean_total_risk = e$risk$summary$mean_total_risk,
    convergence = e$convergence,
    max_disc = e$max_disc,
    selection_risk = e$selection_risk
  ))
  summary <- do.call(rbind, rows)
  summary <- summary[order(summary$lambda), , drop = FALSE]
  rownames(summary) <- NULL

  if (all(is.infinite(summary$selection_risk))) {
    warning(
      "No lambda candidate converged with finite ability-score risk. ",
      "Returning lambda = ", lo, " (human-only when lower bound is 0).",
      call. = FALSE
    )
    best <- evaluate(lo)
  } else {
    best <- evaluate(summary$lambda[which.min(summary$selection_risk)])
  }

  fits  <- lapply(summary$lambda, function(l) evaluate(l)$fit)
  risks <- lapply(summary$lambda, function(l) evaluate(l)$risk)

  list(
    summary = summary,
    best_lambda = best$lambda,
    best_fit = best$fit,
    fits = fits,
    risks = risks,
    method = method
  )
}
