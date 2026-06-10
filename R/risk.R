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
  q_observed  <- summaries[[1]]
  q_predicted <- summaries[[2]]
  q_generated <- summaries[[3]]

  item_pars <- object$item_pars
  n_items   <- nrow(item_pars)
  lambda    <- validate_lambda_vector(object$lambda, n = n_items)
  scalar_lam <- length(lambda) == 1

  H_observed  <- avg_hessian_counts(q_observed$counts,  item_pars)
  H_predicted <- avg_hessian_counts(q_predicted$counts, item_pars)
  H_generated <- avg_hessian_counts(q_generated$counts, item_pars)

  if (scalar_lam) {
    # Standard scalar path
    bread <- H_observed + lambda * (H_generated - H_predicted)
  } else {
    # Vector-lambda path: item j's 2×2 diagonal block uses λ_j (not λ_j²).
    # Build the bread explicitly block-by-block to avoid incorrect sweep scaling.
    H_corr <- H_generated - H_predicted
    bread   <- H_observed
    for (j in seq_len(n_items)) {
      idx <- c(j, n_items + j)
      bread[idx, idx] <- H_observed[idx, idx] + lambda[j] * H_corr[idx, idx]
    }
  }

  S_observed  <- person_scores_2pl(q_observed$resp,  q_observed$weights,  item_pars)
  S_predicted <- person_scores_2pl(q_predicted$resp, q_predicted$weights, item_pars)
  S_generated <- person_scores_2pl(q_generated$resp, q_generated$weights, item_pars)

  if (nrow(S_observed) != nrow(S_predicted)) {
    stop("Observed and paired predicted score matrices must have the same rows.",
         call. = FALSE)
  }

  if (scalar_lam) {
    S_labeled <- S_observed - lambda * S_predicted
    meat <- safe_cov(S_labeled) / nrow(S_labeled) +
      lambda^2 * safe_cov(S_generated) / nrow(S_generated)
  } else {
    # Item j's score columns are scaled by λ_j
    lambda_2j <- rep(lambda, 2L)
    S_labeled <- S_observed - sweep(S_predicted, 2L, lambda_2j, `*`)
    # Meat: Cov(S_h - Λ S_p)/n + Cov(S_g Λ)/N
    S_gen_scaled <- sweep(S_generated, 2L, lambda_2j, `*`)
    meat <- safe_cov(S_labeled) / nrow(S_labeled) +
      safe_cov(S_gen_scaled) / nrow(S_generated)
  }

  bread_inv <- stable_inverse(bread, ridge = ridge)
  Sigma <- bread_inv %*% meat %*% t(bread_inv)
  Sigma <- (Sigma + t(Sigma)) / 2

  nms <- c(paste0("a_", item_pars$item), paste0("d_", item_pars$item))
  dimnames(Sigma) <- list(nms, nms)
  attr(Sigma, "bread") <- bread
  attr(Sigma, "meat")  <- meat
  Sigma
}

louis_missing_info <- function(resp, weights, item_pars) {
  # Computes the Louis (1982) missing-information matrix for the 2PL marginal
  # IRT likelihood.  Used by vcov_mixed_subjects_mml() to form the marginal
  # observed-information bread via A_marg = H_comp - I_miss.
  #
  # For person i with posterior weights w_ik and responses y_i, the
  # complete-data score at quadrature node k is
  #   s_ik^{a_j} = (y_ij - p_j(θ_k)) * θ_k
  #   s_ik^{d_j} = (y_ij - p_j(θ_k))
  #
  # Louis' identity gives the per-person missing information:
  #   M_i = Var_{w_i}(s_ik) = E_{w_i}[s_ik s_ik'] - E_{w_i}[s_ik] E_{w_i}[s_ik]'
  #
  # Averaged across persons: I_miss = (1/n) Σ_i M_i.
  #
  # The decomposition into J×J blocks (aa, ad, dd) is
  #   I_miss_{jl}^{αβ} = sm_αβ[j,l] - S_α[,j]' S_β[,l] / n
  # where sm_αβ is the second-moment matrix and S_α are the marginal scores.
  #
  # All computations are vectorised over persons and quadrature nodes.

  n       <- nrow(resp)
  n_items <- nrow(item_pars)
  theta   <- attr(weights, "theta")
  K       <- length(theta)

  # p[j, k] = P(correct | θ_k; γ)
  eta <- outer(item_pars$a, theta, `*`) +
    matrix(item_pars$d, nrow = n_items, ncol = K)
  p   <- stats::plogis(eta)

  W        <- as.matrix(weights)
  W_theta  <- sweep(W, 2L, theta,   `*`)
  W_theta2 <- sweep(W, 2L, theta^2, `*`)

  # Per-person expected θ^α
  E1 <- drop(W_theta  %*% rep(1, K))
  E2 <- drop(W_theta2 %*% rep(1, K))

  # Per-person E_w[p_l(θ) * θ^α]: n × J
  EP0 <- W        %*% t(p)
  EP1 <- W_theta  %*% t(p)
  EP2 <- W_theta2 %*% t(p)

  Y        <- as.matrix(resp)
  Y[is.na(Y)] <- 0   # treat missing as 0 (score = 0 for missing items)
  A_k <- colMeans(W)           # estimated marginal node weights

  # sm_block computes the J×J second-moment matrix for block (α):
  #   (1/n) Σ_i Σ_k w_ik r_ij(k) r_il(k) θ_k^α
  # where r_ij(k) = y_ij - p_j(θ_k), with NA responses treated as y_ij = 0.
  # NOTE: this approximation is exact for complete data.  For responses with
  # item-level missingness, it is approximately correct when missingness is
  # sparse and uncorrelated with item difficulty.  See vcov_mixed_subjects_mml()
  # for a guard that warns when missingness is present.
  sm_block <- function(EW_alpha, EP_alpha, A_alpha) {
    yy <- crossprod(Y * EW_alpha, Y) / n
    yp <- crossprod(Y, EP_alpha)     / n
    pp <- tcrossprod(sweep(p, 2L, A_alpha, `*`), p)
    yy - yp - t(yp) + pp
  }

  sm_aa <- sm_block(E2,        EP2, A_k * theta^2)
  sm_ad <- sm_block(E1,        EP1, A_k * theta)
  sm_dd <- sm_block(rep(1, n), EP0, A_k)

  # Marginal per-person scores (= E_{w_i}[s_ik]): already in person_scores_2pl
  S   <- person_scores_2pl(resp, weights, item_pars)
  S_a <- S[, seq_len(n_items),           drop = FALSE]
  S_d <- S[, n_items + seq_len(n_items), drop = FALSE]

  # I_miss = second_moment - outer_product_of_marginal_scores
  I_aa <- sm_aa - crossprod(S_a)     / n
  I_ad <- sm_ad - crossprod(S_a, S_d) / n
  I_dd <- sm_dd - crossprod(S_d)     / n

  dim_p <- 2L * n_items
  I_miss <- matrix(0, dim_p, dim_p)
  ia <- seq_len(n_items)
  id <- n_items + seq_len(n_items)

  I_miss[ia, ia] <- I_aa
  I_miss[ia, id] <- I_ad
  I_miss[id, ia] <- t(I_ad)
  I_miss[id, id] <- I_dd

  I_miss
}

#' Marginal-MML sandwich covariance for a mixed-subjects fit
#'
#' Computes the full sandwich covariance for the scalar marginal-MML PPI++
#' estimator from [fit_mixed_subjects_mml()].  The bread uses Louis's (1982)
#' observed marginal-information formula
#'
#' \deqn{A_\lambda^\mathrm{marg} = H_\lambda^\mathrm{comp} - I_\lambda^\mathrm{miss}}
#'
#' rather than the EM/complete-data Hessian used by [vcov_mixed_subjects()].
#' Using the complete-data Hessian as the bread for a marginal-MML estimator
#' would over-state efficiency by ignoring the missing-information correction.
#'
#' The meat uses the standard marginal per-person score vectors (posteriors at
#' the converged parameters), which is identical to [vcov_mixed_subjects()].
#'
#' **When is this function called automatically?** The `vcov()` method for
#' `"mixedsubjects_fit"` objects (see [stats::vcov()]) dispatches here whenever
#' `isTRUE(object$mml) && length(object$lambda) == 1`. For vector-lambda fits, or
#' for frozen expected-count fits, the existing [vcov_mixed_subjects()] is used.
#'
#' @param object A scalar-lambda [fit_mixed_subjects_mml()] fit.
#' @param ridge Ridge regularization for bread inversion.
#' @param ... Unused.
#'
#' @return A \eqn{2J \times 2J} covariance matrix with attributes `bread` and `meat`.
#' @export
#'
#' @seealso [vcov_mixed_subjects()] for the frozen expected-count version. The
#'   internal `louis_missing_info()` helper computes the missing-information
#'   correction.
vcov_mixed_subjects_mml <- function(object, ridge = 1e-8, ...) {
  summaries  <- extract_q_for_vcov(object)
  q_observed  <- summaries[[1]]
  q_predicted <- summaries[[2]]
  q_generated <- summaries[[3]]

  # Guard: louis_missing_info() treats NA responses as y=0, which is exact
  # only for complete data.  Warn when missingness is present so users are
  # aware that the Louis bread is approximate under item-level missing data.
  has_na <- any(
    is.na(q_observed$resp),
    is.na(q_predicted$resp),
    is.na(q_generated$resp)
  )
  if (has_na) {
    warning(
      "vcov_mixed_subjects_mml(): one or more response matrices contain NA. ",
      "louis_missing_info() treats missing responses as y = 0, which gives ",
      "an approximate (not exact) Louis bread under item-level missingness. ",
      "Use vcov_mixed_subjects() for exact results under the frozen EC estimator.",
      call. = FALSE
    )
  }

  item_pars  <- object$item_pars
  n_items    <- nrow(item_pars)
  lambda     <- validate_lambda(object$lambda)

  # Marginal observed-information bread: H_comp - I_miss (Louis formula)
  louis_bread <- function(q) {
    H    <- avg_hessian_counts(q$counts, item_pars)
    Imiss <- louis_missing_info(q$resp, q$weights, item_pars)
    H - Imiss
  }

  bread <- louis_bread(q_observed) +
    lambda * (louis_bread(q_generated) - louis_bread(q_predicted))

  # Marginal per-person score vectors (posteriors at converged parameters)
  S_observed  <- person_scores_2pl(q_observed$resp,  q_observed$weights,  item_pars)
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
  Sigma     <- bread_inv %*% meat %*% t(bread_inv)
  Sigma     <- (Sigma + t(Sigma)) / 2

  nms <- c(paste0("a_", item_pars$item), paste0("d_", item_pars$item))
  dimnames(Sigma) <- list(nms, nms)
  attr(Sigma, "bread") <- bread
  attr(Sigma, "meat")  <- meat
  Sigma
}

#' @export
vcov.mixedsubjects_fit <- function(object, ...) {
  # Dispatch:
  # - Scalar MML fits → Louis-corrected marginal bread (vcov_mixed_subjects_mml)
  # - All other fits → EM Hessian bread (vcov_mixed_subjects)
  #
  # Split-sample fits (from fit_mixed_subjects_split / tune_lambda_ability_risk_crossfit)
  # store pooled expected counts but NOT raw response matrices or posterior
  # weights, so vcov_mixed_subjects will error.  To get ability uncertainty for
  # a split final fit, re-fit with fit_mixed_subjects_mml (scalar lambda, full
  # sample) and call vcov() on that fit.
  if (!is.null(object$split)) {
    stop(
      "vcov() is not yet available for split-sample fits. ",
      "Re-fit with fit_mixed_subjects_mml(lambda = mean(lambda_by_split)) ",
      "on the full sample to obtain a covariance estimate.",
      call. = FALSE
    )
  }
  if (isTRUE(object$mml) && length(object$lambda) == 1L) {
    vcov_mixed_subjects_mml(object, ...)
  } else {
    vcov_mixed_subjects(object, ...)
  }
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

  boundary_tol <- sqrt(.Machine$double.eps)   # ~1.5e-8

  for (i in seq_len(nrow(resp))) {
    observed <- which(!is.na(resp[i, ]))
    if (length(observed) == 0 || !is.finite(theta[i])) {
      next
    }

    # Implicit-differentiation gradient is valid only for interior optima.
    # At a boundary estimate (all-correct or all-incorrect response patterns),
    # theta_hat equals the bound, the score equation does not hold, and the
    # gradient is theoretically undefined.  Set to NA to signal this.
    if (abs(theta[i] - bounds[1]) < boundary_tol ||
        abs(theta[i] - bounds[2]) < boundary_tol) {
      next     # row stays NA — excluded from mean_param_var via na.rm = TRUE
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
      # Use the S3 generic so MML fits get vcov_mixed_subjects_mml (Louis bread)
      # and frozen EC fits get vcov_mixed_subjects (EM bread).
      vcov <- stats::vcov(fit_or_pars)
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
#' Duchi and Zrnic (2023) for the lambda that minimizes the trace of the
#' asymptotic item-parameter covariance matrix `Tr(Sigma_gamma)`.
#'
#' **This is the item-parameter variance objective, not the psychometric
#' scoring objective.** For IRT applications where accurate ability scoring
#' is the goal, use [tune_lambda_ability_risk()] or
#' [tune_lambda_ability_risk_crossfit()] instead. Those functions directly
#' minimize the propagated ability-score risk `E[g' Sigma_gamma g]` — the
#' quantity that matters for test scoring — rather than item-parameter
#' estimation efficiency. `tune_lambda_ppi_score()`
#' is provided as a theoretical diagnostic and to facilitate method validation.
#'
#' The formula uses the **same** human posterior weights for both the human and
#' paired-LLM score vectors. This symmetry is required for the PPI++
#' unbiasedness condition `E[grad_gen] = E[grad_pred]` at the true parameters.
#'
#' @param observed Human response matrix.
#' @param predicted Paired binary LLM responses (0/1) for the same rows as
#'   `observed`. Probabilities are not accepted; sample binary responses first.
#' @param item_pars Item parameters in slope-intercept form at which to
#'   evaluate the score vectors. Typically the human 2PL MLE from [fit_2pl()].
#' @param n_generated Number of generated (unpaired) LLM subjects, used to
#'   compute the ratio `r` (`n / n_generated`).
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
                                        allow_fractional = FALSE)
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
#' sandwich covariance for each, and chooses the lambda that minimizes average
#' propagated ability-score risk on a target response matrix.
#'
#' This function minimizes `E[g' Sigma_gamma g]` — the propagated ability-score
#' risk — which is the appropriate objective for IRT applications where accurate
#' test scoring is the goal. This is **distinct** from [tune_lambda_ppi_score()],
#' which minimizes the trace of the item-parameter covariance matrix
#' `Tr(Sigma_gamma)` (the PPI++ theoretical objective). The two criteria
#' generally yield different lambda values:
#'
#' - `tune_lambda_ability_risk()` asks: which lambda produces the most accurate
#'   ability scores for the target population? Use this for operational scoring.
#' - [tune_lambda_ppi_score()] asks: which lambda minimizes item-parameter
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
#' @param max_discrimination Upper bound on plausible item discrimination. Any
#'   candidate fit whose maximum `|a|` exceeds this value is treated as
#'   degenerate and excluded from selection. This guards against runaway
#'   discrimination fits, which can "converge" with a spuriously low model-based
#'   risk (huge discrimination collapses the item-parameter covariance). The
#'   default of 10 is far above any realistic 2PL discrimination.
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
#'   minimizes the trace of the item-parameter covariance matrix;
#'   [fit_mixed_subjects_mml()] for the marginal-likelihood estimator.
tune_lambda_ability_risk <- function(lambda_grid, observed, predicted, generated,
                                     target_resp = NULL, theta_true = NULL,
                                     n_quad = 31, initial_pars = NULL,
                                     fit_fn = fit_mixed_subjects,
                                     bounds = c(-6, 6),
                                     max_discrimination = 10,
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

    # Use the S3 generic so MML fits get vcov_mixed_subjects_mml (Louis bread).
    Sigma <- tryCatch(stats::vcov(fits[[i]]),
                      error = function(e) NULL)

    if (is.null(Sigma)) {
      risks[[i]] <- list(summary = data.frame(
        mean_param_var     = Inf,
        mean_squared_error = Inf,
        mean_total_risk    = Inf
      ))
    } else {
      risks[[i]] <- ability_risk(
        target_resp, fits[[i]], vcov = Sigma,
        theta_true = theta_true, bounds = bounds
      )
    }

    rows[[i]] <- data.frame(
      lambda             = lambda,
      mean_param_var     = risks[[i]]$summary$mean_param_var,
      mean_squared_error = risks[[i]]$summary$mean_squared_error,
      mean_total_risk    = risks[[i]]$summary$mean_total_risk,
      convergence        = fits[[i]]$convergence,
      max_disc           = max(abs(fits[[i]]$item_pars$a))
    )
  }

  summary <- do.call(rbind, rows)

  # Filter ineligible candidates. A candidate is excluded if it failed to
  # converge, produced non-finite risk, OR has a degenerate (runaway)
  # discrimination estimate. The last guard is essential: a fit whose
  # discrimination blows up reports huge item information, which collapses the
  # model-based covariance and makes the ability risk g'Sigma g spuriously
  # small. Such a fit can "converge" with the lowest apparent risk while being
  # numerically garbage (e.g. a = 1000), so it must be rejected explicitly.
  # Keep the full summary for diagnostics; use selection_risk for argmin.
  summary$selection_risk <- summary$mean_total_risk
  summary$selection_risk[
    !is.finite(summary$selection_risk) |
      summary$convergence != 0 |
      summary$max_disc > max_discrimination
  ] <- Inf

  if (all(is.infinite(summary$selection_risk))) {
    warning(
      "No lambda candidate converged with finite ability-score risk. ",
      "Returning lambda = 0 (human-only estimate).",
      call. = FALSE
    )
    best_idx <- which(summary$lambda == 0)
    best_idx <- if (length(best_idx) > 0L) best_idx[1L] else 1L
  } else {
    best_idx <- which.min(summary$selection_risk)
  }

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
#' labeled rows, then fits a final model with those fold-specific lambda values.
#'
#' @inheritParams tune_lambda_ability_risk
#' @param n_splits Number of sample splits.
#' @param split_id Optional integer split assignment for labeled rows.
#' @param seed Optional seed used when `split_id` is omitted.
#' @param target_mode How `target_resp` is handled in each fold.
#'   `"fixed"` (default): the full `target_resp` is used to evaluate risk in
#'   every fold, suitable when the target population is fixed and independent of
#'   the labeled-data split (e.g. an operational scoring population).
#'   `"row_aligned"`: only the training rows of `target_resp` are used, which
#'   is valid when `target_resp = observed` and fold-matched evaluation is
#'   desired.
#' @param final_fit_fn Function used to produce the final combined-data fit.
#'   Defaults to [fit_mixed_subjects_split()], which accepts a per-fold
#'   `lambda` vector natively.
#'   Pass [fit_mixed_subjects_mml()] to get a scalar marginal-MML final fit: the
#'   fold-specific lambdas are averaged (weighted by fold size) into a single
#'   scalar, avoiding the accidental per-item lambda problem that occurs when a
#'   length-`n_splits` vector is passed directly to `fit_mixed_subjects_mml()`.
#'   Note that mixing MML fold-tuning with a frozen final fit is an
#'   approximation; document this when reporting results.
#'
#' @param fit_fn Fitting function used for each fold's ability-risk tuning
#'   (passed to [tune_lambda_ability_risk()]). Defaults to
#'   [fit_mixed_subjects()] (frozen expected-count). Pass
#'   [fit_mixed_subjects_mml()] for marginal-MML fold tuning.
#' @param tuning_args Named list of extra arguments forwarded only to the
#'   fold-level [tune_lambda_ability_risk()] calls (and through them to
#'   `fit_fn`). For example, `tuning_args = list(slope_upper = 4)`.
#' @param final_args Named list of extra arguments forwarded only to
#'   `final_fit_fn`. For example, `final_args = list(mml_pred_weights = "own")`.
#'   This keeps tuning-specific and final-fit-specific arguments cleanly
#'   separated, avoiding the earlier `...` leakage between the two.
#' @param ... Deprecated; forwarded to `tuning_args` for backward compatibility
#'   with a one-time message. Prefer `tuning_args` / `final_args`.
#'
#' @return A list with fold-specific lambda values, fold tuning objects, and the
#'   final fit.
#' @export
tune_lambda_ability_risk_crossfit <- function(lambda_grid, observed, predicted,
                                              generated, target_resp = NULL,
                                              theta_true = NULL, n_splits = 2,
                                              split_id = NULL, seed = NULL,
                                              n_quad = 31, initial_pars = NULL,
                                              target_mode = c("fixed", "row_aligned"),
                                              fit_fn = fit_mixed_subjects,
                                              final_fit_fn = fit_mixed_subjects_split,
                                              tuning_args = list(),
                                              final_args = list(),
                                              bounds = c(-6, 6),
                                              control = list(maxit = 500), ...) {
  # Backward compatibility: route any legacy ... into tuning_args.
  dots <- list(...)
  if (length(dots) > 0L) {
    message(
      "tune_lambda_ability_risk_crossfit(): passing extra arguments via `...` ",
      "is deprecated. Use `tuning_args` (fold tuning) or `final_args` (final ",
      "fit) instead. Routing `...` into `tuning_args` for now."
    )
    tuning_args <- utils::modifyList(tuning_args, dots)
  }
  if (!is.list(tuning_args) || !is.list(final_args)) {
    stop("tuning_args and final_args must be lists.", call. = FALSE)
  }
  target_mode <- match.arg(target_mode)

  observed <- validate_response_matrix(observed, "observed",
                                       allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted",
                                        allow_fractional = FALSE)
  generated <- validate_response_matrix(generated, "generated",
                                        allow_fractional = FALSE)
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
    # When target_resp is derived from observed, default to row_aligned
    # unless the user explicitly requested fixed.
    if (target_mode == "fixed") target_mode <- "row_aligned"
  } else {
    target_resp <- validate_response_matrix(target_resp, "target_resp",
                                            allow_fractional = TRUE)
    check_same_items(observed, target_resp, "observed", "target_resp")
    if (target_mode == "row_aligned" &&
        nrow(target_resp) != nrow(observed)) {
      stop(
        "target_mode = 'row_aligned' requires nrow(target_resp) == nrow(observed).",
        call. = FALSE
      )
    }
  }

  split_values    <- sort(unique(split_id))
  lambda_by_split <- numeric(length(split_values))
  fold_tuning     <- vector("list", length(split_values))

  for (s in seq_along(split_values)) {
    fold  <- split_values[s]
    train <- split_id != fold

    # Target and theta for this fold
    if (target_mode == "row_aligned") {
      target_train <- target_resp[train, , drop = FALSE]
      theta_train  <- if (is.null(theta_true)) NULL else theta_true[train]
    } else {
      target_train <- target_resp      # fixed target population
      theta_train  <- theta_true       # fixed true thetas (or NULL)
    }

    # Fold tuning receives the core args plus the user's tuning_args only.
    fold_tuning[[s]] <- do.call(tune_lambda_ability_risk, c(
      list(
        lambda_grid  = lambda_grid,
        observed     = observed[train, , drop = FALSE],
        predicted    = predicted[train, , drop = FALSE],
        generated    = generated,
        target_resp  = target_train,
        theta_true   = theta_train,
        n_quad       = n_quad,
        initial_pars = initial_pars,
        fit_fn       = fit_fn,
        bounds       = bounds,
        control      = control
      ),
      tuning_args
    ))
    lambda_by_split[s] <- fold_tuning[[s]]$best_lambda
  }

  # Determine what lambda to pass to the final fit function.
  # fit_mixed_subjects_split() and related functions accept a per-fold vector.
  # fit_mixed_subjects_mml() accepts only a scalar or per-item vector, so
  # passing a length-n_splits vector risks accidental per-item interpretation.
  # When final_fit_fn is not fit_mixed_subjects_split, average the fold lambdas.
  uses_split_fn <- identical(final_fit_fn, fit_mixed_subjects_split)
  if (uses_split_fn) {
    lambda_final  <- lambda_by_split
    extra_split   <- list(split_id = split_id)
  } else {
    # Weighted mean: weight each fold's lambda by the number of labeled rows in
    # that fold's held-out set.
    fold_sizes <- vapply(
      split_values,
      function(v) sum(split_id == v),
      integer(1)
    )
    lambda_final <- sum(fold_sizes * lambda_by_split) / sum(fold_sizes)
    extra_split  <- list()
  }

  # final_fit_fn receives the core args, split metadata, plus the user's
  # final_args only.  Tuning-specific args (fit_fn, tuning_args) never reach it.
  final_fit <- do.call(final_fit_fn, c(
    list(
      observed     = observed,
      predicted    = predicted,
      generated    = generated,
      lambda       = lambda_final,
      n_quad       = n_quad,
      initial_pars = initial_pars,
      control      = control
    ),
    extra_split,
    final_args
  ))

  list(
    lambda_by_split = lambda_by_split,
    lambda_final    = lambda_final,   # scalar mean or fold vector, as passed to final_fit_fn
    split_id        = split_id,
    fold_tuning     = fold_tuning,
    final_fit       = final_fit
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
#' the λ that minimizes the marginal variance of `(a_j, d_j)` independently for
#' each item.
#'
#' **Use case.** When a single global λ is forced to zero because a few items
#' have poor LLM predictions, per-item λ_j allows well-predicted items to still
#' benefit from the LLM data. Pass the returned vector to
#' [fit_mixed_subjects_mml()] as the `lambda` argument.
#'
#' This is a **theoretical diagnostic**: it minimizes item-parameter variance,
#' not ability-score risk. For operational scoring use
#' [tune_lambda_ability_risk_item()] instead.
#'
#' @param observed Human response matrix.
#' @param predicted Paired binary LLM responses (0/1) for the same rows as
#'   `observed`. Probabilities are not accepted; sample binary responses first.
#' @param item_pars Item parameters at which to evaluate the score vectors.
#' @param n_generated Number of generated (unpaired) LLM subjects.
#' @param quadrature Optional quadrature grid.
#' @param n_quad Number of quadrature nodes when `quadrature` is omitted.
#'
#' @return A list with `lambda` (numeric vector of length `n_items`), `item`
#'   (item names), `n`, `n_generated`, and `r` (the ratio `n / n_generated`).
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
                                        allow_fractional = FALSE)
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
#' Finds a per-item vector of lambda values `λ_j ∈ [0, 1]` that minimizes
#' propagated ability-score risk `E[g' Σ_γ g]` using coordinate descent on the
#' items. Each coordinate step selects the `λ_j` in `lambda_grid` that gives
#' the smallest mean ability risk while holding all other `λ_{j'}` fixed.
#'
#' Calls [fit_mixed_subjects_mml()] with a per-item lambda vector at each
#' candidate evaluation. Because the lambda is a vector, that function
#' **switches to its frozen expected-count Q-function path** — posteriors are
#' frozen at `initial_pars`, not recomputed continuously. This is an
#' approximation; see the `@note` below. The resulting lambda vector can be
#' used directly with [fit_mixed_subjects_mml()].
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
#'   `init_lambda = 0.5`) to start the search around a useful operating point.
#'   Starting from all-zeros is not recommended: each single-item improvement is
#'   too small to detect when other items are at zero. A scalar is broadcast to
#'   all items; a vector of length `n_items` sets per-item starting values.
#' @param bounds Bounds passed to [score_theta()].
#' @param max_discrimination Upper bound on plausible item discrimination; any
#'   candidate fit whose maximum `|a|` exceeds it is treated as degenerate and
#'   skipped. See [tune_lambda_ability_risk()] for the rationale. Default 10.
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to [fit_mixed_subjects_mml()].
#'
#' @note **Approximation status.** The coordinate descent fits use the frozen
#'   expected-count Q-function (not the full marginal-MML objective) because the
#'   IRT marginal likelihood integrates over the joint response pattern and does
#'   not decompose item-wise. The approach is approximately correct when
#'   `initial_pars` is close to the converged parameters. Report per-item
#'   results as experimental / approximate.
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
                                          max_discrimination = 10,
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

        # Reject degenerate (runaway discrimination) fits before trusting their
        # model-based risk; see tune_lambda_ability_risk() for rationale.
        if (max(abs(fit_j$item_pars$a)) > max_discrimination) next

        risk_j <- tryCatch({
          Sigma_j <- stats::vcov(fit_j)
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
