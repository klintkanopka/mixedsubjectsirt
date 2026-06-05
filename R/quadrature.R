#' Create a standard-normal Gauss-Hermite quadrature grid
#'
#' `rmutil::gauss.hermite()` returns nodes and weights for integrals of the form
#' `integral f(x) exp(-x^2) dx`. This function rescales those nodes and weights
#' to approximate expectations under a standard normal latent trait
#' distribution.
#'
#' @param n_quad Number of quadrature nodes.
#' @param iterlim Maximum number of Newton-Raphson iterations passed to
#'   [rmutil::gauss.hermite()].
#'
#' @return A data frame with node index, `theta`, `weight`, and backward
#'   compatible aliases `X_k` and `A_k`.
#' @export
#'
#' @examples
#' quad <- make_quadrature(7)
#' sum(quad$weight)
make_quadrature <- function(n_quad = 31, iterlim = 1e5) {
  if (!is.numeric(n_quad) || length(n_quad) != 1 || n_quad < 2) {
    stop("n_quad must be a single integer greater than 1.", call. = FALSE)
  }

  gh <- as.data.frame(rmutil::gauss.hermite(as.integer(n_quad), iterlim = iterlim))
  theta <- sqrt(2) * gh$Points
  weight <- gh$Weights / sqrt(pi)
  o <- order(theta)

  data.frame(
    node = seq_along(theta),
    theta = theta[o],
    weight = weight[o] / sum(weight[o]),
    X_k = theta[o],
    A_k = weight[o] / sum(weight[o])
  )
}

#' Compute posterior quadrature weights for a 2PL model
#'
#' Computes each subject's posterior distribution over a fixed quadrature grid
#' under a 2PL model, using stable log-likelihood calculations. Fractional
#' responses in `[0, 1]` are allowed, which is useful when LLM output is stored
#' as probabilities rather than sampled binary responses.
#'
#' @param resp A response matrix with rows for subjects and columns for items.
#'   Values may be binary, fractional in `[0, 1]`, or `NA`.
#' @param item_pars Item parameters in slope-intercept form. Supply a data frame
#'   or matrix with columns `a`/`a1` and `d`, or a fitted `mirt` model.
#' @param quadrature Optional quadrature data frame with `theta` and `weight`
#'   columns. If omitted, a standard-normal grid is created.
#' @param n_quad Number of quadrature nodes used when `quadrature` is omitted.
#' @param iterlim Maximum number of Newton-Raphson iterations passed to
#'   [rmutil::gauss.hermite()] when `quadrature` is omitted.
#'
#' @return A matrix with one row per subject and one column per quadrature node.
#'   Rows sum to one. Attributes `theta` and `weight` contain the grid.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' W <- posterior_weights_2pl(resp, pars, n_quad = 5)
#' rowSums(W)
posterior_weights_2pl <- function(resp, item_pars, quadrature = NULL,
                                  n_quad = 31, iterlim = 1e5) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  quadrature <- check_quadrature(quadrature, n_quad = n_quad, iterlim = iterlim)
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = ncol(resp),
    item_names = colnames(resp)
  )

  theta <- quadrature$theta
  prior <- quadrature$weight / sum(quadrature$weight)
  n <- nrow(resp)
  k <- length(theta)
  weights <- matrix(NA_real_, nrow = n, ncol = k)

  eta <- outer(theta, item_pars$a, `*`) +
    matrix(item_pars$d, nrow = k, ncol = ncol(resp), byrow = TRUE)
  log_p <- -safe_log1pexp(-eta)
  log_q <- -safe_log1pexp(eta)
  log_prior <- log(pmax(prior, .Machine$double.eps))

  for (i in seq_len(n)) {
    observed <- which(!is.na(resp[i, ]))
    ll <- log_prior

    if (length(observed) > 0) {
      y <- resp[i, observed]
      ll <- ll +
        as.vector((log_p[, observed, drop = FALSE] -
                     log_q[, observed, drop = FALSE]) %*% y) +
        rowSums(log_q[, observed, drop = FALSE])
    }

    weights[i, ] <- exp(ll - logsumexp(ll))
  }

  colnames(weights) <- paste0("node", seq_len(k))
  attr(weights, "theta") <- theta
  attr(weights, "weight") <- prior
  weights
}

#' Summarize response data as expected quadrature counts
#'
#' Converts response data and posterior quadrature weights into Bock-Aitkin style
#' expected counts. For each item and quadrature node, `N` is the expected number
#' of observed responses and `R` is the expected number correct.
#'
#' @param resp A response matrix with rows for subjects and columns for items.
#' @param weights Posterior quadrature weights, usually returned by
#'   [posterior_weights_2pl()].
#'
#' @return A list of class `"mixedsubjects_counts"` containing matrices `N` and
#'   `R`, sample size `n`, quadrature nodes, quadrature weights, and item names.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' W <- posterior_weights_2pl(resp, pars, n_quad = 5)
#' counts <- summarize_expected_counts(resp, W)
#' counts$N
summarize_expected_counts <- function(resp, weights) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  weights <- as.matrix(weights)

  if (nrow(resp) != nrow(weights)) {
    stop("resp and weights must have the same number of rows.", call. = FALSE)
  }
  if (any(!is.finite(weights)) || any(weights < -1e-12)) {
    stop("weights must be finite and non-negative.", call. = FALSE)
  }

  row_total <- rowSums(weights)
  if (any(abs(row_total - 1) > 1e-6)) {
    weights <- weights / row_total
  }

  n_items <- ncol(resp)
  n_nodes <- ncol(weights)
  N <- matrix(0, nrow = n_items, ncol = n_nodes)
  R <- matrix(0, nrow = n_items, ncol = n_nodes)

  for (j in seq_len(n_items)) {
    observed <- !is.na(resp[, j])
    if (any(observed)) {
      W_obs <- weights[observed, , drop = FALSE]
      y <- resp[observed, j]
      N[j, ] <- colSums(W_obs)
      R[j, ] <- colSums(W_obs * matrix(y, nrow = length(y), ncol = n_nodes))
    }
  }

  theta <- attr(weights, "theta")
  prior <- attr(weights, "weight")
  if (is.null(theta)) {
    theta <- seq_len(n_nodes)
  }
  if (is.null(prior)) {
    prior <- rep(1 / n_nodes, n_nodes)
  }

  colnames(N) <- colnames(R) <- paste0("node", seq_len(n_nodes))
  rownames(N) <- rownames(R) <- colnames(resp)

  out <- list(
    N = N,
    R = R,
    n = nrow(resp),
    n_items = n_items,
    n_nodes = n_nodes,
    item_names = colnames(resp),
    theta = as.numeric(theta),
    weight = as.numeric(prior) / sum(as.numeric(prior)),
    node_count = colSums(weights)
  )
  class(out) <- "mixedsubjects_counts"
  out
}

#' Convert responses to quadrature form
#'
#' Fits or accepts a 2PL model, computes posterior quadrature weights for each
#' subject, and returns expected counts for mixed-subjects calibration. This is a
#' lower-level helper; most analyses should call [fit_mixed_subjects()] or
#' [fit_mixed_subjects_split()].
#'
#' @param resp A response matrix with rows for subjects and columns for items.
#' @param N_quad Number of quadrature nodes to compute. Kept for backward
#'   compatibility; prefer `n_quad` in new code.
#' @param eps Retained for backward compatibility. Stable log computations are
#'   used instead of probability clipping.
#' @param iterlim Maximum number of Newton-Raphson iterations passed to
#'   [rmutil::gauss.hermite()].
#' @param irt_pars Optional target item parameters for mean-mean linking. This
#'   argument is kept for backward compatibility with earlier package versions.
#' @param item_pars Optional item parameters. If omitted, a 2PL model is fit to
#'   `resp` using [fit_2pl()].
#' @param quadrature Optional quadrature grid with `theta` and `weight` columns.
#' @param link_method Linking method used when `irt_pars` is supplied.
#' @param ... Additional arguments passed to [fit_2pl()] when `item_pars` is
#'   omitted.
#'
#' @return A list with `quad`, `counts`, `weights`, `irt_pars`, `quadrature`, and
#'   `theta`.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
#' names(q)
posterior_and_log_lik_2pl <- function(resp, item_pars, quadrature) {
  # Computes posterior weights AND per-subject marginal log-likelihoods.
  # Identical to posterior_weights_2pl but additionally saves the log-normaliser
  # log Z_i = log[Σ_k A_k p(Y_i | θ_k; γ)] for each subject.  These
  # log-normalisers form the true IRT marginal log-likelihood and are needed
  # by fit_mixed_subjects_mml() to evaluate the exact PPI++ objective.
  theta     <- quadrature$theta
  prior     <- quadrature$weight / sum(quadrature$weight)
  n         <- nrow(resp)
  k         <- length(theta)
  weights   <- matrix(NA_real_, nrow = n, ncol = k)
  log_norms <- numeric(n)

  eta   <- outer(theta, item_pars$a, `*`) +
    matrix(item_pars$d, nrow = k, ncol = ncol(resp), byrow = TRUE)
  log_p <- -safe_log1pexp(-eta)
  log_q <- -safe_log1pexp(eta)
  log_prior <- log(pmax(prior, .Machine$double.eps))

  for (i in seq_len(n)) {
    obs_j <- which(!is.na(resp[i, ]))
    ll <- log_prior
    if (length(obs_j) > 0) {
      y  <- resp[i, obs_j]
      ll <- ll +
        as.vector((log_p[, obs_j, drop = FALSE] -
                   log_q[, obs_j, drop = FALSE]) %*% y) +
        rowSums(log_q[, obs_j, drop = FALSE])
    }
    log_z      <- logsumexp(ll)
    log_norms[i] <- log_z
    weights[i, ] <- exp(ll - log_z)
  }

  colnames(weights) <- paste0("node", seq_len(k))
  attr(weights, "theta")  <- theta
  attr(weights, "weight") <- prior
  list(weights = weights, log_normalizers = log_norms)
}

mixed_subjects_quadrature <- function(resp, N_quad = 31, eps = 1e-15,
                                      iterlim = 1e5, irt_pars = NULL,
                                      item_pars = NULL, quadrature = NULL,
                                      link_method = "mean_mean", ...) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  quadrature <- check_quadrature(quadrature, n_quad = N_quad, iterlim = iterlim)

  if (is.null(item_pars)) {
    fitted <- fit_2pl(resp, ...)
    item_pars <- fitted$pars
  } else {
    item_pars <- standardize_item_pars(item_pars, n_items = ncol(resp),
                                       item_names = colnames(resp))
  }

  if (!is.null(irt_pars)) {
    item_pars <- link_item_parameters(
      source = item_pars,
      target = irt_pars,
      method = link_method
    )$pars
  }

  build_quadrature_summary(resp, item_pars, quadrature)
}
