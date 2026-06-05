# ---------------------------------------------------------------------------- #
# 1PL (one-parameter logistic) model fitting
#
# The 1PL estimates a single shared discrimination a (equal across all items)
# plus per-item intercepts d_j.  The response probability is
#   P(x_j = 1 | theta) = plogis(a * theta + d_j)
# Parameter vector: c(a_shared, d_1, ..., d_J) — length J+1, not 2J.
#
# All quadrature, posterior, and expected-count machinery from the 2PL is
# reused unchanged.  Only the gradient and parameter-vector bookkeeping differ.
# ---------------------------------------------------------------------------- #


# ---------- Internal helpers ------------------------------------------------

item_pars_from_vector_1pl <- function(par, item_names) {
  n_items <- length(item_names)
  if (length(par) != n_items + 1L) {
    stop("1PL par must have length n_items + 1 (one shared a, then d_1 ... d_J).",
         call. = FALSE)
  }
  out <- data.frame(
    item = item_names,
    a    = par[1],
    d    = par[2:(n_items + 1L)],
    stringsAsFactors = FALSE
  )
  out$b <- -out$d / out$a
  out
}

vector_from_item_pars_1pl <- function(item_pars) {
  item_pars <- standardize_item_pars(item_pars)
  c(item_pars$a[1], item_pars$d)
}

gradient_expected_counts_1pl <- function(counts, item_pars) {
  # Returns a (J+1)-vector: c(grad_a_shared, grad_d_1, ..., grad_d_J).
  # The shared-discrimination gradient is the SUM of per-item a contributions.
  item_pars <- standardize_item_pars(
    item_pars, n_items = counts$n_items, item_names = counts$item_names
  )

  eta <- outer(item_pars$a, counts$theta, `*`) +
    matrix(item_pars$d, nrow = counts$n_items, ncol = counts$n_nodes)
  resid <- counts$N * stats::plogis(eta) - counts$R

  theta_mat <- matrix(counts$theta, nrow = counts$n_items, ncol = counts$n_nodes,
                      byrow = TRUE)

  # Shared a: sum over ALL items × nodes
  grad_a_shared <- sum(resid * theta_mat) / counts$n

  # Per-item d
  grad_d <- rowSums(resid) / counts$n

  c(grad_a_shared, grad_d)
}

marginal_gradient_1pl <- function(resp, item_pars, quadrature) {
  # Gradient of the marginal IRT loss for the 1PL model at item_pars.
  # Identical to gradient_expected_counts_1pl evaluated at current posteriors.
  weights <- posterior_weights_2pl(resp, item_pars, quadrature = quadrature)
  counts  <- summarize_expected_counts(resp, weights)
  gradient_expected_counts_1pl(counts, item_pars)
}

fit_from_counts_1pl <- function(counts_observed, counts_predicted, counts_generated,
                                 initial_pars, lambda,
                                 slope_lower = 1e-4, slope_upper = NULL,
                                 control = list(maxit = 500)) {
  check_counts_compatible(list(counts_observed, counts_predicted, counts_generated))
  lambda <- validate_lambda(lambda)

  item_names <- counts_observed$item_names
  n_items    <- length(item_names)
  initial_pars <- standardize_item_pars(initial_pars, n_items = n_items,
                                         item_names = item_names)

  objective <- function(par) {
    if (!is.finite(par[1]) || par[1] <= 0) return(.Machine$double.xmax)
    ip  <- item_pars_from_vector_1pl(par, item_names)
    val <- loss_expected_counts(counts_observed, ip) +
      lambda * (loss_expected_counts(counts_generated, ip) -
                loss_expected_counts(counts_predicted, ip))
    if (!is.finite(val)) .Machine$double.xmax else val
  }

  gradient <- function(par) {
    if (!is.finite(par[1]) || par[1] <= 0) return(rep(0, n_items + 1L))
    ip <- item_pars_from_vector_1pl(par, item_names)
    g  <- gradient_expected_counts_1pl(counts_observed, ip) +
      lambda * (gradient_expected_counts_1pl(counts_generated, ip) -
                gradient_expected_counts_1pl(counts_predicted, ip))
    ifelse(is.finite(g), g, 0)
  }

  start <- vector_from_item_pars_1pl(initial_pars)

  lower <- c(if (is.null(slope_lower)) -Inf else slope_lower, rep(-Inf, n_items))
  upper <- c(if (is.null(slope_upper)) Inf  else slope_upper, rep( Inf, n_items))
  if (!is.null(slope_lower)) start[1] <- max(start[1], slope_lower)
  if (!is.null(slope_upper)) start[1] <- min(start[1], slope_upper)

  ctrl <- utils::modifyList(list(maxit = 500), control)
  opt  <- stats::optim(par = start, fn = objective, gr = gradient,
                        method = "L-BFGS-B", lower = lower, upper = upper,
                        control = ctrl)

  list(
    item_pars   = item_pars_from_vector_1pl(opt$par, item_names),
    par         = opt$par,
    value       = opt$value,
    convergence = opt$convergence,
    message     = opt$message,
    optimizer   = opt
  )
}


# ---------- Exported fitting functions --------------------------------------

#' Fit a 1PL (one-parameter logistic) model
#'
#' Estimates a shared discrimination parameter `a` (equal across all items)
#' and per-item intercepts `d_j` by maximising the IRT marginal likelihood
#' under a standard-normal ability prior using L-BFGS-B.
#'
#' The response probability is `P(x_j = 1 | theta) = plogis(a * theta + d_j)`.
#' The parameter vector has length `J + 1`: one shared discrimination followed
#' by J per-item intercepts.
#'
#' @param resp Binary response matrix.
#' @param n_quad Number of standard-normal quadrature nodes.
#' @param initial_pars Optional starting item parameters (data frame with `a`
#'   and `d` columns). If omitted, `a = 1` and `d_j = qlogis(p_j)` where
#'   `p_j` is the observed proportion correct for item `j`.
#' @param quadrature Optional quadrature grid.
#' @param slope_lower,slope_upper Bounds on the shared discrimination.
#' @param control Control list passed to [stats::optim()].
#'
#' @return A list with `pars` (item parameter data frame with all `a` equal),
#'   `par` (the raw parameter vector), and optimizer details.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
#' resp <- simulate_2pl(rnorm(60), pars)
#' fit <- fit_1pl(resp, n_quad = 7)
#' fit$pars
fit_1pl <- function(resp, n_quad = 31, initial_pars = NULL, quadrature = NULL,
                     slope_lower = 1e-4, slope_upper = NULL,
                     control = list(maxit = 500)) {
  resp       <- validate_response_matrix(resp, name = "resp", allow_fractional = FALSE)
  quadrature <- check_quadrature(quadrature, n_quad = n_quad)
  item_names <- colnames(resp)
  n_items    <- length(item_names)

  if (is.null(initial_pars)) {
    p_cor <- colMeans(resp, na.rm = TRUE)
    d_init <- stats::qlogis(pmax(0.05, pmin(0.95, p_cor)))
    initial_pars <- data.frame(item = item_names, a = 1, d = d_init,
                                stringsAsFactors = FALSE)
    initial_pars$b <- -initial_pars$d / initial_pars$a
  } else {
    initial_pars <- standardize_item_pars(initial_pars, n_items = n_items,
                                           item_names = item_names)
  }

  .cache     <- new.env(parent = emptyenv())
  .cache$par <- NULL
  .cache$val <- NULL
  .cache$grd <- NULL

  recompute <- function(par) {
    if (!is.null(.cache$par) && isTRUE(all.equal(par, .cache$par, tolerance = 0))) {
      return(invisible(NULL))
    }
    if (!is.finite(par[1]) || par[1] <= 0) {
      .cache$val <- .Machine$double.xmax
      .cache$grd <- rep(0, n_items + 1L)
      .cache$par <- par
      return(invisible(NULL))
    }
    ip  <- item_pars_from_vector_1pl(par, item_names)
    val <- marginal_loss_2pl(resp, ip, quadrature)
    g   <- marginal_gradient_1pl(resp, ip, quadrature)
    .cache$par <- par
    .cache$val <- if (is.finite(val)) val else .Machine$double.xmax
    .cache$grd <- ifelse(is.finite(g), g, 0)
    invisible(NULL)
  }

  objective <- function(par) { recompute(par); .cache$val }
  gradient  <- function(par) { recompute(par); .cache$grd }

  start <- vector_from_item_pars_1pl(initial_pars)
  lower <- c(if (is.null(slope_lower)) -Inf else slope_lower, rep(-Inf, n_items))
  upper <- c(if (is.null(slope_upper)) Inf  else slope_upper, rep( Inf, n_items))
  if (!is.null(slope_lower)) start[1] <- max(start[1], slope_lower)
  if (!is.null(slope_upper)) start[1] <- min(start[1], slope_upper)

  ctrl <- utils::modifyList(list(maxit = 500), control)
  opt  <- stats::optim(par = start, fn = objective, gr = gradient,
                        method = "L-BFGS-B", lower = lower, upper = upper,
                        control = ctrl)

  list(
    pars        = item_pars_from_vector_1pl(opt$par, item_names),
    par         = opt$par,
    value       = opt$value,
    convergence = opt$convergence,
    message     = opt$message,
    model       = NULL
  )
}


#' Fit a mixed-subjects 1PL calibration (frozen expected-count)
#'
#' Analogous to [fit_mixed_subjects()] but estimates a shared discrimination
#' parameter `a` across all items (1PL model). Posterior quadrature weights
#' are frozen at the initial parameter estimates.
#'
#' @inheritParams fit_mixed_subjects
#' @param ... Additional arguments passed to [fit_1pl()] when `initial_pars`
#'   is omitted.
#'
#' @return An object of class `c("mixedsubjects_1pl_fit", "mixedsubjects_fit")`.
#' @export
#'
#' @seealso [fit_mixed_subjects_mml_1pl()] for the marginal-likelihood version;
#'   [fit_mixed_subjects()] for the 2PL version.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
#' observed  <- simulate_2pl(rnorm(40), pars)
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects_1pl(
#'   observed, observed, generated,
#'   lambda = 0.5, initial_pars = pars, n_quad = 7,
#'   control = list(maxit = 50)
#' )
#' fit$item_pars
fit_mixed_subjects_1pl <- function(observed, predicted, generated, lambda = 1,
                                    n_quad = 31, initial_pars = NULL,
                                    quadrature = NULL,
                                    common_predicted_weights = TRUE,
                                    slope_lower = 1e-4, slope_upper = NULL,
                                    control = list(maxit = 500), ...) {
  lambda    <- validate_lambda(lambda)
  observed  <- validate_response_matrix(observed,  "observed",  allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted", allow_fractional = TRUE)
  generated <- validate_response_matrix(generated, "generated", allow_fractional = TRUE)
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")
  if (nrow(observed) != nrow(predicted)) {
    stop("observed and predicted must have the same number of rows.", call. = FALSE)
  }

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  if (is.null(initial_pars)) {
    initial_fit  <- fit_1pl(observed, n_quad = n_quad, ...)
    initial_pars <- initial_fit$pars
  } else {
    initial_fit  <- NULL
    initial_pars <- standardize_item_pars(initial_pars, n_items = ncol(observed),
                                           item_names = colnames(observed))
  }

  q_obs  <- build_quadrature_summary(observed, initial_pars, quadrature)
  w_pred <- if (isTRUE(common_predicted_weights)) q_obs$weights else NULL
  q_pred <- build_quadrature_summary(predicted, initial_pars, quadrature,
                                      weights = w_pred)
  q_gen  <- build_quadrature_summary(generated, initial_pars, quadrature)

  fit <- fit_from_counts_1pl(
    counts_observed  = q_obs$counts,
    counts_predicted = q_pred$counts,
    counts_generated = q_gen$counts,
    initial_pars     = initial_pars,
    lambda           = lambda,
    slope_lower      = slope_lower,
    slope_upper      = slope_upper,
    control          = control
  )

  out <- c(fit, list(
    lambda                   = lambda,
    initial_pars             = initial_pars,
    initial_model            = if (is.null(initial_fit)) NULL else initial_fit$model,
    quadrature               = quadrature,
    q_observed               = q_obs,
    q_predicted              = q_pred,
    q_generated              = q_gen,
    common_predicted_weights = common_predicted_weights,
    mml                      = FALSE,
    model_type               = "1pl",
    paired_missing           = "match_observed",
    split                    = NULL,
    call                     = match.call()
  ))
  class(out) <- c("mixedsubjects_1pl_fit", "mixedsubjects_fit")
  out
}


#' Fit a mixed-subjects 1PL calibration via marginal maximum likelihood
#'
#' Analogous to [fit_mixed_subjects_mml()] but estimates a shared discrimination
#' parameter `a` across all items (1PL model). Posteriors are recomputed at
#' every gradient evaluation — no frozen-posterior gradient asymmetry.
#'
#' Only scalar `lambda` is supported; per-item lambda is not meaningful for
#' the 1PL because the discrimination is shared across items.
#'
#' @inheritParams fit_mixed_subjects_mml
#' @param ... Additional arguments passed to [fit_1pl()] when `initial_pars`
#'   is omitted.
#'
#' @return An object of class `c("mixedsubjects_1pl_fit", "mixedsubjects_fit")`.
#' @export
#'
#' @seealso [fit_mixed_subjects_1pl()] for the frozen expected-count version;
#'   [fit_mixed_subjects_mml()] for the 2PL version.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
#' observed  <- simulate_2pl(rnorm(40), pars)
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects_mml_1pl(
#'   observed, observed, generated,
#'   lambda = 0.5, initial_pars = pars, n_quad = 7,
#'   control = list(maxit = 100)
#' )
#' fit$item_pars
fit_mixed_subjects_mml_1pl <- function(observed, predicted, generated, lambda = 1,
                                        n_quad = 31, initial_pars = NULL,
                                        quadrature = NULL,
                                        mml_pred_weights = c("own", "human"),
                                        slope_lower = 1e-4, slope_upper = NULL,
                                        control = list(maxit = 500), ...) {
  lambda           <- validate_lambda(lambda)
  mml_pred_weights <- match.arg(mml_pred_weights)

  observed  <- validate_response_matrix(observed,  "observed",  allow_fractional = FALSE)
  predicted <- validate_response_matrix(predicted, "predicted", allow_fractional = TRUE)
  generated <- validate_response_matrix(generated, "generated", allow_fractional = TRUE)
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")
  if (nrow(observed) != nrow(predicted)) {
    stop("observed and predicted must have the same number of rows.", call. = FALSE)
  }

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  if (is.null(initial_pars)) {
    initial_fit <- fit_1pl(observed, n_quad = n_quad, ...)
    init_std    <- initial_fit$pars
  } else {
    initial_fit <- NULL
    init_std    <- standardize_item_pars(initial_pars, n_items = ncol(observed),
                                          item_names = colnames(observed))
  }

  item_names <- colnames(observed)
  n_items    <- length(item_names)

  .cache     <- new.env(parent = emptyenv())
  .cache$par <- NULL
  .cache$val <- NULL
  .cache$grd <- NULL

  recompute <- function(par) {
    if (!is.null(.cache$par) &&
        isTRUE(all.equal(par, .cache$par, tolerance = 0))) {
      return(invisible(NULL))
    }
    if (!is.finite(par[1]) || par[1] <= 0) {
      .cache$val <- .Machine$double.xmax
      .cache$grd <- rep(0, n_items + 1L)
      .cache$par <- par
      return(invisible(NULL))
    }

    ip <- item_pars_from_vector_1pl(par, item_names)

    obs_r <- posterior_and_log_lik_2pl(observed,  ip, quadrature)
    gen_r <- posterior_and_log_lik_2pl(generated, ip, quadrature)

    obs_counts <- summarize_expected_counts(observed,  obs_r$weights)
    gen_counts <- summarize_expected_counts(generated, gen_r$weights)

    if (mml_pred_weights == "own") {
      pred_r      <- posterior_and_log_lik_2pl(predicted, ip, quadrature)
      pred_counts <- summarize_expected_counts(predicted, pred_r$weights)
      l_pred <- -mean(pred_r$log_normalizers)
      g_pred <- gradient_expected_counts_1pl(pred_counts, ip)
    } else {
      pred_counts <- summarize_expected_counts(predicted, obs_r$weights)
      l_pred <- loss_expected_counts(pred_counts, ip)
      g_pred <- gradient_expected_counts_1pl(pred_counts, ip)
    }

    val <- -mean(obs_r$log_normalizers) +
      lambda * (-mean(gen_r$log_normalizers) - l_pred)
    g   <- gradient_expected_counts_1pl(obs_counts, ip) +
      lambda * (gradient_expected_counts_1pl(gen_counts, ip) - g_pred)

    .cache$par <- par
    .cache$val <- if (is.finite(val)) val else .Machine$double.xmax
    .cache$grd <- ifelse(is.finite(g), g, 0)
    invisible(NULL)
  }

  objective <- function(par) { recompute(par); .cache$val }
  gradient  <- function(par) { recompute(par); .cache$grd }

  start <- vector_from_item_pars_1pl(init_std)
  lower <- c(if (is.null(slope_lower)) -Inf else slope_lower, rep(-Inf, n_items))
  upper <- c(if (is.null(slope_upper)) Inf  else slope_upper, rep( Inf, n_items))
  if (!is.null(slope_lower)) start[1] <- max(start[1], slope_lower)
  if (!is.null(slope_upper)) start[1] <- min(start[1], slope_upper)

  ctrl <- utils::modifyList(list(maxit = 500), control)
  opt  <- stats::optim(par = start, fn = objective, gr = gradient,
                        method = "L-BFGS-B", lower = lower, upper = upper,
                        control = ctrl)

  conv_pars   <- item_pars_from_vector_1pl(opt$par, item_names)
  q_obs_final <- build_quadrature_summary(observed, conv_pars, quadrature)
  q_pred_final <- if (mml_pred_weights == "own") {
    build_quadrature_summary(predicted, conv_pars, quadrature)
  } else {
    build_quadrature_summary(predicted, conv_pars, quadrature,
                              weights = q_obs_final$weights)
  }
  q_gen_final <- build_quadrature_summary(generated, conv_pars, quadrature)

  out <- list(
    item_pars                = conv_pars,
    par                      = opt$par,
    value                    = opt$value,
    convergence              = opt$convergence,
    message                  = opt$message,
    optimizer                = opt,
    lambda                   = lambda,
    initial_pars             = init_std,
    initial_model            = if (is.null(initial_fit)) NULL else initial_fit$model,
    quadrature               = quadrature,
    q_observed               = q_obs_final,
    q_predicted              = q_pred_final,
    q_generated              = q_gen_final,
    common_predicted_weights = (mml_pred_weights == "human"),
    mml_pred_weights         = mml_pred_weights,
    mml                      = TRUE,
    model_type               = "1pl",
    paired_missing           = NA,
    split                    = NULL,
    call                     = match.call()
  )
  class(out) <- c("mixedsubjects_1pl_fit", "mixedsubjects_fit")
  out
}


#' @export
print.mixedsubjects_1pl_fit <- function(x, ...) {
  cat("mixedsubjectsirt 1PL fit\n")
  cat("  items:      ", nrow(x$item_pars), "\n", sep = "")
  cat("  a (shared): ", signif(x$item_pars$a[1], 4), "\n", sep = "")
  cat("  lambda:     ", x$lambda, "\n", sep = "")
  cat("  loss:       ", signif(x$value, 6), "\n", sep = "")
  cat("  convergence:", x$convergence, "\n")
  if (isTRUE(x$mml)) cat("  estimator:  marginal MML PPI++ (1PL)\n")
  invisible(x)
}
