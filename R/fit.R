#' Mixed-subjects objective function
#'
#' Evaluates the rectified mixed-subjects loss for 2PL item parameters. The
#' parameter vector must contain all discriminations first, followed by all
#' intercepts. The response probability is `plogis(d + a * theta)`.
#'
#' The objective is
#' `L_observed(pars) + lambda * (L_generated(pars) - L_predicted(pars))`.
#' Setting `lambda = 0` gives the human-only expected-count objective.
#'
#' @param pars Numeric vector of item parameters: all discriminations `a`
#'   followed by all intercepts `d`.
#' @param q_observed Quadrature summary for observed human responses, usually
#'   returned by [mixed_subjects_quadrature()].
#' @param q_predicted Quadrature summary for LLM responses/predictions on the
#'   same labeled human subjects.
#' @param q_llm Quadrature summary for generated or unlabeled LLM responses.
#' @param lambda Power-tuning parameter in `[0, 1]`.
#'
#' @return A scalar loss.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
#' mixed_subjects_loss(c(pars$a, pars$d), q, q, q, lambda = 0.5)
mixed_subjects_loss <- function(pars, q_observed, q_predicted, q_llm,
                                lambda = 0) {
  lambda <- validate_lambda(lambda)
  counts_observed <- q_to_counts(q_observed, "q_observed")
  counts_predicted <- q_to_counts(q_predicted, "q_predicted")
  counts_llm <- q_to_counts(q_llm, "q_llm")
  check_counts_compatible(list(counts_observed, counts_predicted, counts_llm))

  item_pars <- item_pars_from_vector(pars, counts_observed$item_names)

  loss_expected_counts(counts_observed, item_pars) +
    lambda * (
      loss_expected_counts(counts_llm, item_pars) -
        loss_expected_counts(counts_predicted, item_pars)
    )
}

#' Fit a mixed-subjects 2PL calibration
#'
#' Fits item parameters using observed human responses, paired LLM
#' responses/predictions for those same subjects, and generated or unlabeled LLM
#' responses. This implements the expected-count objective
#'
#' `L_human + lambda * (L_generated - L_paired_llm)`.
#'
#' By default the paired LLM responses reuse the posterior quadrature weights
#' from the observed human responses. This keeps the paired human and LLM terms
#' on the same latent covariate distribution, which is the closest analog to
#' prediction-powered inference with paired labels.
#'
#' @param observed Human response matrix, with rows for subjects and columns for
#'   items. Values must be binary when `initial_pars` is omitted.
#' @param predicted LLM responses or probabilities for the same rows and items as
#'   `observed`.
#' @param generated Generated or unlabeled LLM responses or probabilities for
#'   the same item columns.
#' @param lambda Power-tuning parameter in `[0, 1]`.
#' @param n_quad Number of standard-normal quadrature nodes.
#' @param initial_pars Optional starting item parameters. If omitted, a 2PL model
#'   is fit to `observed`.
#' @param quadrature Optional quadrature grid with `theta` and `weight` columns.
#' @param common_predicted_weights Logical; if `TRUE`, reuse the observed human
#'   posterior weights for `predicted`.
#' @param paired_missing How to handle missingness when
#'   `common_predicted_weights = TRUE`. The default, `"match_observed"`, requires
#'   `observed` and `predicted` to have the same missingness pattern so the paired
#'   LLM correction is evaluated only where a human label is present. Use
#'   `"allow"` only for explicit sensitivity analyses.
#' @param slope_lower Lower bound for discrimination parameters during
#'   optimization. Use `NULL` for no lower bound.
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to [fit_2pl()] when `initial_pars` is
#'   omitted.
#'
#' @return An object of class `"mixedsubjects_fit"` with fitted `item_pars`,
#'   optimizer details, quadrature summaries, and input settings.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' predicted <- observed
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects(
#'   observed, predicted, generated,
#'   lambda = 0.5, initial_pars = pars, n_quad = 7,
#'   control = list(maxit = 50)
#' )
#' fit$item_pars
fit_mixed_subjects <- function(observed, predicted, generated, lambda = 1,
                               n_quad = 31, initial_pars = NULL,
                               quadrature = NULL,
                               common_predicted_weights = TRUE,
                               paired_missing = c("match_observed", "allow"),
                               slope_lower = 1e-4,
                               slope_upper = NULL,
                               control = list(maxit = 500), ...) {
  lambda <- validate_lambda(lambda)
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
  paired_missing <- check_paired_missingness(
    observed = observed,
    predicted = predicted,
    paired_missing = paired_missing,
    common_predicted_weights = common_predicted_weights
  )

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)
  if (is.null(initial_pars)) {
    initial_fit <- fit_2pl(observed, ...)
    initial_pars <- initial_fit$pars
  } else {
    initial_fit <- NULL
    initial_pars <- standardize_item_pars(
      initial_pars,
      n_items = ncol(observed),
      item_names = colnames(observed)
    )
  }

  q_observed <- build_quadrature_summary(observed, initial_pars, quadrature)

  predicted_weights <- if (isTRUE(common_predicted_weights)) {
    q_observed$weights
  } else {
    NULL
  }
  q_predicted <- build_quadrature_summary(
    predicted,
    initial_pars,
    quadrature,
    weights = predicted_weights
  )
  q_generated <- build_quadrature_summary(generated, initial_pars, quadrature)

  fit <- fit_from_counts(
    counts_observed = q_observed$counts,
    counts_predicted = q_predicted$counts,
    counts_generated = q_generated$counts,
    initial_pars = initial_pars,
    lambda = lambda,
    slope_lower = slope_lower,
    slope_upper = slope_upper,
    control = control
  )

  out <- c(
    fit,
    list(
      lambda = lambda,
      initial_pars = initial_pars,
      initial_model = initial_fit$model,
      quadrature = quadrature,
      q_observed = q_observed,
      q_predicted = q_predicted,
      q_generated = q_generated,
      common_predicted_weights = common_predicted_weights,
      paired_missing = paired_missing,
      split = NULL,
      call = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' Fit from precomputed quadrature summaries
#'
#' Fits the mixed-subjects 2PL objective from quadrature/count summaries rather
#' than raw response matrices. This lower-level interface is useful when the
#' human, paired LLM, and generated LLM summaries have already been linked onto a
#' common scale outside the package.
#'
#' @param q_observed Quadrature summary for observed human responses. Usually
#'   returned by [mixed_subjects_quadrature()], but a raw counts object returned
#'   by [summarize_expected_counts()] is also accepted.
#' @param q_predicted Quadrature summary for paired LLM responses/predictions on
#'   the labeled human rows.
#' @param q_generated Quadrature summary for generated or unlabeled LLM
#'   responses.
#' @param lambda Power-tuning parameter in `[0, 1]`.
#' @param initial_pars Starting item parameters in slope-intercept form. If
#'   omitted, `q_observed$irt_pars` is used when available.
#' @param slope_lower Lower bound for discrimination parameters during
#'   optimization. Use `NULL` for no lower bound.
#' @param control Control list passed to [stats::optim()].
#'
#' @return An object of class `"mixedsubjects_fit"`.
#' @export
#'
#' @examples
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
#' q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
#' fit_mixed_subjects_from_quadrature(q, q, q, lambda = 0.5)$item_pars
fit_mixed_subjects_from_quadrature <- function(q_observed, q_predicted,
                                               q_generated, lambda = 1,
                                               initial_pars = NULL,
                                               slope_lower = 1e-4,
                                               slope_upper = NULL,
                                               control = list(maxit = 500)) {
  counts_observed <- q_to_counts(q_observed, "q_observed")
  counts_predicted <- q_to_counts(q_predicted, "q_predicted")
  counts_generated <- q_to_counts(q_generated, "q_generated")
  check_counts_compatible(list(counts_observed, counts_predicted, counts_generated))

  if (is.null(initial_pars)) {
    if (is.list(q_observed) && !is.null(q_observed$irt_pars)) {
      initial_pars <- q_observed$irt_pars
    } else {
      stop("initial_pars is required when q_observed does not contain irt_pars.",
           call. = FALSE)
    }
  }

  fit <- fit_from_counts(
    counts_observed = counts_observed,
    counts_predicted = counts_predicted,
    counts_generated = counts_generated,
    initial_pars = initial_pars,
    lambda = lambda,
    slope_lower = slope_lower,
    slope_upper = slope_upper,
    control = control
  )

  out <- c(
    fit,
    list(
      lambda = lambda,
      initial_pars = standardize_item_pars(
        initial_pars,
        n_items = counts_observed$n_items,
        item_names = counts_observed$item_names
      ),
      initial_model = NULL,
      quadrature = if (is.list(q_observed)) q_observed$quadrature else NULL,
      q_observed = q_observed,
      q_predicted = q_predicted,
      q_generated = q_generated,
      common_predicted_weights = NA,
      paired_missing = NA,
      split = NULL,
      call = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' Fit a split-sample mixed-subjects 2PL calibration
#'
#' Fits the same objective as [fit_mixed_subjects()], but constructs labeled
#' expected counts with cross-fitted posterior weights. For each split, the
#' initial human 2PL model is fit on the other splits and then used to compute
#' posterior weights for the held-out split. Each human row contributes to the
#' final estimating equation exactly once.
#'
#' Generated LLM counts are computed once per fold and averaged across folds so
#' that the generated sample keeps its original sample-size scale.
#'
#' @param observed Human response matrix, with rows for subjects and columns for
#'   items. Values must be binary when `initial_pars` is omitted.
#' @param predicted LLM responses or probabilities for the same rows and items as
#'   `observed`.
#' @param generated Generated or unlabeled LLM responses or probabilities for
#'   the same item columns.
#' @param lambda Power-tuning parameter in `[0, 1]`. Supply a scalar for a fixed
#'   lambda or a vector with one value per split for a precomputed
#'   cross-fitted-lambda analysis. When a vector is supplied, the generated term
#'   uses the split-size weighted mean lambda.
#' @param n_splits Number of sample splits.
#' @param split_id Optional integer vector assigning each observed row to a
#'   split. If omitted, splits are sampled at random.
#' @param seed Optional random seed used when `split_id` is omitted.
#' @param n_quad Number of standard-normal quadrature nodes.
#' @param initial_pars Optional item parameters to use in every fold instead of
#'   fitting fold-specific human models. This is mainly useful for testing or
#'   sensitivity analyses.
#' @param quadrature Optional quadrature grid with `theta` and `weight` columns.
#' @param common_predicted_weights Logical; if `TRUE`, reuse each held-out
#'   observed posterior weight matrix for its paired LLM responses.
#' @param paired_missing How to handle missingness when
#'   `common_predicted_weights = TRUE`. The default, `"match_observed"`, requires
#'   `observed` and `predicted` to have the same missingness pattern.
#' @param slope_lower Lower bound for discrimination parameters during
#'   optimization. Use `NULL` for no lower bound.
#' @param control Control list passed to [stats::optim()].
#' @param ... Additional arguments passed to [fit_2pl()] when fold-specific
#'   initial models are fit.
#'
#' @return An object of class `"mixedsubjects_fit"` with `split` metadata and
#'   fold-level initial parameters.
#' @export
#'
#' @examples
#' set.seed(2)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' predicted <- observed
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects_split(
#'   observed, predicted, generated,
#'   lambda = 0.5, initial_pars = pars, n_splits = 2,
#'   n_quad = 7, control = list(maxit = 50)
#' )
#' fit$item_pars
fit_mixed_subjects_split <- function(observed, predicted, generated,
                                     lambda = 1, n_splits = 2,
                                     split_id = NULL, seed = NULL,
                                     n_quad = 31, initial_pars = NULL,
                                     quadrature = NULL,
                                     common_predicted_weights = TRUE,
                                     paired_missing = c("match_observed", "allow"),
                                     slope_lower = 1e-4,
                                     slope_upper = NULL,
                                     control = list(maxit = 500), ...) {
  lambda <- validate_lambda_vector(lambda)
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
  paired_missing <- check_paired_missingness(
    observed = observed,
    predicted = predicted,
    paired_missing = paired_missing,
    common_predicted_weights = common_predicted_weights
  )

  if (is.null(split_id)) {
    split_id <- make_split_id(nrow(observed), n_splits, seed = seed)
  } else {
    split_id <- as.integer(split_id)
    if (length(split_id) != nrow(observed)) {
      stop("split_id must have length nrow(observed).", call. = FALSE)
    }
    n_splits <- length(unique(split_id))
  }

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)
  split_values <- sort(unique(split_id))
  lambda <- validate_lambda_vector(lambda, n = length(split_values))
  lambda_by_split <- if (length(lambda) == 1) {
    rep(lambda, length(split_values))
  } else {
    lambda
  }

  observed_counts <- vector("list", length(split_values))
  predicted_counts <- vector("list", length(split_values))
  generated_counts <- vector("list", length(split_values))
  fold_initial_pars <- vector("list", length(split_values))
  fold_models <- vector("list", length(split_values))

  for (s in seq_along(split_values)) {
    fold <- split_values[s]
    held_out <- split_id == fold
    train <- !held_out

    if (is.null(initial_pars)) {
      fold_fit <- fit_2pl(observed[train, , drop = FALSE], ...)
      fold_pars <- fold_fit$pars
      fold_models[[s]] <- fold_fit$model
    } else {
      fold_pars <- standardize_item_pars(
        initial_pars,
        n_items = ncol(observed),
        item_names = colnames(observed)
      )
      fold_models[[s]] <- NULL
    }

    fold_initial_pars[[s]] <- fold_pars

    q_observed <- build_quadrature_summary(
      observed[held_out, , drop = FALSE],
      fold_pars,
      quadrature
    )

    predicted_weights <- if (isTRUE(common_predicted_weights)) {
      q_observed$weights
    } else {
      NULL
    }
    q_predicted <- build_quadrature_summary(
      predicted[held_out, , drop = FALSE],
      fold_pars,
      quadrature,
      weights = predicted_weights
    )
    q_generated <- build_quadrature_summary(generated, fold_pars, quadrature)

    observed_counts[[s]] <- q_observed$counts
    predicted_counts[[s]] <- q_predicted$counts
    generated_counts[[s]] <- q_generated$counts
  }

  counts_observed <- combine_counts(observed_counts, mode = "sum")
  counts_predicted <- combine_counts(predicted_counts, mode = "sum")
  counts_generated <- combine_counts(generated_counts, mode = "mean")
  initial_average <- average_item_pars(fold_initial_pars)

  component_weights <- vapply(predicted_counts, `[[`, numeric(1), "n")
  component_weights <- component_weights / sum(component_weights)

  fit <- fit_from_count_components(
    counts_observed = counts_observed,
    counts_predicted = predicted_counts,
    counts_generated = counts_generated,
    initial_pars = initial_average,
    lambda = lambda_by_split,
    component_weights = component_weights,
    slope_lower = slope_lower,
    slope_upper = slope_upper,
    control = control
  )
  lambda_generated <- fit$lambda_generated
  fit$lambda <- NULL
  fit$lambda_generated <- NULL

  q_observed <- list(
    quad = counts_to_quad(counts_observed),
    counts = counts_observed,
    weights = NULL,
    irt_pars = initial_average,
    quadrature = quadrature,
    theta = quadrature$theta
  )
  q_predicted <- list(
    quad = counts_to_quad(counts_predicted),
    counts = counts_predicted,
    weights = NULL,
    irt_pars = initial_average,
    quadrature = quadrature,
    theta = quadrature$theta
  )
  q_generated <- list(
    quad = counts_to_quad(counts_generated),
    counts = counts_generated,
    weights = NULL,
    irt_pars = initial_average,
    quadrature = quadrature,
    theta = quadrature$theta
  )

  out <- c(
    fit,
    list(
      lambda = lambda,
      lambda_by_split = lambda_by_split,
      lambda_generated = lambda_generated,
      initial_pars = initial_average,
      initial_model = NULL,
      quadrature = quadrature,
      q_observed = q_observed,
      q_predicted = q_predicted,
      q_generated = q_generated,
      common_predicted_weights = common_predicted_weights,
      paired_missing = paired_missing,
      split = list(
        n_splits = n_splits,
        split_id = split_id,
        lambda_by_split = lambda_by_split,
        fold_initial_pars = fold_initial_pars,
        fold_models = fold_models
      ),
      call = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' Diagnose lambda values over a grid
#'
#' Fits [fit_mixed_subjects()] or [fit_mixed_subjects_split()] over a set of
#' candidate lambda values. The returned summary reports the fitted
#' mixed-subjects objective and the observed human expected-count loss for each
#' candidate. This is a sensitivity diagnostic, not a statistically valid PPI++
#' tuning rule.
#'
#' @param lambda_grid Numeric vector of lambda values in `[0, 1]`.
#' @param observed,predicted,generated Response matrices passed to
#'   [fit_mixed_subjects()].
#' @param split Logical; if `TRUE`, call [fit_mixed_subjects_split()].
#' @param ... Additional arguments passed to the selected fitting function.
#'
#' @return A list with `summary`, `lowest_observed_loss_lambda`, and all fitted
#'   model objects.
#' @export
#'
#' @examples
#' set.seed(3)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(30), pars)
#' predicted <- observed
#' generated <- simulate_2pl(rnorm(80), pars)
#' tuned <- diagnose_lambda_grid(
#'   c(0, 0.5),
#'   observed, predicted, generated,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$summary
diagnose_lambda_grid <- function(lambda_grid, observed, predicted, generated,
                                 split = FALSE, ...) {
  if (!is.numeric(lambda_grid) || length(lambda_grid) == 0) {
    stop("lambda_grid must be a non-empty numeric vector.", call. = FALSE)
  }
  lambda_grid <- sort(unique(lambda_grid))
  vapply(lambda_grid, validate_lambda, numeric(1))

  fits <- vector("list", length(lambda_grid))
  rows <- vector("list", length(lambda_grid))

  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    fits[[i]] <- if (isTRUE(split)) {
      fit_mixed_subjects_split(
        observed = observed,
        predicted = predicted,
        generated = generated,
        lambda = lambda,
        ...
      )
    } else {
      fit_mixed_subjects(
        observed = observed,
        predicted = predicted,
        generated = generated,
        lambda = lambda,
        ...
      )
    }

    observed_loss <- loss_expected_counts(
      fits[[i]]$q_observed$counts,
      fits[[i]]$item_pars
    )
    rows[[i]] <- data.frame(
      lambda = lambda,
      mixed_loss = fits[[i]]$value,
      observed_loss = observed_loss,
      convergence = fits[[i]]$convergence
    )
  }

  summary <- do.call(rbind, rows)
  list(
    summary = summary,
    lowest_observed_loss_lambda = summary$lambda[which.min(summary$observed_loss)],
    fits = fits
  )
}

#' Tune lambda over a grid
#'
#' @description
#' `tune_lambda_grid()` is retained as a backward-compatible wrapper for
#' [diagnose_lambda_grid()]. It emits a warning because the grid output is a
#' sensitivity diagnostic, not a statistically valid PPI++ tuning rule.
#'
#' @inheritParams diagnose_lambda_grid
#'
#' @return A list returned by [diagnose_lambda_grid()]. For backward
#'   compatibility, the list also contains `best_lambda_by_observed_loss`.
#' @export
tune_lambda_grid <- function(lambda_grid, observed, predicted, generated,
                             split = FALSE, ...) {
  warning(
    "tune_lambda_grid() is a diagnostic sensitivity helper, not a valid ",
    "PPI++ tuning rule. Prefer diagnose_lambda_grid() or a study-specific ",
    "cross-fitted tuning procedure.",
    call. = FALSE
  )

  out <- diagnose_lambda_grid(
    lambda_grid = lambda_grid,
    observed = observed,
    predicted = predicted,
    generated = generated,
    split = split,
    ...
  )
  out$best_lambda_by_observed_loss <- out$lowest_observed_loss_lambda
  out
}

#' Fit a mixed-subjects 2PL calibration with iterative EM
#'
#' Extends [fit_mixed_subjects()] by iterating the E-step and M-step until
#' convergence rather than fixing posterior quadrature weights at the initial
#' parameter estimates. At every iteration the posterior weights for all three
#' datasets (observed, predicted, generated) are recomputed using the same
#' current item parameters. This keeps the posteriors internally consistent and
#' avoids the asymmetry between `L_pred` and `L_gen` that arises when frozen
#' human-MLE weights are applied to LLM data with different item parameters.
#'
#' **Note on lambda selection.** This function accepts a fixed `lambda`. For
#' psychometric applications where accurate ability scoring is the goal, select
#' `lambda` with [tune_lambda_ability()] rather than [tune_lambda_ppi_score()].
#' The PPI++ score objective minimises the trace of the item-parameter
#' covariance matrix; [tune_lambda_ability()] minimises the propagated
#' ability-score risk `g' Sigma g`, which is the quantity that matters for
#' downstream test scoring.
#'
#' @param observed Human response matrix.
#' @param predicted LLM responses or probabilities for the same rows as
#'   `observed`.
#' @param generated Generated or unlabeled LLM responses.
#' @param lambda Power-tuning parameter in `[0, 1]`.
#' @param n_quad Number of standard-normal quadrature nodes.
#' @param initial_pars Optional starting item parameters. If omitted, a 2PL
#'   model is fit to `observed`.
#' @param quadrature Optional quadrature grid.
#' @param common_predicted_weights Logical; if `TRUE`, reuse observed posterior
#'   weights for `predicted` at each iteration.
#' @param paired_missing Missingness check passed to [fit_mixed_subjects()].
#' @param slope_lower Lower bound on discrimination parameters.
#' @param slope_upper Upper bound on discrimination parameters. **Strongly
#'   recommended** when `lambda > 0` — the iterative EM updates posteriors at
#'   each step, and without an upper bound the gradient asymmetry between
#'   `L_pred` and `L_gen` can compound across iterations, driving
#'   discrimination estimates to extreme values. A typical choice is
#'   `slope_upper = 4` or `slope_upper = 6`.
#' @param tol Convergence tolerance: maximum absolute change in any parameter
#'   across an EM iteration.
#' @param em_maxit Maximum number of EM iterations.
#' @param control Control list passed to [stats::optim()] for the M-step.
#' @param ... Additional arguments passed to [fit_2pl()] when `initial_pars`
#'   is omitted.
#'
#' @return An object of class `"mixedsubjects_fit"` with the standard fields
#'   plus `em_iterations` (number of EM cycles completed) and `em_converged`
#'   (logical).
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(40), pars)
#' predicted <- observed
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects_iterative(
#'   observed, predicted, generated,
#'   lambda = 0.5, initial_pars = pars, n_quad = 7,
#'   control = list(maxit = 50), em_maxit = 5
#' )
#' fit$item_pars
fit_mixed_subjects_iterative <- function(observed, predicted, generated,
                                          lambda = 1,
                                          n_quad = 31,
                                          initial_pars = NULL,
                                          quadrature = NULL,
                                          common_predicted_weights = TRUE,
                                          paired_missing = c("match_observed", "allow"),
                                          slope_lower = 1e-4,
                                          slope_upper = NULL,
                                          tol = 1e-4,
                                          em_maxit = 30,
                                          control = list(maxit = 200),
                                          ...) {
  lambda <- validate_lambda(lambda)
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
  paired_missing <- check_paired_missingness(
    observed = observed,
    predicted = predicted,
    paired_missing = paired_missing,
    common_predicted_weights = common_predicted_weights
  )

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  if (is.null(initial_pars)) {
    initial_fit <- fit_2pl(observed, ...)
    current_pars <- initial_fit$pars
  } else {
    initial_fit <- NULL
    current_pars <- standardize_item_pars(
      initial_pars,
      n_items = ncol(observed),
      item_names = colnames(observed)
    )
  }
  stored_initial_pars <- current_pars

  if (is.null(slope_upper) && lambda > 0) {
    warning(
      "fit_mixed_subjects_iterative() with lambda > 0 and no slope_upper ",
      "can diverge to extreme discrimination values. Setting slope_upper ",
      "(e.g. slope_upper = 6) is strongly recommended.",
      call. = FALSE
    )
  }

  em_iter <- 0L
  converged <- FALSE

  for (em_iter in seq_len(em_maxit)) {
    # E-step: use the SAME current_pars for all three datasets to keep
    # posteriors internally consistent across the three terms.
    q_obs <- build_quadrature_summary(observed, current_pars, quadrature)

    predicted_weights <- if (isTRUE(common_predicted_weights)) {
      q_obs$weights
    } else {
      NULL
    }
    q_pred <- build_quadrature_summary(
      predicted, current_pars, quadrature, weights = predicted_weights
    )
    q_gen <- build_quadrature_summary(generated, current_pars, quadrature)

    # M-step
    fit <- fit_from_counts(
      counts_observed  = q_obs$counts,
      counts_predicted = q_pred$counts,
      counts_generated = q_gen$counts,
      initial_pars     = current_pars,
      lambda           = lambda,
      slope_lower      = slope_lower,
      slope_upper      = slope_upper,
      control          = control
    )

    new_pars <- fit$item_pars
    delta    <- max(abs(new_pars$a - current_pars$a),
                    abs(new_pars$d - current_pars$d))
    current_pars <- new_pars

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  # Final E-step with converged parameters (for vcov and ability risk)
  q_obs_final <- build_quadrature_summary(observed, current_pars, quadrature)
  pred_w_final <- if (isTRUE(common_predicted_weights)) q_obs_final$weights else NULL
  q_pred_final <- build_quadrature_summary(
    predicted, current_pars, quadrature, weights = pred_w_final
  )
  q_gen_final <- build_quadrature_summary(generated, current_pars, quadrature)

  out <- c(
    fit,
    list(
      lambda                   = lambda,
      initial_pars             = stored_initial_pars,
      initial_model            = if (is.null(initial_fit)) NULL else initial_fit$model,
      quadrature               = quadrature,
      q_observed               = q_obs_final,
      q_predicted              = q_pred_final,
      q_generated              = q_gen_final,
      common_predicted_weights = common_predicted_weights,
      paired_missing           = paired_missing,
      split                    = NULL,
      em_iterations            = em_iter,
      em_converged             = converged,
      call                     = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' @export
print.mixedsubjects_fit <- function(x, ...) {
  cat("mixedsubjectsirt 2PL fit\n")
  cat("  items:      ", nrow(x$item_pars), "\n", sep = "")
  cat("  lambda:     ", x$lambda, "\n", sep = "")
  if (!is.null(x$lambda_by_split) && length(unique(x$lambda_by_split)) > 1) {
    cat("  lambda bar: ", signif(x$lambda_generated, 6), "\n", sep = "")
  }
  cat("  loss:       ", signif(x$value, 6), "\n", sep = "")
  cat("  convergence:", x$convergence, "\n")
  if (!is.null(x$split)) {
    cat("  splits:     ", x$split$n_splits, "\n", sep = "")
  }
  if (!is.null(x$em_iterations)) {
    cat("  EM iters:   ", x$em_iterations,
        if (isTRUE(x$em_converged)) " (converged)" else " (max reached)",
        "\n", sep = "")
  }
  invisible(x)
}
