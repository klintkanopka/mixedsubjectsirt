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
mixed_subjects_loss <- function(
  pars,
  q_observed,
  q_predicted,
  q_llm,
  lambda = 0
) {
  lambda <- validate_lambda(lambda)
  counts_observed <- q_to_counts(q_observed, "q_observed")
  counts_predicted <- q_to_counts(q_predicted, "q_predicted")
  counts_llm <- q_to_counts(q_llm, "q_llm")
  check_counts_compatible(list(counts_observed, counts_predicted, counts_llm))

  item_pars <- item_pars_from_vector(pars, counts_observed$item_names)

  loss_expected_counts(counts_observed, item_pars) +
    lambda *
      (loss_expected_counts(counts_llm, item_pars) -
        loss_expected_counts(counts_predicted, item_pars))
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
#' @param predicted Binary LLM responses (0/1) for the same rows and items as
#'   `observed`. Probabilities are not accepted: fractional values are not a
#'   valid likelihood input for the marginal IRT objective and break the PPI
#'   correction, so sample binary responses from any probabilities first (e.g.
#'   `rbinom`).
#' @param generated Binary generated or unlabeled LLM responses (0/1) for the
#'   same item columns. Probabilities are not accepted (see `predicted`).
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
#' @param slope_upper Upper bound for discrimination parameters during
#'   optimization. Use `NULL` (the default) for no upper bound. Setting a finite
#'   bound (e.g. 4) can stabilize the frozen expected-count fit when the LLM
#'   parameters differ substantially from the human parameters.
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
fit_mixed_subjects <- function(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  common_predicted_weights = TRUE,
  paired_missing = c("match_observed", "allow"),
  slope_lower = 1e-4,
  slope_upper = NULL,
  control = list(maxit = 500),
  ...
) {
  lambda <- validate_lambda(lambda)
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
  generated <- validate_response_matrix(
    generated,
    "generated",
    allow_fractional = FALSE
  )
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")

  if (nrow(observed) != nrow(predicted)) {
    stop(
      "observed and predicted must have the same number of rows.",
      call. = FALSE
    )
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
#' @param slope_upper Upper bound for discrimination parameters during
#'   optimization. Use `NULL` (the default) for no upper bound.
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
fit_mixed_subjects_from_quadrature <- function(
  q_observed,
  q_predicted,
  q_generated,
  lambda = 1,
  initial_pars = NULL,
  slope_lower = 1e-4,
  slope_upper = NULL,
  control = list(maxit = 500)
) {
  counts_observed <- q_to_counts(q_observed, "q_observed")
  counts_predicted <- q_to_counts(q_predicted, "q_predicted")
  counts_generated <- q_to_counts(q_generated, "q_generated")
  check_counts_compatible(list(
    counts_observed,
    counts_predicted,
    counts_generated
  ))

  if (is.null(initial_pars)) {
    if (is.list(q_observed) && !is.null(q_observed$irt_pars)) {
      initial_pars <- q_observed$irt_pars
    } else {
      stop(
        "initial_pars is required when q_observed does not contain irt_pars.",
        call. = FALSE
      )
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
#' @param predicted Binary LLM responses (0/1) for the same rows and items as
#'   `observed`. Probabilities are not accepted: fractional values are not a
#'   valid likelihood input for the marginal IRT objective and break the PPI
#'   correction, so sample binary responses from any probabilities first (e.g.
#'   `rbinom`).
#' @param generated Binary generated or unlabeled LLM responses (0/1) for the
#'   same item columns. Probabilities are not accepted (see `predicted`).
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
#' @param slope_upper Upper bound for discrimination parameters during
#'   optimization. Use `NULL` (the default) for no upper bound.
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
fit_mixed_subjects_split <- function(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_splits = 2,
  split_id = NULL,
  seed = NULL,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  common_predicted_weights = TRUE,
  paired_missing = c("match_observed", "allow"),
  slope_lower = 1e-4,
  slope_upper = NULL,
  control = list(maxit = 500),
  ...
) {
  lambda <- validate_lambda_vector(lambda)
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
  generated <- validate_response_matrix(
    generated,
    "generated",
    allow_fractional = FALSE
  )
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")

  if (nrow(observed) != nrow(predicted)) {
    stop(
      "observed and predicted must have the same number of rows.",
      call. = FALSE
    )
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
#' candidate. This is a sensitivity diagnostic, not a valid tuning rule.
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
diagnose_lambda_grid <- function(
  lambda_grid,
  observed,
  predicted,
  generated,
  split = FALSE,
  ...
) {
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
    lowest_observed_loss_lambda = summary$lambda[which.min(
      summary$observed_loss
    )],
    fits = fits
  )
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
#' `lambda` with [tune_lambda_ability_risk()] rather than [tune_lambda_ppi_score()].
#' The PPI++ score objective minimizes the trace of the item-parameter
#' covariance matrix; [tune_lambda_ability_risk()] minimizes the propagated
#' ability-score risk `g' Sigma g`, which is the quantity that matters for
#' downstream test scoring.
#'
#' @inheritParams fit_mixed_subjects
#' @param slope_upper Upper bound on discrimination parameters. **Strongly
#'   recommended** when `lambda > 0` — the iterative EM updates posteriors at
#'   each step, and without an upper bound the gradient asymmetry between
#'   `L_pred` and `L_gen` can compound across iterations, driving
#'   discrimination estimates to extreme values. A typical choice is
#'   `slope_upper = 4` or `slope_upper = 6`.
#' @param tol Convergence tolerance: maximum absolute change in any parameter
#'   across an EM iteration.
#' @param em_maxit Maximum number of EM iterations.
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
fit_mixed_subjects_iterative <- function(
  observed,
  predicted,
  generated,
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
  ...
) {
  lambda <- validate_lambda(lambda)
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
  generated <- validate_response_matrix(
    generated,
    "generated",
    allow_fractional = FALSE
  )
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")

  if (nrow(observed) != nrow(predicted)) {
    stop(
      "observed and predicted must have the same number of rows.",
      call. = FALSE
    )
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
      predicted,
      current_pars,
      quadrature,
      weights = predicted_weights
    )
    q_gen <- build_quadrature_summary(generated, current_pars, quadrature)

    # M-step
    fit <- fit_from_counts(
      counts_observed = q_obs$counts,
      counts_predicted = q_pred$counts,
      counts_generated = q_gen$counts,
      initial_pars = current_pars,
      lambda = lambda,
      slope_lower = slope_lower,
      slope_upper = slope_upper,
      control = control
    )

    new_pars <- fit$item_pars
    delta <- max(
      abs(new_pars$a - current_pars$a),
      abs(new_pars$d - current_pars$d)
    )
    current_pars <- new_pars

    if (delta < tol) {
      converged <- TRUE
      break
    }
  }

  # Final E-step with converged parameters (for vcov and ability risk)
  q_obs_final <- build_quadrature_summary(observed, current_pars, quadrature)
  pred_w_final <- if (isTRUE(common_predicted_weights)) {
    q_obs_final$weights
  } else {
    NULL
  }
  q_pred_final <- build_quadrature_summary(
    predicted,
    current_pars,
    quadrature,
    weights = pred_w_final
  )
  q_gen_final <- build_quadrature_summary(generated, current_pars, quadrature)

  out <- c(
    fit,
    list(
      lambda = lambda,
      initial_pars = stored_initial_pars,
      initial_model = if (is.null(initial_fit)) NULL else initial_fit$model,
      quadrature = quadrature,
      q_observed = q_obs_final,
      q_predicted = q_pred_final,
      q_generated = q_gen_final,
      common_predicted_weights = common_predicted_weights,
      paired_missing = paired_missing,
      split = NULL,
      em_iterations = em_iter,
      em_converged = converged,
      call = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' Fit a mixed-subjects 2PL calibration via marginal maximum likelihood
#'
#' Estimates item parameters using the true IRT marginal likelihood for all
#' three loss terms.  Unlike [fit_mixed_subjects()], which freezes posterior
#' quadrature weights at the initial parameter estimates before optimizing,
#' this function recomputes posterior weights at every gradient evaluation.
#' This eliminates the gradient asymmetry that causes [fit_mixed_subjects()] to
#' converge to false minima at inflated discrimination values when LLM item
#' parameters differ from human parameters.
#'
#' **Why it matters for lambda selection.**  With the frozen expected-count
#' implementation, the gradient of `L_pred` uses concentrated human posteriors
#' while `L_gen` uses diffuse LLM posteriors, making
#' `grad(L_pred) >> grad(L_gen)` and systematically pushing discriminations
#' upward at any `lambda > 0`.  In the marginal-MML formulation all three terms
#' use their own current-parameter posteriors, so the asymmetry is absent at the
#' true optimum. As a result [tune_lambda_ability_risk()] selects `lambda > 0`
#' whenever the LLM predictions are genuinely informative (e.g. `predicted =
#' observed`), rather than collapsing to `lambda = 0` for all misaligned LLMs.
#'
#' **`mml_pred_weights`.**
#' \describe{
#'   \item{`"own"` (default)}{L_pred uses posteriors computed from the
#'     *predicted* response matrix at the current parameter values.  All three
#'     terms are true marginal likelihoods; objective and gradient are
#'     internally consistent.  Recommended for most applications and required
#'     for [vcov_mixed_subjects_mml()] to produce the fully correct
#'     Louis-formula bread.}
#'   \item{`"human"`}{L_pred uses posteriors computed from the *observed*
#'     (human) response matrix, frozen at `initial_pars`.  This is a
#'     **fixed-nuisance Q-function**: the predicted term is treated as a frozen
#'     expected-count lower bound rather than a true marginal likelihood.
#'     Objective and gradient are mutually consistent (both use the same frozen
#'     posteriors) so L-BFGS-B converges correctly.  Useful when strong
#'     ability-level pairing is needed.  Note that [vcov_mixed_subjects_mml()]
#'     applies Louis' formula to the stored fixed posteriors, which is
#'     approximately correct when `initial_pars` ≈ `conv_pars`.}
#' }
#'
#' **Per-item lambda (vector `lambda`).**  When `lambda` is a length-`n_items`
#' vector rather than a scalar, `fit_mixed_subjects_mml` switches to a
#' **frozen Q-function** objective: expected-count counts are computed once from
#' `initial_pars` and held fixed during L-BFGS-B, with item `j`'s counts
#' weighted by `lambda[j]`.  This is a consistent (objective, gradient) pair
#' but is *not* the full marginal-MML objective — it is a frozen expected-count
#' approximation analogous to [fit_mixed_subjects()].  Per-item lambda values
#' obtained from [tune_lambda_ability_risk_item()] assign `lambda_j ≈ 0` to
#' items where the LLM correction is harmful, containing the frozen-posterior
#' gradient asymmetry.  Document per-item lambda results as approximate.
#'
#' @inheritParams fit_mixed_subjects
#' @param mml_pred_weights How to compute posteriors for the paired `predicted`
#'   term.  `"own"` uses posteriors from the predicted responses; `"human"`
#'   uses posteriors from the observed human responses.  See Details.
#' @param slope_upper Upper bound on discrimination parameters. Unlike
#'   [fit_mixed_subjects()], this function should not require capping for
#'   well-posed problems because the true marginal objective has no false
#'   minimum at large discrimination.
#'
#' @return An object of class `"mixedsubjects_fit"` with the same structure as
#'   [fit_mixed_subjects()].  For **scalar** lambda fits, the quadrature
#'   summaries store posteriors at the converged parameters, and
#'   `stats::vcov()` dispatches automatically to
#'   [vcov_mixed_subjects_mml()] to compute the Louis-corrected marginal
#'   sandwich covariance.  Calling [vcov_mixed_subjects()] directly bypasses
#'   the Louis correction.  For **vector** lambda fits, the summaries store
#'   the frozen posteriors used during optimization, and `stats::vcov()`
#'   dispatches to [vcov_mixed_subjects()] (EM bread) for consistency with the
#'   frozen Q-function objective.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed  <- simulate_2pl(rnorm(40), pars)
#' generated <- simulate_2pl(rnorm(100), pars)
#' fit <- fit_mixed_subjects_mml(
#'   observed, observed, generated,
#'   lambda = 0.5, initial_pars = pars, n_quad = 7,
#'   control = list(maxit = 100)
#' )
#' fit$item_pars
fit_mixed_subjects_mml <- function(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  mml_pred_weights = c("own", "human"),
  slope_lower = 1e-4,
  slope_upper = NULL,
  control = list(maxit = 500),
  ...
) {
  # Defer final lambda validation to after n_items is known (line below)
  lambda <- validate_lambda_vector(lambda) # basic range/type check
  mml_pred_weights <- match.arg(mml_pred_weights)

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
  generated <- validate_response_matrix(
    generated,
    "generated",
    allow_fractional = FALSE
  )
  check_same_items(observed, predicted, "observed", "predicted")
  check_same_items(observed, generated, "observed", "generated")

  if (nrow(observed) != nrow(predicted)) {
    stop(
      "observed and predicted must have the same number of rows.",
      call. = FALSE
    )
  }

  quadrature <- check_quadrature(quadrature, n_quad = n_quad)

  if (is.null(initial_pars)) {
    initial_fit <- fit_2pl(observed, ...)
    init_std <- initial_fit$pars
  } else {
    initial_fit <- NULL
    init_std <- standardize_item_pars(
      initial_pars,
      n_items = ncol(observed),
      item_names = colnames(observed)
    )
  }

  item_names <- colnames(observed)
  n_items <- length(item_names)

  # Validate and expand lambda: scalar or length-n_items vector.
  # Scalar path (Phase 2): true marginal likelihoods; posteriors recomputed at
  #   every gradient evaluation. Fully consistent objective and gradient.
  # Vector path (Phase 3): frozen expected-count objective with per-item λ_j.
  #   Posteriors frozen at initial_pars (same as fit_mixed_subjects). The frozen
  #   Q-function and its gradient are fully consistent for L-BFGS-B.
  #   Per-item λ_j from tune_lambda_ability_risk_item assigns λ_j ≈ 0 to items
  #   where the LLM correction is harmful, containing the gradient asymmetry.
  lambda <- validate_lambda_vector(lambda, n = n_items)
  scalar_lambda <- length(lambda) == 1
  lambda_vec <- if (scalar_lambda) rep(lambda, n_items) else lambda

  # For the vector path, precompute frozen expected counts from initial_pars.
  if (!scalar_lambda) {
    obs_w_frozen <- posterior_weights_2pl(observed, init_std, quadrature)
    gen_w_frozen <- posterior_weights_2pl(generated, init_std, quadrature)
    obs_cts_frozen <- summarize_expected_counts(observed, obs_w_frozen)
    pred_cts_frozen <- summarize_expected_counts(predicted, obs_w_frozen)
    gen_cts_frozen <- summarize_expected_counts(generated, gen_w_frozen)
  }

  # For mml_pred_weights = "human": precompute FIXED-NUISANCE human posteriors
  # from init_std and freeze them throughout optimization. This makes the
  # objective and gradient fully consistent (both use the same frozen Q-function
  # for L_pred). The old behavior — recomputing human posteriors at every
  # candidate ip — was inconsistent: the objective and gradient differed by the
  # chain-rule term through the posteriors, causing L-BFGS-B line-search issues.
  if (scalar_lambda && mml_pred_weights == "human") {
    pred_w_fixed <- posterior_weights_2pl(observed, init_std, quadrature)
    pred_cts_fixed <- summarize_expected_counts(predicted, pred_w_fixed)
  }

  # Shared cache: reuse posteriors when fn and gr are called at the same par.
  .cache <- new.env(parent = emptyenv())
  .cache$par <- NULL
  .cache$val <- NULL
  .cache$grd <- NULL

  recompute <- function(par) {
    if (
      !is.null(.cache$par) &&
        isTRUE(all.equal(par, .cache$par, tolerance = 0))
    ) {
      return(invisible(NULL))
    }

    ip <- item_pars_from_vector(par, item_names)

    if (scalar_lambda) {
      # --- Scalar path: true marginal likelihoods (Phase 2) ---
      obs_r <- posterior_and_log_lik_2pl(observed, ip, quadrature)
      gen_r <- posterior_and_log_lik_2pl(generated, ip, quadrature)

      obs_counts <- summarize_expected_counts(observed, obs_r$weights)
      gen_counts <- summarize_expected_counts(generated, gen_r$weights)

      if (mml_pred_weights == "own") {
        # Pure MML: L_pred uses its own posteriors (consistent marginal objective)
        pred_r <- posterior_and_log_lik_2pl(predicted, ip, quadrature)
        pred_counts <- summarize_expected_counts(predicted, pred_r$weights)
        l_pred <- -mean(pred_r$log_normalizers)
      } else {
        # Fixed-nuisance Q-function: L_pred uses posteriors frozen at init_std.
        # pred_cts_fixed is precomputed outside this closure and never updated.
        # Both objective and gradient use these fixed counts → fully consistent.
        pred_counts <- pred_cts_fixed
        l_pred <- loss_expected_counts(pred_cts_fixed, ip)
      }

      .cache$val <- -mean(obs_r$log_normalizers) +
        lambda * (-mean(gen_r$log_normalizers) - l_pred)
      .cache$grd <- gradient_expected_counts(obs_counts, ip) +
        lambda *
          (gradient_expected_counts(gen_counts, ip) -
            gradient_expected_counts(pred_counts, ip))
    } else {
      # --- Vector-lambda path: frozen Q-function with per-item λ_j (Phase 3) ---
      # Objective: Σ_j [ L_obs_j^EC + λ_j*(L_gen_j^EC - L_pred_j^EC) ]
      # Gradient : Σ_j [ ∂L_obs_j/∂(a_j,d_j) + λ_j*(∂L_gen_j - ∂L_pred_j) ]
      # Both use frozen counts — fully consistent pair for L-BFGS-B.
      obs_ig <- item_loss_and_grad(obs_cts_frozen, ip)
      gen_ig <- item_loss_and_grad(gen_cts_frozen, ip)
      pred_ig <- item_loss_and_grad(pred_cts_frozen, ip)

      .cache$val <- sum(obs_ig$loss) +
        sum(lambda_vec * (gen_ig$loss - pred_ig$loss))
      .cache$grd <- c(
        obs_ig$grad_a + lambda_vec * (gen_ig$grad_a - pred_ig$grad_a),
        obs_ig$grad_d + lambda_vec * (gen_ig$grad_d - pred_ig$grad_d)
      )
    }

    .cache$par <- par
    invisible(NULL)
  }

  objective <- function(par) {
    disc <- par[seq_len(n_items)]
    if (any(!is.finite(disc)) || any(disc <= 0)) {
      return(.Machine$double.xmax)
    }
    recompute(par)
    val <- .cache$val
    if (!is.finite(val)) .Machine$double.xmax else val
  }

  gradient <- function(par) {
    disc <- par[seq_len(n_items)]
    if (any(!is.finite(disc)) || any(disc <= 0)) {
      return(rep(0, 2L * n_items))
    }
    recompute(par)
    g <- .cache$grd
    ifelse(is.finite(g), g, 0)
  }

  start <- vector_from_item_pars(init_std)

  if (is.null(slope_lower)) {
    lower <- rep(-Inf, length(start))
  } else {
    lower <- c(rep(slope_lower, n_items), rep(-Inf, n_items))
    start[seq_len(n_items)] <- pmax(start[seq_len(n_items)], slope_lower)
  }
  if (is.null(slope_upper)) {
    upper <- rep(Inf, length(start))
  } else {
    upper <- c(rep(slope_upper, n_items), rep(Inf, n_items))
    start[seq_len(n_items)] <- pmin(start[seq_len(n_items)], slope_upper)
  }

  ctrl <- utils::modifyList(list(maxit = 500), control)
  opt <- stats::optim(
    par = start,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = ctrl
  )

  conv_pars <- item_pars_from_vector(opt$par, item_names)

  # Build final quadrature summaries.
  # - Scalar MML path: posteriors at the CONVERGED parameters so that
  #   vcov_mixed_subjects_mml() applies Louis' formula at the MLE.
  # - Vector-lambda path: posteriors frozen at INIT_STD (same as those used
  #   during optimization) so that vcov_mixed_subjects() is consistent with
  #   the frozen Q-function that was actually optimized.  Storing converged
  #   posteriors here would mean vcov uses a different estimating equation
  #   from the one that determined conv_pars, making the sandwich invalid.
  if (scalar_lambda) {
    q_obs_final <- build_quadrature_summary(observed, conv_pars, quadrature)

    q_pred_final <- if (mml_pred_weights == "own") {
      build_quadrature_summary(predicted, conv_pars, quadrature)
    } else {
      # "human": fixed-nuisance posteriors from init_std
      build_quadrature_summary(
        predicted,
        conv_pars,
        quadrature,
        weights = pred_w_fixed
      )
    }

    q_gen_final <- build_quadrature_summary(generated, conv_pars, quadrature)
  } else {
    # Vector-lambda: store the frozen weights/counts that were used in fitting
    q_obs_final <- build_quadrature_summary(
      observed,
      init_std,
      quadrature,
      weights = obs_w_frozen
    )
    q_pred_final <- build_quadrature_summary(
      predicted,
      init_std,
      quadrature,
      weights = obs_w_frozen
    )
    q_gen_final <- build_quadrature_summary(
      generated,
      init_std,
      quadrature,
      weights = gen_w_frozen
    )
  }

  out <- c(
    list(
      item_pars = conv_pars,
      par = opt$par,
      value = opt$value,
      convergence = opt$convergence,
      message = opt$message,
      optimizer = opt
    ),
    list(
      lambda = if (scalar_lambda) lambda else lambda_vec,
      initial_pars = init_std,
      initial_model = if (is.null(initial_fit)) NULL else initial_fit$model,
      quadrature = quadrature,
      q_observed = q_obs_final,
      q_predicted = q_pred_final,
      q_generated = q_gen_final,
      common_predicted_weights = (mml_pred_weights == "human"),
      mml_pred_weights = mml_pred_weights,
      mml = TRUE,
      paired_missing = NA,
      split = NULL,
      call = match.call()
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
    cat(
      "  EM iters:   ",
      x$em_iterations,
      if (isTRUE(x$em_converged)) " (converged)" else " (max reached)",
      "\n",
      sep = ""
    )
  }
  if (isTRUE(x$mml)) {
    cat("  estimator:  marginal MML PPI++\n")
  }
  invisible(x)
}
