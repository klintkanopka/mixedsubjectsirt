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
                               slope_lower = 1e-4,
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
#' @param lambda Power-tuning parameter in `[0, 1]`.
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
                                     slope_lower = 1e-4,
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

  fit <- fit_from_counts(
    counts_observed = counts_observed,
    counts_predicted = counts_predicted,
    counts_generated = counts_generated,
    initial_pars = initial_average,
    lambda = lambda,
    slope_lower = slope_lower,
    control = control
  )

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
      initial_pars = initial_average,
      initial_model = NULL,
      quadrature = quadrature,
      q_observed = q_observed,
      q_predicted = q_predicted,
      q_generated = q_generated,
      common_predicted_weights = common_predicted_weights,
      split = list(
        n_splits = n_splits,
        split_id = split_id,
        fold_initial_pars = fold_initial_pars,
        fold_models = fold_models
      ),
      call = match.call()
    )
  )
  class(out) <- "mixedsubjects_fit"
  out
}

#' Tune lambda over a grid
#'
#' Fits [fit_mixed_subjects()] or [fit_mixed_subjects_split()] over a set of
#' candidate lambda values. The returned summary reports the fitted
#' mixed-subjects objective and the observed human expected-count loss for each
#' candidate. These diagnostics are not a replacement for a study-specific
#' validation or bootstrap procedure, but they are useful for sensitivity checks.
#'
#' @param lambda_grid Numeric vector of lambda values in `[0, 1]`.
#' @param observed,predicted,generated Response matrices passed to
#'   [fit_mixed_subjects()].
#' @param split Logical; if `TRUE`, call [fit_mixed_subjects_split()].
#' @param ... Additional arguments passed to the selected fitting function.
#'
#' @return A list with `summary`, `best_lambda_by_observed_loss`, and all fitted
#'   model objects.
#' @export
#'
#' @examples
#' set.seed(3)
#' pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
#' observed <- simulate_2pl(rnorm(30), pars)
#' predicted <- observed
#' generated <- simulate_2pl(rnorm(80), pars)
#' tuned <- tune_lambda_grid(
#'   c(0, 0.5),
#'   observed, predicted, generated,
#'   initial_pars = pars, n_quad = 5, control = list(maxit = 30)
#' )
#' tuned$summary
tune_lambda_grid <- function(lambda_grid, observed, predicted, generated,
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
    best_lambda_by_observed_loss = summary$lambda[which.min(summary$observed_loss)],
    fits = fits
  )
}

#' @export
print.mixedsubjects_fit <- function(x, ...) {
  cat("mixedsubjectsirt 2PL fit\n")
  cat("  items:      ", nrow(x$item_pars), "\n", sep = "")
  cat("  lambda:     ", x$lambda, "\n", sep = "")
  cat("  loss:       ", signif(x$value, 6), "\n", sep = "")
  cat("  convergence:", x$convergence, "\n")
  if (!is.null(x$split)) {
    cat("  splits:     ", x$split$n_splits, "\n", sep = "")
  }
  invisible(x)
}
