#' Fit a unidimensional 2PL IRT model
#'
#' Estimates per-item discriminations `a_j` and intercepts `d_j` by maximizing
#' the IRT marginal likelihood under a standard-normal ability prior using
#' L-BFGS-B. The response probability is `plogis(d + a * theta)`, where `a` is
#' the discrimination and `d` is the intercept. Difficulty is returned as
#' `b = -d / a`.
#'
#' This is the 2PL counterpart of [fit_1pl()] and uses the same marginal
#' likelihood, quadrature, and gradient machinery as
#' [fit_mixed_subjects_mml()] at `lambda = 0`, so a human-only baseline fit and
#' the mixed-subjects estimator are optimized on a common objective.
#'
#' @param resp A numeric item response matrix with rows for subjects and columns
#'   for items. Values must be binary `0`/`1`; `NA` is allowed.
#' @param n_quad Number of standard-normal quadrature nodes.
#' @param initial_pars Optional starting item parameters (matrix or data frame
#'   with `a`/`a1` and `d` columns). If omitted, `a_j = 1` and
#'   `d_j = qlogis(p_j)` where `p_j` is the observed proportion correct for
#'   item `j`.
#' @param quadrature Optional quadrature grid with `theta` and `weight` columns.
#' @param slope_lower,slope_upper Bounds on the discriminations. `NULL` leaves
#'   the corresponding side unbounded.
#' @param control Control list passed to [stats::optim()].
#'
#' @return A list with `pars`, a data frame containing `item`, `a`, `d`, and
#'   `b`; `par`, the raw parameter vector; optimizer details (`value`,
#'   `convergence`, `message`); and `model`, which is `NULL` and retained for
#'   backward compatibility.
#' @export
#'
#' @seealso [fit_1pl()] for the shared-discrimination version.
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9, 1.1, 0.8), d = c(0, 0.5, -0.5, 0.2, -0.3))
#' resp <- simulate_2pl(rnorm(500), pars)
#' fit <- fit_2pl(resp)
#' fit$pars
fit_2pl <- function(
  resp,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  slope_lower = 1e-4,
  slope_upper = NULL,
  control = list(maxit = 500)
) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = FALSE)
  quadrature <- check_quadrature(quadrature, n_quad = n_quad)
  item_names <- colnames(resp)
  n_items <- length(item_names)

  if (is.null(initial_pars)) {
    p_cor <- colMeans(resp, na.rm = TRUE)
    d_init <- stats::qlogis(pmax(0.05, pmin(0.95, p_cor)))
    initial_pars <- data.frame(
      item = item_names,
      a = 1,
      d = d_init,
      stringsAsFactors = FALSE
    )
    initial_pars$b <- -initial_pars$d / initial_pars$a
  } else {
    initial_pars <- standardize_item_pars(
      initial_pars,
      n_items = n_items,
      item_names = item_names
    )
  }

  .cache <- new.env(parent = emptyenv())
  .cache$par <- NULL
  .cache$val <- NULL
  .cache$grd <- NULL

  recompute <- function(par) {
    if (
      !is.null(.cache$par) && isTRUE(all.equal(par, .cache$par, tolerance = 0))
    ) {
      return(invisible(NULL))
    }
    disc <- par[seq_len(n_items)]
    if (any(!is.finite(disc)) || any(disc <= 0)) {
      .cache$val <- .Machine$double.xmax
      .cache$grd <- rep(0, 2L * n_items)
      .cache$par <- par
      return(invisible(NULL))
    }
    ip <- item_pars_from_vector(par, item_names)
    val <- marginal_loss_2pl(resp, ip, quadrature)
    g <- marginal_gradient_2pl(resp, ip, quadrature)
    .cache$par <- par
    .cache$val <- if (is.finite(val)) val else .Machine$double.xmax
    .cache$grd <- ifelse(is.finite(g), g, 0)
    invisible(NULL)
  }

  objective <- function(par) {
    recompute(par)
    .cache$val
  }
  gradient <- function(par) {
    recompute(par)
    .cache$grd
  }

  start <- vector_from_item_pars(initial_pars)
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

  list(
    pars = item_pars_from_vector(opt$par, item_names),
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    model = NULL
  )
}

#' Link item parameters onto a target scale
#'
#' Applies mean-mean linking to express source item parameters on the scale of a
#' target calibration. Both parameter sets must be in slope-intercept form for
#' the model `plogis(d + a * theta)`.
#'
#' If `theta_target = A * theta_source + B`, then source parameters transform as
#' `a_target = a_source / A` and `b_target = A * b_source + B`, with
#' `d_target = -a_target * b_target`. Mean-mean linking chooses `A` and `B` so
#' that the transformed source parameters match the target mean discrimination
#' and mean difficulty.
#'
#' @param source Item parameters to transform. A matrix or data frame with
#'   columns `a`/`a1` and `d`, or a fitted `mirt` model.
#' @param target Item parameters defining the target scale. Uses the same
#'   accepted formats as `source`.
#' @param method Linking method. Currently `"mean_mean"` and `"none"` are
#'   supported.
#'
#' @return A list with transformed `pars`, linking constants `A` and `B`, and
#'   the selected `method`.
#' @export
#'
#' @examples
#' source <- data.frame(a = c(0.8, 1.2), d = c(-0.2, 0.5))
#' target <- data.frame(a = c(1.0, 1.5), d = c(-0.1, 0.4))
#' link_item_parameters(source, target)$pars
link_item_parameters <- function(source, target,
                                 method = c("mean_mean", "none")) {
  method <- match.arg(method)
  source <- standardize_item_pars(source, name = "source")
  target <- standardize_item_pars(target, name = "target")

  if (!identical(source$item, target$item)) {
    if (all(source$item %in% target$item)) {
      target <- target[match(source$item, target$item), , drop = FALSE]
    } else {
      stop("source and target must contain the same items.", call. = FALSE)
    }
  }

  if (method == "none") {
    return(list(pars = source, A = 1, B = 0, method = method))
  }

  A <- mean(source$a) / mean(target$a)
  B <- mean(target$b) - A * mean(source$b)

  out <- source
  out$a <- source$a / A
  out$b <- A * source$b + B
  out$d <- -out$a * out$b

  list(pars = out, A = A, B = B, method = method)
}
