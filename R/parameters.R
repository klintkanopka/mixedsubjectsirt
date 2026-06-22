#' Fit a unidimensional 2PL IRT model
#'
#' Fits a two-parameter logistic model with `mirt` and returns item parameters in
#' slope-intercept form. The response probability is
#' `plogis(d + a * theta)`, where `a` is the discrimination and `d` is the
#' intercept. Difficulty is returned as `b = -d / a`.
#'
#' @param resp A numeric item response matrix with rows for subjects and columns
#'   for items. Values must be binary `0`/`1`; `NA` is allowed.
#' @param technical A list passed to the `technical` argument of
#'   [mirt::mirt()].
#' @param verbose Logical; passed to [mirt::mirt()].
#' @param ... Additional arguments passed to [mirt::mirt()].
#'
#' @return A list with `pars`, a data frame containing `item`, `a`, `d`, and
#'   `b`, and `model`, the fitted `mirt` model.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2, 0.9, 1.1, 0.8), d = c(0, 0.5, -0.5, 0.2, -0.3))
#' resp <- simulate_2pl(rnorm(500), pars)
#' fit <- fit_2pl(resp)
#' fit$pars
fit_2pl <- function(resp, technical = list(NCYCLES = 1000),
                    verbose = FALSE, ...) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = FALSE)

  model <- mirt::mirt(
    data = as.data.frame(resp),
    model = 1,
    itemtype = "2PL",
    technical = technical,
    verbose = verbose,
    ...
  )

  list(
    pars = standardize_item_pars(model, n_items = ncol(resp),
                                 item_names = colnames(resp)),
    model = model
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
