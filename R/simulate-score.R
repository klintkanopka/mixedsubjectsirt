#' Simulate 2PL item responses
#'
#' Generates binary item responses from the model `plogis(d + a * theta)`.
#'
#' @param theta Numeric vector of latent trait values.
#' @param item_pars Item parameters in slope-intercept form. Supply a data frame
#'   or matrix with columns `a`/`a1` and `d`, or a fitted `mirt` model.
#'
#' @return A binary response matrix with one row per value of `theta` and one
#'   column per item.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' simulate_2pl(rnorm(5), pars)
simulate_2pl <- function(theta, item_pars) {
  theta <- as.numeric(theta)
  if (length(theta) == 0 || any(!is.finite(theta))) {
    stop("theta must be a non-empty finite numeric vector.", call. = FALSE)
  }

  item_pars <- standardize_item_pars(item_pars)
  eta <- outer(theta, item_pars$a, `*`) +
    matrix(item_pars$d, nrow = length(theta), ncol = nrow(item_pars),
           byrow = TRUE)
  p <- stats::plogis(eta)
  out <- matrix(
    stats::rbinom(length(p), size = 1, prob = as.vector(p)),
    nrow = length(theta),
    ncol = nrow(item_pars)
  )
  colnames(out) <- item_pars$item
  out
}

#' Estimate ability scores from a 2PL calibration
#'
#' Computes bounded maximum-likelihood ability estimates for response patterns
#' under fixed item parameters. This is a scoring helper for inspecting fitted
#' calibrations; it does not account for uncertainty in the item parameters.
#'
#' @param resp Response matrix with rows for subjects and columns for items.
#' @param item_pars Item parameters in slope-intercept form. Supply a data frame
#'   or matrix with columns `a`/`a1` and `d`, or a fitted `mirt` model.
#' @param bounds Numeric vector of length two giving the optimization interval
#'   for theta.
#'
#' @return A numeric vector of ability estimates.
#' @export
#'
#' @examples
#' set.seed(1)
#' pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
#' resp <- simulate_2pl(rnorm(5), pars)
#' score_theta(resp, pars)
score_theta <- function(resp, item_pars, bounds = c(-6, 6)) {
  resp <- validate_response_matrix(resp, name = "resp", allow_fractional = TRUE)
  item_pars <- standardize_item_pars(item_pars, n_items = ncol(resp),
                                     item_names = colnames(resp))

  if (!is.numeric(bounds) || length(bounds) != 2 ||
      any(!is.finite(bounds)) || bounds[1] >= bounds[2]) {
    stop("bounds must be a finite increasing numeric vector of length two.",
         call. = FALSE)
  }

  theta_hat <- rep(NA_real_, nrow(resp))

  for (i in seq_len(nrow(resp))) {
    observed <- which(!is.na(resp[i, ]))
    if (length(observed) == 0) {
      next
    }

    objective <- function(theta) {
      eta <- item_pars$d[observed] + item_pars$a[observed] * theta
      sum(safe_log1pexp(eta) - resp[i, observed] * eta)
    }

    theta_hat[i] <- stats::optimize(objective, interval = bounds)$minimum
  }

  theta_hat
}
