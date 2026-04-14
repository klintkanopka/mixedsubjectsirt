#' Converts an item response matrix into quadrature form and estimates a 2PL IRT model with parameters in slope-intercept form
#'
#' @param resp An item response matrix in wide form
#' @param N_quad The number of quadrature points to compute. Higher numbers can induce numerical errors
#' @param eps A tolerance for the minimum values
#' @param iterlim Maximum number of Newton-Raphson iterations passed to `rmutil::gauss.hermite()`
#' @param irt_pars IRT parameters from human calibration for rescaling parameters estimated from predicted responses
#' @return A list object with two components. `$quad` contains a dataframe with the expected sample sizes and number of correct responses at each quadrature point. `$irt_pars` contains the parameter estimates from the fitted IRT model that generated the expected counts
mixed_subjects_quadrature <- function(
  resp,
  N_quad = 10,
  eps = 1e-15,
  iterlim = 1e5,
  irt_pars = NULL
) {
  require(rmutil)
  require(dplyr)
  require(mirt)

  N_person <- nrow(resp)
  N_item <- ncol(resp)

  m <- mirt::mirt(
    resp,
    model = 1,
    itemtype = '2PL',
    technical = list(NCYCLES = 1e4)
  )

  pars <- data.frame(coef(m, simplify = TRUE)$items)

  if (!is.null(irt_pars)) {
    # rescale based on mean-mean equating using irt_pars
    irt_pars$b <- -1 * irt_pars$d / irt_pars$a1
    pars$b <- -1 * pars$d / pars$a1
    A <- mean(irt_pars$a1) / mean(pars$a1)
    B <- mean(pars$b) - A * mean(irt_pars$b)
    pars$a1 <- irt_pars$a1 / A
    pars$b <- A * irt_pars$b + B
    pars$d <- -1 * pars$a1 * pars$b
  }

  quad <- data.frame(rmutil::gauss.hermite(N_quad, iterlim)) |>
    select(X_k = Points, A_k = Weights)

  q_t_mat <- matrix(quad$X_k, byrow = FALSE, nrow = N_quad, ncol = N_item)
  q_a_mat <- matrix(pars$a1, byrow = TRUE, nrow = N_quad, ncol = N_item)
  q_b_mat <- matrix(pars$d, byrow = TRUE, nrow = N_quad, ncol = N_item)

  q_p_mat <- plogis(q_a_mat * q_t_mat + q_b_mat)

  q_p_mat[q_p_mat > (1 - eps)] <- 1 - eps
  q_p_mat[q_p_mat < eps] <- eps

  seqs <- resp |>
    dplyr::group_by(dplyr::pick(dplyr::everything())) |>
    dplyr::summarize(r = n(), .groups = 'drop')

  X <- seqs |>
    dplyr::select(-r) |>
    as.matrix()

  L <- exp(log(q_p_mat) %*% t(X) + log(1 - q_p_mat) %*% t(1 - X))

  P_tilde <- as.vector(t(L) %*% quad$A_k)
  N_k <- colSums(t(L * quad$A_k) * seqs$r / P_tilde)

  quad$N_k <- N_k

  for (i in 1:ncol(X)) {
    name <- paste0('p_', i)
    quad[[name]] <- colSums(t(L * quad$A_k) * X[, i] * seqs$r / P_tilde) /
      N_k
  }

  out <- list(quad = quad, irt_pars = pars)

  return(out)
}


#' Loss function for parameter estimation. Generally not called directly by users.
#'
#' @param pars A vector of item parameters in slope-intercept form. Passed with discriminations first, location second, in item order.
#' @param q_observed A quadrature object constructed from running mixed_subjects_quadrature on observed item responses
#' @param q_predicted A quadrature object constructed from running mixed_subjects_quadrature on predicted item responses for the human calibration sample
#' @param q_observed A quadrature object constructed from running mixed_subjects_quadrature on LLM-generated item responses
#' @param lambda The power tuning parameter. Takes values between 0 and 1. When lambda is zero, disregards predicted data completely
#' @return Loss for the mixed subjects estimator at the values of `pars`
mixed_subjects_loss <- function(
  pars,
  q_observed,
  q_predicted,
  q_llm,
  lambda = 0
) {
  # unpack parameter estimates

  n_items <- ncol(q_observed$quad) - 3
  a <- pars[1:n_items]
  d <- pars[(n_items + 1):(2 * n_items)]

  loss <- function(q, a, d, n_items) {
    q <- q$quad
    l <- 0

    for (i in 1:n_items) {
      l <- l +
        sum(
          q$N_k *
            (-1 *
              q[[paste0('p_', i)]] *
              (a[i] * q$X_k + d[i]) +
              log(1 + exp(a[i] * q$X_k + d[i])))
        ) /
          sum(q$N_k)
    }
    return(l)
  }
  total_loss <- loss(q_observed, a, d, n_items) +
    lambda * (loss(q_llm, a, d, n_items) - loss(q_predicted, a, d, n_items))

  return(total_loss)
}
