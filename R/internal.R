safe_log1pexp <- function(x) {
  out <- numeric(length(x))
  positive <- x > 0
  out[positive] <- x[positive] + log1p(exp(-x[positive]))
  out[!positive] <- log1p(exp(x[!positive]))
  dim(out) <- dim(x)
  out
}

logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

validate_response_matrix <- function(x, name = "resp", allow_fractional = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop(name, " must be a matrix or data frame.", call. = FALSE)
  }
  if (!is.numeric(x) && !is.integer(x) && !is.logical(x)) {
    stop(name, " must contain numeric, integer, or logical responses.", call. = FALSE)
  }
  if (length(dim(x)) != 2 || nrow(x) == 0 || ncol(x) == 0) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }

  x <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x),
              dimnames = dimnames(x))

  values <- x[!is.na(x)]
  if (length(values) > 0 && any(!is.finite(values))) {
    stop(name, " contains non-finite response values.", call. = FALSE)
  }
  if (length(values) > 0 && (any(values < -1e-12) || any(values > 1 + 1e-12))) {
    stop(name, " must contain values in [0, 1] or NA.", call. = FALSE)
  }
  if (!allow_fractional && length(values) > 0) {
    is_binary <- abs(values - round(values)) < 1e-12
    if (!all(is_binary)) {
      stop(name, " must contain binary 0/1 responses when fitting a 2PL model.",
           call. = FALSE)
    }
  }

  x[x < 0 & !is.na(x)] <- 0
  x[x > 1 & !is.na(x)] <- 1

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("Item", seq_len(ncol(x)))
  }

  x
}

check_same_items <- function(x, y, x_name = "x", y_name = "y") {
  if (ncol(x) != ncol(y)) {
    stop(x_name, " and ", y_name, " must have the same number of items.",
         call. = FALSE)
  }

  x_names <- colnames(x)
  y_names <- colnames(y)
  if (!is.null(x_names) && !is.null(y_names) && !identical(x_names, y_names)) {
    stop(x_name, " and ", y_name, " must have matching item column names.",
         call. = FALSE)
  }

  invisible(TRUE)
}

standardize_item_pars <- function(item_pars, n_items = NULL, item_names = NULL,
                                  name = "item_pars") {
  if (inherits(item_pars, "SingleGroupClass")) {
    item_pars <- as.data.frame(mirt::coef(item_pars, simplify = TRUE)$items)
  } else if (is.matrix(item_pars) || is.data.frame(item_pars)) {
    item_pars <- as.data.frame(item_pars)
  } else {
    stop(name, " must be a matrix, data frame, or mirt model.", call. = FALSE)
  }

  nm <- names(item_pars)
  slope_col <- match(TRUE, nm %in% c("a", "a1", "slope", "discrimination"))
  intercept_col <- match(TRUE, nm %in% c("d", "intercept", "gamma0"))

  if (is.na(slope_col) || is.na(intercept_col)) {
    if (ncol(item_pars) >= 2) {
      slope_col <- 1
      intercept_col <- 2
    } else {
      stop(name, " must contain slope and intercept columns.", call. = FALSE)
    }
  }

  item_values <- if ("item" %in% nm) {
    as.character(item_pars[["item"]])
  } else {
    rownames(item_pars)
  }

  out <- data.frame(
    item = item_values,
    a = as.numeric(item_pars[[slope_col]]),
    d = as.numeric(item_pars[[intercept_col]]),
    stringsAsFactors = FALSE
  )

  if (is.null(out$item) || any(out$item == "")) {
    out$item <- paste0("Item", seq_len(nrow(out)))
  }
  if (!is.null(item_names)) {
    out$item <- item_names
  }

  if (!is.null(n_items) && nrow(out) != n_items) {
    stop(name, " must contain one row per item.", call. = FALSE)
  }
  if (any(!is.finite(out$a)) || any(!is.finite(out$d))) {
    stop(name, " contains non-finite item parameters.", call. = FALSE)
  }
  if (any(abs(out$a) < .Machine$double.eps)) {
    stop(name, " contains zero discrimination parameters.", call. = FALSE)
  }

  out$b <- -out$d / out$a
  out
}

vector_from_item_pars <- function(item_pars) {
  item_pars <- standardize_item_pars(item_pars)
  c(item_pars$a, item_pars$d)
}

item_pars_from_vector <- function(pars, item_names) {
  n_items <- length(item_names)
  if (length(pars) != 2 * n_items) {
    stop("pars must contain all discriminations followed by all intercepts.",
         call. = FALSE)
  }

  out <- data.frame(
    item = item_names,
    a = as.numeric(pars[seq_len(n_items)]),
    d = as.numeric(pars[n_items + seq_len(n_items)]),
    stringsAsFactors = FALSE
  )
  out$b <- -out$d / out$a
  out
}

check_quadrature <- function(quadrature = NULL, n_quad = 31, iterlim = 1e5) {
  if (is.null(quadrature)) {
    return(make_quadrature(n_quad = n_quad, iterlim = iterlim))
  }

  if (!is.data.frame(quadrature)) {
    quadrature <- as.data.frame(quadrature)
  }

  theta_col <- if ("theta" %in% names(quadrature)) "theta" else "X_k"
  weight_col <- if ("weight" %in% names(quadrature)) "weight" else "A_k"
  if (!theta_col %in% names(quadrature) || !weight_col %in% names(quadrature)) {
    stop("quadrature must contain theta/weight or X_k/A_k columns.",
         call. = FALSE)
  }

  theta <- as.numeric(quadrature[[theta_col]])
  weight <- as.numeric(quadrature[[weight_col]])
  if (length(theta) == 0 || length(theta) != length(weight)) {
    stop("quadrature nodes and weights must have the same non-zero length.",
         call. = FALSE)
  }
  if (any(!is.finite(theta)) || any(!is.finite(weight)) || any(weight < 0)) {
    stop("quadrature nodes and weights must be finite, with non-negative weights.",
         call. = FALSE)
  }

  o <- order(theta)
  weight <- weight[o] / sum(weight[o])
  theta <- theta[o]
  data.frame(
    node = seq_along(theta),
    theta = theta,
    weight = weight,
    X_k = theta,
    A_k = weight
  )
}

counts_to_quad <- function(counts) {
  quad <- data.frame(
    X_k = counts$theta,
    A_k = counts$weight,
    N_k = counts$node_count
  )

  for (j in seq_len(counts$n_items)) {
    p <- counts$R[j, ] / counts$N[j, ]
    p[!is.finite(p)] <- NA_real_
    quad[[paste0("p_", j)]] <- p
  }

  quad
}

build_quadrature_summary <- function(resp, item_pars, quadrature, weights = NULL) {
  if (is.null(weights)) {
    weights <- posterior_weights_2pl(resp, item_pars, quadrature = quadrature)
  }

  counts <- summarize_expected_counts(resp, weights)
  list(
    quad = counts_to_quad(counts),
    counts = counts,
    weights = weights,
    irt_pars = standardize_item_pars(item_pars, n_items = ncol(resp),
                                     item_names = colnames(resp)),
    quadrature = quadrature,
    theta = quadrature$theta
  )
}

loss_expected_counts <- function(counts, item_pars) {
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = counts$n_items,
    item_names = counts$item_names
  )

  eta <- outer(item_pars$a, counts$theta, `*`) +
    matrix(item_pars$d, nrow = counts$n_items, ncol = counts$n_nodes)
  sum(counts$N * safe_log1pexp(eta) - counts$R * eta) / counts$n
}

gradient_expected_counts <- function(counts, item_pars) {
  item_pars <- standardize_item_pars(
    item_pars,
    n_items = counts$n_items,
    item_names = counts$item_names
  )

  eta <- outer(item_pars$a, counts$theta, `*`) +
    matrix(item_pars$d, nrow = counts$n_items, ncol = counts$n_nodes)
  resid <- counts$N * stats::plogis(eta) - counts$R

  grad_a <- rowSums(resid * matrix(counts$theta,
                                   nrow = counts$n_items,
                                   ncol = counts$n_nodes,
                                   byrow = TRUE)) / counts$n
  grad_d <- rowSums(resid) / counts$n
  c(grad_a, grad_d)
}

q_to_counts <- function(q, name) {
  if (is.list(q) && !is.null(q$counts)) {
    return(q$counts)
  }
  if (is.list(q) && all(c("N", "R", "theta", "n") %in% names(q))) {
    return(q)
  }
  stop(name, " must be a quadrature summary returned by mixed_subjects_quadrature().",
       call. = FALSE)
}

check_counts_compatible <- function(counts) {
  theta <- counts[[1]]$theta
  item_names <- counts[[1]]$item_names

  for (i in seq_along(counts)) {
    if (!isTRUE(all.equal(theta, counts[[i]]$theta, tolerance = 1e-10)) ||
        !identical(item_names, counts[[i]]$item_names)) {
      stop("Quadrature summaries must use the same items and quadrature nodes.",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}

combine_counts <- function(counts, mode = c("sum", "mean")) {
  mode <- match.arg(mode)
  check_counts_compatible(counts)

  divisor <- if (mode == "mean") length(counts) else 1
  n <- if (mode == "sum") sum(vapply(counts, `[[`, numeric(1), "n")) else counts[[1]]$n

  out <- counts[[1]]
  out$N <- Reduce(`+`, lapply(counts, `[[`, "N")) / divisor
  out$R <- Reduce(`+`, lapply(counts, `[[`, "R")) / divisor
  out$node_count <- Reduce(`+`, lapply(counts, `[[`, "node_count")) / divisor
  out$n <- n
  out
}

average_item_pars <- function(pars_list) {
  item_names <- pars_list[[1]]$item
  for (pars in pars_list) {
    if (!identical(item_names, pars$item)) {
      stop("Cannot average item parameters with different item names.",
           call. = FALSE)
    }
  }

  a <- rowMeans(do.call(cbind, lapply(pars_list, `[[`, "a")))
  d <- rowMeans(do.call(cbind, lapply(pars_list, `[[`, "d")))
  out <- data.frame(item = item_names, a = a, d = d, stringsAsFactors = FALSE)
  out$b <- -out$d / out$a
  out
}

validate_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda)) {
    stop("lambda must be a single finite number.", call. = FALSE)
  }
  if (lambda < 0 || lambda > 1) {
    stop("lambda must be between 0 and 1.", call. = FALSE)
  }
  lambda
}

fit_from_counts <- function(counts_observed, counts_predicted, counts_generated,
                            initial_pars, lambda, slope_lower = 1e-4,
                            control = list(maxit = 500)) {
  check_counts_compatible(list(counts_observed, counts_predicted, counts_generated))
  lambda <- validate_lambda(lambda)

  item_names <- counts_observed$item_names
  initial_pars <- standardize_item_pars(
    initial_pars,
    n_items = length(item_names),
    item_names = item_names
  )

  objective <- function(par) {
    item_pars <- item_pars_from_vector(par, item_names)
    loss_expected_counts(counts_observed, item_pars) +
      lambda * (
        loss_expected_counts(counts_generated, item_pars) -
          loss_expected_counts(counts_predicted, item_pars)
      )
  }

  gradient <- function(par) {
    item_pars <- item_pars_from_vector(par, item_names)
    gradient_expected_counts(counts_observed, item_pars) +
      lambda * (
        gradient_expected_counts(counts_generated, item_pars) -
          gradient_expected_counts(counts_predicted, item_pars)
      )
  }

  start <- vector_from_item_pars(initial_pars)
  if (is.null(slope_lower)) {
    lower <- rep(-Inf, length(start))
  } else {
    lower <- c(rep(slope_lower, length(item_names)), rep(-Inf, length(item_names)))
    start[seq_along(item_names)] <- pmax(start[seq_along(item_names)], slope_lower)
  }

  control <- utils::modifyList(list(maxit = 500), control)
  opt <- stats::optim(
    par = start,
    fn = objective,
    gr = gradient,
    method = "L-BFGS-B",
    lower = lower,
    upper = rep(Inf, length(start)),
    control = control
  )

  list(
    item_pars = item_pars_from_vector(opt$par, item_names),
    par = opt$par,
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    optimizer = opt
  )
}

make_split_id <- function(n, n_splits, seed = NULL) {
  if (n_splits < 2 || n_splits > n) {
    stop("n_splits must be at least 2 and no larger than nrow(observed).",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (old_seed_exists) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  sample(rep(seq_len(n_splits), length.out = n))
}
