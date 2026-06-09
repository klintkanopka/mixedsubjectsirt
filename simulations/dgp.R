# ---------------------------------------------------------------------------- #
# Shared data-generating process and predictor regimes for simulation
# validation of the mixed-subjects MML estimator.
#
# This file is sourced by the run_*.R scripts.  It is NOT part of the installed
# package; it depends on the package being loadable via devtools::load_all() or
# library(mixedsubjectsirt).
# ---------------------------------------------------------------------------- #

# Default true 2PL item parameters used across all scenarios.
default_true_pars <- function(n_items = 8) {
  pars <- data.frame(
    item = paste0("Item", seq_len(n_items)),
    a    = seq(0.8, 1.6, length.out = n_items),
    d    = seq(-1.1, 1.1, length.out = n_items),
    stringsAsFactors = FALSE
  )
  pars$b <- -pars$d / pars$a
  pars
}

# Oracle conditional-mean ("probability") predictions: P(x=1 | theta; pars).
# Returns a fractional response matrix in [0, 1].
conditional_mean_pred <- function(theta, pars) {
  eta <- outer(theta, pars$a, `*`) +
    matrix(pars$d, nrow = length(theta), ncol = nrow(pars), byrow = TRUE)
  out <- stats::plogis(eta)
  colnames(out) <- pars$item
  out
}

# "Scrambled" item parameters: near-random discrimination, shuffled difficulty.
# Used for the independent-noise regime (R4).
scrambled_pars <- function(true_pars) {
  n_items <- nrow(true_pars)
  sc <- true_pars
  sc$a <- pmax(0.05, abs(stats::rnorm(n_items, 0, 0.1)))
  sc$d <- stats::rnorm(n_items, 0, 2)         # difficulty unrelated to true
  sc$b <- -sc$d / sc$a
  sc
}

# LLM-shift item parameters: ~10% attenuated discrimination, +0.25 logit shift.
# Used for the realistic LLM regime (R5).
llm_shift_pars <- function(true_pars) {
  n_items <- nrow(true_pars)
  llm <- true_pars
  llm$a <- pmax(0.4, 0.9 * true_pars$a + stats::rnorm(n_items, 0, 0.05))
  llm$d <- true_pars$d + 0.25 + stats::rnorm(n_items, 0, 0.15)
  llm$b <- -llm$d / llm$a
  llm
}

# ---------------------------------------------------------------------------- #
# Generate one replication's worth of data for a given predictor regime.
#
# Returns a list with:
#   observed   : n_human x n_items binary human responses
#   predicted   : paired predictions for the same rows (binary or fractional)
#   generated   : n_generated x n_items LLM responses
#   theta_human : true abilities for the human sample
#   true_pars   : the true item parameters
#   regime      : the regime label
#
# Regimes:
#   R1 "perfect"          predicted = observed                       (F = Y)
#   R2 "conditional_mean" predicted = P(x=1|theta; true)             (oracle prob)
#   R3 "same_dgp"         predicted = fresh draw from true DGP
#   R4 "independent"      predicted = draw from scrambled parameters
#   R5 "llm_shift"        predicted = draw from attenuated/shifted params
# ---------------------------------------------------------------------------- #
generate_regime <- function(regime,
                            n_human     = 400,
                            n_generated = 1200,
                            true_pars   = default_true_pars()) {
  theta_human <- stats::rnorm(n_human)
  observed    <- simulate_2pl(theta_human, true_pars)

  if (regime == "perfect") {
    predicted <- observed
    gen_pars  <- true_pars

  } else if (regime == "conditional_mean") {
    predicted <- conditional_mean_pred(theta_human, true_pars)
    gen_pars  <- true_pars

  } else if (regime == "same_dgp") {
    predicted <- simulate_2pl(theta_human, true_pars)   # fresh Bernoulli draw
    gen_pars  <- true_pars

  } else if (regime == "independent") {
    sc        <- scrambled_pars(true_pars)
    predicted <- simulate_2pl(theta_human, sc)
    gen_pars  <- sc

  } else if (regime == "llm_shift") {
    llm       <- llm_shift_pars(true_pars)
    predicted <- simulate_2pl(theta_human, llm)
    gen_pars  <- llm

  } else {
    stop("Unknown regime: ", regime, call. = FALSE)
  }

  # Generated LLM sample is drawn from the same parameters that produced the
  # paired predictions, from an independent ability sample.
  generated <- simulate_2pl(stats::rnorm(n_generated), gen_pars)
  colnames(generated) <- true_pars$item

  list(
    observed    = observed,
    predicted   = predicted,
    generated   = generated,
    theta_human = theta_human,
    true_pars   = true_pars,
    gen_pars    = gen_pars,
    regime      = regime
  )
}

# All regime labels, ordered from best to worst predictor quality.
all_regimes <- function() {
  c("perfect", "conditional_mean", "same_dgp", "independent", "llm_shift")
}

# Human-readable regime descriptions for tables/figures.
regime_labels <- function() {
  c(
    perfect          = "R1 perfect (F=Y)",
    conditional_mean = "R2 conditional mean",
    same_dgp         = "R3 same-DGP draw",
    independent      = "R4 independent noise",
    llm_shift        = "R5 LLM shift"
  )
}

# Default lambda grid used by the tuning scenarios.
default_lambda_grid <- function() seq(0, 1, by = 0.1)

# ---------------------------------------------------------------------------- #
# Parallel execution helpers.
#
# Each scenario builds a flat list of tasks (one per regime x rep), assigns a
# deterministic seed to each task, and maps run_one() over the tasks with
# par_map().  Because the seed lives on the task (not on the loop position),
# results are IDENTICAL whether run serially or in parallel, on any core count.
# ---------------------------------------------------------------------------- #

# Parse the number of cores from the command line. By convention the cores
# argument is the SECOND positional argument (the first is n_reps).
parse_cores <- function(args, pos = 2L) {
  if (length(args) >= pos) {
    max(1L, as.integer(args[pos]))
  } else {
    1L
  }
}

# Map FUN over X, forking across `cores` workers when cores > 1.  Forking
# (mclapply) is unavailable on Windows, so we fall back to lapply there and
# whenever cores == 1.  mc.preschedule = FALSE distributes tasks dynamically,
# which matters because per-task cost varies a lot across regimes.
par_map <- function(X, FUN, cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(X, FUN)
  }
}

# Build the flat task list for a scenario: one entry per (regime, rep) with a
# deterministic seed derived from a base seed, the regime index, and the rep.
make_tasks <- function(regimes, n_reps, base_seed = 20260605L) {
  tasks <- list()
  for (ri in seq_along(regimes)) {
    for (r in seq_len(n_reps)) {
      tasks[[length(tasks) + 1L]] <- list(
        regime = regimes[ri],
        rep    = r,
        # Distinct, deterministic seed per task; regime index is spaced out so
        # regimes never share a seed even for large n_reps.
        seed   = base_seed + ri * 1000000L + r
      )
    }
  }
  tasks
}

# Drop tasks that returned NULL or a parallel try-error, then row-bind the rest.
collect_rows <- function(results) {
  ok <- Filter(function(x) is.data.frame(x), results)
  if (length(ok) == 0L) {
    stop("All simulation tasks failed; nothing to summarise.", call. = FALSE)
  }
  do.call(rbind, ok)
}
