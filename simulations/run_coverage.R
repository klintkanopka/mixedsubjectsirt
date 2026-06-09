# ---------------------------------------------------------------------------- #
# Scenario set 2: Louis-corrected SE coverage (the key validation of the Louis
# marginal-information bread implemented in Round 2).
#
# We compare empirical coverage of nominal 90% / 95% Wald intervals for the item
# parameters using two covariance estimates of the SAME scalar-MML fit:
#   * Louis-corrected marginal sandwich  : stats::vcov(fit)  (dispatch)
#   * EM complete-data Hessian bread     : vcov_mixed_subjects(fit) (bypass)
#
# Expected finding: the EM bread UNDER-covers (intervals too narrow, coverage
# below nominal) because it ignores missing information about theta; the Louis
# bread restores coverage to roughly nominal.
#
# Consistency note. Coverage of the TRUE DGP item parameters is only meaningful
# when the combined estimator is consistent for `true_pars`.  This holds for:
#   R1 perfect (F=Y)   : paired and generated terms both derive from true pars
#   R3 same_dgp        : paired and generated both drawn from true pars
# It does NOT hold for R5 llm_shift at lambda > 0 (generated comes from shifted
# parameters, so the combined estimator is biased).  We therefore restrict the
# coverage study to R1 and R3, at a fixed lambda, where the estimand is the true
# parameter vector.
#
# Usage:
#   Rscript simulations/run_coverage.R [n_reps] [cores]
# Output:
#   simulations/results/coverage.rds
# ---------------------------------------------------------------------------- #

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(mixedsubjectsirt)
  }
})
source("simulations/dgp.R")

args   <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 200L
cores  <- parse_cores(args)

n_human     <- 400
n_generated <- 1200
n_quad      <- 11
lambda_fix  <- 0.5                  # fixed lambda; estimand = true_pars
regimes     <- c("perfect", "same_dgp")
z90         <- stats::qnorm(0.95)   # two-sided 90%
z95         <- stats::qnorm(0.975)  # two-sided 95%

# Extract item-parameter point estimates in the (a_1..a_J, d_1..d_J) order
# matching the vcov dimnames.
param_vector <- function(fit) {
  c(fit$item_pars$a, fit$item_pars$d)
}
true_vector <- function(true_pars) {
  c(true_pars$a, true_pars$d)
}

run_one <- function(task) {
  set.seed(task$seed)
  regime <- task$regime
  dat <- generate_regime(regime, n_human = n_human, n_generated = n_generated)
  human_start <- fit_2pl(dat$observed, technical = list(NCYCLES = 500))$pars

  fit <- tryCatch(
    fit_mixed_subjects_mml(
      observed = dat$observed, predicted = dat$predicted,
      generated = dat$generated, lambda = lambda_fix,
      initial_pars = human_start, n_quad = n_quad,
      control = list(maxit = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || fit$convergence != 0) return(NULL)

  est  <- param_vector(fit)
  true <- true_vector(dat$true_pars)

  Sig_louis <- tryCatch(stats::vcov(fit),            error = function(e) NULL)
  Sig_em    <- tryCatch(vcov_mixed_subjects(fit),    error = function(e) NULL)
  if (is.null(Sig_louis) || is.null(Sig_em)) return(NULL)

  se_louis <- sqrt(pmax(diag(Sig_louis), 0))
  se_em    <- sqrt(pmax(diag(Sig_em),    0))

  in_ci <- function(se, z) as.integer(abs(est - true) <= z * se)

  # One row per parameter
  n_par <- length(est)
  data.frame(
    regime    = regime,
    par_idx   = seq_len(n_par),
    par_type  = rep(c("a", "d"), each = n_par / 2),
    cov_louis_90 = in_ci(se_louis, z90),
    cov_louis_95 = in_ci(se_louis, z95),
    cov_em_90    = in_ci(se_em,    z90),
    cov_em_95    = in_ci(se_em,    z95),
    se_louis     = se_louis,
    se_em        = se_em,
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Scenario 2: Louis SE coverage, %d reps x %d regimes (lambda=%.2f, %d cores)",
  n_reps, length(regimes), lambda_fix, cores))

tasks   <- make_tasks(regimes, n_reps)
results <- collect_rows(par_map(tasks, run_one, cores = cores))
n_par   <- max(results$par_idx)            # parameters per rep (= 2 * n_items)
for (rg in regimes) {
  message(sprintf("  %s: %d/%d reps usable",
                  rg, sum(results$regime == rg) / n_par, n_reps))
}

# Coverage summary per regime (averaged over reps and parameters)
summ <- do.call(rbind, lapply(regimes, function(rg) {
  sub <- results[results$regime == rg, ]
  data.frame(
    regime        = rg,
    label         = regime_labels()[[rg]],
    louis_cov_90  = round(mean(sub$cov_louis_90), 3),
    louis_cov_95  = round(mean(sub$cov_louis_95), 3),
    em_cov_90     = round(mean(sub$cov_em_90),    3),
    em_cov_95     = round(mean(sub$cov_em_95),    3),
    mean_se_ratio = round(mean(sub$se_louis / sub$se_em), 3),  # >1 means Louis wider
    stringsAsFactors = FALSE
  )
}))

dir.create("simulations/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(reps = results, summary = summ,
       settings = list(n_reps = n_reps, cores = cores, n_human = n_human,
                       n_generated = n_generated, n_quad = n_quad,
                       lambda_fix = lambda_fix, regimes = regimes)),
  "simulations/results/coverage.rds"
)

print(summ, row.names = FALSE)
message("Nominal targets: 0.90 and 0.95. EM bread should under-cover; ",
        "Louis should be close to nominal.")
message("Saved: simulations/results/coverage.rds")
