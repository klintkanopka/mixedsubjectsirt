# ---------------------------------------------------------------------------- #
# Scenario set 5: OPERATIONAL coverage of the Louis-corrected intervals.
#
# Scenario 2 (run_coverage.R) validates the covariance FORMULA at a FIXED lambda:
# vcov_mixed_subjects_mml() is the sandwich for the estimator gamma-hat(lambda) at
# a *known* lambda, and there it covers nominally. But the lambda you actually
# report is TUNED from the data, and vcov() treats lambda as fixed -- it does not
# propagate the uncertainty in lambda-hat. Using it at the tuned lambda is a
# POST-SELECTION inference problem and may under-cover.
#
# This scenario compares Wald-interval coverage of the true item parameters at
# three operating points, per replication:
#   (1) fixed lambda = 0.5                       -- the formula check (= Scenario 2)
#   (2) same-data tuned lambda-hat + vcov(best_fit)
#       (tune_lambda_ability_risk, optimize)     -- the naive workflow
#   (3) cross-fit tuned lambda-hat + vcov(final_fit)
#       (tune_lambda_ability_risk_crossfit, 2-fold, MML scalar-mean final)
#
# Expectation: (1) nominal; (2) mildly anti-conservative (post-selection); (3)
# back toward nominal, since cross-fitting tunes lambda on held-out folds. If (2)
# and (3) are both ~nominal, the post-selection penalty is negligible here (a flat
# risk surface where gamma-hat is insensitive to lambda) -- also worth knowing.
#
# Uses the same per-task seeds as run_coverage.R, so operating point (1) reproduces
# the Scenario 2 Louis coverage exactly.
#
# Usage:
#   Rscript simulations/run_coverage_tuned.R [n_reps] [cores]
# Output:
#   simulations/results/coverage_tuned.rds
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
lambda_fix  <- 0.5
regimes     <- all_regimes()
grid        <- default_lambda_grid()        # bounds the optimizer's search to [0, 1]
z90         <- stats::qnorm(0.95)
z95         <- stats::qnorm(0.975)

param_vector <- function(fit) c(fit$item_pars$a, fit$item_pars$d)
true_vector  <- function(true_pars) c(true_pars$a, true_pars$d)

# Per-parameter Wald coverage of `true` from a fitted object's Louis vcov.
# Returns NULL if the fit failed / did not converge / vcov is unavailable.
cover_of <- function(fit, true) {
  if (is.null(fit) || isTRUE(fit$convergence != 0)) return(NULL)
  Sig <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(Sig)) return(NULL)
  se  <- sqrt(pmax(diag(Sig), 0))
  est <- param_vector(fit)
  list(c90 = as.integer(abs(est - true) <= z90 * se),
       c95 = as.integer(abs(est - true) <= z95 * se),
       err = est - true)   # signed error -> item-parameter bias (no-free-lunch check)
}

run_one <- function(task) {
  set.seed(task$seed)
  regime <- task$regime
  dat <- generate_regime(regime, n_human = n_human, n_generated = n_generated)
  hs   <- fit_2pl(dat$observed, technical = list(NCYCLES = 500))$pars
  true <- true_vector(dat$true_pars)
  J2   <- length(true)
  na2  <- rep(NA_integer_, J2)
  na2r <- rep(NA_real_, J2)

  # (1) fixed lambda
  f_fix <- tryCatch(fit_mixed_subjects_mml(
    dat$observed, dat$predicted, dat$generated, lambda = lambda_fix,
    initial_pars = hs, n_quad = n_quad, control = list(maxit = 300)),
    error = function(e) NULL)
  cov_fix <- cover_of(f_fix, true)

  # (2) same-data tuned (optimize) -> vcov at the tuned best_fit
  sd <- tryCatch(tune_lambda_ability_risk(
    grid, dat$observed, dat$predicted, dat$generated, target_resp = dat$observed,
    initial_pars = hs, fit_fn = fit_mixed_subjects_mml, n_quad = n_quad,
    control = list(maxit = 300)), error = function(e) NULL)
  cov_sd <- if (is.null(sd)) NULL else cover_of(sd$best_fit, true)
  lam_sd <- if (is.null(sd)) NA_real_ else sd$best_lambda

  # (3) cross-fit tuned (2-fold) -> vcov at the scalar-mean MML final fit
  xf <- tryCatch(tune_lambda_ability_risk_crossfit(
    grid, dat$observed, dat$predicted, dat$generated, initial_pars = hs,
    fit_fn = fit_mixed_subjects_mml, final_fit_fn = fit_mixed_subjects_mml,
    n_splits = 2L, n_quad = n_quad, control = list(maxit = 300)),
    error = function(e) NULL)
  cov_xf <- if (is.null(xf) || is.null(xf$final_fit)) NULL else cover_of(xf$final_fit, true)
  lam_xf <- if (is.null(xf)) NA_real_ else xf$lambda_final

  data.frame(
    regime   = regime,
    par_idx  = seq_len(J2),
    par_type = rep(c("a", "d"), each = J2 / 2),
    fix_90 = if (is.null(cov_fix)) na2 else cov_fix$c90,
    fix_95 = if (is.null(cov_fix)) na2 else cov_fix$c95,
    sd_90  = if (is.null(cov_sd))  na2 else cov_sd$c90,
    sd_95  = if (is.null(cov_sd))  na2 else cov_sd$c95,
    xf_90  = if (is.null(cov_xf))  na2 else cov_xf$c90,
    xf_95  = if (is.null(cov_xf))  na2 else cov_xf$c95,
    err_fix = if (is.null(cov_fix)) na2r else cov_fix$err,
    err_sd  = if (is.null(cov_sd))  na2r else cov_sd$err,
    err_xf  = if (is.null(cov_xf))  na2r else cov_xf$err,
    lambda_sd = lam_sd,
    lambda_xf = lam_xf,
    stringsAsFactors = FALSE
  )
}

message(sprintf(
  "Scenario 5: operational coverage (fixed / same-data / cross-fit), %d reps x %d regimes (%d cores)",
  n_reps, length(regimes), cores))

tasks   <- make_tasks(regimes, n_reps)
results <- collect_rows(par_map(tasks, run_one, cores = cores))

mc <- function(x) round(mean(x, na.rm = TRUE), 3)
summ <- do.call(rbind, lapply(regimes, function(rg) {
  s <- results[results$regime == rg, ]
  # Item-parameter bias by parameter type (the no-free-lunch check): mean signed
  # error over all `a` (or `d`) parameters and replications, per operating point.
  bias <- function(col, ptype) round(mean(s[[col]][s$par_type == ptype], na.rm = TRUE), 4)
  data.frame(
    regime    = rg,
    label     = regime_labels()[[rg]],
    lambda_sd = mc(s$lambda_sd),
    lambda_xf = mc(s$lambda_xf),
    fixed_90  = mc(s$fix_90), fixed_95 = mc(s$fix_95),
    samedata_90 = mc(s$sd_90), samedata_95 = mc(s$sd_95),
    crossfit_90 = mc(s$xf_90), crossfit_95 = mc(s$xf_95),
    bias_a_fixed = bias("err_fix", "a"), bias_a_samedata = bias("err_sd", "a"),
    bias_a_crossfit = bias("err_xf", "a"),
    bias_d_fixed = bias("err_fix", "d"), bias_d_samedata = bias("err_sd", "d"),
    bias_d_crossfit = bias("err_xf", "d"),
    stringsAsFactors = FALSE
  )
}))

dir.create("simulations/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(reps = results, summary = summ,
       settings = list(n_reps = n_reps, cores = cores, n_human = n_human,
                       n_generated = n_generated, n_quad = n_quad,
                       lambda_fix = lambda_fix, regimes = regimes)),
  "simulations/results/coverage_tuned.rds"
)

print(summ, row.names = FALSE)
message("Nominal targets 0.90 / 0.95. fixed = formula check; same-data may ",
        "under-cover (post-selection); cross-fit should recover.")
