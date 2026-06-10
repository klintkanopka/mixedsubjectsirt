# ---------------------------------------------------------------------------- #
# Scenario set 4: cross-fitted vs non-cross-fitted ability-risk tuning.
#
# Cross-fitting was a central design point of the package: it tunes lambda on
# held-out folds to avoid using a fold's own labels to tune the lambda applied
# to that fold's paired correction. This scenario reports whether cross-fitting
# changes the selected lambda, item-parameter bias, item-parameter SE coverage,
# and held-out ability-score RMSE relative to the non-cross-fitted tuner.
#
#   Non-crossfit : tune_lambda_ability_risk(fit_fn = fit_mixed_subjects_mml)
#   Crossfit     : tune_lambda_ability_risk_crossfit(
#                    fit_fn = fit_mixed_subjects_mml,        # per-fold tuning
#                    final_fit_fn = fit_mixed_subjects_mml)  # scalar-mean final
#
# We run the informative regimes (R1, R2, R4) where the tuner selects a positive
# lambda; R3 collapses to ~0 under both methods and adds little contrast.
#
# Usage:
#   Rscript simulations/run_crossfit.R [n_reps] [cores]
# Output:
#   simulations/results/crossfit.rds
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
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 50L
cores  <- parse_cores(args)

n_human     <- 400
n_generated <- 1200
n_score     <- 1000
n_quad      <- 11
lambda_grid <- default_lambda_grid()
regimes     <- c("perfect", "same_dgp", "llm_shift")

rmse <- function(a, b) sqrt(mean((a - b)^2))

scoring_rmse <- function(item_pars, score_resp, theta_true) {
  theta_hat <- score_theta(score_resp, item_pars, bounds = c(-6, 6))
  ok <- is.finite(theta_hat) & abs(theta_hat) < 6
  rmse(theta_hat[ok], theta_true[ok])
}

# Fraction of item parameters whose true value lies in the 95% Wald interval.
cover95 <- function(fit, true_pars) {
  Sig <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(Sig)) return(NA_real_)
  est <- c(fit$item_pars$a, fit$item_pars$d)
  tru <- c(true_pars$a, true_pars$d)
  se  <- sqrt(pmax(diag(Sig), 0))
  mean(abs(est - tru) <= stats::qnorm(0.975) * se)
}

run_one <- function(task) {
  set.seed(task$seed)
  regime <- task$regime
  dat <- generate_regime(regime, n_human = n_human, n_generated = n_generated)
  human_start <- fit_2pl(dat$observed, technical = list(NCYCLES = 500))$pars

  theta_score <- stats::rnorm(n_score)
  score_resp  <- simulate_2pl(theta_score, dat$true_pars)

  # --- Non-cross-fitted ---
  nocf <- tryCatch(
    tune_lambda_ability_risk(
      lambda_grid  = lambda_grid,
      observed     = dat$observed, predicted = dat$predicted,
      generated    = dat$generated, target_resp = dat$observed,
      initial_pars = human_start, fit_fn = fit_mixed_subjects_mml,
      n_quad = n_quad, control = list(maxit = 200)),
    error = function(e) NULL
  )

  # --- Cross-fitted (2 folds), MML per-fold tuning + scalar-mean MML final ---
  split_id <- rep(1:2, length.out = nrow(dat$observed))
  cf <- tryCatch(
    tune_lambda_ability_risk_crossfit(
      lambda_grid  = lambda_grid,
      observed     = dat$observed, predicted = dat$predicted,
      generated    = dat$generated, split_id = split_id,
      initial_pars = human_start,
      fit_fn       = fit_mixed_subjects_mml,
      final_fit_fn = fit_mixed_subjects_mml,
      n_quad = n_quad, control = list(maxit = 200)),
    error = function(e) NULL
  )

  if (is.null(nocf) || is.null(cf) || is.null(cf$final_fit)) return(NULL)

  fit_nocf <- nocf$best_fit
  fit_cf   <- cf$final_fit

  data.frame(
    regime      = regime,
    rep         = task$rep,
    lambda_nocf = nocf$best_lambda,
    lambda_cf   = cf$lambda_final,
    bias_a_nocf = mean(fit_nocf$item_pars$a - dat$true_pars$a),
    bias_a_cf   = mean(fit_cf$item_pars$a   - dat$true_pars$a),
    rmse_nocf   = scoring_rmse(fit_nocf$item_pars, score_resp, theta_score),
    rmse_cf     = scoring_rmse(fit_cf$item_pars,   score_resp, theta_score),
    cover_nocf  = cover95(fit_nocf, dat$true_pars),
    cover_cf    = cover95(fit_cf,   dat$true_pars),
    stringsAsFactors = FALSE
  )
}

message(sprintf("Scenario 4: crossfit comparison, %d reps x %d regimes (%d cores)",
                n_reps, length(regimes), cores))

tasks   <- make_tasks(regimes, n_reps)
results <- collect_rows(par_map(tasks, run_one, cores = cores))

summ <- do.call(rbind, lapply(regimes, function(rg) {
  sub <- results[results$regime == rg, ]
  data.frame(
    regime        = rg,
    label         = regime_labels()[[rg]],
    lambda_nocf   = round(mean(sub$lambda_nocf), 3),
    lambda_cf     = round(mean(sub$lambda_cf),   3),
    bias_a_nocf   = round(mean(sub$bias_a_nocf), 4),
    bias_a_cf     = round(mean(sub$bias_a_cf),   4),
    rmse_nocf     = round(mean(sub$rmse_nocf),   4),
    rmse_cf       = round(mean(sub$rmse_cf),     4),
    cover_nocf    = round(mean(sub$cover_nocf, na.rm = TRUE), 3),
    cover_cf      = round(mean(sub$cover_cf,   na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}))

dir.create("simulations/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(reps = results, summary = summ,
       settings = list(n_reps = n_reps, cores = cores, n_human = n_human,
                       n_generated = n_generated, n_score = n_score,
                       n_quad = n_quad, lambda_grid = lambda_grid,
                       regimes = regimes)),
  "simulations/results/crossfit.rds"
)

print(summ, row.names = FALSE)
message("Compare lambda, bias_a, RMSE, and 95% coverage between the ",
        "non-cross-fitted and cross-fitted tuners.")
message("Saved: simulations/results/crossfit.rds")
