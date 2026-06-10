# ---------------------------------------------------------------------------- #
# Scenario set 3: downstream ability-score payoff and the no-harm property.
#
# The method's value proposition: lower ability-score RMSE when the predictor is
# informative, and no WORSE than human-only on average when it is not. We measure
# scoring RMSE on a held-out scoring sample with known theta, comparing three
# item calibrations:
#   (a) human-only      : MML at lambda = 0
#   (b) MML tuned        : MML at the ability-risk-tuned lambda
#   (c) MML full PPI     : MML at lambda = 1
#
# Expected behaviour (100 reps):
#   R1 (perfect, lambda~0.75)  : tuned RMSE slightly < human-only
#   R2, R4 (lambda~0.10)       : tuned RMSE marginally < human-only
#   R3 (lambda~0)              : tuned falls back to human-only (no harm)
# Effect sizes are small because an 8-item test is measurement-limited; the
# larger payoff is in item-parameter precision (see run_coverage.R). The summary
# reports the paired SE / 95% CI of mean_delta and prop_improve so small but
# stable effects can be told apart from simulation noise.
#
# We also record item-parameter bias at the tuned lambda to confirm the MML
# estimator is (approximately) unbiased at the selected lambda.
#
# Usage:
#   Rscript simulations/run_downstream.R [n_reps] [cores]
# Output:
#   simulations/results/downstream.rds
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
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
cores  <- parse_cores(args)

n_human     <- 400
n_generated <- 1200
n_score     <- 1000          # held-out scoring sample size
n_quad      <- 11
lambda_grid <- default_lambda_grid()
regimes     <- all_regimes()

rmse <- function(a, b) sqrt(mean((a - b)^2))

# Score a held-out sample with calibrated item parameters and return ability RMSE
scoring_rmse <- function(item_pars, score_resp, theta_true) {
  theta_hat <- score_theta(score_resp, item_pars, bounds = c(-6, 6))
  ok <- is.finite(theta_hat) & abs(theta_hat) < 6   # drop boundary scores
  rmse(theta_hat[ok], theta_true[ok])
}

run_one <- function(task) {
  set.seed(task$seed)
  regime <- task$regime
  rep_id <- task$rep
  dat <- generate_regime(regime, n_human = n_human, n_generated = n_generated)
  human_start <- fit_2pl(dat$observed, technical = list(NCYCLES = 500))$pars

  # Held-out scoring sample (human responses from the TRUE DGP)
  theta_score <- stats::rnorm(n_score)
  score_resp  <- simulate_2pl(theta_score, dat$true_pars)

  fit_at <- function(lam) {
    tryCatch(
      fit_mixed_subjects_mml(
        observed = dat$observed, predicted = dat$predicted,
        generated = dat$generated, lambda = lam,
        initial_pars = human_start, n_quad = n_quad,
        control = list(maxit = 300)
      ),
      error = function(e) NULL
    )
  }

  # (a) human-only lambda = 0
  fit0 <- fit_at(0)
  # (b) tuned lambda
  tuned <- tryCatch(
    tune_lambda_ability_risk(
      lambda_grid  = lambda_grid,
      observed     = dat$observed, predicted = dat$predicted,
      generated    = dat$generated, target_resp = dat$observed,
      initial_pars = human_start, fit_fn = fit_mixed_subjects_mml,
      n_quad = n_quad, control = list(maxit = 200)
    ),
    error = function(e) NULL
  )
  # (c) full PPI lambda = 1
  fit1 <- fit_at(1)

  if (is.null(fit0) || is.null(tuned) || is.null(fit1)) return(NULL)
  fit_tuned <- tuned$best_fit
  lam_tuned <- tuned$best_lambda

  # Item-parameter bias at tuned lambda
  bias_a <- mean(fit_tuned$item_pars$a - dat$true_pars$a)
  bias_d <- mean(fit_tuned$item_pars$d - dat$true_pars$d)

  data.frame(
    regime     = regime,
    rep        = rep_id,
    lambda_tuned = lam_tuned,
    rmse_human = scoring_rmse(fit0$item_pars,      score_resp, theta_score),
    rmse_tuned = scoring_rmse(fit_tuned$item_pars, score_resp, theta_score),
    rmse_full  = scoring_rmse(fit1$item_pars,      score_resp, theta_score),
    bias_a     = bias_a,
    bias_d     = bias_d,
    stringsAsFactors = FALSE
  )
}

message(sprintf("Scenario 3: downstream payoff, %d reps x %d regimes (%d cores)",
                n_reps, length(regimes), cores))

tasks   <- make_tasks(regimes, n_reps)
results <- collect_rows(par_map(tasks, run_one, cores = cores))

# Per-regime summary with Monte Carlo uncertainty on the mean delta. The
# downstream RMSE gains are small, so we report the paired standard error of the
# mean delta and a 95% CI (mean_delta is a paired difference across reps), plus
# prop_improve with a binomial 95% CI. This lets readers separate a stable small
# effect from simulation noise and prevents reading "no harm" as a per-rep
# guarantee.
summ <- do.call(rbind, lapply(regimes, function(rg) {
  sub   <- results[results$regime == rg, ]
  delta <- sub$rmse_tuned - sub$rmse_human   # tuned - human (negative = better)
  m     <- length(delta)
  se    <- stats::sd(delta) / sqrt(m)        # paired SE of the mean delta
  z     <- stats::qnorm(0.975)
  p_imp <- mean(delta < 0)
  p_se  <- sqrt(p_imp * (1 - p_imp) / m)     # binomial SE for prop_improve

  data.frame(
    regime        = rg,
    label         = regime_labels()[[rg]],
    mean_lambda   = round(mean(sub$lambda_tuned), 3),
    rmse_human    = round(mean(sub$rmse_human), 4),
    rmse_tuned    = round(mean(sub$rmse_tuned), 4),
    rmse_full     = round(mean(sub$rmse_full),  4),
    mean_delta    = round(mean(delta), 4),                 # tuned - human
    delta_se      = round(se, 4),                          # paired SE of mean
    delta_lo      = round(mean(delta) - z * se, 4),        # 95% CI lower
    delta_hi      = round(mean(delta) + z * se, 4),        # 95% CI upper
    prop_improve  = round(p_imp, 3),                       # fraction tuned helps
    prop_improve_se = round(p_se, 3),
    bias_a        = round(mean(sub$bias_a), 4),
    bias_d        = round(mean(sub$bias_d), 4),
    stringsAsFactors = FALSE
  )
}))

dir.create("simulations/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(reps = results, summary = summ,
       settings = list(n_reps = n_reps, cores = cores, n_human = n_human,
                       n_generated = n_generated, n_score = n_score,
                       n_quad = n_quad, lambda_grid = lambda_grid)),
  "simulations/results/downstream.rds"
)

print(summ, row.names = FALSE)
message("Interpretation: mean_delta < 0 means tuned MML improves scoring; ",
        "for poor predictors it should be ~0 (no harm).")
message("Saved: simulations/results/downstream.rds")
