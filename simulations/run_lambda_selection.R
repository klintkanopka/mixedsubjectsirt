# ---------------------------------------------------------------------------- #
# Scenario set 1: lambda-selection behaviour across predictor regimes.
#
# For each regime, fit with fit_mixed_subjects_mml over a lambda grid, tune by
# ability-score risk, and record the selected lambda. The validation claim is
# that lambda tracks the score-level usefulness of the paired pseudo-responses.
#
# Observed ordering (100 reps):
#   lambda(R1) = 0.75  (= N/(n+N), perfect paired predictor)
#   lambda(R2) ~ lambda(R4) ~ 0.10  (fresh real responses / shifted LLM)
#   lambda(R3) ~ 0                   (independent noise)
#
# We report the full selected-lambda distribution (not just the mean): R2 is the
# subtle regime (a fresh same-DGP draw) and deserves a histogram.
#
# Usage:
#   Rscript simulations/run_lambda_selection.R [n_reps] [cores]
# Output:
#   simulations/results/lambda_selection.rds
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
n_quad      <- 11
lambda_grid <- default_lambda_grid()
regimes     <- all_regimes()

run_one <- function(task) {
  set.seed(task$seed)
  regime <- task$regime
  rep_id <- task$rep
  dat <- generate_regime(regime, n_human = n_human, n_generated = n_generated)

  human_start <- fit_2pl(dat$observed, technical = list(NCYCLES = 500))$pars

  # PPI++ score lambda (theoretical diagnostic)
  ppi <- tryCatch(
    tune_lambda_ppi_score(dat$observed, dat$predicted, human_start,
                          n_generated = n_generated, n_quad = n_quad)$lambda,
    error = function(e) NA_real_
  )

  # Ability-risk lambda with the MML estimator (the recommended workflow)
  tuned <- tryCatch(
    tune_lambda_ability_risk(
      lambda_grid  = lambda_grid,
      observed     = dat$observed,
      predicted    = dat$predicted,
      generated    = dat$generated,
      target_resp  = dat$observed,
      initial_pars = human_start,
      fit_fn       = fit_mixed_subjects_mml,
      n_quad       = n_quad,
      control      = list(maxit = 200)
    ),
    error = function(e) NULL
  )

  data.frame(
    regime       = regime,
    rep          = rep_id,
    lambda_ppi   = ppi,
    lambda_risk  = if (is.null(tuned)) NA_real_ else tuned$best_lambda,
    stringsAsFactors = FALSE
  )
}

message(sprintf("Scenario 1: lambda selection, %d reps x %d regimes (%d cores)",
                n_reps, length(regimes), cores))

tasks   <- make_tasks(regimes, n_reps)
results <- collect_rows(par_map(tasks, run_one, cores = cores))

# Summary table: mean and quantiles of selected lambda per regime
summ <- do.call(rbind, lapply(regimes, function(rg) {
  sub <- results[results$regime == rg, ]
  data.frame(
    regime      = rg,
    label       = regime_labels()[[rg]],
    mean_risk   = round(mean(sub$lambda_risk, na.rm = TRUE), 3),
    median_risk = round(stats::median(sub$lambda_risk, na.rm = TRUE), 3),
    q10_risk    = round(stats::quantile(sub$lambda_risk, 0.10, na.rm = TRUE), 3),
    q90_risk    = round(stats::quantile(sub$lambda_risk, 0.90, na.rm = TRUE), 3),
    mean_ppi    = round(mean(sub$lambda_ppi, na.rm = TRUE), 3),
    prop_zero   = round(mean(sub$lambda_risk == 0, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}))

dir.create("simulations/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(reps = results, summary = summ,
       settings = list(n_reps = n_reps, cores = cores, n_human = n_human,
                       n_generated = n_generated, n_quad = n_quad,
                       lambda_grid = lambda_grid)),
  "simulations/results/lambda_selection.rds"
)

print(summ, row.names = FALSE)
message("Saved: simulations/results/lambda_selection.rds")
