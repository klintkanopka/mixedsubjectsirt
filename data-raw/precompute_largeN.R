# ---------------------------------------------------------------------------- #
# Precompute results for the "weakly-informative LLM" vignette.
#
# Shows that the mixed-subjects (PPI) estimator is unbiased for the true human
# parameters at EVERY lambda, while lambda only moves EFFICIENCY -- and that
# naive pooling of human + LLM responses is biased. With a weakly-informative,
# biased (shifted) LLM, n = 500 human and N = 100000 unlabeled LLM responses.
#
#   Rscript data-raw/precompute_largeN.R [n_reps] [cores] [N]
#   -> vignettes/largeN_results.rds   (a few KB; build-ignored)
#
# Outputs:
#   $curve   : per-lambda bias(a), bias(d) (empirical, Monte Carlo) and risk
#              (the model-based ability-score risk the tuner minimizes), with SEs.
#   $naive/$human : item-parameter bias of the pooled-naive and human-only fits.
#   $opt     : lambda chosen by DIRECT 1-D optimization (optimize over [0,1]) on a
#              single dataset, plus the grid argmin for comparison.
# Each MML fit at N = 100000 takes several seconds, so this runs offline, not at
# vignette-knit time.
# ---------------------------------------------------------------------------- #

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(mixedsubjectsirt)
  }
})
source("simulations/dgp.R")   # par_map()

args   <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 16L
cores  <- if (length(args) >= 2) as.integer(args[2]) else 6L
N      <- if (length(args) >= 3) as.integer(args[3]) else 100000L
n      <- 500L
J      <- 8L

true_pars <- data.frame(
  item = paste0("Item", seq_len(J)),
  a    = seq(0.8, 1.6, length.out = J),
  d    = seq(-1.1, 1.1, length.out = J)
)
true_pars$b <- -true_pars$d / true_pars$a

llm   <- true_pars
llm$a <- 0.9 * true_pars$a          # ~10% attenuated discriminations
llm$d <- true_pars$d + 0.25         # +0.25 intercept shift
llm$b <- -llm$d / llm$a

lam_grid <- seq(0, 1, by = 0.1)     # fine grid -- only to SAMPLE the surface

gen_one <- function(seed) {
  set.seed(seed)
  theta     <- stats::rnorm(n)
  observed  <- simulate_2pl(theta, true_pars)
  predicted <- simulate_2pl(theta, llm)
  generated <- simulate_2pl(stats::rnorm(N), llm)
  list(observed = observed, predicted = predicted, generated = generated)
}

# Model-based ability-score risk of a fit (the quantity the tuner minimizes).
risk_of <- function(fit, target) {
  Sigma <- stats::vcov(fit)
  ability_risk(target, fit, vcov = Sigma)$summary$mean_param_var
}

fit_lambda <- function(d, lam, init) {
  fit_mixed_subjects_mml(d$observed, d$predicted, d$generated, lambda = lam,
                         initial_pars = init, n_quad = 11, control = list(maxit = 200))
}

one_rep <- function(seed) {
  d     <- gen_one(seed)
  human <- fit_2pl(d$observed, technical = list(NCYCLES = 300))$pars
  naive <- fit_2pl(rbind(d$observed, d$generated), technical = list(NCYCLES = 300))$pars

  curve <- do.call(rbind, lapply(lam_grid, function(lam) {
    f <- fit_lambda(d, lam, human)
    data.frame(lambda = lam,
               bias_a = mean(f$item_pars$a - true_pars$a),
               bias_d = mean(f$item_pars$d - true_pars$d),
               risk   = risk_of(f, d$observed))
  }))
  curve$rep <- seed
  # Per-dataset optimal lambda by DIRECT 1-D optimization (continuous) -- the same
  # quantity optimize() / tune_lambda_ability_risk() target on real data, NOT the
  # grid argmin. (The grid above is only used to draw the average risk surface.)
  opt_rep <- stats::optimize(
    function(lam) risk_of(fit_lambda(d, lam, human), d$observed),
    interval = c(0, 1), tol = 0.01)$minimum
  list(curve = curve, opt_rep = opt_rep,
       naive_a = mean(naive$a - true_pars$a), naive_d = mean(naive$d - true_pars$d),
       human_a = mean(human$a - true_pars$a), human_d = mean(human$d - true_pars$d))
}

message(sprintf("Curve MC: %d reps, N=%d, n=%d, %d cores, %d-point grid",
                n_reps, N, n, cores, length(lam_grid)))
res  <- par_map(as.list(seq_len(n_reps)), one_rep, cores = cores)
res  <- res[!vapply(res, is.null, logical(1))]
m    <- length(res)

allc <- do.call(rbind, lapply(res, `[[`, "curve"))
agg  <- function(col) tapply(allc[[col]], allc$lambda, mean)
se   <- function(col) tapply(allc[[col]], allc$lambda, function(x) stats::sd(x) / sqrt(length(x)))
curve <- data.frame(
  lambda    = sort(unique(allc$lambda)),
  bias_a    = as.numeric(agg("bias_a")), bias_a_se = as.numeric(se("bias_a")),
  bias_d    = as.numeric(agg("bias_d")), bias_d_se = as.numeric(se("bias_d")),
  risk      = as.numeric(agg("risk")),   risk_se   = as.numeric(se("risk"))
)

naive_bias <- c(a = mean(vapply(res, `[[`, numeric(1), "naive_a")),
                d = mean(vapply(res, `[[`, numeric(1), "naive_d")))
human_bias <- c(a = mean(vapply(res, `[[`, numeric(1), "human_a")),
                d = mean(vapply(res, `[[`, numeric(1), "human_d")))

# Per-dataset optimal lambda: each replication's own risk-minimizing lambda, found
# by DIRECT 1-D optimization (continuous). Its spread shows how the optimum varies
# around the population optimum (the minimum of the AVERAGED curve).
per_rep_opt <- vapply(res, `[[`, numeric(1), "opt_rep")
opt_lambda  <- per_rep_opt[1]                        # rep 1's continuous optimum
grid_argmin <- curve$lambda[which.min(curve$risk)]   # argmin of the averaged MC curve

saveRDS(
  list(curve = curve, naive_bias = naive_bias, human_bias = human_bias,
       opt_lambda = opt_lambda,
       grid_argmin = grid_argmin, per_rep_opt = per_rep_opt, lam_grid = lam_grid,
       n = n, N = N, n_reps = m, true_pars = true_pars, llm = llm),
  "vignettes/largeN_results.rds"
)

cat(sprintf("\nN=%d, n=%d, %d reps. Direct optimize() lambda=%.3f (grid argmin %.2f).\n",
            N, n, m, opt_lambda, grid_argmin))
cat(sprintf("per-dataset optimal lambda: mean=%.2f, sd=%.2f, range=[%.1f, %.1f]\n",
            mean(per_rep_opt), stats::sd(per_rep_opt),
            min(per_rep_opt), max(per_rep_opt)))
cat("  distribution:", paste(sort(per_rep_opt), collapse = " "), "\n")
cat("naive bias (a,d):", round(naive_bias, 3),
    " | human bias (a,d):", round(human_bias, 3), "\n")
print(round(curve, 4), row.names = FALSE)
