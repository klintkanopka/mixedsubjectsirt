# ---------------------------------------------------------------------------- #
# Render summary figures and tables from the saved simulation results.
#
# Reads:
#   simulations/results/lambda_selection.rds
#   simulations/results/coverage.rds
#   simulations/results/downstream.rds
# Writes:
#   simulations/figures/*.png  (if ggplot2 is available)
#   prints summary tables to the console
#
# Usage:
#   Rscript simulations/figures.R
# ---------------------------------------------------------------------------- #

source("simulations/dgp.R")

has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (has_ggplot) library(ggplot2)
dir.create("simulations/figures", showWarnings = FALSE, recursive = TRUE)

read_if <- function(path) if (file.exists(path)) readRDS(path) else NULL

lam <- read_if("simulations/results/lambda_selection.rds")
cov <- read_if("simulations/results/coverage.rds")
dwn <- read_if("simulations/results/downstream.rds")

regime_order <- all_regimes()
lab          <- regime_labels()

# ----- Scenario 1: lambda selection -----
if (!is.null(lam)) {
  message("\n=== Scenario 1: lambda selection ===")
  print(lam$summary, row.names = FALSE)

  if (has_ggplot) {
    df <- lam$reps
    df$label <- factor(lab[df$regime], levels = lab[regime_order])
    p1 <- ggplot(df, aes(label, lambda_risk)) +
      geom_boxplot(outlier.size = 0.6, fill = "#9ecae1") +
      labs(x = NULL, y = "Selected lambda (ability risk)",
           title = "Lambda selection by predictor regime",
           subtitle = "Ordering should decrease left to right") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave("simulations/figures/lambda_selection.png", p1,
           width = 7, height = 4.2, dpi = 150)
    message("  wrote simulations/figures/lambda_selection.png")
  }
}

# ----- Scenario 2: coverage -----
if (!is.null(cov)) {
  message("\n=== Scenario 2: Louis SE coverage ===")
  print(cov$summary, row.names = FALSE)

  if (has_ggplot) {
    s  <- cov$summary
    df <- rbind(
      data.frame(regime = s$regime, bread = "Louis",
                 cov90 = s$louis_cov_90, cov95 = s$louis_cov_95),
      data.frame(regime = s$regime, bread = "EM",
                 cov90 = s$em_cov_90,    cov95 = s$em_cov_95)
    )
    df$label <- factor(lab[df$regime], levels = lab[regime_order])
    long <- rbind(
      data.frame(label = df$label, bread = df$bread, nominal = "90%", cov = df$cov90),
      data.frame(label = df$label, bread = df$bread, nominal = "95%", cov = df$cov95)
    )
    p2 <- ggplot(long, aes(label, cov, fill = bread)) +
      geom_col(position = position_dodge()) +
      geom_hline(data = data.frame(nominal = c("90%", "95%"),
                                   y = c(0.90, 0.95)),
                 aes(yintercept = y), linetype = "dashed") +
      facet_wrap(~nominal) +
      scale_fill_manual(values = c(Louis = "#31a354", EM = "#de2d26")) +
      labs(x = NULL, y = "Empirical coverage", fill = "Bread",
           title = "Item-parameter CI coverage: Louis vs EM bread",
           subtitle = "Dashed = nominal; EM should under-cover") +
      theme_minimal(base_size = 11)
    ggsave("simulations/figures/coverage.png", p2,
           width = 7, height = 4.2, dpi = 150)
    message("  wrote simulations/figures/coverage.png")
  }
}

# ----- Scenario 3: downstream payoff -----
if (!is.null(dwn)) {
  message("\n=== Scenario 3: downstream ability-score payoff ===")
  print(dwn$summary, row.names = FALSE)

  if (has_ggplot) {
    df <- dwn$reps
    df$delta <- df$rmse_tuned - df$rmse_human
    df$label <- factor(lab[df$regime], levels = lab[regime_order])
    p3 <- ggplot(df, aes(label, delta)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_boxplot(outlier.size = 0.6, fill = "#a1d99b") +
      labs(x = NULL, y = "RMSE(tuned) - RMSE(human-only)",
           title = "Downstream scoring payoff by regime",
           subtitle = "< 0 = improvement; ~0 for poor predictors (no harm)") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave("simulations/figures/downstream.png", p3,
           width = 7, height = 4.2, dpi = 150)
    message("  wrote simulations/figures/downstream.png")
  }
}

if (!has_ggplot) {
  message("\n(ggplot2 not installed; printed tables only, no figures written.)")
}
