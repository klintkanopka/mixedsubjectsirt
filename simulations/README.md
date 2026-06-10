# Simulation validation for the mixed-subjects MML estimator

These standalone scripts validate the scalar 2PL marginal-MML (MML) PPI++
workflow. They are **not** part of the installed package and are excluded from
`R CMD check`; they depend on the package being loadable via
`devtools::load_all()` (run from the package root) or `library(mixedsubjectsirt)`.

## What is validated

Four scenario sets, each over the four predictor regimes:

| Regime | `predicted` construction | Predictor quality |
|--------|--------------------------|-------------------|
| R1 perfect (F=Y) | `predicted = observed` | best possible |
| R2 same-DGP draw | fresh `simulate_2pl(θ, true)` | modest |
| R3 independent noise | binary draw from scrambled parameters | none |
| R4 LLM shift | binary draw from attenuated/shifted `a`, `d` | biased but informative |

All predictors are **binary 0/1 responses**: the package disallows probability
(fractional) predictions for `predicted`/`generated`. Earlier rounds also included
"conditional mean" regimes that fed fractional predictions; those have been
removed, and the four binary regimes are numbered contiguously R1–R4.

1. **`run_lambda_selection.R` — λ selection.** For each regime, tune λ by
   ability-score risk with `fit_mixed_subjects_mml` and record the selected λ.
   Validation claim: λ tracks the score-level usefulness of the paired
   pseudo-responses — R1 near `N/(n+N)`, R3 near 0, R2/R4 small positive.

2. **`run_coverage.R` — Louis SE coverage (the key check).** At a fixed λ = 0.5,
   compute Wald intervals from `stats::vcov(fit)` (Louis-corrected) and from
   `vcov_mixed_subjects(fit)` (EM bread, bypass), for all four regimes. The EM
   bread under-covers; the Louis correction restores coverage toward nominal.
   Crucially, the estimator is **consistent for the true human parameters in
   every regime** — even when the LLM is useless (R3) or biased (R4) — because
   the PPI correction is mean-zero whenever the paired and generated
   pseudo-responses share the same distribution and ability spread (see the
   consistency note in the script). R3 and R4 covering nominally is the flagship
   demonstration that PPI corrects biased LLM outputs.

3. **`run_downstream.R` — downstream payoff + no-harm.** Score a held-out sample
   with known θ. Compare ability-score RMSE for human-only (λ=0), MML at tuned λ,
   and MML at λ=1. Validation claim: in these simulations tuned MML does not
   increase average held-out RMSE relative to human-only, and improves it where
   the predictor is informative. The summary reports `prop_improve` and the
   paired SE / 95% CI of `mean_delta` so small effects can be distinguished from
   simulation noise.

4. **`run_crossfit.R` — cross-fit vs non-cross-fit.** Compares
   `tune_lambda_ability_risk` against `tune_lambda_ability_risk_crossfit` (MML
   per-fold tuning, scalar-mean MML final fit) on selected λ, item-parameter
   bias, 95% Wald coverage, and held-out RMSE, for R1/R2/R4.

## Running

From the package root. Each `run_*.R` script accepts two positional command-line
arguments: `[n_reps] [cores]`.

```sh
# Quick smoke test (a few reps each — NOT inferential)
Rscript simulations/run_lambda_selection.R 5
Rscript simulations/run_coverage.R 5
Rscript simulations/run_downstream.R 5

# Full runs (recommended rep counts), parallel across reps on 8 cores
Rscript simulations/run_lambda_selection.R 100 8
Rscript simulations/run_coverage.R 200 8
Rscript simulations/run_downstream.R 100 8
Rscript simulations/run_crossfit.R 50 8

# Render summary tables + figures from saved results
Rscript simulations/figures.R
```

Results are written to `simulations/results/*.rds`; figures (if `ggplot2` is
installed) to `simulations/figures/*.png`.

### Parallelism

The second argument is the number of CPU cores. With `cores > 1`, reps are
distributed across forked workers via `parallel::mclapply` (Unix/macOS only;
Windows falls back to serial). Speed-up is near-linear in cores because the reps
are independent.

**Reproducibility is independent of the core count.** Each (regime, rep) task
carries a deterministic seed derived from a fixed base seed, the regime index,
and the rep number, and `run_one()` sets that seed before generating data. A run
with `cores = 1` and a run with `cores = 8` produce identical results; only wall
time differs.

The downstream scenario is the most expensive (each rep does an 11-point MML tune
plus scoring), so it benefits most from cores. Serial single-core estimates are
roughly: coverage ~1 h, λ-selection ~15 h, downstream ~24 h; on 8 cores the suite
is a few hours.

> **Smoke tests are not inferential.** Coverage in particular requires ≥200 reps
> to stabilise; with a handful of reps the absolute numbers are pure noise.
> Useful signals visible even at low reps: the Louis/EM SE ratio (> 1) and the
> ordering of selected λ across regimes.

## Files

| File | Purpose |
|------|---------|
| `dgp.R` | True item parameters; the four predictor regimes; `generate_regime()` |
| `run_lambda_selection.R` | Scenario set 1 (λ selection) |
| `run_coverage.R` | Scenario set 2 (Louis SE coverage, all 4 regimes) |
| `run_downstream.R` | Scenario set 3 (downstream payoff + MC uncertainty) |
| `run_crossfit.R` | Scenario set 4 (cross-fit vs non-cross-fit) |
| `figures.R` | Render tables/figures from saved `.rds` |
| `results/` | Saved result objects (created on run) |
| `figures/` | Saved figures (created on run) |

## Notes and known subtleties

- **Why probability predictions are disallowed.** Earlier rounds included a
  "conditional mean" regime whose `predicted` was the *fractional* probability
  `plogis(aθ + d)`. This exposed a genuine flaw, now fixed at the API level: a
  fractional value is not a coherent likelihood term for the marginal IRT
  objective (the response enters inside a log-sum over quadrature). Mixing a
  fractional paired stream with a binary generated stream makes
  `E[∇L_gen] ≠ E[∇L_pred]`, so the PPI correction no longer cancels and the
  estimator is biased at λ > 0; worse, at λ = 0.5 the objective is *unbounded* in
  discrimination (`a → ∞`), so the fit diverges and ~92% of reps were dropped. The
  package now requires **binary** `predicted`/`generated` (sample from any
  probabilities first), which makes this failure impossible by construction. Those
  conditional-mean regimes were removed accordingly.

- **Coverage consistency.** The PPI correction `λ(L_gen − L_pred)` is mean-zero
  whenever the paired and generated pseudo-responses are drawn from the **same**
  distribution with the **same** ability spread, so the human term `L_obs`
  anchors the estimand to the *true* parameters. This holds for every (binary)
  regime — R1, R2, R3, R4 — so the estimator is consistent for the truth even
  when the LLM is useless (R3) or biased (R4). R3 and R4 covering nominally is the
  flagship demonstration that PPI corrects biased LLM outputs.

- **Reproducibility.** Each (regime, rep) task carries a deterministic seed
  derived from a fixed base seed (`20260605`), the regime index, and the rep
  number. Results are therefore identical across serial and parallel runs and
  across core counts. Because each rep has its own independent seed, the first
  `k` reps of a 200-rep run are identical to a `k`-rep run, so increasing the rep
  count only adds reps rather than changing existing ones.
