# Simulation validation for the mixed-subjects MML estimator

These standalone scripts validate the scalar 2PL marginal-MML (MML) PPI++
workflow. They are **not** part of the installed package and are excluded from
`R CMD check`; they depend on the package being loadable via
`devtools::load_all()` (run from the package root) or `library(mixedsubjectsirt)`.

## What is validated

Three scenario sets, each over five predictor regimes:

| Regime | `predicted` construction | Predictor quality |
|--------|--------------------------|-------------------|
| R1 perfect (F=Y) | `predicted = observed` | best possible |
| R2 conditional mean | `plogis(aθ + d)` oracle probabilities | high (in theory) |
| R3 same-DGP draw | fresh `simulate_2pl(θ, true)` | modest |
| R4 independent noise | draw from scrambled parameters | none |
| R5 LLM shift | attenuated `a`, shifted `d` | low |

1. **`run_lambda_selection.R` — λ selection.** For each regime, tune λ by
   ability-score risk with `fit_mixed_subjects_mml` and record the selected λ.
   Validation claim: the selected λ decreases with predictor quality, with R1
   near `N/(n+N)` and R4/R5 near 0.

2. **`run_coverage.R` — Louis SE coverage (the key check).** At a fixed λ where
   the estimator is consistent for the true parameters (R1 and R3), compute Wald
   intervals from `stats::vcov(fit)` (Louis-corrected) and from
   `vcov_mixed_subjects(fit)` (EM bread, bypass). Validation claim: the EM bread
   under-covers; the Louis correction restores coverage toward nominal. This
   directly tests the Louis marginal-information bread.

3. **`run_downstream.R` — downstream payoff + no-harm.** Score a held-out sample
   with known θ. Compare ability-score RMSE for human-only (λ=0), MML at tuned λ,
   and MML at λ=1. Validation claims: tuned MML beats human-only when the
   predictor is informative (R1/R2), and does no worse than human-only when it is
   not (R4/R5).

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
| `dgp.R` | True item parameters; the five predictor regimes; `generate_regime()` |
| `run_lambda_selection.R` | Scenario set 1 |
| `run_coverage.R` | Scenario set 2 (Louis SE coverage) |
| `run_downstream.R` | Scenario set 3 (downstream payoff) |
| `figures.R` | Render tables/figures from saved `.rds` |
| `results/` | Saved result objects (created on run) |
| `figures/` | Saved figures (created on run) |

## Notes and known subtleties

- **R2 (conditional mean) may select λ ≈ 0** despite being an "oracle" predictor.
  This reproduces the Phase 0 diagnostic finding: at the human-MLE evaluation
  point the gradient cross-covariance between the human score and a deterministic
  conditional-mean prediction can be near zero or slightly negative, so the
  control-variate benefit does not materialise in the person-level score
  formulation. The full λ-selection run quantifies how often this happens.

- **Coverage is restricted to R1 and R3** because only there is the combined
  estimator consistent for the true parameters at λ > 0. Under R5 the generated
  sample comes from shifted parameters, so the estimator is biased at λ > 0 and
  coverage of the true parameters is not the right target.

- **Reproducibility.** Each (regime, rep) task carries a deterministic seed
  derived from a fixed base seed (`20260605`), the regime index, and the rep
  number. Results are therefore identical across serial and parallel runs and
  across core counts. Because each rep has its own independent seed, the first
  `k` reps of a 200-rep run are identical to a `k`-rep run, so increasing the rep
  count only adds reps rather than changing existing ones.
