# Simulation Validation of the Mixed-Subjects MML Estimator

This vignette presents a simulation study validating the scalar 2PL
marginal maximum-likelihood (MML) PPI++ estimator
(`fit_mixed_subjects_mml`) and its ability-risk tuning
(`tune_lambda_ability_risk`). The full reproducible harness lives in the
package’s `simulations/` directory; the tables below report results from
runs with the rep counts noted in each section.

The study answers three questions:

1.  **Does λ-selection track predictor quality?** (Scenario 1)
2.  **Are the Louis-corrected standard errors honest?** (Scenario 2 —
    the key check)
3.  **Does the method improve downstream scoring without ever harming
    it?** (Scenario 3)

## Design

Human responses come from a true 8-item 2PL model (`a ∈ [0.8, 1.6]`,
`d ∈ [-1.1, 1.1]`) with `n = 400` human subjects and `N = 1200`
generated subjects, abilities drawn from a standard normal. Five
predictor regimes vary the quality of the paired LLM predictions `F`:

| Regime | Predicted F | Predictor quality |
|:---|:---|:---|
| R1 perfect (F=Y) | predicted = observed (exact) | perfect |
| R2 conditional mean | true probabilities p(theta) (fractional) | oracle (but see below) |
| R3 same-DGP draw | fresh draw from the true model | modest |
| R4 independent noise | draw from scrambled item parameters | none |
| R5 LLM shift | draw from attenuated/shifted parameters | low |

Five predictor regimes, ordered best to worst. {.table}

The control parameter λ ∈ \[0, 1\] sets how strongly the generated LLM
sample is used after subtracting the paired LLM correction; λ = 0 is the
human-only calibration.

## Scenario 1: λ selection (100 reps)

For each regime we tune λ by ability-score risk and record the selected
value. We also report the theoretical PPI++-score λ for reference.

| Regime               | mean lambda | median | prop lambda=0 | PPI score lambda |
|:---------------------|------------:|-------:|--------------:|-----------------:|
| R1 perfect (F=Y)     |       0.750 |   0.75 |          0.00 |            0.750 |
| R2 conditional mean  |       0.028 |   0.00 |          0.96 |            0.016 |
| R3 same-DGP draw     |       0.108 |   0.10 |          0.26 |            0.004 |
| R4 independent noise |       0.032 |   0.00 |          0.74 |            0.000 |
| R5 LLM shift         |       0.102 |   0.10 |          0.31 |            0.002 |

Selected lambda by regime (ability-risk tuning, 100 reps). {.table}

**Findings.**

- **R1 perfect → λ = 0.750**, exactly `N / (n + N) = 1200 / 1600`, the
  theoretical maximum for a perfect paired predictor. The PPI++-score λ
  agrees to the decimal.
- **R4 independent noise → λ ≈ 0** (74% of reps select 0): a useless
  predictor is correctly rejected.
- **R3 same-DGP and R5 LLM shift → λ ≈ 0.10**: fresh real responses (R3)
  and a correlated-but-biased LLM (R5) each carry a little signal, and
  the criterion uses a little of it.
- **R2 conditional mean → λ ≈ 0** (96% of reps select 0). This is the
  one counterintuitive result and is discussed next.

### Why the “oracle” predictor (R2) is declined

R2 feeds the method the exact true probability of a correct response for
each person — seemingly the best possible predictor. Yet λ is almost
always 0. The reason: PPI++ extracts information from the **residuals**
between actual answers and the model’s expectation. A real 0/1 response
carries a residual; the conditional-mean prediction *is* the
expectation, so its residual is zero and it carries nothing to correct
with. A *noisier* predictor that emits actual responses (R3) is
therefore more useful than a *perfect* one that emits smooth
probabilities (R2). The practical lesson for applications: to benefit
from synthetic/LLM data, have it produce sampled responses, not averaged
probabilities.

## Scenario 2: Louis-corrected SE coverage (200 reps) — the key check

We fix λ = 0.5 in the two regimes where the combined estimator is
consistent for the true parameters (R1 and R3), and compare empirical
coverage of Wald intervals built from two covariance estimates of the
*same* fit: the Louis-corrected marginal sandwich
([`vcov()`](https://rdrr.io/r/stats/vcov.html) dispatch) and the EM
complete-data Hessian bread
([`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md),
a deliberate bypass).

| Regime           | Louis 90% | Louis 95% | EM 90% | EM 95% | SE ratio (Louis/EM) |
|:-----------------|----------:|----------:|-------:|-------:|--------------------:|
| R1 perfect (F=Y) |     0.909 |     0.955 |  0.713 |  0.787 |                1.63 |
| R3 same-DGP draw |     0.916 |     0.957 |  0.727 |  0.797 |                1.62 |

Item-parameter CI coverage (200 reps). Nominal targets 0.90 and 0.95.
{.table}

**The Louis correction restores nominal coverage; the EM bread
under-covers.** Louis intervals attain 0.909/0.955 and 0.916/0.957
against nominal 0.90/0.95, while the EM bread covers only ~0.71/0.79 —
it understates uncertainty because it ignores the missing information
about each subject’s latent ability. The Louis standard errors are ~1.6×
wider, which is the difference between an honest interval and an
overconfident one.

Coverage broken down by parameter type confirms the correction works
uniformly:

| Parameter        | Louis 90% | Louis 95% | EM 90% |
|:-----------------|----------:|----------:|-------:|
| discrimination a |      0.91 |      0.96 |   0.60 |
| intercept d      |      0.91 |      0.95 |   0.84 |

Coverage by parameter type (averaged over R1 and R3). {.table}

### A bug this scenario caught

An earlier version of this scenario showed near-zero coverage for
discrimination while intercept coverage was nominal. That pattern — a
uniform downward bias in all slopes with unbiased intercepts — is the
signature of a latent-scale error, and it traced to
[`make_quadrature()`](http://klintkanopka.com/mixedsubjectsirt/reference/make_quadrature.md)
building an N(0, 2) ability grid instead of N(0, 1) (a spurious √2
rescaling of standard-normal Gauss–Hermite nodes). Every discrimination
in the package had been biased downward by ~1/√2 ≈ 0.71. The bug was
invisible to relative checks (λ-selection), to gradient
finite-difference checks (loss and gradient were mutually consistent),
and to the standard-error calibration itself (the SEs correctly
described the *biased* estimator). Only coverage against ground truth
exposed it. After the fix, discrimination coverage rose from 0.04 to
0.91. This is the strongest argument for ground-truth simulation
validation: it catches errors that internal consistency checks cannot.

## Scenario 3: downstream payoff and no-harm (100 reps)

We score a held-out sample of 1000 subjects with known ability and
compare ability-score RMSE for the human-only calibration (λ = 0) and
the tuned MML calibration.

| Regime | mean lambda | RMSE human | RMSE tuned | delta (tuned-human) | bias a |
|:---|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.750 | 1.4961 | 1.4921 | -0.0040 | 0.011 |
| R2 conditional mean | 0.000 | 1.4958 | 1.4958 | 0.0000 | 0.019 |
| R3 same-DGP draw | 0.108 | 1.5020 | 1.5010 | -0.0010 | 0.038 |
| R4 independent noise | 0.032 | 1.5064 | 1.5060 | -0.0004 | 0.005 |
| R5 LLM shift | 0.102 | 1.4994 | 1.4981 | -0.0013 | 0.031 |

Downstream ability-score RMSE (100 reps). {.table}

**Findings.**

- **No harm anywhere.** The tuned RMSE is never worse than human-only
  (delta ≤ 0 in every regime). When the predictor is useless (R2, R4)
  the method falls back to the human-only estimate; when it helps (R1,
  R3, R5) it improves, modestly.
- **Discrimination is unbiased** at the selected λ (\|bias_a\| ≤ 0.04
  everywhere), confirming the quadrature fix and that the tuner is not
  selecting degenerate fits.
- **The RMSE gains are small** because an 8-item test is
  measurement-limited: ability scoring is dominated by the irreducible
  measurement error of a short test, so even a perfect predictor (R1)
  moves the scoring RMSE only from 1.4961 to 1.4921. The method’s larger
  payoff is in **item-parameter precision** — visible in the Scenario 2
  standard errors — which matters most when calibration uncertainty
  propagates into operational use.

### A robustness guard

Building this scenario surfaced a failure mode: with the fractional R2
predictor at large λ, the MML objective has a “converged” minimum at
runaway discrimination (a ≈ 1000). Such a fit reports huge item
information, which collapses its covariance and makes the model-based
ability risk spuriously small — so the tuner could select it even though
it scores worse on held-out data. The ability-risk tuners now reject any
candidate whose maximum discrimination exceeds `max_discrimination`
(default 10); a 2PL discrimination above that is never a real
calibration. With the guard, all R2 reps correctly select λ = 0.

## Summary

| Question | Result |
|----|----|
| Does λ track predictor quality? | Yes: perfect → N/(n+N), noise → 0, intermediate → small positive |
| Are the Louis SEs honest? | Yes: nominal coverage for both `a` and `d`; the EM bread under-covers |
| Does the method ever harm scoring? | No: tuned RMSE ≤ human-only in every regime |

The headline practical guidance: use `fit_mixed_subjects_mml` with
`tune_lambda_ability_risk`, rely on
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (Louis-corrected) for
inference, have synthetic data emit sampled responses rather than
probabilities, and treat the conditional-mean and runaway-discrimination
behaviours as the two documented edge cases.

## Reproducing these results

From the package root:

``` r

# rep counts used above; second argument is the number of cores
# Rscript simulations/run_lambda_selection.R 100 8
# Rscript simulations/run_coverage.R 200 8
# Rscript simulations/run_downstream.R 100 8
# Rscript simulations/figures.R
```

See `simulations/README.md` for the full design, the deterministic
per-task seeding (results are identical serially or in parallel), and
interpretation notes.
