# Simulation Validation of the Mixed-Subjects MML Estimator

This vignette presents a simulation study validating the scalar 2PL
marginal maximum-likelihood (MML) PPI++ estimator
(`fit_mixed_subjects_mml`) and its ability-risk tuning
(`tune_lambda_ability_risk`). The full reproducible harness lives in the
package’s `simulations/` directory; the tables below load the most
recent run from `simulations/results/` when available and otherwise show
cached values. Rep counts: λ-selection 100, coverage 200, downstream
100, cross-fit 100.

The study answers five questions:

1.  **Does λ-selection track predictor quality?** (Scenario 1)
2.  **Are the Louis-corrected standard errors honest, including under a
    biased LLM?** (Scenario 2 — the key check)
3.  **Does the method improve downstream scoring without harming it on
    average?** (Scenario 3)
4.  **Does cross-fitting change the conclusions?** (Scenario 4)
5.  **Is coverage still valid at the *tuned* λ — the no-free-lunch
    question?** (Scenario 5)

## Design

Human responses come from a true 8-item 2PL model (`a ∈ [0.8, 1.6]`,
`d ∈ [-1.1, 1.1]`) with `n = 400` human subjects and `N = 1200`
generated subjects, abilities drawn from a standard normal. Four
predictor regimes vary the quality of the paired LLM predictions `F`:

| Regime | Predicted F | Role |
|:---|:---|:---|
| R1 perfect (F=Y) | predicted = observed (exact) | perfect predictor |
| R2 same-DGP draw | fresh binary draw from the true model | modest real signal |
| R3 independent noise | binary draw from scrambled item parameters | useless LLM |
| R4 LLM shift | binary draw from attenuated/shifted parameters | biased but informative LLM |

Four predictor regimes (all binary). {.table}

All predictors are **binary 0/1 responses** — the package does not
accept probability (fractional) predictions for `predicted`/`generated`
(see the note at the end of Scenario 1). The control parameter λ ∈ \[0,
1\] sets how strongly the generated LLM sample is used after subtracting
the paired LLM correction; λ = 0 is the human-only calibration.

## Scenario 1: λ selection

For each regime we tune λ by ability-score risk and record the selected
value; we also report the theoretical PPI++-score λ for reference.

(λ is selected by
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md),
which optimizes the risk directly, so the selected values are
*continuous*. The means match a grid run to the third decimal and the
regime ordering is identical, but the medians are off-grid and
`prop_zero` — the fraction selecting *exactly* 0 — is lower than a grid
gives, because the optimizer resolves the shallow *positive* optimum
that a grid rounds down to 0.)

| label                | mean_risk | median_risk | prop_zero | mean_ppi |
|:---------------------|----------:|------------:|----------:|---------:|
| R1 perfect (F=Y)     |     0.750 |       0.750 |      0.00 |    0.750 |
| R2 same-DGP draw     |     0.119 |       0.108 |      0.07 |    0.004 |
| R3 independent noise |     0.063 |       0.040 |      0.35 |    0.000 |
| R4 LLM shift         |     0.105 |       0.104 |      0.16 |    0.002 |

Selected lambda by regime (ability-risk tuning). {.table}

**Findings.**

- **R1 perfect → λ = 0.750**, exactly `N / (n + N) = 1200 / 1600`, the
  theoretical maximum for a perfect paired predictor. The PPI++-score λ
  agrees to the decimal.
- **R3 independent noise → λ ≈ 0.06** (median 0.04; about a third of
  reps select exactly 0): a useless predictor is correctly down-weighted
  to near zero. The small positive median is the optimizer resolving a
  shallow, noisy risk minimum — benign, since the downstream scoring is
  unaffected (Scenario 3).
- **R2 same-DGP and R4 LLM shift → λ ≈ 0.11**: fresh real responses (R2)
  and a correlated-but-biased LLM (R4) each carry a little score-level
  signal, and the criterion uses a little of it.

### Predictions must be sampled responses, not probabilities

The package requires **binary** `predicted`/`generated` and rejects
probability inputs. This is not a convenience restriction — it removes a
genuine failure mode. A fractional value such as a conditional mean
`p(θ)` is not a coherent likelihood term for the marginal IRT objective:
the response vector enters inside a log-sum over quadrature points, so
feeding probabilities breaks the identity `E[∇L_gen] = E[∇L_pred]` that
makes the PPI correction mean-zero. The estimator then becomes biased at
λ \> 0, and at moderate λ the objective is even *unbounded* in
discrimination (`a → ∞`), so the fit diverges. An earlier version of
this study included a “conditional-mean” regime that exhibited exactly
this pathology; it has been removed along with probability support. The
practical rule: if you have LLM-derived probabilities, **sample** binary
responses from them (e.g. `rbinom`) before calibrating.

## Scenario 2: Louis-corrected SE coverage — the key check

We fix λ = 0.5 and compare empirical coverage of Wald intervals built
from two covariance estimates of the *same* fit: the Louis-corrected
marginal sandwich ([`vcov()`](https://rdrr.io/r/stats/vcov.html)
dispatch) and the EM complete-data Hessian bread
([`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md),
a deliberate bypass).

**Consistency for the true parameters — even under a biased LLM.** The
PPI correction `λ(L_gen − L_pred)` is *mean-zero* whenever the paired
and generated pseudo-responses are drawn from the same distribution with
the same ability spread:

``` math
E[\Psi(\gamma^\star)] = \underbrace{E[\nabla L_{\text{obs}}]}_{=0}
 + \lambda\,(\,\underbrace{E[\nabla L_{\text{gen}}] - E[\nabla L_{\text{pred}}]}_{=0}\,) = 0.
```

So the human term anchors the estimand to the **true** parameters. This
holds for every (binary) regime — R1, R2, R3 (useless LLM) and R4
(biased LLM) — so the estimator stays consistent for the truth even when
the LLM is wrong. **R3 and R4 covering nominally is the flagship
demonstration that PPI corrects biased LLM outputs.**

| label | louis_cov_90 | louis_cov_95 | em_cov_90 | em_cov_95 | mean_se_ratio |
|:---|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.909 | 0.955 | 0.713 | 0.787 | 1.626 |
| R2 same-DGP draw | 0.916 | 0.957 | 0.727 | 0.797 | 1.616 |
| R3 independent noise | 0.914 | 0.961 | 0.721 | 0.795 | 1.651 |
| R4 LLM shift | 0.905 | 0.954 | 0.726 | 0.792 | 1.622 |

Item-parameter CI coverage, all four regimes (200 reps). Nominal targets
0.90 and 0.95. {.table style="width:100%;"}

**The Louis correction restores nominal coverage; the EM bread
under-covers.** Across all four regimes — including the useless LLM (R3)
and the biased LLM (R4) — Louis intervals attain ~0.91/0.96 against
nominal 0.90/0.95, while the EM bread covers only ~0.71/0.79; it
understates uncertainty because it ignores the missing information about
each subject’s latent ability. The Louis standard errors are ~1.6×
wider, the difference between an honest interval and an overconfident
one. That R3 and R4 cover nominally — with all 200 reps usable — is the
concrete demonstration that PPI keeps inference valid even when the LLM
is wrong.

### A bug this scenario caught

An earlier version showed near-zero coverage for discrimination while
intercept coverage was nominal — a uniform downward bias in all slopes
with unbiased intercepts, the signature of a latent-scale error. It
traced to
[`make_quadrature()`](http://klintkanopka.com/mixedsubjectsirt/reference/make_quadrature.md)
building an N(0, 2) ability grid instead of N(0, 1) (a spurious √2
rescaling of standard-normal Gauss–Hermite nodes); every discrimination
had been biased downward by ~1/√2 ≈ 0.71. The bug was invisible to
relative checks (λ-selection), to gradient finite-difference checks
(loss and gradient were mutually consistent), and to the SE calibration
itself (the SEs correctly described the *biased* estimator). Only
coverage against ground truth exposed it; after the fix, discrimination
coverage rose from 0.04 to 0.91. This is the strongest argument for
ground-truth simulation validation: it catches errors that
internal-consistency checks cannot.

## Scenario 3: downstream payoff and no-harm

We score a held-out sample of 1000 subjects with known ability and
compare ability-score RMSE for the human-only calibration (λ = 0) and
the tuned MML calibration. `mean_delta = RMSE(tuned) − RMSE(human)`;
negative is an improvement. We report the paired 95% CI of `mean_delta`
and `prop_improve` (the fraction of replications where tuned beats
human-only) so a small stable effect can be told apart from noise.

| label | mean_lambda | rmse_human | rmse_tuned | mean_delta | delta_lo | delta_hi | prop_improve | bias_a |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.750 | 1.4961 | 1.4923 | -0.0038 | -0.0073 | -0.0003 | 0.48 | 0.010 |
| R2 same-DGP draw | 0.119 | 1.5023 | 1.5005 | -0.0019 | -0.0027 | -0.0011 | 0.59 | 0.025 |
| R3 independent noise | 0.063 | 1.5037 | 1.5030 | -0.0006 | -0.0009 | -0.0003 | 0.50 | 0.035 |
| R4 LLM shift | 0.105 | 1.5064 | 1.5046 | -0.0018 | -0.0030 | -0.0007 | 0.50 | 0.008 |

Downstream ability-score RMSE. {.table}

**Findings.**

- **No average harm.** In these simulations the tuned RMSE did not
  increase average held-out RMSE relative to human-only in any regime —
  every `mean_delta` is negative and its paired 95% CI excludes zero.
  This is evidence *for these simulated conditions*, not a general
  guarantee: the non-asymptotic PPI++ literature shows that power-tuned
  estimators can underperform classical estimation when pseudo-label
  correlation is very weak. The `prop_improve` column keeps the claim
  honest at the per-replication level — only about half of individual
  replications improve (0.48–0.59), so “no harm” is an average
  statement; the gains are consistent in aggregate but small relative to
  rep-level noise.
- **Discrimination is unbiased** at the selected λ (\|bias_a\| ≤ 0.04
  everywhere), confirming the quadrature fix and that the tuner is not
  selecting degenerate fits.
- **The RMSE gains are small** because an 8-item test is
  measurement-limited: ability scoring is dominated by the irreducible
  measurement error of a short test. The method’s larger payoff is in
  **item-parameter precision** — the Scenario 2 standard errors — which
  matters most when calibration uncertainty propagates into operational
  use.

### A robustness guard

The ability-risk tuners carry a defensive guard: any candidate fit whose
maximum discrimination exceeds `max_discrimination` (default 10) is
excluded from selection. A runaway-discrimination fit reports huge item
information, which collapses its covariance and makes the model-based
ability risk spuriously small, so an unguarded tuner could select it
even though it scores worse on held-out data. A 2PL discrimination above
10 is never a real calibration. With binary inputs the operational
regimes never approach this bound, but the guard keeps the tuner safe
against pathological optimizer behavior.

## Scenario 4: cross-fitted vs non-cross-fitted tuning

Cross-fitting tunes λ on held-out folds so a fold’s own labels are not
used to tune the λ applied to its paired correction. We compare the
non-cross-fitted tuner against `tune_lambda_ability_risk_crossfit` (MML
per-fold tuning, scalar-mean MML final fit) on selected λ,
item-parameter bias, 95% Wald coverage, and held-out RMSE.

| label | lambda_nocf | lambda_cf | bias_a_nocf | bias_a_cf | rmse_nocf | rmse_cf | cover_nocf | cover_cf |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.750 | 0.858 | 0.0104 | 0.0100 | 1.4923 | 1.4927 | 0.954 | 0.949 |
| R2 same-DGP draw | 0.119 | 0.135 | 0.0250 | 0.0258 | 1.5005 | 1.5005 | 0.954 | 0.956 |
| R4 LLM shift | 0.100 | 0.115 | 0.0374 | 0.0372 | 1.5027 | 1.5029 | 0.953 | 0.956 |

Cross-fitted vs non-cross-fitted tuning (100 reps, R1/R2/R4). {.table
style="width:100%;"}

**Cross-fitting does not change the conclusions.** Item-parameter bias,
held-out RMSE, and 95% Wald coverage are essentially identical between
the non-cross-fitted and cross-fitted tuners (differences in the fourth
decimal). The one systematic difference is the selected λ: cross-fitting
picks a slightly *higher* value (e.g. R1 0.75 → 0.86). This is a
sample-size effect of the split, not a better λ — these runs use 2
folds, so each fold tunes on only `n_train = n/2 = 200` labeled
subjects, and the PPI-optimal λ grows as labeled data shrinks (it tracks
`N/(N + n_train) = 1200/1400 ≈ 0.86`). The fold λ’s are tuned for
`n_train` but applied to the full sample, so they run a little high;
more folds shrink the gap. The genuine reason to cross-fit is not the λ
value but that λ is selected without each fold’s own labels, so reported
uncertainty is not optimistic about a λ tuned on the same data. Because
the downstream quantities are unchanged, the cheaper non-cross-fitted
tuner is an adequate default; prefer cross-fitting when the selected λ
or its uncertainty is itself reported.

## Scenario 5: operational coverage at the *tuned* λ

Scenario 2 validated the covariance *formula* at a **fixed** λ. But the
λ you report is **tuned from the data**, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) treats λ as known — it
does not propagate the uncertainty in λ̂. There is no free lunch:
choosing λ on the same data used to estimate the item parameters induces
a **finite-sample bias** in those estimates (and, in principle,
anti-conservative intervals). **Split-sample (cross-fit) λ estimation
removes that bias.** This scenario measures the size of the effect by
comparing Wald coverage of the *true* item parameters at three operating
points: fixed λ = 0.5, the same-data tuned λ̂ with `vcov(best_fit)`, and
the cross-fit tuned λ̂ with `vcov(final_fit)`. It reuses the Scenario 2
seeds, so the fixed-λ column reproduces Scenario 2 exactly.

| label                | lambda_sd | lambda_xf | fixed_95 | samedata_95 | crossfit_95 |
|:---------------------|----------:|----------:|---------:|------------:|------------:|
| R1 perfect (F=Y)     |     0.751 |     0.859 |    0.955 |       0.955 |       0.953 |
| R2 same-DGP draw     |     0.116 |     0.135 |    0.957 |       0.956 |       0.956 |
| R3 independent noise |     0.053 |     0.074 |    0.961 |       0.958 |       0.958 |
| R4 LLM shift         |     0.103 |     0.130 |    0.954 |       0.958 |       0.958 |

95% Wald coverage of the true item parameters at the fixed,
same-data-tuned, and cross-fit-tuned λ (200 reps). lambda_sd / lambda_xf
are the mean selected λ for each tuner. {.table}

And the matching item-parameter (discrimination) bias — the
no-free-lunch term read off directly:

| label                | bias_a_fixed | bias_a_samedata | bias_a_crossfit |
|:---------------------|-------------:|----------------:|----------------:|
| R1 perfect (F=Y)     |       0.0097 |          0.0063 |          0.0051 |
| R2 same-DGP draw     |       0.0280 |          0.0224 |          0.0222 |
| R3 independent noise |       0.0365 |          0.0286 |          0.0289 |
| R4 LLM shift         |       0.0207 |          0.0159 |          0.0166 |

Mean discrimination bias E\[a-hat - a\] at each operating point.
Same-data vs cross-fit differ by \<= 0.001 except R1. {.table}

**Findings.**

- **The fixed-λ column reproduces Scenario 2** (0.955 / 0.957 / 0.961 /
  0.954) — the seed reuse is exact, confirming the harness.
- **Same-data tuning keeps coverage at nominal here** (95% coverage
  0.955–0.958), and the no-free-lunch finite-sample bias it warns about
  is **small**: same-data vs cross-fit discrimination bias differ by ≤
  0.001 in every regime except R1.
- **The one place the bias term is visible is R1** (the perfect
  predictor, where λ̂ is large and the optimum is sharp): cross-fitting
  trims `bias_a` from **0.0063 → 0.0051**. That is the no-free-lunch
  bias, removed by split-sample tuning — present and in the expected
  direction, but tiny.
- **Both tuned estimators are *less* biased than the fixed λ = 0.5**
  (e.g. R3 0.0365 → ≈ 0.029): λ = 0.5 is simply a poor operating point,
  and tuning toward the optimum reduces bias more than same-data tuning
  could ever add. (Intercept bias is smaller still, ≤ 0.014.)
- **Cross-fitting gives the same coverage** (0.953–0.958) and removes
  the finite-sample tuning bias *by construction*; it selects a higher λ̂
  (the small-per-fold inflation of Scenario 4) but, since the estimator
  is unbiased at any λ, coverage is unchanged.

So in this design the cheap **same-data tuner is adequate** — its
reported intervals are valid and its bias is within Monte Carlo noise of
the cross-fit estimator. **Cross-fitting is the principled default**: it
is *guaranteed* free of the finite-sample tuning bias (visible at R1),
and it costs nothing here. Prefer it when the tuning is more aggressive
(a finer or wider λ search), the predictor is flexible, or the selected
λ̂ and its uncertainty are themselves reported.

## Summary

| Question | Result |
|----|----|
| Does λ track predictor quality? | Yes: perfect → N/(n+N), noise → 0, intermediate → small positive |
| Are the Louis SEs honest? | Yes: nominal coverage, including under a biased LLM (R3/R4); EM bread under-covers |
| Does the method harm scoring? | No average harm in these simulations; per-rep improvement is about half of reps |
| Does cross-fitting change conclusions? | No: bias/RMSE/coverage unchanged; cross-fit selects a slightly higher λ |
| Is coverage valid at the tuned λ? | Yes here: same-data ~nominal (the no-free-lunch finite-sample bias is small); cross-fit removes it by construction |

Practical guidance: use `fit_mixed_subjects_mml` with
`tune_lambda_ability_risk`, rely on
[`vcov()`](https://rdrr.io/r/stats/vcov.html) (Louis-corrected) for
inference, and have synthetic/LLM data emit **sampled binary
responses**, not probabilities — the package enforces this, and it is
what keeps the PPI correction well-behaved.

## Reproducing these results

From the package root (second argument is the number of cores):

``` r

# Rscript simulations/run_lambda_selection.R 100 8
# Rscript simulations/run_coverage.R 200 8
# Rscript simulations/run_downstream.R 100 8
# Rscript simulations/run_crossfit.R 50 8
# Rscript simulations/figures.R
```

See `simulations/README.md` for the full design, the deterministic
per-task seeding (results are identical serially or in parallel), and
interpretation notes.
