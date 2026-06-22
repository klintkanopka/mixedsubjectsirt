# Simulation Validation of the Mixed-Subjects MML Estimator

This vignette presents a simulation study validating the scalar 2PL
marginal maximum-likelihood (MML) PPI++ estimator
(`fit_mixed_subjects_mml`) and its ability-risk tuning
(`tune_lambda_ability_risk`). The full reproducible harness lives in the
package’s `simulations/` directory; the tables below load the most
recent run from `simulations/results/` when available and otherwise show
cached values. Rep counts: $`\lambda`$-selection 100, coverage 200,
downstream 100, cross-fit 100.

The study answers five questions:

1.  Does $`\lambda`$-selection track predictor quality?
2.  Do standard errors achieve appropriate coverage?
3.  Does the method improve downstream scoring?
4.  What is the role of cross-fitting?
5.  Is coverage valid at the tuned $`\lambda`$?

## Design

Human responses come from a true 8-item 2PL model (`a ∈ [0.8, 1.6]`,
`d ∈ [-1.1, 1.1]`) with `n = 400` human subjects and `N = 1200`
generated subjects, abilities drawn from a standard normal. Four
predictor regimes vary the quality of the paired LLM predictions, `F`:

| Regime | Predicted F | Role |
|:---|:---|:---|
| R1: perfect prediction (F=Y) | predicted = observed (exact) | perfect predictor |
| R2: same-DGP draw | fresh binary draw from the true model | modest real signal |
| R3: independent noise | binary draw from scrambled item parameters | uninformative LLM |
| R4: LLM shift | binary draw from attenuated/shifted parameters | biased but informative LLM |

Four predictor regimes (all binar responses). {.table}

All responses (observed and predicted) are binary
$`Y_{ij}, F_{ij} \in \{0,1\}`$. Note that the package does not currently
accept probability predictions for `predicted`/`generated` (see the note
at the end of Scenario 1) or polytomous responses. The control parameter
$`\lambda \in [0,1]`$ sets how strongly the generated LLM sample is used
after subtracting the paired LLM correction; $`\lambda = 0`$ is a
human-only calibration.

## Does $`\lambda`$-selection track predictor quality?

For each regime we tune $`\lambda`$ by ability-score risk and record the
selected value; we also report the theoretical PPI++-score $`\lambda`$
for reference.

Here $`\lambda`$ is selected by
[`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md),
which optimizes the risk directly. We do not bother using the cross-fit
estimator, as we only are concerned with $`\lambda`$ magnitude.

| label                 | mean_risk | median_risk | prop_zero | mean_ppi |
|:----------------------|----------:|------------:|----------:|---------:|
| R1: perfect (F=Y)     |     0.750 |       0.750 |      0.00 |    0.750 |
| R2: same-DGP draw     |     0.119 |       0.108 |      0.07 |    0.004 |
| R3: independent noise |     0.063 |       0.040 |      0.35 |    0.000 |
| R4: LLM shift         |     0.105 |       0.104 |      0.16 |    0.002 |

Selected lambda by regime (ability-risk tuning). {.table}

**Takeaways:**

- R1 (perfect predictor) produces $`\lambda = 0.750`$, exactly
  `N / (n + N) = 1200 / 1600`, the theoretical maximum for a perfect
  paired predictor at these values of $`n`$ and $`N`$. The PPI++
  estimated $`\lambda`$ from minimizing
  $`\text{Tr}\big [ \Sigma_\gamma \big]`$ agrees.
- R3 (independent noise) produces $`\lambda ≈ 0.06`$ (median 0.04; about
  a third of reps select exactly 0). The uninformative predictor is
  correctly down-weighted to near zero. The small positive median is the
  optimizer resolving a shallow, noisy risk minimum.
- R2 (same-DGP) and R4 (LLM shift) produce $`\lambda \approx 0.11`$.
  Fresh real responses (R2) and a correlated-but-biased LLM (R4) each
  carry some score-level signal, and the criterion makes some use of it.

**Note that predictions must be sampled responses, not probabilities.**
We require `predicted`/`generated` responses and reject probability
inputs. This is not a convenience restriction. The response vector
enters inside a log-sum over quadrature points, so feeding probabilities
breaks the identity
$`\mathbb{E}\big[\nabla L_{gen}\big] = \mathbb{E}\big[\nabla L_{pred}\big]`$
that makes the PPI loss correction mean-zero. The estimator then becomes
biased at $`\lambda > 0`$, and at moderate $`\lambda`$ the objective is
unbounded for discrimination parameters, causing the fit to diverge.
Practically: if you have model-derived probabilities, sample responses
from them (e.g. `rbinom`) before calibrating.

## Do standard errors achieve appropriate coverage?

We fix $`\lambda`$ = 0.5 and compare empirical coverage of Wald
intervals built from two covariance estimates of the same fit: the
Louis-corrected marginal sandwich
([`vcov()`](https://rdrr.io/r/stats/vcov.html) dispatch) and the EM
complete-data Hessian bread
([`vcov_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)).

Note that the PPI correction `$\lambda$(L_gen − L_pred)` is mean-zero
whenever the paired and generated pseudo-responses are drawn from the
same distribution with the same ability spread, so the human term
anchors the estimand to the true parameters. This holds for every
simulation condition: R1, R2, R3 (useless LLM) and R4 (biased LLM). As
such, the estimator stays consistent for the truth even when the LLM is
uninformative or biased.

| label | louis_cov_90 | louis_cov_95 | em_cov_90 | em_cov_95 | mean_se_ratio |
|:---|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.909 | 0.955 | 0.713 | 0.787 | 1.626 |
| R2 same-DGP draw | 0.916 | 0.957 | 0.727 | 0.797 | 1.616 |
| R3 independent noise | 0.914 | 0.961 | 0.721 | 0.795 | 1.651 |
| R4 LLM shift | 0.905 | 0.954 | 0.726 | 0.792 | 1.622 |

Item-parameter CI coverage, all four regimes (200 reps). Nominal targets
0.90 and 0.95. {.table style="width:100%;"}

See that the Louis correction restores nominal coverage while the EM
bread in the sandwich covariance under-covers. Across all four
conditions, including the useless LLM (R3) and the biased LLM (R4),
Louis intervals attain ~0.91/0.96 against nominal 0.90/0.95, while the
EM bread covers only ~0.71/0.79; it understates uncertainty because it
ignores the missing information about each subject’s latent ability.

## Does the method improve downstream scoring?

We score a held-out sample of 1000 subjects with known ability and
compare ability-score RMSE for the human-only calibration ($`\lambda`$ =
0) and the tuned MML calibration.
`mean_delta = RMSE(tuned) − RMSE(human)`; negative is an improvement. We
report the paired 95% CI of `mean_delta` and `prop_improve` (the
fraction of replications where tuned beats human-only) so a small stable
effect can be told apart from noise.

| label | mean_lambda | rmse_human | rmse_tuned | mean_delta | delta_lo | delta_hi | prop_improve | bias_a |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1: perfect (F=Y) | 0.750 | 1.4961 | 1.4923 | -0.0038 | -0.0073 | -0.0003 | 0.48 | 0.010 |
| R2: same-DGP draw | 0.119 | 1.5023 | 1.5005 | -0.0019 | -0.0027 | -0.0011 | 0.59 | 0.025 |
| R3: independent noise | 0.063 | 1.5037 | 1.5030 | -0.0006 | -0.0009 | -0.0003 | 0.50 | 0.035 |
| R4: LLM shift | 0.105 | 1.5064 | 1.5046 | -0.0018 | -0.0030 | -0.0007 | 0.50 | 0.008 |

Downstream ability-score RMSE. {.table}

**Takeaways.**

- In these simulations the tuned RMSE did not increase average held-out
  RMSE relative to human-only in any regime. Every `mean_delta` is
  negative and its paired 95% CI excludes zero.
- Discrimination is unbiased at the selected $`\lambda`$
  ($`|\text{bias}_a| \leq 0.04`$ everywhere), confirming that the tuner
  is not selecting degenerate fits.
- RMSE gains are small because an 8-item test is measurement-limited:
  ability scoring is dominated by the irreducible measurement error of a
  short test. We do still see an improvement in item parameter
  precision, which contributes to reduced calibration uncertainty.

## What is the role of cross-fitting?

Cross-fitting tunes $`\lambda`$ on held-out folds so a fold’s own labels
are not used to tune the $`\lambda`$ applied to its paired correction.
We compare the non-cross-fitted tuner against
`tune_lambda_ability_risk_crossfit` (MML per-fold tuning, scalar-mean
MML final fit) on selected $`\lambda`$, item-parameter bias, coverage,
and held-out RMSE.

| label | lambda_nocf | lambda_cf | bias_a_nocf | bias_a_cf | rmse_nocf | rmse_cf | cover_nocf | cover_cf |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| R1 perfect (F=Y) | 0.750 | 0.858 | 0.0104 | 0.0100 | 1.4923 | 1.4927 | 0.954 | 0.949 |
| R2 same-DGP draw | 0.119 | 0.135 | 0.0250 | 0.0258 | 1.5005 | 1.5005 | 0.954 | 0.956 |
| R4 LLM shift | 0.100 | 0.115 | 0.0374 | 0.0372 | 1.5027 | 1.5029 | 0.953 | 0.956 |

Cross-fitted vs non-cross-fitted tuning (100 reps, R1/R2/R4). {.table
style="width:100%;"}

Cross-fitting does not change observed behavior in these simulations,
but guards against finite-sample bias. Item-parameter bias, held-out
RMSE, and 95% Wald coverage are essentially identical between the
non-cross-fitted and cross-fitted tuners (differences in the fourth
decimal). Because the downstream quantities are unchanged, the cheaper
non-cross-fitted tuner carries value for intermediate testing and
exploration, but prefer cross-fitting for final model fitting.

## Is coverage valid at the tuned $`\lambda`$?

Scenario 2 validated the covariance *formula* at a fixed $`\lambda`$,
but the $`\lambda`$ we use is estimated from data, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) treats $`\lambda`$ as
fixed. As such, it does not propagate the uncertainty in
$`\hat{\lambda}`$ into $`\hat{\Sigma}_\gamma`$. Additionally, choosing
$`\lambda`$ on the same data used to estimate the item parameters
induces finite-sample bias in those estimates (and, in principle,
anti-conservative intervals). Split-sample (cross-fit) $`\lambda`$
estimation removes that bias.

This scenario measures the size of the effect by comparing coverage of
the true item parameters at three operating points: fixed $`\lambda`$ =
0.5, the same-data tuned $`\hat{\lambda}`$ with `vcov(best_fit)`, and
the cross-fit tuned $`\hat{\lambda}`$ with `vcov(final_fit)`. It reuses
the Scenario 2 seeds, so the fixed-$`\lambda`$ column reproduces
Scenario 2 exactly.

| label                | lambda_sd | lambda_xf | fixed_95 | samedata_95 | crossfit_95 |
|:---------------------|----------:|----------:|---------:|------------:|------------:|
| R1 perfect (F=Y)     |     0.751 |     0.859 |    0.955 |       0.955 |       0.953 |
| R2 same-DGP draw     |     0.116 |     0.135 |    0.957 |       0.956 |       0.956 |
| R3 independent noise |     0.053 |     0.074 |    0.961 |       0.958 |       0.958 |
| R4 LLM shift         |     0.103 |     0.130 |    0.954 |       0.958 |       0.958 |

95% coverage rates of the true item parameters at the fixed,
same-data-tuned, and cross-fit-tuned λ (200 reps). lambda_sd / lambda_xf
are the mean selected λ for each tuner. {.table}

And the matching item-parameter (discrimination) bias:

| label                | bias_a_fixed | bias_a_samedata | bias_a_crossfit |
|:---------------------|-------------:|----------------:|----------------:|
| R1 perfect (F=Y)     |       0.0097 |          0.0063 |          0.0051 |
| R2 same-DGP draw     |       0.0280 |          0.0224 |          0.0222 |
| R3 independent noise |       0.0365 |          0.0286 |          0.0289 |
| R4 LLM shift         |       0.0207 |          0.0159 |          0.0166 |

Mean discrimination bias E\[a-hat - a\] at each operating point.
Same-data vs cross-fit differ by \<= 0.001 except R1. {.table}

**Takeaways.**

- Same-data tuning produces nominal coverage (95% coverage 0.955–0.958),
  and the finite-sample bias is small: same-data vs. cross-fit
  discrimination bias differ by $`\leq 0.001`$ in every regime except
  R1.
- Cross-fitting removes bias in the perfect predictor (R1), but this
  bias is small.
- Cross-fitting gives the same coverage (0.953–0.958) while also
  guarding against finite-sample tuning bias.

Cross-fitting is the principled default, despite the same-data tuning
performing well. It is guaranteed free of the finite-sample tuning bias.

## Summary

| Question | Result |
|----|----|
| Does $`\lambda`$-selection track predictor quality? | Yes: better predictors select higher $`\lambda`$ |
| Do standard errors achieve appropriate coverage? | Yes: nominal coverage, including under a biased predictor (R3/R4) |
| Does the method improve downstream scoring? | No average harm demonstrated and per-rep improvement in about half of reps |
| What is the role of cross-fitting? | Cross-fitting removes finite sample bias, despite reporting similar bias/RMSE/coverage |
| Is coverage valid at the tuned $`\lambda`$? | Yes |

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
