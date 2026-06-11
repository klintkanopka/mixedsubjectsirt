# Calibrating with a Weakly-Informative, Biased LLM

This vignette treats the regime prediction-powered inference is built
for: a **small human sample** (here `n = 500`) alongside a **large
synthetic / LLM sample** (`N = 100000`). The LLM is **biased** — its
item parameters are systematically off — and only **weakly informative**
about the human responses.

The whole story is one distinction, made visual below:

> The mixed-subjects (PPI) estimator is **asymptotically unbiased for
> the true human parameters at *every* λ.** λ is an **efficiency** knob,
> not a bias knob. Bias is a property of the *estimation structure*, not
> of λ. A *naive* fit that pools the human and LLM responses has no such
> protection: the `n = 500` humans are outvoted by the `N = 100000` LLM
> rows, and the estimate inherits the LLM’s bias.

*All numbers are precomputed (`data-raw/precompute_largeN.R`): n = 500,
N = 100000, 16 Monte Carlo replications.*

## The setup

Human responses come from a true 8-item 2PL model (`a ∈ [0.8, 1.6]`,
`d ∈ [-1.1, 1.1]`). The LLM is a **shifted** version — slopes attenuated
by 10%, intercepts shifted up by 0.25 — so its responses are plausible
but biased:

``` r

true_pars <- data.frame(item = paste0("Item", 1:8),
                        a = seq(0.8, 1.6, length.out = 8),
                        d = seq(-1.1, 1.1, length.out = 8))
llm   <- true_pars
llm$a <- 0.9 * true_pars$a       # ~10% attenuated discriminations
llm$d <- true_pars$d + 0.25      # +0.25 intercept shift

theta     <- rnorm(500)
observed  <- simulate_2pl(theta, true_pars)             # n = 500 human
predicted <- simulate_2pl(theta, llm)                   # paired LLM (same people)
generated <- simulate_2pl(rnorm(100000), llm)           # N = 100000 unlabeled LLM
```

## Naive pooling inherits the bias

The obvious move is to pool everything and fit one model:

``` r

naive <- fit_2pl(rbind(observed, generated))   # 500 human + 100000 LLM rows
```

The 500 humans *are* in the fit, but against 100,000 LLM rows their
information is washed out, and the estimate is dragged onto the LLM’s
shifted parameters:

Averaged over the replications, the naive estimator’s item-parameter
bias is **-0.119** in the slopes and **+0.248** in the intercepts —
essentially the LLM’s shift (−0.1·a, +0.25). Because `N = 100000`, that
wrong answer is estimated *very precisely* (a tiny standard error); more
LLM data only sharpens the bias.

## λ moves efficiency, not bias

The mixed-subjects estimator minimizes

``` math
L_o^{\mathrm{marg}}(\gamma) \;+\; \lambda\bigl[L_g^{\mathrm{marg}}(\gamma) - L_p^{\mathrm{marg}}(\gamma)\bigr].
```

At the true parameters the human score is mean-zero and the paired
correction `L_g − L_p` is *also* mean-zero (generated and paired LLM
responses share a distribution and ability spread), so the estimating
equation is mean-zero **for every λ**. Unbiasedness comes from this
*structure*, not from λ.

To see it directly, fit the estimator across a fine grid of λ and track
two things: the **item-parameter bias** (Monte Carlo mean of
`estimate − truth`) and the **model-based ability-score risk**
`E[g'Σ_γ(λ) g]` — the quantity the tuner actually minimizes.

![Item-parameter bias of the mixed-subjects estimator is flat near zero
across all lambda, far from the naive pooled bias shown as dashed
reference
lines.](weakly-informative-llm_files/figure-html/bias-curve-1.png)

The mixed-subjects bias sits on zero across the entire range of λ (the
shaded band is ±2 Monte Carlo SE); the dashed red lines mark the naive
pooled bias, off the chart. **No λ buys or removes bias.** What λ does
change is efficiency:

![Model-based ability-score risk as a function of lambda, with a shallow
minimum near the optimized
lambda.](weakly-informative-llm_files/figure-html/risk-curve-1.png)

For this weakly-informative LLM the averaged risk curve is **shallow**
and rises for larger λ: leaning on a poorly-correlated predictor *adds*
item-parameter variance. Its minimum sits near λ ≈ 0.1, only about 2%
below the λ = 0 (human-only) risk — almost no efficiency to be gained.
Because the curve is so flat, each individual dataset’s optimum scatters
around this value (the red ticks); see the next section. Every point on
the curve is unbiased.

## Choosing λ: the tuner optimizes directly

The curve above was sampled on a grid only to *draw the surface* — in
fact with `tune_lambda_ability_risk(..., method = "grid")`. To
**choose** an operating λ you do not need a grid at all. By default
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
selects λ by **direct 1-D optimization** of the risk over `[0, 1]`
([`stats::optimize()`](https://rdrr.io/r/stats/optimize.html)),
returning a continuous value with no grid rounding and no evaluations
wasted far from the optimum:

``` r

# Direct optimization is the default (method = "optimize").
tuned <- tune_lambda_ability_risk(
  observed = observed, predicted = predicted, generated = generated,
  target_resp = observed, initial_pars = human_start$pars,
  fit_fn = fit_mixed_subjects_mml, n_quad = 11
)
tuned$best_lambda            # continuous lambda

# Pass method = "grid" (and a lambda_grid) to scan instead -- how the curve above
# was drawn. lambda_grid otherwise just bounds the optimizer's search.
```

The optimizer returns the minimizer of *this dataset’s* risk surface
continuously — here λ = **0.27**. That is a genuine per-dataset optimum,
not noise: every dataset has its own (noisy) risk surface, so its
optimal λ varies around the population value. Across the 16 replications
the per-dataset optimum averaged 0.14 and ranged \[0.0, 0.3\],
scattering around the minimum of the *averaged* curve (≈ 0.1). (These
are not the same point — the minimum of the average risk is not the
average of the per-dataset minima.) The scatter is wide here because the
surface is shallow; an informative LLM sharpens it. `method = "grid"`
would snap each of these to the nearest node; the default optimizer
resolves them continuously. Either way the estimate is unbiased.

(The 2-fold cross-fitted tuner,
[`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md),
lands at the same place: at `N ≫ n` the cross-fit λ-inflation
`N/(N + n/2)` vs `N/(N + n)` is negligible, so cross-fitting changes
neither the bias nor the selected λ here.)

## Takeaways

1.  **Bias is set by the estimation structure, not by λ.** The
    mixed-subjects estimator is unbiased for the true human parameters
    at every λ; pooling lets a large biased LLM sample outvote the human
    anchor and inherits its bias.
2.  **λ tunes efficiency, and the objective decides which.**
    `tune_lambda_ppi_score` minimizes the item-parameter trace
    `Tr(Σ_γ)`; `tune_lambda_ability_risk` minimizes the propagated
    ability-recovery risk `E[g'Σ_γ g]`. For a weakly-informative LLM
    both curves are shallow — little efficiency to gain — but neither
    sacrifices unbiasedness.
3.  **You don’t need a grid to tune.**
    [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
    selects λ by direct 1-D optimization by default; a grid
    (`method = "grid"`) is just a convenient way to *see* the whole risk
    surface, as above.

## Reproducing

`data-raw/precompute_largeN.R` runs the Monte Carlo over the λ grid and
the direct optimization, and writes the cached results
(`Rscript data-raw/precompute_largeN.R [n_reps] [cores] [N]`). At
`N = 100000` each fit takes several seconds, so it is run once offline
rather than during vignette knitting; pass a larger `N` to confirm the
picture is unchanged.
