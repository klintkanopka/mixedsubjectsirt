# Choosing Lambda in Mixed-Subjects IRT

The mixed-subjects estimator calibrates a 2PL model by minimizing a
PPI++-style combined loss over observed human responses, paired LLM
responses, and additional generated LLM responses.

This vignette explains how to choose $`\lambda`$ and which functions to
use. The short answer:

| Task | Function |
|----|----|
| Theoretical diagnostic | [`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md) |
| **Recommended practical tuning** | **`tune_lambda_ability_risk(..., fit_fn = fit_mixed_subjects_mml)`** |
| Frozen EC (fast approximation) | [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md) (default `fit_fn`) |
| Cross-fitted tuning | [`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md) |
| Per-item tuning (experimental) | [`tune_lambda_ability_risk_item()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_item.md) |

## Two objectives, two estimators

**Why there are two lambda objectives:**

- [`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
  minimizes $`\text{Tr}(\Sigma_\gamma)`$ — the trace of the
  item-parameter covariance matrix (PPI++ Proposition 2). This is a
  *theoretical diagnostic*: it measures whether the LLM prediction
  reduces item-parameter variance. It answers “is the predictor
  statistically useful?”

- [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  minimizes $`E[g'\Sigma_\gamma g]`$ — the propagated ability-score
  risk. This is the *practical criterion*: it answers “does using the
  LLM improve downstream test scoring?” Use this for operational
  calibration.

**Why there are two estimators:**

The original
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
uses frozen expected-count posteriors. This creates a gradient asymmetry
when the LLM item parameters differ from human parameters,
systematically inflating discriminations and driving
`tune_lambda_ability_risk` to select $`\lambda = 0`$ even when the LLM
is genuinely informative.

[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
recomputes posteriors at every gradient evaluation, eliminating this
asymmetry. Its sandwich covariance (`vcov_mixed_subjects_mml`, called
automatically via [`vcov()`](https://rdrr.io/r/stats/vcov.html)) uses
Louis’ (1982) observed-information correction for the bread, giving
honest uncertainty estimates. Use the MML estimator unless speed is a
binding constraint.

## Example data

``` r

library(mixedsubjectsirt)

set.seed(2027)

n_human     <- 120
n_generated <- 350
n_items     <- 5

human_pars <- data.frame(
  item = paste0("Item", seq_len(n_items)),
  a    = c(0.8, 1.0, 1.2, 1.35, 1.55),
  d    = c(-0.8, -0.3, 0.1, 0.45, 0.9)
)
human_pars$b <- -human_pars$d / human_pars$a

llm_pars   <- human_pars
llm_pars$a <- pmax(0.35, 0.9 * human_pars$a)
llm_pars$d <- human_pars$d + 0.25
llm_pars$b <- -llm_pars$d / llm_pars$a

theta_human <- rnorm(n_human)
observed    <- simulate_2pl(theta_human, human_pars)
predicted   <- simulate_2pl(theta_human, llm_pars)
generated   <- simulate_2pl(rnorm(n_generated), llm_pars)

lambda_grid <- c(0, 0.25, 0.5, 0.75, 1)
```

The examples use `human_pars` as `initial_pars` for speed.

## Step 1: PPI++ score diagnostic

[`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
estimates the PPI++ Proposition 2 lambda using the same human posterior
weights for both human and paired-LLM score vectors. F = Y (identical
predictions) gives exactly $`N/(n+N)`$.

``` r

ppi_score <- tune_lambda_ppi_score(
  observed    = observed,
  predicted   = predicted,
  item_pars   = human_pars,
  n_generated = nrow(generated),
  n_quad      = 7
)

cat("PPI++ score lambda:", round(ppi_score$lambda, 3),
    " r =", round(ppi_score$r, 3),
    " N/(n+N) =", round(1 / (1 + ppi_score$r), 3), "\n")
#> PPI++ score lambda: 0  r = 0.343  N/(n+N) = 0.745
```

A value near zero means the paired LLM responses do not systematically
reduce gradient variance in the person-level score formulation. This is
expected for stochastic binary LLM responses. For F = Y the formula
recovers $`N/(n+N)`$.

## Step 2: Recommended — ability-risk tuning with MML estimator

The key result from the [Linking and Gradient
Asymmetry](http://klintkanopka.com/mixedsubjectsirt/articles/linking-comparison.md)
vignette: when the LLM parameters differ from human parameters, the
frozen expected-count estimator drives $`\lambda \to 0`$ due to an
artificial gradient asymmetry. The MML estimator removes this asymmetry.
With a good predictor (F = Y), `tune_lambda_ability_risk` with MML
correctly selects $`\lambda > 0`$.

``` r

ability_tuned_mml <- tune_lambda_ability_risk(
  lambda_grid  = lambda_grid,
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  target_resp  = observed,
  initial_pars = human_pars,
  fit_fn       = fit_mixed_subjects_mml,     # marginal MML estimator
  n_quad       = 7,
  control      = list(maxit = 100)
)

ability_tuned_mml$summary[, c("lambda", "mean_param_var", "mean_total_risk",
                               "convergence", "selection_risk")]
#>   lambda mean_param_var mean_total_risk convergence selection_risk
#> 1   0.00      0.4341733       0.4341733           0      0.4341733
#> 2   0.25      0.4402754       0.4402754           0      0.4402754
#> 3   0.50      0.6286167       0.6286167           0      0.6286167
#> 4   0.75      0.6339994       0.6339994           0            Inf
#> 5   1.00      1.3329374       1.3329374           0            Inf
ability_tuned_mml$best_lambda
#> [1] 0
```

`selection_risk` is `Inf` for any candidate with non-zero convergence
code or non-finite risk, protecting selection from numerical failures.

The selected calibration carries a correctly-sized Louis-corrected
covariance:

``` r

Sigma_mml <- vcov(ability_tuned_mml$best_fit)  # dispatches to vcov_mixed_subjects_mml
dim(Sigma_mml)
#> [1] 10 10
```

## Step 3: Ability-score risk (inspect components)

``` r

risk <- ability_risk(
  resp        = observed,
  fit_or_pars = ability_tuned_mml$best_fit,
  vcov        = Sigma_mml
)
risk$summary
#>   mean_param_var mean_squared_error mean_total_risk
#> 1      0.4341733                 NA       0.4341733
```

`mean_squared_error` is `NA` because `theta_true` was not supplied; it
is only computed in simulation studies. `mean_param_var` is the
propagated item-parameter uncertainty — the criterion that
`tune_lambda_ability_risk` minimizes.

In simulation studies, pass `theta_true = theta_human` to also include
squared ability-estimation error:

``` r

ability_tuned_truth <- tune_lambda_ability_risk(
  lambda_grid  = lambda_grid,
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  target_resp  = observed,
  theta_true   = theta_human,
  initial_pars = human_pars,
  fit_fn       = fit_mixed_subjects_mml,
  n_quad       = 7,
  control      = list(maxit = 100)
)

ability_tuned_truth$summary[, c("lambda", "mean_param_var",
                                 "mean_squared_error", "mean_total_risk")]
#>   lambda mean_param_var mean_squared_error mean_total_risk
#> 1   0.00      0.4341733           4.689478        5.123652
#> 2   0.25      0.4402754           4.676972        5.117247
#> 3   0.50      0.6286167           4.650725        5.279342
#> 4   0.75      0.6339994           4.575805        5.209805
#> 5   1.00      1.3329374           4.595487        5.928425
```

## Step 4: Cross-fitted ability-risk tuning

[`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md)
estimates $`\lambda`$ per fold on training data, then fits a final
split-sample model. Two important arguments:

- `target_mode = "fixed"` (default): the full `target_resp` is used for
  every fold’s risk evaluation, which is correct when the target is an
  operational scoring population independent of the calibration sample.
- `target_mode = "row_aligned"`: subsets `target_resp` to training rows,
  valid when `target_resp = observed`.

``` r

split_id <- rep(1:2, length.out = nrow(observed))

crossfit_tuned <- tune_lambda_ability_risk_crossfit(
  lambda_grid  = c(0, 0.5, 1),
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  split_id     = split_id,
  initial_pars = human_pars,
  n_quad       = 7,
  control      = list(maxit = 100)
)

crossfit_tuned$lambda_by_split
#> [1] 0 0
```

## Step 5: Frozen expected-count estimator (fast approximation)

The default `fit_fn = fit_mixed_subjects` uses the older frozen
expected-count estimator. This is faster but produces inflated
discriminations when the LLM parameters differ from human parameters,
driving $`\lambda \to 0`$ even when the LLM is informative.

``` r

ability_tuned_ec <- tune_lambda_ability_risk(
  lambda_grid  = lambda_grid,
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  target_resp  = observed,
  initial_pars = human_pars,
  n_quad       = 7,
  slope_upper  = 4,           # required to prevent divergence
  control      = list(maxit = 100)
)

ability_tuned_ec$best_lambda
#> [1] 0
```

`slope_upper = 4` is required to prevent discriminations from diverging.
The MML estimator does not need this cap because it has no false minimum
at large discrimination values.

## Choosing a procedure

| Procedure | Objective | When to use |
|----|----|----|
| [`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md) | $`\text{Tr}(\Sigma_\gamma)`$ | Method diagnostics only |
| `tune_lambda_ability_risk(..., fit_fn = fit_mixed_subjects_mml)` | $`E[g'\Sigma_\gamma g]`$, Louis bread | **Recommended default** |
| [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md) (default `fit_fn`) | Same risk, EM bread | Fast approximation; requires `slope_upper` |
| [`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md) | Out-of-sample risk | Final inferential analyses |
| [`tune_lambda_ability_risk_item()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_item.md) | Per-item risk (approx.) | Experimental; some items poor predictors |

The target population matters. `target_resp = observed` tunes for the
observed human response patterns. In operational scoring, use a larger
target matrix representing the actual scoring population.
