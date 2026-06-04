# Choosing Lambda in Mixed-Subjects IRT

The mixed-subjects estimator fits a 2PL calibration by minimizing

`L_human + lambda * (L_generated - L_paired_llm)`

The tuning parameter `lambda` controls how strongly generated LLM
responses are used after subtracting the paired LLM correction. Setting
`lambda = 0` gives the human-only expected-count calibration. Larger
values use more information from the generated response sample, but they
can also amplify LLM bias when the paired correction is weak.

This vignette distinguishes four package workflows:

- [`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md):
  a sensitivity diagnostic. It is not a valid final tuning rule.
- `tune_lambda_ppi_score()`: returns the PPI++ Proposition 2 plug-in
  estimate, the $`\lambda`$ that minimises the *trace of the
  item-parameter covariance matrix* $`\text{Tr}(\Sigma_\gamma)`$. **This
  is a theoretical diagnostic.** For psychometric applications, use
  [`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
  instead.
- [`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md):
  chooses a lambda by minimizing propagated *ability-score risk*
  $`\mathbb{E}[g' \Sigma_\gamma g]`$ on a target scoring population.
  **This is the recommended practical criterion** for IRT applications
  where accurate test scoring is the goal.
- [`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md):
  estimates lambda out of sample by fold, then fits the final
  split-sample estimator with fold-specific lambdas.

**The two lambda objectives.** The PPI++ score objective minimises
item-parameter estimation efficiency ($`\text{Tr}(\Sigma_\gamma)`$). The
ability-risk objective minimises downstream scoring accuracy
($`\mathbb{E}[g' \Sigma_\gamma g]`$). These are different quantities and
generally select different $`\lambda`$ values. Because the practical
goal in psychometrics is accurate ability scoring, not item-parameter
precision,
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
is the recommended default.

## Example data

The example uses simulated human responses, paired LLM responses for the
same human rows, and additional generated LLM responses. In an applied
project, replace these matrices with your study data.

``` r

library(mixedsubjectsirt)

set.seed(2027)

n_human <- 120
n_generated <- 350
n_items <- 5

human_pars <- data.frame(
  item = paste0("Item", seq_len(n_items)),
  a = c(0.8, 1.0, 1.2, 1.35, 1.55),
  d = c(-0.8, -0.3, 0.1, 0.45, 0.9)
)
human_pars$b <- -human_pars$d / human_pars$a

llm_pars <- human_pars
llm_pars$a <- pmax(0.35, 0.9 * human_pars$a)
llm_pars$d <- human_pars$d + 0.25
llm_pars$b <- -llm_pars$d / llm_pars$a

theta_human <- rnorm(n_human)
observed <- simulate_2pl(theta_human, human_pars)
predicted <- simulate_2pl(theta_human, llm_pars)
generated <- simulate_2pl(rnorm(n_generated), llm_pars)

lambda_grid <- c(0, 0.25, 0.5, 0.75, 1)
```

The examples below use `human_pars` as `initial_pars` so the vignette
runs quickly. In a real package workflow, a typical starting point is
`initial_pars = fit_2pl(observed)$pars`.

## Fixed Lambda

When `lambda` is fixed by design, call
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
directly.

``` r

fixed_fit <- fit_mixed_subjects(
  observed = observed,
  predicted = predicted,
  generated = generated,
  lambda = 0.5,
  initial_pars = human_pars,
  n_quad = 7,
  slope_upper = 4,
  control = list(maxit = 100)
)

fixed_fit$item_pars
#>    item        a          d          b
#> 1 Item1 1.136801 -1.5492749  1.3628377
#> 2 Item2 2.357070 -1.1258387  0.4776433
#> 3 Item3 1.485432 -0.1796543  0.1209441
#> 4 Item4 1.945808  0.6330096 -0.3253196
#> 5 Item5 4.000000  1.2612756 -0.3153189
```

Use this route when a simulation study, design analysis, or
preregistered plan already specifies `lambda`.

## Sensitivity Diagnostics

[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
fits candidate lambda values and reports the mixed objective and the
observed-human expected-count loss. This is useful for seeing how
sensitive the fitted calibration is to lambda.

The value in `lowest_observed_loss_lambda` is not a final PPI++ tuning
rule. It is the lambda with the smallest observed-human expected-count
loss after fitting. That criterion can overfit the labeled sample and is
not the same as minimizing downstream ability-score uncertainty.

When the LLM parameters differ from human parameters (here: attenuated
discrimination, shifted intercepts), the correction gradient
$`\nabla(L_\text{gen} - L_\text{pred})`$ systematically pushes item
discriminations upward at any $`\lambda > 0`$. Without an upper bound on
discriminations this leads to convergence code 52 (L-BFGS-B line-search
failure) at large $`\lambda`$. The `slope_upper = 4` argument caps
discriminations so the optimizer converges at the bound instead of
diverging.

``` r

diagnostic <- diagnose_lambda_grid(
  lambda_grid = lambda_grid,
  observed = observed,
  predicted = predicted,
  generated = generated,
  initial_pars = human_pars,
  n_quad = 7,
  slope_upper = 4,
  control = list(maxit = 100)
)

diagnostic$summary
#>   lambda mixed_loss observed_loss convergence
#> 1   0.00  2.5828751      2.582875           0
#> 2   0.25  2.3030718      2.640973           0
#> 3   0.50  1.8506861      2.919215           0
#> 4   0.75  1.1908046      3.451707           0
#> 5   1.00  0.3898802      3.772615           0
diagnostic$lowest_observed_loss_lambda
#> [1] 0
```

## PPI++ Score Tuning (Theoretical Diagnostic)

`tune_lambda_ppi_score()` implements the closed-form plug-in estimator
from Proposition 2 of Angelopoulos, Duchi and Zrnic (2023). It returns
the $`\lambda`$ that minimises the trace of the item-parameter
asymptotic covariance matrix $`\text{Tr}(\Sigma_\gamma)`$.

This objective measures **item-parameter estimation efficiency** and is
suitable for method validation and theoretical analysis. For
psychometric applications where accurate test scoring is the goal, use
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
instead — it directly minimises the propagated ability-score risk
$`\mathbb{E}[g' \Sigma_\gamma g]`$.

The function uses the **same** human posterior weights for both the
human and paired-LLM score vectors. This symmetry satisfies the PPI++
unbiasedness condition.

``` r

ppi_score <- tune_lambda_ppi_score(
  observed    = observed,
  predicted   = predicted,
  item_pars   = human_pars,
  n_generated = nrow(generated),
  n_quad      = 7
)

cat("PPI++ score lambda (minimises Tr(Sigma_gamma)):", round(ppi_score$lambda, 3), "\n")
#> PPI++ score lambda (minimises Tr(Sigma_gamma)): 0
cat("  r = n/N =", round(ppi_score$r, 3),
    "  => N/(n+N) upper bound =", round(1/(1 + ppi_score$r), 3), "\n")
#>   r = n/N = 0.343   => N/(n+N) upper bound = 0.745
```

The PPI++ lambda is typically lower than the upper bound $`N/(n+N)`$
because the paired predictions are stochastic (independent draws from
the LLM, not exact copies of the human responses). For exact predictions
(`predicted = observed`), the formula recovers
$`\lambda^* \approx N/(n+N)`$.

## Ability-Risk Tuning

For scoring applications,
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
targets the quantity

`mean(g_i' Sigma_lambda g_i)`

where `Sigma_lambda` is the full sandwich covariance matrix of the item
parameters and `g_i` is the gradient of the estimated ability score for
response pattern `i` with respect to those item parameters.

``` r

ability_tuned <- tune_lambda_ability(
  lambda_grid = lambda_grid,
  observed = observed,
  predicted = predicted,
  generated = generated,
  target_resp = observed,
  initial_pars = human_pars,
  n_quad = 7,
  slope_upper = 4,
  control = list(maxit = 100)
)

ability_tuned$summary
#>   lambda mean_param_var mean_squared_error mean_total_risk convergence
#> 1   0.00      0.0765121                 NA       0.0765121           0
#> 2   0.25      0.1137183                 NA       0.1137183           0
#> 3   0.50      0.4609819                 NA       0.4609819           0
#> 4   0.75      2.0858900                 NA       2.0858900           0
#> 5   1.00      4.3075417                 NA       4.3075417           0
ability_tuned$best_lambda
#> [1] 0
```

For this scenario the criterion selects $`\lambda = 0`$: gradient
asymmetry (LLM discrimination is ~10% attenuated) means any
$`\lambda > 0`$ inflates item parameters, increasing scoring
uncertainty. This is the correct behavior — the human-only estimate is
optimal when the LLM model is substantially misaligned.
`slope_upper = 4` prevents convergence failures at large $`\lambda`$ by
capping discriminations at 4.

The selected calibration is available as `ability_tuned$best_fit`.

``` r

ability_tuned$best_fit$item_pars
#>    item         a           d           b
#> 1 Item1 0.7361200 -1.18833985  1.61432899
#> 2 Item2 1.0936166 -0.78032833  0.71353010
#> 3 Item3 0.9518657  0.08528796 -0.08960084
#> 4 Item4 1.0296401  0.66922356 -0.64995873
#> 5 Item5 1.4362114  1.04603021 -0.72832607
```

You can also inspect the covariance and ability-risk components
directly.

``` r

Sigma <- vcov(ability_tuned$best_fit)
dim(Sigma)
#> [1] 10 10

risk <- ability_risk(
  resp = observed,
  fit_or_pars = ability_tuned$best_fit,
  vcov = Sigma
)

risk$summary
#>   mean_param_var mean_squared_error mean_total_risk
#> 1      0.0765121                 NA       0.0765121
```

`mean_squared_error` is `NA` because no `theta_true` was supplied. This
column is only populated in simulation studies where true ability is
known. `mean_param_var` is the propagated item-parameter uncertainty and
is the relevant quantity for calibration studies with real data.

In simulation studies, pass `theta_true` to include squared
ability-estimation error in addition to propagated item-parameter
uncertainty.

``` r

ability_tuned_with_truth <- tune_lambda_ability(
  lambda_grid = lambda_grid,
  observed = observed,
  predicted = predicted,
  generated = generated,
  target_resp = observed,
  theta_true = theta_human,
  initial_pars = human_pars,
  n_quad = 7,
  slope_upper = 4,
  control = list(maxit = 100)
)

ability_tuned_with_truth$summary
#>   lambda mean_param_var mean_squared_error mean_total_risk convergence
#> 1   0.00      0.0765121           4.779969        4.856481           0
#> 2   0.25      0.1137183           4.617626        4.731344           0
#> 3   0.50      0.4609819           4.565434        5.026416           0
#> 4   0.75      2.0858900           4.577775        6.663665           0
#> 5   1.00      4.3075417           4.593427        8.900968           0
```

For real data, `theta_true` is normally unavailable. The default risk is
then the propagated calibration component, `g_i' Sigma_lambda g_i`.

## Cross-Fitted Ability-Risk Tuning

[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
uses the same labeled responses to tune and fit. That is useful for
prototyping, but final finite-sample analyses should avoid using a
fold’s labels to tune the lambda applied to that fold’s paired
correction.

[`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md)
estimates lambda on training folds and then fits the final split-sample
estimator with the resulting fold-specific lambdas.

``` r

split_id <- rep(1:2, length.out = nrow(observed))

crossfit_tuned <- tune_lambda_ability_crossfit(
  lambda_grid = c(0, 0.5, 1),
  observed = observed,
  predicted = predicted,
  generated = generated,
  split_id = split_id,
  initial_pars = human_pars,
  n_quad = 7,
  slope_upper = 4,
  control = list(maxit = 100)
)

crossfit_tuned$lambda_by_split
#> [1] 0 0
crossfit_tuned$final_fit
#> mixedsubjectsirt 2PL fit
#>   items:      5
#>   lambda:     00
#>   loss:       2.58288
#>   convergence: 0 
#>   splits:     2
```

Cross-fitting also selects $`\lambda = 0`$ for both folds for the same
reason as `tune_lambda_ability`: gradient asymmetry from LLM-human
misalignment means the correction term hurts estimation at any positive
$`\lambda`$, and this shows up on training folds as well as the full
sample.

The final fit stores both the fold-specific lambdas and their split-size
weighted mean for the generated term.

``` r

crossfit_tuned$final_fit$split$lambda_by_split
#> [1] 0 0
crossfit_tuned$final_fit$lambda_generated
#> [1] 0
```

## Choosing a Procedure

| Procedure | Objective | When to use |
|----|----|----|
| [`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md) | Sensitivity diagnostic | Exploratory analysis only — not a valid tuning rule |
| `tune_lambda_ppi_score()` | Minimises $`\text{Tr}(\Sigma_\gamma)`$ (item-parameter variance) | Method validation; theoretical benchmarking |
| [`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md) | Minimises $`\mathbb{E}[g' \Sigma_\gamma g]`$ (ability-score risk) | **Recommended default** for psychometric applications |
| [`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md) | Same as above, cross-fitted | Final analyses requiring out-of-sample lambda estimation |

Use a fixed lambda when the value is determined by design or by a
simulation study. Use
[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
to understand sensitivity, not to make a final inferential choice. Use
`tune_lambda_ppi_score()` to inspect the theoretical PPI++ optimum as a
method diagnostic — note that this minimises item-parameter variance,
not scoring accuracy. **Use
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
when the practical target is ability scoring** and you want the lambda
that minimizes propagated scoring risk on a chosen target population.
Use
[`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md)
for the same target when finite-sample adaptivity matters and you want
lambda estimated out of sample.

The target population matters for ability-risk tuning.
`target_resp = observed` tunes for the observed human response patterns.
In operational scoring, you may prefer a larger target matrix
representing the population that will actually be scored.
