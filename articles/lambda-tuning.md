# Choosing Lambda in Mixed-Subjects IRT

The mixed-subjects estimator fits a 2PL calibration by minimizing

`L_human + lambda * (L_generated - L_paired_llm)`

The tuning parameter `lambda` controls how strongly generated LLM
responses are used after subtracting the paired LLM correction. Setting
`lambda = 0` gives the human-only expected-count calibration. Larger
values use more information from the generated response sample, but they
can also amplify LLM bias when the paired correction is weak.

This vignette distinguishes three package workflows:

- [`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md):
  a sensitivity diagnostic. It is not a valid final tuning rule.
- [`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md):
  chooses a lambda by minimizing propagated ability-score risk on a
  target scoring population.
- [`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md):
  estimates lambda out of sample by fold, then fits the final
  split-sample estimator with fold-specific lambdas.

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
  control = list(maxit = 100)
)

fixed_fit$item_pars
#>    item            a             d             b
#> 1 Item1 1.408295e+37 -9.683049e+36  0.6875726503
#> 2 Item2 2.255703e+38 -1.994754e+37  0.0884316076
#> 3 Item3 2.622970e+36  3.615787e+34 -0.0137850861
#> 4 Item4 4.408536e+37  3.190817e+35 -0.0072378151
#> 5 Item5 9.703142e+38  4.395959e+35 -0.0004530449
```

Use this route when a simulation study, design analysis, or
preregistered plan already specifies `lambda`.

## Sensitivity Diagnostics

[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
fits candidate lambda values and reports the mixed objective and the
observed-human expected-count loss. This is useful for seeing how
sensitive the fitted calibration is to lambda.

``` r

diagnostic <- diagnose_lambda_grid(
  lambda_grid = lambda_grid,
  observed = observed,
  predicted = predicted,
  generated = generated,
  initial_pars = human_pars,
  n_quad = 7,
  control = list(maxit = 100)
)

diagnostic$summary
#>   lambda    mixed_loss observed_loss convergence
#> 1   0.00  2.582875e+00  2.582875e+00           0
#> 2   0.25  2.303072e+00  2.640973e+00           0
#> 3   0.50 -2.135974e+37  1.251187e+38           1
#> 4   0.75 -1.103359e+28  2.255674e+28           1
#> 5   1.00 -3.073320e+36  3.409178e+36           1
diagnostic$lowest_observed_loss_lambda
#> [1] 0
```

The value in `lowest_observed_loss_lambda` is not a final PPI++ tuning
rule. It is the lambda with the smallest observed-human expected-count
loss after fitting. That criterion can overfit the labeled sample and is
not the same as minimizing downstream ability-score uncertainty.

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
  control = list(maxit = 100)
)

ability_tuned$summary
#>   lambda mean_param_var mean_squared_error mean_total_risk convergence
#> 1   0.00      0.0765121                 NA       0.0765121           0
#> 2   0.25      0.1137183                 NA       0.1137183           0
#> 3   0.50            NaN                 NA             NaN           1
#> 4   0.75            NaN                 NA             NaN           1
#> 5   1.00            NaN                 NA             NaN           1
ability_tuned$best_lambda
#> [1] 0
```

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
  control = list(maxit = 100)
)

ability_tuned_with_truth$summary
#>   lambda mean_param_var mean_squared_error mean_total_risk convergence
#> 1   0.00      0.0765121           4.779969        4.856481           0
#> 2   0.25      0.1137183           4.617626        4.731344           0
#> 3   0.50            NaN           1.359333             NaN           1
#> 4   0.75            NaN           1.412640             NaN           1
#> 5   1.00            NaN           1.145156             NaN           1
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

The final fit stores both the fold-specific lambdas and their split-size
weighted mean for the generated term.

``` r

crossfit_tuned$final_fit$split$lambda_by_split
#> [1] 0 0
crossfit_tuned$final_fit$lambda_generated
#> [1] 0
```

## Choosing a Procedure

Use a fixed lambda when the value is determined by design or by a
simulation study. Use
[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
to understand sensitivity, not to make a final inferential choice. Use
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
when the practical target is ability scoring and you want the lambda
that minimizes propagated scoring risk on a chosen target population.
Use
[`tune_lambda_ability_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_crossfit.md)
for the same target when finite-sample adaptivity matters and you want
lambda estimated out of sample.

The target population matters. `target_resp = observed` tunes for the
observed human response patterns. In operational scoring, you may prefer
a larger target matrix representing the population that will actually be
scored.
