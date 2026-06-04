# Mixed-Subjects IRT Calibration

This vignette shows the basic mixed-subjects workflow for a
unidimensional 2PL model. The package expects three response matrices
with the same item columns:

- `observed`: human responses.
- `predicted`: LLM responses or probabilities for the same human rows.
- `generated`: additional generated or unlabeled LLM responses.

The fitted objective is

`L_human + lambda * (L_generated - L_paired_llm)`

Setting `lambda = 0` gives the human-only expected-count calibration.

## Simulate example data

This example uses simulated LLM responses so the vignette can run
without any external API calls.

``` r

library(mixedsubjectsirt)
library(ggplot2)

set.seed(2026)

n_human <- 400
n_generated <- 1200
n_items <- 8

true_pars <- data.frame(
  item = paste0("Item", seq_len(n_items)),
  a = seq(0.8, 1.6, length.out = n_items),
  d = seq(-1.1, 1.1, length.out = n_items)
)
true_pars$b <- -true_pars$d / true_pars$a

theta_human <- rnorm(n_human)
observed <- simulate_2pl(theta_human, true_pars)

llm_pars <- true_pars
llm_pars$a <- pmax(0.4, 0.9 * true_pars$a + rnorm(n_items, 0, 0.05))
llm_pars$d <- true_pars$d + 0.25 + rnorm(n_items, 0, 0.15)
llm_pars$b <- -llm_pars$d / llm_pars$a

predicted <- simulate_2pl(theta_human, llm_pars)
generated <- simulate_2pl(rnorm(n_generated, mean = 0.1, sd = 1.05), llm_pars)
```

## Workflow without split samples

The no-split workflow uses all human responses to construct the human
posterior weights and uses those same posterior weights for the paired
LLM responses. That default is controlled by
`common_predicted_weights = TRUE`.

``` r

human_start <- fit_2pl(
  observed,
  technical = list(NCYCLES = 500)
)

human_only <- fit_mixed_subjects(
  observed = observed,
  predicted = predicted,
  generated = generated,
  lambda = 0,
  initial_pars = human_start$pars,
  n_quad = 11
)

mixed_no_split <- fit_mixed_subjects(
  observed = observed,
  predicted = predicted,
  generated = generated,
  lambda = 0.5,
  initial_pars = human_start$pars,
  n_quad = 11,
  slope_upper = 4   # cap discriminations; gradient asymmetry drives them up at lambda > 0
)

mixed_no_split
#> mixedsubjectsirt 2PL fit
#>   items:      8
#>   lambda:     0.5
#>   loss:       3.46725
#>   convergence: 0
```

The fitted item parameters are returned in slope-intercept form.

``` r

head(mixed_no_split$item_pars)
#>    item         a           d            b
#> 1 Item1 0.8347858 -1.26460926  1.514890762
#> 2 Item2 1.5692739 -1.16988198  0.745492554
#> 3 Item3 1.2331265 -0.20684323  0.167738861
#> 4 Item4 1.3177043 -0.49184243  0.373257052
#> 5 Item5 4.0000000 -0.02148907  0.005372268
#> 6 Item6 2.1663434  0.57442082 -0.265156871
```

## Workflow with split samples

The split workflow cross-fits the posterior weights used for the labeled
human and paired LLM terms. For each split, the initial human 2PL model
is estimated on the other splits, then used to summarize the held-out
split. Every human row still contributes to the final estimating
equation once.

``` r

mixed_split <- fit_mixed_subjects_split(
  observed = observed,
  predicted = predicted,
  generated = generated,
  lambda = 0.5,
  n_splits = 2,
  seed = 2026,
  n_quad = 11,
  slope_upper = 4,
  technical = list(NCYCLES = 500)
)

mixed_split
#> mixedsubjectsirt 2PL fit
#>   items:      8
#>   lambda:     0.5
#>   loss:       3.485
#>   convergence: 0 
#>   splits:     2
```

## Compare the calibrations

In a real calibration study the true item parameters are unknown. In
this simulation they are useful for checking that the workflow is
behaving as expected.

``` r

estimates <- rbind(
  data.frame(estimator = "human only", human_only$item_pars),
  data.frame(estimator = "mixed no split", mixed_no_split$item_pars),
  data.frame(estimator = "mixed split", mixed_split$item_pars)
)
estimates$true_b <- true_pars$b[match(estimates$item, true_pars$item)]
estimates$true_a <- true_pars$a[match(estimates$item, true_pars$item)]

ggplot(estimates, aes(true_b, b, color = estimator)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  geom_point(size = 2) +
  labs(x = "True difficulty", y = "Estimated difficulty", color = NULL) +
  theme_minimal()
```

![](mixed-subjects-workflow_files/figure-html/compare-1.png)

The helper
[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
can be used for quick sensitivity checks over candidate lambda values.
The summary reports the optimized mixed objective and the human
expected-count loss, but final lambda selection should be tied to the
study’s inferential target or evaluated with a
bootstrap/cross-validation plan.

When the LLM has attenuated discrimination relative to humans (as here),
the combined gradient
$`\nabla L_\text{obs} + \lambda\nabla(L_\text{gen} - L_\text{pred})`$ is
systematically negative for high-discrimination items — it pushes $`a`$
upward at any $`\lambda > 0`$. Without an upper bound, this drives
discriminations to extreme values and the optimizer reports convergence
code 52 (L-BFGS-B line-search failure). The `slope_upper` argument caps
discriminations during optimization. With it, the optimizer converges at
the bound, ability risk is high for capped items, and
`tune_lambda_ability` correctly selects $`\lambda \approx 0`$ —
recovering the human-only estimate.

``` r

tuned <- diagnose_lambda_grid(
  lambda_grid = c(0, 0.25, 0.5, 0.75),
  observed = observed,
  predicted = predicted,
  generated = generated,
  initial_pars = human_start$pars,
  n_quad = 11,
  slope_upper = 4,
  control = list(maxit = 200)
)

tuned$summary
#>   lambda mixed_loss observed_loss convergence
#> 1   0.00   4.081765      4.081765           0
#> 2   0.25   3.825618      4.117991           0
#> 3   0.50   3.467245      4.395663           0
#> 4   0.75   2.973209      4.573353           0
```

## Tune lambda by ability risk

For a target scoring population,
[`tune_lambda_ability()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability.md)
evaluates candidate calibrations using propagated ability-score risk.
For each candidate lambda it fits item parameters, estimates the full
sandwich covariance matrix for the item parameters, computes
ability-score gradients, and summarizes the average `g' Sigma g` risk.

For this scenario the ability-risk criterion selects $`\lambda = 0`$:
the LLM’s attenuated discrimination causes gradient asymmetry that
inflates item parameters at any $`\lambda > 0`$, increasing scoring
uncertainty. This is the expected and correct behavior — the method
recovers the human-only estimate when the LLM model is substantially
misaligned.

``` r

ability_tuned <- tune_lambda_ability(
  lambda_grid = c(0, 0.25, 0.5, 0.75),
  observed = observed,
  predicted = predicted,
  generated = generated,
  target_resp = observed,
  initial_pars = human_start$pars,
  n_quad = 11,
  slope_upper = 4,
  control = list(maxit = 200)
)

ability_tuned$summary
#>   lambda mean_param_var mean_squared_error mean_total_risk convergence
#> 1   0.00     0.01539735                 NA      0.01539735           0
#> 2   0.25     0.02186879                 NA      0.02186879           0
#> 3   0.50     0.07231729                 NA      0.07231729           0
#> 4   0.75     0.13123481                 NA      0.13123481           0
ability_tuned$best_lambda
#> [1] 0
```
