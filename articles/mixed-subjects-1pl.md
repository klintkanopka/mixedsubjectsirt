# Mixed-Subjects 1PL Calibration

The 1PL (one-parameter logistic) model estimates a single shared
discrimination $`a`$ across all items together with per-item intercepts
$`d_j`$:

``` math
P(x_j = 1 \mid \theta) = \text{logistic}(a\,\theta + d_j)
```

The parameter vector has length $`J+1`$ rather than $`2J`$. The package
provides exact analogues of the 2PL mixed-subjects functions for the 1PL
case.

**When to prefer 1PL over 2PL:** - Ability-focused tests where the items
are designed to be equally discriminating. - Tests built from a single
item pool with homogeneous item characteristics. - When the 2PL
discrimination estimates are very noisy (small $`n`$).

> **Note on vcov.**
> [`vcov_mixed_subjects_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_1pl.md)
> currently uses the EM complete-data Hessian (not Louis’
> marginal-information correction). The uncertainty estimates are
> slightly over-precise. A Louis-corrected 1PL bread is planned for a
> future release.

## Simulate a 1PL test

``` r

library(mixedsubjectsirt)
library(ggplot2)

set.seed(2026)

n_human     <- 400
n_generated <- 1200
n_items     <- 8

# True 1PL: shared discrimination a = 1.2, varying difficulties
true_1pl <- data.frame(
  item = paste0("Item", seq_len(n_items)),
  a    = 1.2,
  d    = seq(-1.1, 1.1, length.out = n_items)
)
true_1pl$b <- -true_1pl$d / true_1pl$a

theta_human <- rnorm(n_human)
observed    <- simulate_2pl(theta_human, true_1pl)

# LLM: same 1PL structure, small intercept shift
llm_1pl   <- true_1pl
llm_1pl$d <- true_1pl$d + 0.25
llm_1pl$b <- -llm_1pl$d / llm_1pl$a

predicted <- simulate_2pl(theta_human, llm_1pl)
generated <- simulate_2pl(rnorm(n_generated), llm_1pl)
```

## Step 1: Fit the 1PL baseline

[`fit_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
estimates $`a`$ and $`d_1, \ldots, d_J`$ by maximizing the IRT marginal
likelihood under a standard-normal ability prior.

``` r

fit1 <- fit_1pl(observed, n_quad = 15)
cat("Shared a:", round(fit1$pars$a[1], 3), " (true:", true_1pl$a[1], ")\n")
#> Shared a: 1.287  (true: 1.2 )
cat("Convergence:", fit1$convergence, "\n\n")
#> Convergence: 0
fit1$pars
#>    item        a          d           b
#> 1 Item1 1.286555 -1.1603124  0.90187542
#> 2 Item2 1.286555 -0.8618391  0.66988122
#> 3 Item3 1.286555 -0.1642846  0.12769343
#> 4 Item4 1.286555 -0.4177364  0.32469371
#> 5 Item5 1.286555  0.1133618 -0.08811268
#> 6 Item6 1.286555  0.5983198 -0.46505573
#> 7 Item7 1.286555  0.8240620 -0.64051824
#> 8 Item8 1.286555  1.3131182 -1.02064664
```

All items in the output have the same `a` value, confirming the 1PL
constraint.

## Step 2: Fit mixed-subjects MML (1PL)

[`fit_mixed_subjects_mml_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md)
uses the true marginal likelihood with a 1PL-specific gradient: the
shared discrimination gradient accumulates contributions from all $`J`$
items, while each intercept has its own gradient.

``` r

fit_mml_1pl <- fit_mixed_subjects_mml_1pl(
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  lambda       = 0.5,
  initial_pars = fit1$pars,
  n_quad       = 15,
  control      = list(maxit = 300)
)

print(fit_mml_1pl)
#> mixedsubjectsirt 1PL fit
#>   items:      8
#>   a (shared): 1.327
#>   lambda:     0.5
#>   loss:       4.72331
#>   convergence: 0 
#>   estimator:  marginal MML PPI++ (1PL)
fit_mml_1pl$item_pars
#>    item        a          d          b
#> 1 Item1 1.327115 -1.2025380  0.9061297
#> 2 Item2 1.327115 -0.8728234  0.6576850
#> 3 Item3 1.327115 -0.1832574  0.1380871
#> 4 Item4 1.327115 -0.4426497  0.3335429
#> 5 Item5 1.327115  0.2317895 -0.1746567
#> 6 Item6 1.327115  0.6473713 -0.4878036
#> 7 Item7 1.327115  0.8889582 -0.6698428
#> 8 Item8 1.327115  1.5297000 -1.1526510
```

## Step 3: Correct covariance — $`(J+1) \times (J+1)`$ sandwich

[`vcov()`](https://rdrr.io/r/stats/vcov.html) dispatches to
[`vcov_mixed_subjects_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_1pl.md)
for 1PL fits, returning a $`(J+1) \times (J+1)`$ matrix with `a_shared`
and per-item `d_j` as rows/columns.

``` r

Sigma_1pl <- vcov(fit_mml_1pl)
dim(Sigma_1pl)
#> [1] 9 9
rownames(Sigma_1pl)
#> [1] "a_shared" "d_Item1"  "d_Item2"  "d_Item3"  "d_Item4"  "d_Item5"  "d_Item6" 
#> [8] "d_Item7"  "d_Item8"
```

## Step 4: Ability-score risk and lambda tuning

[`tune_lambda_ability_risk_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_1pl.md)
uses the 1PL-parameterized gradient
$`\partial\hat\theta / \partial (a_\text{shared}, d_1, \ldots, d_J)`$
for the ability-score risk. The chain rule gives
$`\partial\hat\theta / \partial a_\text{shared} = \sum_j \partial\hat\theta / \partial a_j`$.

``` r

tuned_1pl <- tune_lambda_ability_risk_1pl(
  lambda_grid  = seq(0, 1, by = 0.2),
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  initial_pars = fit1$pars,
  n_quad       = 15,
  control      = list(maxit = 300)
)

tuned_1pl$summary[, c("lambda", "mean_param_var", "mean_total_risk",
                       "convergence")]
#>   lambda mean_param_var mean_total_risk convergence
#> 1    0.0    0.003293615     0.003293615           0
#> 2    0.2    0.002997535     0.002997535           0
#> 3    0.4    0.003032710     0.003032710           0
#> 4    0.6    0.003395853     0.003395853           0
#> 5    0.8    0.004086590     0.004086590           0
#> 6    1.0    0.005108417     0.005108417           0
tuned_1pl$best_lambda
#> [1] 0.2
```

## Step 5: Verify — F = Y gives lambda \> 0

With `predicted = observed` (perfect paired predictor), the ability-risk
criterion should select a positive lambda.

``` r

tuned_fy <- tune_lambda_ability_risk_1pl(
  lambda_grid  = seq(0, 1, by = 0.2),
  observed     = observed,
  predicted    = observed,     # F = Y
  generated    = simulate_2pl(rnorm(n_generated), true_1pl),
  initial_pars = fit1$pars,
  n_quad       = 15,
  control      = list(maxit = 300)
)

cat("F=Y best lambda:", tuned_fy$best_lambda,
    " (theory: N/(n+N) =", round(n_generated / (n_human + n_generated), 3), ")\n")
#> F=Y best lambda: 0.8  (theory: N/(n+N) = 0.75 )
tuned_fy$summary[, c("lambda", "mean_param_var")]
#>   lambda mean_param_var
#> 1    0.0   0.0032936154
#> 2    0.2   0.0022043602
#> 3    0.4   0.0014197985
#> 4    0.6   0.0009634564
#> 5    0.8   0.0008625318
#> 6    1.0   0.0011476931
```

## Compare 1PL and 2PL

On a well-specified 1PL test, how do the 1PL and 2PL estimators compare?

``` r

fit_2pl_mml <- fit_mixed_subjects_mml(
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  lambda       = tuned_1pl$best_lambda,
  initial_pars = fit_2pl(observed, technical = list(NCYCLES = 500))$pars,
  n_quad       = 15,
  control      = list(maxit = 300)
)

rmse <- function(x, y) sqrt(mean((x - y)^2))
cat("1PL RMSE(a):", round(rmse(tuned_1pl$best_fit$item_pars$a, true_1pl$a), 4), "\n")
#> 1PL RMSE(a): 0.1021
cat("2PL RMSE(a):", round(rmse(fit_2pl_mml$item_pars$a, true_1pl$a), 4), "\n")
#> 2PL RMSE(a): 0.292

# Difficulty recovery
cat("1PL RMSE(d):", round(rmse(tuned_1pl$best_fit$item_pars$d, true_1pl$d), 4), "\n")
#> 1PL RMSE(d): 0.1901
cat("2PL RMSE(d):", round(rmse(fit_2pl_mml$item_pars$d, true_1pl$d), 4), "\n")
#> 2PL RMSE(d): 0.196
```

The 1PL uses fewer parameters ($`J+1`$ vs $`2J`$), which can give lower
RMSE on a test generated from a true 1PL DGP — especially for small
$`n`$.

## Ability-score risk: 1PL vs 2PL parameterization

The 1PL ability-score risk is smaller in the $`(J+1)`$-parameter space
because the shared $`a`$ concentrates all discrimination information in
a single parameter.

``` r

Sigma_2pl <- vcov(fit_2pl_mml)  # 2J × 2J Louis-corrected

risk_1pl <- ability_risk_1pl(observed, tuned_1pl$best_fit)
risk_2pl <- ability_risk(observed, fit_2pl_mml, vcov = Sigma_2pl)

cat("1PL mean param_var:", round(risk_1pl$summary$mean_param_var, 5), "\n")
#> 1PL mean param_var: 0.003
cat("2PL mean param_var:", round(risk_2pl$summary$mean_param_var, 5), "\n")
#> 2PL mean param_var: 0.04243
```
