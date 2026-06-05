# Per-Item Lambda (Experimental)

> **Experimental.** Per-item lambda uses a frozen expected-count
> approximation, not the full marginal-MML objective. The IRT marginal
> likelihood integrates the *joint* response pattern and does not
> decompose item-wise, so a theoretically principled per-item marginal
> objective does not yet exist in this package. Results from per-item
> tuning should be treated as approximate and validated against the
> scalar-lambda MML baseline. See the documentation of
> `tune_lambda_ability_risk_item()` for details.

## Why per-item lambda?

A single global $`\lambda`$ is fragile: if even one item has
poorly-correlated LLM predictions, the ability-risk criterion may force
the global optimum to $`\lambda = 0`$, preventing all items from
benefiting. Per-item $`\lambda_j`$ allows each item to draw on the LLM
data at its own optimal level.

Consider a test with eight items where: - Items 1–4 are straightforward
factual questions — the LLM predicts these well. - Items 5–8 require
nuanced reasoning — the LLM is nearly random on these.

With a scalar $`\lambda`$, the four poor items push the optimum toward
0. With per-item $`\lambda_j`$, items 1–4 get
$`\lambda_j \approx 0.5`$–0.7 while items 5–8 get
$`\lambda_j \approx 0`$.

## Simulate a heterogeneous test

``` r

library(mixedsubjectsirt)
library(ggplot2)

set.seed(2026)

n_human    <- 400
n_generated <- 1200
n_items    <- 8
n_good     <- 4   # items where LLM predicts well

true_pars <- data.frame(
  item = paste0("Item", seq_len(n_items)),
  a    = seq(0.8, 1.6, length.out = n_items),
  d    = seq(-1.1, 1.1, length.out = n_items)
)
true_pars$b <- -true_pars$d / true_pars$a

theta_human <- rnorm(n_human)
observed    <- simulate_2pl(theta_human, true_pars)

# LLM: good for items 1–4 (same DGP), poor for 5–8 (random noise)
llm_pars_good <- true_pars
llm_pars_poor <- true_pars
llm_pars_poor$a <- pmax(0.05, rnorm(n_items, 0, 0.1))  # near-random
llm_pars_poor$d <- rnorm(n_items, 0, 2)
llm_pars_poor$b <- -llm_pars_poor$d / llm_pars_poor$a

# Build predicted (same subjects as human)
predicted <- observed   # F = Y for first 4 items
predicted[, 5:8] <- simulate_2pl(theta_human, llm_pars_poor)[, 5:8]

# Build generated
generated_good <- simulate_2pl(rnorm(n_generated), true_pars)
generated_poor <- simulate_2pl(rnorm(n_generated), llm_pars_poor)
generated <- cbind(generated_good[, 1:4], generated_poor[, 5:8])
colnames(generated) <- true_pars$item
```

The first four items have perfect paired predictions (F = Y); items 5–8
have near-random LLM predictions.

## Step 1: Fit 2PL baseline and get global scalar lambda

``` r

human_pars <- fit_2pl(observed, technical = list(NCYCLES = 500))$pars

global_tuned <- tune_lambda_ability_risk(
  lambda_grid  = seq(0, 1, by = 0.1),
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  initial_pars = human_pars,
  fit_fn       = fit_mixed_subjects_mml,
  n_quad       = 11,
  control      = list(maxit = 200)
)

cat("Global scalar best lambda:", global_tuned$best_lambda, "\n")
#> Global scalar best lambda: 0.2
```

The global scalar is forced to a compromise — the four poor items
constrain it to a value smaller than what items 1–4 could support.

## Step 2: PPI++ score per item (fast diagnostic)

`tune_lambda_ppi_score_item()` applies the Proposition 2 formula
independently per item using the 2×2 diagonal block of $`H^{-1}`$ and
the item-level sub-vectors of the score matrices. This is fast (no
fitting required) and shows which items are well-predicted.

``` r

ppi_item <- tune_lambda_ppi_score_item(
  observed    = observed,
  predicted   = predicted,
  item_pars   = human_pars,
  n_generated = n_generated,
  n_quad      = 11
)

cat("Per-item PPI++ lambda:\n")
#> Per-item PPI++ lambda:
print(data.frame(item = ppi_item$item, lambda = round(ppi_item$lambda, 3)))
#>    item lambda
#> 1 Item1   0.75
#> 2 Item2   0.75
#> 3 Item3   0.75
#> 4 Item4   0.75
#> 5 Item5   0.00
#> 6 Item6   0.00
#> 7 Item7   0.00
#> 8 Item8   0.00
cat("N/(n+N) upper bound:", round(n_generated / (n_human + n_generated), 3), "\n")
#> N/(n+N) upper bound: 0.75
```

Items 1–4 (F = Y) should show $`\lambda_j \approx N/(n+N) = 0.75`$;
items 5–8 (random LLM) should show $`\lambda_j \approx 0`$.

## Step 3: Per-item ability-risk tuning

`tune_lambda_ability_risk_item()` uses coordinate descent: for each item
$`j`$, it finds the $`\lambda_j`$ in `lambda_grid` that minimises
ability-score risk while holding all other $`\lambda_{j'}`$ fixed. Each
evaluation fits with the **frozen expected-count Q-function** (not the
full marginal-MML objective) because the IRT marginal likelihood does
not decompose item-wise. Starting from the global scalar optimum (not
from all-zeros) is essential — see the note below.

``` r

item_tuned <- tune_lambda_ability_risk_item(
  lambda_grid  = seq(0, 1, by = 0.25),
  observed     = observed,
  predicted    = predicted,
  generated    = generated,
  initial_pars = human_pars,
  init_lambda  = global_tuned$best_lambda,   # start from global best
  n_quad       = 11,
  n_pass       = 1,
  control      = list(maxit = 200)
)

cat("Per-item ability-risk lambda:\n")
#> Per-item ability-risk lambda:
print(data.frame(item = item_tuned$item, lambda = round(item_tuned$lambda, 3)))
#>    item lambda
#> 1 Item1   0.75
#> 2 Item2   0.75
#> 3 Item3   0.75
#> 4 Item4   0.75
#> 5 Item5   0.00
#> 6 Item6   0.00
#> 7 Item7   0.00
#> 8 Item8   0.00
```

Items 1–4 should receive positive $`\lambda_j`$ (good predictor); items
5–8 should be near zero (poor predictor).

## Step 4: Compare scalar vs. per-item parameter recovery

``` r

fit_scalar   <- global_tuned$best_fit
fit_per_item <- item_tuned$final_fit

rmse <- function(x, y) sqrt(mean((x - y)^2))

comparison <- data.frame(
  item    = true_pars$item,
  true_a  = round(true_pars$a, 3),
  human_a = round(human_pars$a, 3),
  scalar_a = round(fit_scalar$item_pars$a, 3),
  item_a   = if (is.null(fit_per_item)) NA_real_ else
    round(fit_per_item$item_pars$a, 3)
)
knitr::kable(comparison, row.names = FALSE,
  caption = "Discrimination recovery: scalar lambda vs. per-item lambda")
```

| item  | true_a | human_a | scalar_a | item_a |
|:------|-------:|--------:|---------:|-------:|
| Item1 |  0.800 |   0.727 |    0.523 |  0.633 |
| Item2 |  0.914 |   1.256 |    0.856 |  0.863 |
| Item3 |  1.029 |   1.117 |    0.793 |  0.835 |
| Item4 |  1.143 |   1.091 |    0.754 |  0.902 |
| Item5 |  1.257 |   1.758 |    1.187 |  1.586 |
| Item6 |  1.371 |   1.610 |    1.113 |  1.449 |
| Item7 |  1.486 |   1.369 |    1.002 |  1.187 |
| Item8 |  1.600 |   1.925 |    1.453 |  1.683 |

Discrimination recovery: scalar lambda vs. per-item lambda {.table}

``` r


cat("RMSE(a) human-only:   ",
    round(rmse(human_pars$a, true_pars$a), 4), "\n")
#> RMSE(a) human-only:    0.2645
cat("RMSE(a) scalar MML:   ",
    round(rmse(fit_scalar$item_pars$a, true_pars$a), 4), "\n")
#> RMSE(a) scalar MML:    0.277
if (!is.null(fit_per_item)) {
  cat("RMSE(a) per-item MML: ",
      round(rmse(fit_per_item$item_pars$a, true_pars$a), 4), "\n")
}
#> RMSE(a) per-item MML:  0.205
```

## Important note on initialisation

Starting coordinate descent from all-zeros is not recommended. When all
other items are at $`\lambda_j = 0`$, each single-item improvement is
diluted across the full ability-risk criterion, making improvements hard
to detect. The recommended workflow is:

1.  Fit the scalar global optimum with
    `tune_lambda_ability_risk(..., fit_fn = fit_mixed_subjects_mml)`.
2.  Pass that global value as `init_lambda` to
    `tune_lambda_ability_risk_item()`.
3.  The coordinate descent then identifies which items should deviate
    from the global default — raising $`\lambda_j`$ for well-predicted
    items and lowering it for poorly-predicted ones.

## Approximation caveat

The per-item lambda coordinate descent uses the **frozen Q-function**
(not the full marginal-MML objective) for each candidate evaluation.
This is necessary because the IRT marginal likelihood integrates the
joint response pattern and does not decompose item-wise. The
approximation is good when `initial_pars` is close to the converged
parameters. For final reporting, always:

- Compare the per-item fit against the scalar MML baseline.
- Report per-item results as approximate / experimental.
- Use `vcov(item_tuned$final_fit)` for uncertainty (which applies
  `vcov_mixed_subjects` with the vector-lambda bread and meat scaling).
