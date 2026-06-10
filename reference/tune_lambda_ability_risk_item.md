# Per-item ability-risk lambda tuning via coordinate descent

Finds a per-item vector of lambda values `λ_j ∈ [0, 1]` that minimizes
propagated ability-score risk `E[g' Σ_γ g]` using coordinate descent on
the items. Each coordinate step selects the `λ_j` in `lambda_grid` that
gives the smallest mean ability risk while holding all other `λ_{j'}`
fixed.

## Usage

``` r
tune_lambda_ability_risk_item(
  lambda_grid,
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
  n_quad = 31,
  initial_pars = NULL,
  n_pass = 1,
  init_lambda = 0,
  bounds = c(-6, 6),
  max_discrimination = 10,
  control = list(maxit = 300),
  ...
)
```

## Arguments

- lambda_grid:

  Numeric vector of candidate λ values in `[0, 1]` to try for each item
  independently.

- observed, predicted, generated:

  Response matrices passed to
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md).

- target_resp:

  Target scoring population. If omitted, `observed` is used.

- theta_true:

  Optional true theta values, used to add squared scoring error to the
  risk.

- n_quad:

  Number of quadrature nodes.

- initial_pars:

  Optional starting item parameters.

- n_pass:

  Number of coordinate-descent passes (default 1).

- init_lambda:

  Starting lambda vector for coordinate descent. Supply the global
  scalar optimum from
  [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  (e.g. `init_lambda = 0.5`) to start the search around a useful
  operating point. Starting from all-zeros is not recommended: each
  single-item improvement is too small to detect when other items are at
  zero. A scalar is broadcast to all items; a vector of length `n_items`
  sets per-item starting values.

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- max_discrimination:

  Upper bound on plausible item discrimination; any candidate fit whose
  maximum `|a|` exceeds it is treated as degenerate and skipped. See
  [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  for the rationale. Default 10.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md).

## Value

A list with `lambda` (per-item vector), `item` (item names), `n_pass`,
and `final_fit` (the
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
fit at the selected lambda).

## Details

Calls
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
with a per-item lambda vector at each candidate evaluation. Because the
lambda is a vector, that function **switches to its frozen
expected-count Q-function path** — posteriors are frozen at
`initial_pars`, not recomputed continuously. This is an approximation;
see the `@note` below. The resulting lambda vector can be used directly
with
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md).

**Computational cost.** Each pass evaluates
`n_items × length(lambda_grid)` fits. For `n_items = 8` and a 5-point
grid this is 40 fits per pass. Use `n_pass = 1` (the default) for a
single greedy sweep, which is usually sufficient.

## Note

**Approximation status.** The coordinate descent fits use the frozen
expected-count Q-function (not the full marginal-MML objective) because
the IRT marginal likelihood integrates over the joint response pattern
and does not decompose item-wise. The approach is approximately correct
when `initial_pars` is close to the converged parameters. Report
per-item results as experimental / approximate.

## See also

[`tune_lambda_ppi_score_item()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score_item.md)
for the faster PPI++-score version;
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
for the global scalar version.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed  <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
tuned <- tune_lambda_ability_risk_item(
  c(0, 0.5), observed, observed, generated,
  initial_pars = pars, n_quad = 5, control = list(maxit = 30)
)
tuned$lambda
#> [1] 0.5 0.5 0.5
```
