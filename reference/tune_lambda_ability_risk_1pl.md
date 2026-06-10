# Tune lambda by downstream ability-score risk for a 1PL model

Grid-searches lambda values using
[`fit_mixed_subjects_mml_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md)
(by default) and selects the value minimizing `E[g' Sigma_1pl g]` — the
propagated ability-score risk in the 1PL parameterization.

## Usage

``` r
tune_lambda_ability_risk_1pl(
  lambda_grid,
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
  n_quad = 31,
  initial_pars = NULL,
  fit_fn = fit_mixed_subjects_mml_1pl,
  bounds = c(-6, 6),
  max_discrimination = 10,
  control = list(maxit = 500),
  ...
)
```

## Arguments

- lambda_grid:

  Numeric vector of candidate lambda values in `[0, 1]`.

- observed, predicted, generated:

  Response matrices passed to
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).

- target_resp:

  Response matrix defining the target scoring population. If omitted,
  `observed` is used.

- theta_true:

  Optional true theta values for `target_resp`, used in simulation
  studies to add squared scoring error to the risk. When omitted,
  `mean_squared_error` in the summary is `NA`; only `mean_param_var` is
  computed.

- n_quad:

  Number of quadrature nodes.

- initial_pars:

  Optional starting item parameters.

- fit_fn:

  Fitting function. Defaults to
  [`fit_mixed_subjects_mml_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md).

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- max_discrimination:

  Upper bound on plausible item discrimination. Any candidate fit whose
  maximum `|a|` exceeds this value is treated as degenerate and excluded
  from selection. This guards against runaway discrimination fits, which
  can "converge" with a spuriously low model-based risk (huge
  discrimination collapses the item-parameter covariance). The default
  of 10 is far above any realistic 2PL discrimination.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to `fit_fn`.

## Value

A list with `summary`, `best_lambda`, `best_fit`, `fits`, `risks`.

## Details

Passes `fit_fn` to allow switching between the frozen expected-count
estimator
([`fit_mixed_subjects_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_1pl.md))
and the marginal-MML estimator
([`fit_mixed_subjects_mml_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md)).

## See also

[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
for the 2PL version;
[`tune_lambda_ppi_score_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score_1pl.md)
for the PPI++ score diagnostic.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
obs  <- simulate_2pl(rnorm(40), pars)
gen  <- simulate_2pl(rnorm(100), pars)
tuned <- tune_lambda_ability_risk_1pl(
  c(0, 0.5), obs, obs, gen,
  initial_pars = pars, n_quad = 5, control = list(maxit = 30)
)
tuned$best_lambda
#> [1] 0.5
```
