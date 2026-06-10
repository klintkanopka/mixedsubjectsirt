# Tune lambda by downstream ability-score risk

Fits candidate mixed-subjects calibrations, estimates the item-parameter
sandwich covariance for each, and chooses the lambda that minimizes
average propagated ability-score risk on a target response matrix.

## Usage

``` r
tune_lambda_ability_risk(
  lambda_grid,
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
  n_quad = 31,
  initial_pars = NULL,
  fit_fn = fit_mixed_subjects,
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

  Fitting function to use. Defaults to
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
  (frozen expected-count). Pass
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  to use the marginal-likelihood PPI++ estimator, which eliminates
  gradient asymmetry and should select `lambda > 0` for genuinely
  informative predictors.

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

A list with `summary`, `best_lambda`, `best_fit`, and all fitted
candidate objects.

## Details

This function minimizes `E[g' Sigma_gamma g]` — the propagated
ability-score risk — which is the appropriate objective for IRT
applications where accurate test scoring is the goal. This is
**distinct** from
[`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md),
which minimizes the trace of the item-parameter covariance matrix
`Tr(Sigma_gamma)` (the PPI++ theoretical objective). The two criteria
generally yield different lambda values:

- `tune_lambda_ability_risk()` asks: which lambda produces the most
  accurate ability scores for the target population? Use this for
  operational scoring.

- [`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
  asks: which lambda minimizes item-parameter estimation variance? Use
  this for method validation and diagnostics.

Diagnostic note: if `tune_lambda_ability_risk()` selects `lambda = 0`
for a misaligned LLM (one whose item parameters differ from the human
calibration), this is the correct mathematical outcome under the current
fixed-posterior expected-count implementation. The frozen posteriors
create a gradient asymmetry that inflates item parameters at any
`lambda > 0`, increasing ability risk. This is not a bug in the risk
function; it is a property of the estimating equations. See
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
for a marginal-likelihood implementation that removes this asymmetry.

## See also

[`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
for the PPI++ theoretical lambda that minimizes the trace of the
item-parameter covariance matrix;
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
for the marginal-likelihood estimator.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
tuned <- tune_lambda_ability_risk(
  c(0, 0.5), observed, observed, generated,
  initial_pars = pars, n_quad = 5, control = list(maxit = 30)
)
tuned$best_lambda
#> [1] 0.5
```
