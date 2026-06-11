# Changelog

## mixedsubjectsirt 0.0.1 (development version)

### Lambda tuning

- [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  and
  [`tune_lambda_ability_risk_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_1pl.md)
  now select lambda by **direct 1-D optimization**
  ([`stats::optimize()`](https://rdrr.io/r/stats/optimize.html)) by
  default, returning a continuous lambda with no grid rounding. Pass
  `method = "grid"` to recover the previous behaviour (evaluate every
  value of `lambda_grid` and take the argmin), which is still useful for
  inspecting the risk surface. For `method = "optimize"`, `lambda_grid`
  only sets the search range via `range(lambda_grid)`. The cross-fitted
  tuner
  ([`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md))
  inherits the new default through its per-fold calls. The
  runaway-discrimination guard and the `lambda = 0` (human-only)
  fallback are unchanged.

### Inputs

- `predicted` and `generated` must now be **binary 0/1 responses** in
  all fitting and PPI-score functions; probability (fractional) inputs
  are rejected with a message to sample from them first. Fractional
  values are not a valid likelihood input for the marginal IRT objective
  and break the PPI correction. The low-level quadrature utilities still
  accept fractional input, and now document that the high-level fitters
  do not.

### Robustness

- [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  now wraps each candidate fit in
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html), so a single
  failed fit (from bad starting values, aggressive bounds, or unusual
  response patterns) is treated as an ineligible (infinite-risk)
  candidate instead of aborting the whole tuning run — matching the 1PL
  tuner.

### Documentation

- Softened the README’s headline claims to match the finite-sample
  validation results, and added a “What should I use?”
  function-selection table.
- Added a `R-CMD-check` GitHub Actions workflow.
