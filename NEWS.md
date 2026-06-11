# mixedsubjectsirt 0.0.1 (development version)

## Lambda tuning

* `tune_lambda_ability_risk()` and `tune_lambda_ability_risk_1pl()` now select
  lambda by **direct 1-D optimization** (`stats::optimize()`) by default, returning
  a continuous lambda with no grid rounding. Pass `method = "grid"` to recover the
  previous behaviour (evaluate every value of `lambda_grid` and take the argmin),
  which is still useful for inspecting the risk surface. For `method = "optimize"`,
  `lambda_grid` only sets the search range via `range(lambda_grid)`. The
  cross-fitted tuner (`tune_lambda_ability_risk_crossfit()`) inherits the new
  default through its per-fold calls. The runaway-discrimination guard and the
  `lambda = 0` (human-only) fallback are unchanged.

## Inputs

* `predicted` and `generated` must now be **binary 0/1 responses** in all fitting
  and PPI-score functions; probability (fractional) inputs are rejected with a
  message to sample from them first. Fractional values are not a valid likelihood
  input for the marginal IRT objective and break the PPI correction.
