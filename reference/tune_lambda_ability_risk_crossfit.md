# Cross-fit ability-score-risk lambda tuning

Estimates lambda separately for each held-out split using only the
remaining labeled rows, then fits a final model with those fold-specific
lambda values.

## Usage

``` r
tune_lambda_ability_risk_crossfit(
  lambda_grid,
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
  n_splits = 2,
  split_id = NULL,
  seed = NULL,
  n_quad = 31,
  initial_pars = NULL,
  target_mode = c("fixed", "row_aligned"),
  fit_fn = fit_mixed_subjects,
  final_fit_fn = fit_mixed_subjects_split,
  tuning_args = list(),
  final_args = list(),
  bounds = c(-6, 6),
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

- n_splits:

  Number of sample splits.

- split_id:

  Optional integer split assignment for labeled rows.

- seed:

  Optional seed used when `split_id` is omitted.

- n_quad:

  Number of quadrature nodes.

- initial_pars:

  Optional starting item parameters.

- target_mode:

  How `target_resp` is handled in each fold. `"fixed"` (default): the
  full `target_resp` is used to evaluate risk in every fold, suitable
  when the target population is fixed and independent of the
  labeled-data split (e.g. an operational scoring population).
  `"row_aligned"`: only the training rows of `target_resp` are used,
  which is valid when `target_resp = observed` and fold-matched
  evaluation is desired.

- fit_fn:

  Fitting function used for each fold's ability-risk tuning (passed to
  [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)).
  Defaults to
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
  (frozen expected-count). Pass
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  for marginal-MML fold tuning.

- final_fit_fn:

  Function used to produce the final combined-data fit. Defaults to
  [`fit_mixed_subjects_split()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md),
  which accepts a per-fold `lambda` vector natively. Pass
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  to get a scalar marginal-MML final fit: the fold-specific lambdas are
  averaged (weighted by fold size) into a single scalar, avoiding the
  accidental per-item lambda problem that occurs when a
  length-`n_splits` vector is passed directly to
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md).
  Note that mixing MML fold-tuning with a frozen final fit is an
  approximation; document this when reporting results.

- tuning_args:

  Named list of extra arguments forwarded only to the fold-level
  [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  calls (and through them to `fit_fn`). For example,
  `tuning_args = list(slope_upper = 4)`.

- final_args:

  Named list of extra arguments forwarded only to `final_fit_fn`. For
  example, `final_args = list(mml_pred_weights = "own")`. This keeps
  tuning-specific and final-fit-specific arguments cleanly separated,
  avoiding the earlier `...` leakage between the two.

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Deprecated; forwarded to `tuning_args` for backward compatibility with
  a one-time message. Prefer `tuning_args` / `final_args`.

## Value

A list with fold-specific lambda values, fold tuning objects, and the
final fit.
