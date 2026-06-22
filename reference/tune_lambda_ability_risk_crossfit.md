# Cross-fit ability-score-risk lambda tuning

Estimates lambda separately for each held-out split using only the
remaining labeled rows, then fits a final model. By default
(`final_fit_fn = fit_mixed_subjects_mml`) the fold lambdas are averaged
(weighted by fold size) into a single scalar and the full sample is
refit; pass `final_fit_fn = fit_mixed_subjects_split` to instead fit
each fold's rows with its own out-of-fold lambda.

## Usage

``` r
tune_lambda_ability_risk_crossfit(
  lambda_grid = seq(0, 1, by = 0.1),
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
  fit_fn = fit_mixed_subjects_mml,
  final_fit_fn = fit_mixed_subjects_mml,
  tuning_args = list(),
  final_args = list(),
  bounds = c(-6, 6),
  control = list(maxit = 500),
  ...
)
```

## Arguments

- lambda_grid:

  Numeric vector of candidate lambda values in `[0, 1]`. For
  `method = "grid"` these are the evaluated candidates; for
  `method = "optimize"` only `range(lambda_grid)` matters and bounds the
  search (e.g. `lambda_grid = c(0, 0.8)` caps lambda at 0.8). Defaults
  to `seq(0, 1, by = 0.1)`.

- observed, predicted, generated:

  Response matrices passed to
  [`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).

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
  [`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)).
  Defaults to
  [`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  (marginal MML, recommended). The frozen expected-count estimator
  [`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
  is still available but discouraged.

- final_fit_fn:

  Function used to produce the final combined-data fit. Defaults to
  [`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md),
  giving a scalar marginal-MML final fit: the fold-specific lambdas are
  averaged (weighted by fold size) into a single scalar and the full
  sample is refit. Pass
  [`fit_mixed_subjects_split()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md)
  to instead keep the per-fold `lambda` vector and fit each fold's rows
  with its own out-of-fold lambda — the textbook cross-fit decoupling,
  but it uses the discouraged frozen expected-count split estimator.

- tuning_args:

  Named list of extra arguments forwarded only to the fold-level
  [`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  calls (and through them to `fit_fn`). For example,
  `tuning_args = list(slope_upper = 4)`.

- final_args:

  Named list of extra arguments forwarded only to `final_fit_fn`. For
  example, `final_args = list(mml_pred_weights = "own")`. This keeps
  tuning-specific and final-fit-specific arguments cleanly separated,
  avoiding the earlier `...` leakage between the two.

- bounds:

  Bounds passed to
  [`score_theta()`](https://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Deprecated; forwarded to `tuning_args` for backward compatibility with
  a one-time message. Prefer `tuning_args` / `final_args`.

## Value

A list with fold-specific lambda values, fold tuning objects, and the
final fit.
