# Cross-fit ability-risk lambda tuning

Estimates lambda separately for each held-out split using only the
remaining labeled rows, then fits the final split-sample mixed-subjects
estimator with those fold-specific lambda values.

## Usage

``` r
tune_lambda_ability_crossfit(
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
  studies to add squared scoring error to the risk.

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

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).

## Value

A list with fold-specific lambda values, fold tuning objects, and the
final split-sample fit.
