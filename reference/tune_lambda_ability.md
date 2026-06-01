# Tune lambda by downstream ability risk

Fits candidate mixed-subjects calibrations, estimates the item-parameter
sandwich covariance for each, and chooses the lambda minimizing average
propagated ability risk on a target response matrix.

## Usage

``` r
tune_lambda_ability(
  lambda_grid,
  observed,
  predicted,
  generated,
  target_resp = NULL,
  theta_true = NULL,
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

A list with `summary`, `best_lambda`, `best_fit`, and all fitted
candidate objects.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
tuned <- tune_lambda_ability(
  c(0, 0.5), observed, observed, generated,
  initial_pars = pars, n_quad = 5, control = list(maxit = 30)
)
tuned$best_lambda
#> [1] 0.5
```
