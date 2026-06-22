# Propagated ability risk from item-parameter uncertainty

Computes `g_i' Sigma g_i` for each response pattern, where `g_i` is the
gradient of the ability estimate with respect to item parameters. If
`theta_true` is supplied, the returned total risk also includes squared
ability estimation error.

## Usage

``` r
ability_risk(
  resp,
  fit_or_pars,
  vcov = NULL,
  theta_true = NULL,
  bounds = c(-6, 6)
)
```

## Arguments

- resp:

  Target response matrix.

- fit_or_pars:

  A `"mixedsubjects_fit"` object or item-parameter data frame.

- vcov:

  Optional covariance matrix. Required when `fit_or_pars` is item
  parameters rather than a fitted mixed-subjects object.

- theta_true:

  Optional true theta values for simulation studies.

- bounds:

  Bounds passed to
  [`score_theta()`](https://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

## Value

A list with `summary` and per-pattern `details`.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- simulate_2pl(rnorm(30), pars)
Sigma <- diag(0.01, 4)
ability_risk(resp, pars, vcov = Sigma)$summary
#>   mean_param_var mean_squared_error mean_total_risk
#> 1     0.09088749                 NA      0.09088749
```
