# Propagated ability risk for a 1PL fit

Computes `g_i' Sigma_1pl g_i` for each response pattern, where `g_i` is
the `(J+1)`-dimensional gradient of the ability estimate with respect to
`(a_shared, d_1, ..., d_J)` and `Sigma_1pl` is the sandwich covariance
from
[`vcov_mixed_subjects_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_1pl.md).

## Usage

``` r
ability_risk_1pl(
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

  A `"mixedsubjects_1pl_fit"` object or item-parameter data frame.

- vcov:

  Optional `(J+1) × (J+1)` covariance matrix. Required when
  `fit_or_pars` is not a fitted object.

- theta_true:

  Optional true theta values for simulation studies.

- bounds:

  Bounds passed to
  [`score_theta()`](https://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

## Value

A list with `summary` and per-pattern `details`, the same structure as
[`ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/ability_risk.md).
