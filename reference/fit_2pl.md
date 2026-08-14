# Fit a unidimensional 2PL IRT model

Estimates per-item discriminations `a_j` and intercepts `d_j` by
maximizing the IRT marginal likelihood under a standard-normal ability
prior using L-BFGS-B. The response probability is
`plogis(d + a * theta)`, where `a` is the discrimination and `d` is the
intercept. Difficulty is returned as `b = -d / a`.

## Usage

``` r
fit_2pl(
  resp,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  slope_lower = 1e-04,
  slope_upper = NULL,
  control = list(maxit = 500)
)
```

## Arguments

- resp:

  A numeric item response matrix with rows for subjects and columns for
  items. Values must be binary `0`/`1`; `NA` is allowed.

- n_quad:

  Number of standard-normal quadrature nodes.

- initial_pars:

  Optional starting item parameters (matrix or data frame with `a`/`a1`
  and `d` columns). If omitted, `a_j = 1` and `d_j = qlogis(p_j)` where
  `p_j` is the observed proportion correct for item `j`.

- quadrature:

  Optional quadrature grid with `theta` and `weight` columns.

- slope_lower, slope_upper:

  Bounds on the discriminations. `NULL` leaves the corresponding side
  unbounded.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `pars`, a data frame containing `item`, `a`, `d`, and `b`;
`par`, the raw parameter vector; optimizer details (`value`,
`convergence`, `message`); and `model`, which is `NULL` and retained for
backward compatibility.

## Details

This is the 2PL counterpart of
[`fit_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
and uses the same marginal likelihood, quadrature, and gradient
machinery as
[`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
at `lambda = 0`, so a human-only baseline fit and the mixed-subjects
estimator are optimized on a common objective.

## See also

[`fit_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
for the shared-discrimination version.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9, 1.1, 0.8), d = c(0, 0.5, -0.5, 0.2, -0.3))
resp <- simulate_2pl(rnorm(500), pars)
fit <- fit_2pl(resp)
fit$pars
#>   item         a           d           b
#> 1    1 1.4084028 -0.04161168  0.02954530
#> 2    2 0.7733562  0.59261605 -0.76629121
#> 3    3 0.6796751 -0.40146365  0.59066992
#> 4    4 1.2116623  0.10508500 -0.08672796
#> 5    5 0.6915083 -0.32087732  0.46402525
```
