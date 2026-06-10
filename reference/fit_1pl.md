# Fit a 1PL (one-parameter logistic) model

Estimates a shared discrimination parameter `a` (equal across all items)
and per-item intercepts `d_j` by maximizing the IRT marginal likelihood
under a standard-normal ability prior using L-BFGS-B.

## Usage

``` r
fit_1pl(
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

  Binary response matrix.

- n_quad:

  Number of standard-normal quadrature nodes.

- initial_pars:

  Optional starting item parameters (data frame with `a` and `d`
  columns). If omitted, `a = 1` and `d_j = qlogis(p_j)` where `p_j` is
  the observed proportion correct for item `j`.

- quadrature:

  Optional quadrature grid.

- slope_lower, slope_upper:

  Bounds on the shared discrimination.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

A list with `pars` (item parameter data frame with all `a` equal), `par`
(the raw parameter vector), and optimizer details.

## Details

The response probability is
`P(x_j = 1 | theta) = plogis(a * theta + d_j)`. The parameter vector has
length `J + 1`: one shared discrimination followed by J per-item
intercepts.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
resp <- simulate_2pl(rnorm(60), pars)
fit <- fit_1pl(resp, n_quad = 7)
fit$pars
#>   item         a           d          b
#> 1    1 0.1013721 -0.40647830  4.0097660
#> 2    2 0.1013721  0.06686741 -0.6596236
#> 3    3 0.1013721  0.84933748 -8.3784165
```
