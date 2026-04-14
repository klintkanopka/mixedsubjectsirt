# Converts an item response matrix into quadrature form and estimates a 2PL IRT model with parameters in slope-intercept form

Converts an item response matrix into quadrature form and estimates a
2PL IRT model with parameters in slope-intercept form

## Usage

``` r
mixed_subjects_quadrature(
  resp,
  N_quad = 10,
  eps = 1e-15,
  iterlim = 1e+05,
  irt_pars = NULL
)
```

## Arguments

- resp:

  An item response matrix in wide form

- N_quad:

  The number of quadrature points to compute. Higher numbers can induce
  numerical errors

- eps:

  A tolerance for the minimum values

- iterlim:

  Maximum number of Newton-Raphson iterations passed to
  [`rmutil::gauss.hermite()`](https://rdrr.io/pkg/rmutil/man/gauss.hermite.html)

- irt_pars:

  IRT parameters from human calibration for rescaling parameters
  estimated from predicted responses

## Value

A list object with two components. `$quad` contains a dataframe with the
expected sample sizes and number of correct responses at each quadrature
point. `$irt_pars` contains the parameter estimates from the fitted IRT
model that generated the expected counts
