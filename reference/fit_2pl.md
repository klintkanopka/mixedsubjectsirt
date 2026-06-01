# Fit a unidimensional 2PL IRT model

Fits a two-parameter logistic model with `mirt` and returns item
parameters in slope-intercept form. The response probability is
`plogis(d + a * theta)`, where `a` is the discrimination and `d` is the
intercept. Difficulty is returned as `b = -d / a`.

## Usage

``` r
fit_2pl(resp, technical = list(NCYCLES = 1000), verbose = FALSE, ...)
```

## Arguments

- resp:

  A numeric item response matrix with rows for subjects and columns for
  items. Values must be binary `0`/`1`; `NA` is allowed.

- technical:

  A list passed to the `technical` argument of
  [`mirt::mirt()`](https://philchalmers.github.io/mirt/reference/mirt.html).

- verbose:

  Logical; passed to
  [`mirt::mirt()`](https://philchalmers.github.io/mirt/reference/mirt.html).

- ...:

  Additional arguments passed to
  [`mirt::mirt()`](https://philchalmers.github.io/mirt/reference/mirt.html).

## Value

A list with `pars`, a data frame containing `item`, `a`, `d`, and `b`,
and `model`, the fitted `mirt` model.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fit_2pl(response_matrix)
fit$pars
} # }
```
