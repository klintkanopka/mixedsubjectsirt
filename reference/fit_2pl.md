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
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9, 1.1, 0.8), d = c(0, 0.5, -0.5, 0.2, -0.3))
resp <- simulate_2pl(rnorm(500), pars)
fit <- fit_2pl(resp)
fit$pars
#>   item         a           d           b
#> 1    1 1.4074559 -0.04157374  0.02953822
#> 2    2 0.7732065  0.59257775 -0.76638999
#> 3    3 0.6796210 -0.40148209  0.59074406
#> 4    4 1.2125743  0.10512718 -0.08669752
#> 5    5 0.6915828 -0.32088333  0.46398399
```
