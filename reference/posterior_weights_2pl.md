# Compute posterior quadrature weights for a 2PL model

Computes each subject's posterior distribution over a fixed quadrature
grid under a 2PL model, using stable log-likelihood calculations.
Fractional responses in `[0, 1]` are allowed at this low level, which is
useful when LLM output is stored as probabilities rather than sampled
binary responses.

## Usage

``` r
posterior_weights_2pl(
  resp,
  item_pars,
  quadrature = NULL,
  n_quad = 31,
  iterlim = 1e+05
)
```

## Arguments

- resp:

  A response matrix with rows for subjects and columns for items. Values
  may be binary, fractional in `[0, 1]`, or `NA`.

- item_pars:

  Item parameters in slope-intercept form. Supply a data frame or matrix
  with columns `a`/`a1` and `d`, or a fitted `mirt` model.

- quadrature:

  Optional quadrature data frame with `theta` and `weight` columns. If
  omitted, a standard-normal grid is created.

- n_quad:

  Number of quadrature nodes used when `quadrature` is omitted.

- iterlim:

  Maximum number of Newton-Raphson iterations passed to
  [`rmutil::gauss.hermite()`](https://rdrr.io/pkg/rmutil/man/gauss.hermite.html)
  when `quadrature` is omitted.

## Value

A matrix with one row per subject and one column per quadrature node.
Rows sum to one. Attributes `theta` and `weight` contain the grid.

## Details

Note: the high-level mixed-subjects fitting functions
([`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
and relatives) require **binary** `predicted` and `generated`;
fractional input is supported only in these low-level quadrature
utilities. If you have LLM-derived probabilities, sample binary
responses from them (e.g. with
[`stats::rbinom()`](https://rdrr.io/r/stats/Binomial.html)) before
calibrating.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
W <- posterior_weights_2pl(resp, pars, n_quad = 5)
rowSums(W)
#> [1] 1 1
```
