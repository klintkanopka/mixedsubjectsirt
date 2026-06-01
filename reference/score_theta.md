# Estimate ability scores from a 2PL calibration

Computes bounded maximum-likelihood ability estimates for response
patterns under fixed item parameters. This is a scoring helper for
inspecting fitted calibrations; it does not account for uncertainty in
the item parameters.

## Usage

``` r
score_theta(resp, item_pars, bounds = c(-6, 6))
```

## Arguments

- resp:

  Response matrix with rows for subjects and columns for items.

- item_pars:

  Item parameters in slope-intercept form. Supply a data frame or matrix
  with columns `a`/`a1` and `d`, or a fitted `mirt` model.

- bounds:

  Numeric vector of length two giving the optimization interval for
  theta.

## Value

A numeric vector of ability estimates.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- simulate_2pl(rnorm(5), pars)
score_theta(resp, pars)
#> [1] -5.9999390  5.9999241  0.4122068  5.9999241  0.4122068
```
