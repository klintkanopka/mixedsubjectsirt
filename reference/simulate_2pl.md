# Simulate 2PL item responses

Generates binary item responses from the model `plogis(d + a * theta)`.

## Usage

``` r
simulate_2pl(theta, item_pars)
```

## Arguments

- theta:

  Numeric vector of latent trait values.

- item_pars:

  Item parameters in slope-intercept form. Supply a data frame or matrix
  with columns `a`/`a1` and `d`, or a fitted `mirt` model.

## Value

A binary response matrix with one row per value of `theta` and one
column per item.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
simulate_2pl(rnorm(5), pars)
#>      1 2
#> [1,] 0 0
#> [2,] 1 1
#> [3,] 0 1
#> [4,] 1 1
#> [5,] 0 1
```
