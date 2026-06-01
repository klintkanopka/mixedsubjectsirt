# Gradient of ML ability scores with respect to item parameters

Computes the implicit derivative of bounded maximum-likelihood ability
scores with respect to 2PL item parameters. The column order is all
discriminations followed by all intercepts.

## Usage

``` r
ability_gradient(resp, item_pars, theta = NULL, bounds = c(-6, 6), eps = 1e-10)
```

## Arguments

- resp:

  Response matrix with rows for subjects and columns for items.

- item_pars:

  Item parameters in slope-intercept form, or a `"mixedsubjects_fit"`
  object.

- theta:

  Optional precomputed ability estimates. If omitted,
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md)
  is used.

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md)
  when `theta` is omitted.

- eps:

  Tolerance used to mark near-zero test information as undefined.

## Value

A matrix with one row per response pattern and one column per item
parameter.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
ability_gradient(resp, pars)
#>         a_Item1    a_Item2    d_Item1    d_Item2
#> [1,]  0.7734458 -0.7106064 -0.4193319 -0.4838901
#> [2,] -1.1679932  0.6298063 -0.3996755 -0.5002704
```
