# Summarize response data as expected quadrature counts

Converts response data and posterior quadrature weights into Bock-Aitkin
style expected counts. For each item and quadrature node, `N` is the
expected number of observed responses and `R` is the expected number
correct.

## Usage

``` r
summarize_expected_counts(resp, weights)
```

## Arguments

- resp:

  A response matrix with rows for subjects and columns for items.

- weights:

  Posterior quadrature weights, usually returned by
  [`posterior_weights_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/posterior_weights_2pl.md).

## Value

A list of class `"mixedsubjects_counts"` containing matrices `N` and
`R`, sample size `n`, quadrature nodes, quadrature weights, and item
names.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
W <- posterior_weights_2pl(resp, pars, n_quad = 5)
counts <- summarize_expected_counts(resp, W)
counts$N
#>             node1     node2    node3     node4       node5
#> Item1 0.003717598 0.2836301 1.312301 0.3944443 0.005906825
#> Item2 0.003717598 0.2836301 1.312301 0.3944443 0.005906825
```
