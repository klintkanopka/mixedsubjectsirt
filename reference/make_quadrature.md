# Create a standard-normal Gauss-Hermite quadrature grid

[`rmutil::gauss.hermite()`](https://rdrr.io/pkg/rmutil/man/gauss.hermite.html)
returns nodes and weights for integrals of the form
`integral f(x) exp(-x^2) dx`. This function rescales those nodes and
weights to approximate expectations under a standard normal latent trait
distribution.

## Usage

``` r
make_quadrature(n_quad = 31, iterlim = 1e+05)
```

## Arguments

- n_quad:

  Number of quadrature nodes.

- iterlim:

  Maximum number of Newton-Raphson iterations passed to
  [`rmutil::gauss.hermite()`](https://rdrr.io/pkg/rmutil/man/gauss.hermite.html).

## Value

A data frame with node index, `theta`, `weight`, and backward compatible
aliases `X_k` and `A_k`.

## Examples

``` r
quad <- make_quadrature(7)
sum(quad$weight)
#> [1] 1
```
