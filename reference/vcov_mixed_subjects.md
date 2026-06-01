# Sandwich covariance for a mixed-subjects fit

Estimates the full sandwich covariance matrix for item parameters from
the fixed-posterior expected-count estimating equations. The parameter
order is all discriminations followed by all intercepts, matching
`fit$par`.

## Usage

``` r
vcov_mixed_subjects(object, ridge = 1e-08, ...)
```

## Arguments

- object:

  A fitted object returned by
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
  or
  [`fit_mixed_subjects_from_quadrature()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_from_quadrature.md)
  with response matrices and posterior weights available in its
  quadrature summaries.

- ridge:

  Small ridge value used when inverting the Hessian.

- ...:

  Unused; included for method compatibility.

## Value

A covariance matrix with attributes `bread` and `meat`.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
fit <- fit_mixed_subjects(
  observed, observed, simulate_2pl(rnorm(80), pars),
  lambda = 0.5, initial_pars = pars, n_quad = 7
)
dim(vcov_mixed_subjects(fit))
#> [1] 6 6
```
