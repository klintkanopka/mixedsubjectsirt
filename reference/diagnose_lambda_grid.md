# Diagnose lambda values over a grid

Fits
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
or
[`fit_mixed_subjects_split()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md)
over a set of candidate lambda values. The returned summary reports the
fitted mixed-subjects objective and the observed human expected-count
loss for each candidate. This is a sensitivity diagnostic, not a
statistically valid PPI++ tuning rule.

## Usage

``` r
diagnose_lambda_grid(
  lambda_grid,
  observed,
  predicted,
  generated,
  split = FALSE,
  ...
)
```

## Arguments

- lambda_grid:

  Numeric vector of lambda values in `[0, 1]`.

- observed, predicted, generated:

  Response matrices passed to
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).

- split:

  Logical; if `TRUE`, call
  [`fit_mixed_subjects_split()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md).

- ...:

  Additional arguments passed to the selected fitting function.

## Value

A list with `summary`, `lowest_observed_loss_lambda`, and all fitted
model objects.

## Examples

``` r
set.seed(3)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(30), pars)
predicted <- observed
generated <- simulate_2pl(rnorm(80), pars)
tuned <- diagnose_lambda_grid(
  c(0, 0.5),
  observed, predicted, generated,
  initial_pars = pars, n_quad = 5, control = list(maxit = 30)
)
tuned$summary
#>   lambda mixed_loss observed_loss convergence
#> 1    0.0   1.654606      1.654606           0
#> 2    0.5   1.793915      1.684653           0
```
