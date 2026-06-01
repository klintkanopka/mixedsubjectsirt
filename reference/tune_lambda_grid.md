# Tune lambda over a grid

`tune_lambda_grid()` is retained as a backward-compatible wrapper for
[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md).
It emits a warning because the grid output is a sensitivity diagnostic,
not a statistically valid PPI++ tuning rule.

## Usage

``` r
tune_lambda_grid(
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

A list returned by
[`diagnose_lambda_grid()`](http://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md).
For backward compatibility, the list also contains
`best_lambda_by_observed_loss`.
