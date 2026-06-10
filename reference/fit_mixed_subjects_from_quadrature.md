# Fit from precomputed quadrature summaries

Fits the mixed-subjects 2PL objective from quadrature/count summaries
rather than raw response matrices. This lower-level interface is useful
when the human, paired LLM, and generated LLM summaries have already
been linked onto a common scale outside the package.

## Usage

``` r
fit_mixed_subjects_from_quadrature(
  q_observed,
  q_predicted,
  q_generated,
  lambda = 1,
  initial_pars = NULL,
  slope_lower = 1e-04,
  slope_upper = NULL,
  control = list(maxit = 500)
)
```

## Arguments

- q_observed:

  Quadrature summary for observed human responses. Usually returned by
  [`mixed_subjects_quadrature()`](http://klintkanopka.com/mixedsubjectsirt/reference/mixed_subjects_quadrature.md),
  but a raw counts object returned by
  [`summarize_expected_counts()`](http://klintkanopka.com/mixedsubjectsirt/reference/summarize_expected_counts.md)
  is also accepted.

- q_predicted:

  Quadrature summary for paired LLM responses/predictions on the labeled
  human rows.

- q_generated:

  Quadrature summary for generated or unlabeled LLM responses.

- lambda:

  Power-tuning parameter in `[0, 1]`.

- initial_pars:

  Starting item parameters in slope-intercept form. If omitted,
  `q_observed$irt_pars` is used when available.

- slope_lower:

  Lower bound for discrimination parameters during optimization. Use
  `NULL` for no lower bound.

- slope_upper:

  Upper bound for discrimination parameters during optimization. Use
  `NULL` (the default) for no upper bound.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

## Value

An object of class `"mixedsubjects_fit"`.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
fit_mixed_subjects_from_quadrature(q, q, q, lambda = 0.5)$item_pars
#>    item         a             d           b
#> 1 Item1 0.0001000  2.446106e-06 -0.02446106
#> 2 Item2 0.1999931 -1.564259e-02  0.07821564
```
