# Fit a mixed-subjects 2PL calibration

Fits item parameters using observed human responses, paired LLM
responses/predictions for those same subjects, and generated or
unlabeled LLM responses. This implements the expected-count objective

## Usage

``` r
fit_mixed_subjects(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  common_predicted_weights = TRUE,
  paired_missing = c("match_observed", "allow"),
  slope_lower = 1e-04,
  control = list(maxit = 500),
  ...
)
```

## Arguments

- observed:

  Human response matrix, with rows for subjects and columns for items.
  Values must be binary when `initial_pars` is omitted.

- predicted:

  LLM responses or probabilities for the same rows and items as
  `observed`.

- generated:

  Generated or unlabeled LLM responses or probabilities for the same
  item columns.

- lambda:

  Power-tuning parameter in `[0, 1]`.

- n_quad:

  Number of standard-normal quadrature nodes.

- initial_pars:

  Optional starting item parameters. If omitted, a 2PL model is fit to
  `observed`.

- quadrature:

  Optional quadrature grid with `theta` and `weight` columns.

- common_predicted_weights:

  Logical; if `TRUE`, reuse the observed human posterior weights for
  `predicted`.

- paired_missing:

  How to handle missingness when `common_predicted_weights = TRUE`. The
  default, `"match_observed"`, requires `observed` and `predicted` to
  have the same missingness pattern so the paired LLM correction is
  evaluated only where a human label is present. Use `"allow"` only for
  explicit sensitivity analyses.

- slope_lower:

  Lower bound for discrimination parameters during optimization. Use
  `NULL` for no lower bound.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  when `initial_pars` is omitted.

## Value

An object of class `"mixedsubjects_fit"` with fitted `item_pars`,
optimizer details, quadrature summaries, and input settings.

## Details

`L_human + lambda * (L_generated - L_paired_llm)`.

By default the paired LLM responses reuse the posterior quadrature
weights from the observed human responses. This keeps the paired human
and LLM terms on the same latent covariate distribution, which is the
closest analog to prediction-powered inference with paired labels.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
predicted <- observed
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects(
  observed, predicted, generated,
  lambda = 0.5, initial_pars = pars, n_quad = 7,
  control = list(maxit = 50)
)
fit$item_pars
#>   item         a           d           b
#> 1    1 0.9345154 -0.05509868  0.05895963
#> 2    2 1.0991610 -0.84500017  0.76876835
#> 3    3 0.7051842  0.24631388 -0.34929016
```
