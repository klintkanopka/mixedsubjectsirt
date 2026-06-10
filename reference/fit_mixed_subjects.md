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
  slope_upper = NULL,
  control = list(maxit = 500),
  ...
)
```

## Arguments

- observed:

  Human response matrix, with rows for subjects and columns for items.
  Values must be binary when `initial_pars` is omitted.

- predicted:

  Binary LLM responses (0/1) for the same rows and items as `observed`.
  Probabilities are not accepted: fractional values are not a valid
  likelihood input for the marginal IRT objective and break the PPI
  correction, so sample binary responses from any probabilities first
  (e.g. `rbinom`).

- generated:

  Binary generated or unlabeled LLM responses (0/1) for the same item
  columns. Probabilities are not accepted (see `predicted`).

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

- slope_upper:

  Upper bound for discrimination parameters during optimization. Use
  `NULL` (the default) for no upper bound. Setting a finite bound
  (e.g. 4) can stabilize the frozen expected-count fit when the LLM
  parameters differ substantially from the human parameters.

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
#> 1    1 1.0503170 -0.08682871  0.08266905
#> 2    2 1.2202867 -0.82420815  0.67542171
#> 3    3 0.8047358  0.20886422 -0.25954384
```
