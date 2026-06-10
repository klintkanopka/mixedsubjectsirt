# Fit a split-sample mixed-subjects 2PL calibration

Fits the same objective as
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md),
but constructs labeled expected counts with cross-fitted posterior
weights. For each split, the initial human 2PL model is fit on the other
splits and then used to compute posterior weights for the held-out
split. Each human row contributes to the final estimating equation
exactly once.

## Usage

``` r
fit_mixed_subjects_split(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_splits = 2,
  split_id = NULL,
  seed = NULL,
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

  Power-tuning parameter in `[0, 1]`. Supply a scalar for a fixed lambda
  or a vector with one value per split for a precomputed
  cross-fitted-lambda analysis. When a vector is supplied, the generated
  term uses the split-size weighted mean lambda.

- n_splits:

  Number of sample splits.

- split_id:

  Optional integer vector assigning each observed row to a split. If
  omitted, splits are sampled at random.

- seed:

  Optional random seed used when `split_id` is omitted.

- n_quad:

  Number of standard-normal quadrature nodes.

- initial_pars:

  Optional item parameters to use in every fold instead of fitting
  fold-specific human models. This is mainly useful for testing or
  sensitivity analyses.

- quadrature:

  Optional quadrature grid with `theta` and `weight` columns.

- common_predicted_weights:

  Logical; if `TRUE`, reuse each held-out observed posterior weight
  matrix for its paired LLM responses.

- paired_missing:

  How to handle missingness when `common_predicted_weights = TRUE`. The
  default, `"match_observed"`, requires `observed` and `predicted` to
  have the same missingness pattern.

- slope_lower:

  Lower bound for discrimination parameters during optimization. Use
  `NULL` for no lower bound.

- slope_upper:

  Upper bound for discrimination parameters during optimization. Use
  `NULL` (the default) for no upper bound.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  when fold-specific initial models are fit.

## Value

An object of class `"mixedsubjects_fit"` with `split` metadata and
fold-level initial parameters.

## Details

Generated LLM counts are computed once per fold and averaged across
folds so that the generated sample keeps its original sample-size scale.

## Examples

``` r
set.seed(2)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
predicted <- observed
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects_split(
  observed, predicted, generated,
  lambda = 0.5, initial_pars = pars, n_splits = 2,
  n_quad = 7, control = list(maxit = 50)
)
fit$item_pars
#>   item        a          d          b
#> 1    1 0.965167  0.1019467 -0.1056260
#> 2    2 1.232267 -0.5651595  0.4586339
#> 3    3 1.004847  0.4347639 -0.4326667
```
