# Fit a mixed-subjects 1PL calibration (frozen expected-count)

Analogous to
[`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
but estimates a shared discrimination parameter `a` across all items
(1PL model). Posterior quadrature weights are frozen at the initial
parameter estimates.

## Usage

``` r
fit_mixed_subjects_1pl(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  common_predicted_weights = TRUE,
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
  [`fit_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
  when `initial_pars` is omitted.

## Value

An object of class `c("mixedsubjects_1pl_fit", "mixedsubjects_fit")`.

## See also

[`fit_mixed_subjects_mml_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md)
for the marginal-likelihood version;
[`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
for the 2PL version.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
observed  <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects_1pl(
  observed, observed, generated,
  lambda = 0.5, initial_pars = pars, n_quad = 7,
  control = list(maxit = 50)
)
fit$item_pars
#>   item        a           d           b
#> 1    1 1.049967 -0.38433720  0.36604710
#> 2    2 1.049967 -0.04797307  0.04569009
#> 3    3 1.049967  0.37394087 -0.35614552
```
