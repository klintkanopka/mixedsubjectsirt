# Fit a mixed-subjects 1PL calibration via marginal maximum likelihood

Analogous to
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
but estimates a shared discrimination parameter `a` across all items
(1PL model). Posteriors are recomputed at every gradient evaluation — no
frozen-posterior gradient asymmetry.

## Usage

``` r
fit_mixed_subjects_mml_1pl(
  observed,
  predicted,
  generated,
  lambda = 1,
  n_quad = 31,
  initial_pars = NULL,
  quadrature = NULL,
  mml_pred_weights = c("own", "human"),
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

- mml_pred_weights:

  How to compute posteriors for the paired `predicted` term. `"own"`
  uses posteriors from the predicted responses; `"human"` uses
  posteriors from the observed human responses. See Details.

- slope_lower:

  Lower bound for discrimination parameters during optimization. Use
  `NULL` for no lower bound.

- slope_upper:

  Upper bound on discrimination parameters. Unlike
  [`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md),
  this function should not require capping for well-posed problems
  because the true marginal objective has no false minimum at large
  discrimination.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
  when `initial_pars` is omitted.

## Value

An object of class `c("mixedsubjects_1pl_fit", "mixedsubjects_fit")`.

## Details

Only scalar `lambda` is supported; per-item lambda is not meaningful for
the 1PL because the discrimination is shared across items.

## See also

[`fit_mixed_subjects_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_1pl.md)
for the frozen expected-count version;
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
for the 2PL version.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
observed  <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects_mml_1pl(
  observed, observed, generated,
  lambda = 0.5, initial_pars = pars, n_quad = 7,
  control = list(maxit = 100)
)
fit$item_pars
#>   item        a           d           b
#> 1    1 1.152006 -0.40800463  0.35416883
#> 2    2 1.152006 -0.06273155  0.05445419
#> 3    3 1.152006  0.37033652 -0.32147099
```
