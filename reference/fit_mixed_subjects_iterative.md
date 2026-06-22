# Fit a mixed-subjects 2PL calibration with iterative EM

Extends
[`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
by iterating the E-step and M-step until convergence rather than fixing
posterior quadrature weights at the initial parameter estimates. At
every iteration the posterior weights for all three datasets (observed,
predicted, generated) are recomputed using the same current item
parameters. This keeps the posteriors internally consistent and avoids
the asymmetry between `L_pred` and `L_gen` that arises when frozen
human-MLE weights are applied to LLM data with different item
parameters.

## Usage

``` r
fit_mixed_subjects_iterative(
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
  tol = 1e-04,
  em_maxit = 30,
  control = list(maxit = 200),
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

  Upper bound on discrimination parameters. **Strongly recommended**
  when `lambda > 0` — the iterative EM updates posteriors at each step,
  and without an upper bound the gradient asymmetry between `L_pred` and
  `L_gen` can compound across iterations, driving discrimination
  estimates to extreme values. A typical choice is `slope_upper = 4` or
  `slope_upper = 6`.

- tol:

  Convergence tolerance: maximum absolute change in any parameter across
  an EM iteration.

- em_maxit:

  Maximum number of EM iterations.

- control:

  Control list passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`fit_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  when `initial_pars` is omitted.

## Value

An object of class `"mixedsubjects_fit"` with the standard fields plus
`em_iterations` (number of EM cycles completed) and `em_converged`
(logical).

## Details

**Note on lambda selection.** This function accepts a fixed `lambda`.
For psychometric applications where accurate ability scoring is the
goal, select `lambda` with
[`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
rather than
[`tune_lambda_ppi_score()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md).
The PPI++ score objective minimizes the trace of the item-parameter
covariance matrix;
[`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
minimizes the propagated ability-score risk `g' Sigma g`, which is the
quantity that matters for downstream test scoring.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
predicted <- observed
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects_iterative(
  observed, predicted, generated,
  lambda = 0.5, initial_pars = pars, n_quad = 7,
  control = list(maxit = 50), em_maxit = 5
)
#> Warning: fit_mixed_subjects_iterative() with lambda > 0 and no slope_upper can diverge to extreme discrimination values. Setting slope_upper (e.g. slope_upper = 6) is strongly recommended.
fit$item_pars
#>   item         a          d          b
#> 1    1 1.2132598 -0.2037054  0.1678992
#> 2    2 1.3220410 -0.9830630  0.7435949
#> 3    3 0.5810652  0.1202592 -0.2069634
```
