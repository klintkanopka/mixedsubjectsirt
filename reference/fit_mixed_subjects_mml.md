# Fit a mixed-subjects 2PL calibration via marginal maximum likelihood

Estimates item parameters using the true IRT marginal likelihood for all
three loss terms. Unlike
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md),
which freezes posterior quadrature weights at the initial parameter
estimates before optimizing, this function recomputes posterior weights
at every gradient evaluation. This eliminates the gradient asymmetry
that causes
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
to converge to false minima at inflated discrimination values when LLM
item parameters differ from human parameters.

## Usage

``` r
fit_mixed_subjects_mml(
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
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  when `initial_pars` is omitted.

## Value

An object of class `"mixedsubjects_fit"` with the same structure as
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).
For **scalar** lambda fits, the quadrature summaries store posteriors at
the converged parameters, and
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) dispatches
automatically to
[`vcov_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md)
to compute the Louis-corrected marginal sandwich covariance. Calling
[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)
directly bypasses the Louis correction. For **vector** lambda fits, the
summaries store the frozen posteriors used during optimization, and
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html) dispatches to
[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)
(EM bread) for consistency with the frozen Q-function objective.

## Details

**Why it matters for lambda selection.** With the frozen expected-count
implementation, the gradient of `L_pred` uses concentrated human
posteriors while `L_gen` uses diffuse LLM posteriors, making
`grad(L_pred) >> grad(L_gen)` and systematically pushing discriminations
upward at any `lambda > 0`. In the marginal-MML formulation all three
terms use their own current-parameter posteriors, so the asymmetry is
absent at the true optimum. As a result
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
selects `lambda > 0` whenever the LLM predictions are genuinely
informative (e.g. `predicted = observed`), rather than collapsing to
`lambda = 0` for all misaligned LLMs.

**`mml_pred_weights`.**

- `"own"` (default):

  L_pred uses posteriors computed from the *predicted* response matrix
  at the current parameter values. All three terms are true marginal
  likelihoods; objective and gradient are internally consistent.
  Recommended for most applications and required for
  [`vcov_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md)
  to produce the fully correct Louis-formula bread.

- `"human"`:

  L_pred uses posteriors computed from the *observed* (human) response
  matrix, frozen at `initial_pars`. This is a **fixed-nuisance
  Q-function**: the predicted term is treated as a frozen expected-count
  lower bound rather than a true marginal likelihood. Objective and
  gradient are mutually consistent (both use the same frozen posteriors)
  so L-BFGS-B converges correctly. Useful when strong ability-level
  pairing is needed. Note that
  [`vcov_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md)
  applies Louis' formula to the stored fixed posteriors, which is
  approximately correct when `initial_pars` ≈ `conv_pars`.

**Per-item lambda (vector `lambda`).** When `lambda` is a
length-`n_items` vector rather than a scalar, `fit_mixed_subjects_mml`
switches to a **frozen Q-function** objective: expected-count counts are
computed once from `initial_pars` and held fixed during L-BFGS-B, with
item `j`'s counts weighted by `lambda[j]`. This is a consistent
(objective, gradient) pair but is *not* the full marginal-MML objective
— it is a frozen expected-count approximation analogous to
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md).
Per-item lambda values obtained from
[`tune_lambda_ability_risk_item()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_item.md)
assign `lambda_j ≈ 0` to items where the LLM correction is harmful,
containing the frozen-posterior gradient asymmetry. Document per-item
lambda results as approximate.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed  <- simulate_2pl(rnorm(40), pars)
generated <- simulate_2pl(rnorm(100), pars)
fit <- fit_mixed_subjects_mml(
  observed, observed, generated,
  lambda = 0.5, initial_pars = pars, n_quad = 7,
  control = list(maxit = 100)
)
fit$item_pars
#>   item        a          d          b
#> 1    1 1.486806 -0.2247580  0.1511684
#> 2    2 1.294408 -0.9772844  0.7550047
#> 3    3 0.480298  0.1161510 -0.2418311
```
