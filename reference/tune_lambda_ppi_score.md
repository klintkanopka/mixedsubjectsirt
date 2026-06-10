# Plug-in PPI++ optimal tuning parameter

Implements the closed-form estimator from Proposition 2 of Angelopoulos,
Duchi and Zrnic (2023) for the lambda that minimizes the trace of the
asymptotic item-parameter covariance matrix `Tr(Sigma_gamma)`.

## Usage

``` r
tune_lambda_ppi_score(
  observed,
  predicted,
  item_pars,
  n_generated,
  quadrature = NULL,
  n_quad = 31
)
```

## Arguments

- observed:

  Human response matrix.

- predicted:

  Paired binary LLM responses (0/1) for the same rows as `observed`.
  Probabilities are not accepted; sample binary responses first.

- item_pars:

  Item parameters in slope-intercept form at which to evaluate the score
  vectors. Typically the human 2PL MLE from
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md).

- n_generated:

  Number of generated (unpaired) LLM subjects, used to compute the ratio
  `r` (`n / n_generated`).

- quadrature:

  Optional quadrature grid. If omitted, a standard-normal grid with
  `n_quad` nodes is created.

- n_quad:

  Number of quadrature nodes when `quadrature` is omitted.

## Value

A list with elements `lambda` (the plug-in estimate, clipped to \[0,
1\]), `n`, `n_generated`, `r`, and the intermediate matrices `C_hf`
(cross-covariance of human and paired-LLM score vectors) and `V_f`
(variance of paired-LLM score vectors).

## Details

**This is the item-parameter variance objective, not the psychometric
scoring objective.** For IRT applications where accurate ability scoring
is the goal, use
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
or
[`tune_lambda_ability_risk_crossfit()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md)
instead. Those functions directly minimize the propagated ability-score
risk `E[g' Sigma_gamma g]` — the quantity that matters for test scoring
— rather than item-parameter estimation efficiency.
`tune_lambda_ppi_score()` is provided as a theoretical diagnostic and to
facilitate method validation.

The formula uses the **same** human posterior weights for both the human
and paired-LLM score vectors. This symmetry is required for the PPI++
unbiasedness condition `E[grad_gen] = E[grad_pred]` at the true
parameters.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
predicted <- observed
tune_lambda_ppi_score(observed, predicted, pars, n_generated = 100, n_quad = 7)$lambda
#> [1] 0.7142857
```
