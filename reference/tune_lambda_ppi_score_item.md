# Per-item PPI++ optimal tuning parameters

Applies the PPI++ Proposition 2 plug-in formula independently for each
item, producing a vector of item-specific lambda values `lambda_j` in
`[0, 1]`.

## Usage

``` r
tune_lambda_ppi_score_item(
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

  Item parameters at which to evaluate the score vectors.

- n_generated:

  Number of generated (unpaired) LLM subjects.

- quadrature:

  Optional quadrature grid.

- n_quad:

  Number of quadrature nodes when `quadrature` is omitted.

## Value

A list with `lambda` (numeric vector of length `n_items`), `item` (item
names), `n`, `n_generated`, and `r` (the ratio `n / n_generated`).

## Details

The global
[`tune_lambda_ppi_score()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
uses the full parameter covariance matrix `Tr(Sigma_gamma)` as the
objective. This function instead applies the same formula using only the
2x2 diagonal block of the inverse Hessian for item `j`, and the 2D
sub-vectors of the human and paired-LLM score vectors. The result is the
lambda that minimizes the marginal variance of `(a_j, d_j)`
independently for each item.

**Use case.** When a single global lambda is forced to zero because a
few items have poor LLM predictions, per-item lambda_j allows
well-predicted items to still benefit from the LLM data. Pass the
returned vector to
[`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
as the `lambda` argument.

This is a **theoretical diagnostic**: it minimizes item-parameter
variance, not ability-score risk. For operational scoring use
[`tune_lambda_ability_risk_item()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_item.md)
instead.

## See also

[`tune_lambda_ppi_score()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
for the global version;
[`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
to fit with a per-item lambda vector.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = c(1, 1.2, 0.9), d = c(0, -0.5, 0.3))
observed <- simulate_2pl(rnorm(40), pars)
tune_lambda_ppi_score_item(observed, observed, pars, n_generated = 100, n_quad = 7)$lambda
#> [1] 0.7142857 0.7142857 0.7142857
```
