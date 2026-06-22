# Plug-in PPI++ optimal tuning parameter for a 1PL model

Applies the PPI++ Proposition 2 formula using `(J+1)`-dimensional score
vectors for the 1PL parameterization `(a_shared, d_1, ..., d_J)`.

## Usage

``` r
tune_lambda_ppi_score_1pl(
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
  [`fit_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md).

- n_generated:

  Number of generated (unpaired) LLM subjects, used to compute the ratio
  `r` (`n / n_generated`).

- quadrature:

  Optional quadrature grid. If omitted, a standard-normal grid with
  `n_quad` nodes is created.

- n_quad:

  Number of quadrature nodes when `quadrature` is omitted.

## Value

A list with `lambda`, `n`, `n_generated`, `r`, `C_hf`, `V_f`.

## Details

This is the **item-parameter variance** objective — it minimizes
`Tr(Sigma_1pl)`. For practical scoring applications use
[`tune_lambda_ability_risk_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_1pl.md)
instead.

## Examples

``` r
set.seed(1)
pars <- data.frame(a = 1, d = c(-0.5, 0, 0.5))
obs  <- simulate_2pl(rnorm(40), pars)
tune_lambda_ppi_score_1pl(obs, obs, pars, n_generated = 100, n_quad = 7)$lambda
#> [1] 0.7142857
```
