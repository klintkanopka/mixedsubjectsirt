# Gradient of ML ability scores w.r.t. 1PL item parameters

Computes the implicit derivative of bounded maximum-likelihood ability
scores with respect to the 1PL parameters `(a_shared, d_1, ..., d_J)`.

## Usage

``` r
ability_gradient_1pl(
  resp,
  item_pars,
  theta = NULL,
  bounds = c(-6, 6),
  eps = 1e-10
)
```

## Arguments

- resp:

  Response matrix.

- item_pars:

  Item parameters with all `a` equal (1PL), or a
  `"mixedsubjects_1pl_fit"` object.

- theta:

  Optional precomputed ability estimates.

- bounds:

  Bounds passed to
  [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md).

- eps:

  Tolerance for near-zero test information.

## Value

A matrix with one row per response pattern and `J + 1` columns
(`a_shared`, then one column per item's `d_j`).

## Details

The gradient for the shared discrimination is the sum of the per-item
discrimination gradients: `da_shared = sum_j da_j` (chain rule via the
constraint `a_j = a_shared`).
