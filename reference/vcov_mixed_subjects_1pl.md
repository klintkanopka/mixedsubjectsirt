# Sandwich covariance for a 1PL mixed-subjects fit

Estimates the `(J+1) × (J+1)` sandwich covariance matrix for the shared
discrimination and per-item intercepts of a 1PL mixed-subjects
calibration.

## Usage

``` r
vcov_mixed_subjects_1pl(object, ridge = 1e-08, ...)
```

## Arguments

- object:

  A `"mixedsubjects_1pl_fit"` object from
  [`fit_mixed_subjects_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_1pl.md)
  or
  [`fit_mixed_subjects_mml_1pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md).

- ridge:

  Ridge regularization for Hessian inversion.

- ...:

  Unused.

## Value

A `(J+1) × (J+1)` covariance matrix. Row/column names are `"a_shared"`
and `"d_Item1"`, `"d_Item2"`, etc.

## Note

**Bread approximation.** The bread uses `avg_hessian_counts_1pl()`, the
EM complete-data Hessian for the 1PL model, rather than the Louis (1982)
marginal observed-information correction implemented for 2PL in
[`vcov_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md).
The EM bread over-states efficiency by ignoring missing information
about θ. A Louis-corrected 1PL bread is planned for a future release.
