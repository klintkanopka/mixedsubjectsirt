# Marginal-MML sandwich covariance for a mixed-subjects fit

Computes the full sandwich covariance for the scalar marginal-MML PPI++
estimator from
[`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md).
The bread uses Louis's (1982) observed marginal-information formula

## Usage

``` r
vcov_mixed_subjects_mml(object, ridge = 1e-08, ...)
```

## Arguments

- object:

  A scalar-lambda
  [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  fit.

- ridge:

  Ridge regularization for bread inversion.

- ...:

  Unused.

## Value

A \\2J \times 2J\\ covariance matrix with attributes `bread` and `meat`.

## Details

\$\$A\_\lambda^\mathrm{marg} = H\_\lambda^\mathrm{comp} -
I\_\lambda^\mathrm{miss}\$\$

rather than the EM/complete-data Hessian used by
[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md).
Using the complete-data Hessian as the bread for a marginal-MML
estimator would over-state efficiency by ignoring the
missing-information correction.

The meat uses the standard marginal per-person score vectors (posteriors
at the converged parameters), which is identical to
[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md).

**When is this function called automatically?** The
[`vcov()`](https://rdrr.io/r/stats/vcov.html) method for
`"mixedsubjects_fit"` objects (see
[`stats::vcov()`](https://rdrr.io/r/stats/vcov.html)) dispatches here
whenever `isTRUE(object$mml) && length(object$lambda) == 1`. For
vector-lambda fits, or for frozen expected-count fits, the existing
[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)
is used.

## See also

[`vcov_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)
for the frozen expected-count version. The internal
`louis_missing_info()` helper computes the missing-information
correction.
