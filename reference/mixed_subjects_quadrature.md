# Convert responses to quadrature form

Fits or accepts a 2PL model, computes posterior quadrature weights for
each subject, and returns expected counts for mixed-subjects
calibration. This is a lower-level helper; most analyses should call
[`fit_mixed_subjects()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
or
[`fit_mixed_subjects_split()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md).

## Usage

``` r
mixed_subjects_quadrature(
  resp,
  N_quad = 31,
  eps = 1e-15,
  iterlim = 1e+05,
  irt_pars = NULL,
  item_pars = NULL,
  quadrature = NULL,
  link_method = "mean_mean",
  ...
)
```

## Arguments

- resp:

  A response matrix with rows for subjects and columns for items.

- N_quad:

  Number of quadrature nodes to compute. Kept for backward
  compatibility; prefer `n_quad` in new code.

- eps:

  Retained for backward compatibility. Stable log computations are used
  instead of probability clipping.

- iterlim:

  Maximum number of Newton-Raphson iterations passed to
  [`rmutil::gauss.hermite()`](https://rdrr.io/pkg/rmutil/man/gauss.hermite.html).

- irt_pars:

  Optional target item parameters for mean-mean linking. This argument
  is kept for backward compatibility with earlier package versions.

- item_pars:

  Optional item parameters. If omitted, a 2PL model is fit to `resp`
  using
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md).

- quadrature:

  Optional quadrature grid with `theta` and `weight` columns.

- link_method:

  Linking method used when `irt_pars` is supplied.

- ...:

  Additional arguments passed to
  [`fit_2pl()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  when `item_pars` is omitted.

## Value

A list with `quad`, `counts`, `weights`, `irt_pars`, `quadrature`, and
`theta`.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
names(q)
#> [1] "quad"       "counts"     "weights"    "resp"       "irt_pars"  
#> [6] "quadrature" "theta"     
```
