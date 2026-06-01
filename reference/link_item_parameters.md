# Link item parameters onto a target scale

Applies mean-mean linking to express source item parameters on the scale
of a target calibration. Both parameter sets must be in slope-intercept
form for the model `plogis(d + a * theta)`.

## Usage

``` r
link_item_parameters(source, target, method = c("mean_mean", "none"))
```

## Arguments

- source:

  Item parameters to transform. A matrix or data frame with columns
  `a`/`a1` and `d`, or a fitted `mirt` model.

- target:

  Item parameters defining the target scale. Uses the same accepted
  formats as `source`.

- method:

  Linking method. Currently `"mean_mean"` and `"none"` are supported.

## Value

A list with transformed `pars`, linking constants `A` and `B`, and the
selected `method`.

## Details

If `theta_target = A * theta_source + B`, then source parameters
transform as `a_target = a_source / A` and
`b_target = A * b_source + B`, with `d_target = -a_target * b_target`.
Mean-mean linking chooses `A` and `B` so that the transformed source
parameters match the target mean discrimination and mean difficulty.

## Examples

``` r
source <- data.frame(a = c(0.8, 1.2), d = c(-0.2, 0.5))
target <- data.frame(a = c(1.0, 1.5), d = c(-0.1, 0.4))
link_item_parameters(source, target)$pars
#>   item   a          d          b
#> 1    1 1.0 -0.1833333  0.1833333
#> 2    2 1.5  0.5250000 -0.3500000
```
