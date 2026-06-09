# Mixed-subjects objective function

Evaluates the rectified mixed-subjects loss for 2PL item parameters. The
parameter vector must contain all discriminations first, followed by all
intercepts. The response probability is `plogis(d + a * theta)`.

## Usage

``` r
mixed_subjects_loss(pars, q_observed, q_predicted, q_llm, lambda = 0)
```

## Arguments

- pars:

  Numeric vector of item parameters: all discriminations `a` followed by
  all intercepts `d`.

- q_observed:

  Quadrature summary for observed human responses, usually returned by
  [`mixed_subjects_quadrature()`](http://klintkanopka.com/mixedsubjectsirt/reference/mixed_subjects_quadrature.md).

- q_predicted:

  Quadrature summary for LLM responses/predictions on the same labeled
  human subjects.

- q_llm:

  Quadrature summary for generated or unlabeled LLM responses.

- lambda:

  Power-tuning parameter in `[0, 1]`.

## Value

A scalar loss.

## Details

The objective is
`L_observed(pars) + lambda * (L_generated(pars) - L_predicted(pars))`.
Setting `lambda = 0` gives the human-only expected-count objective.

## Examples

``` r
pars <- data.frame(a = c(1, 1.2), d = c(0, -0.5))
resp <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
q <- mixed_subjects_quadrature(resp, item_pars = pars, N_quad = 5)
mixed_subjects_loss(c(pars$a, pars$d), q, q, q, lambda = 0.5)
#> [1] 1.57952
```
