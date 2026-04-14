# Loss function for parameter estimation. Generally not called directly by users.

Loss function for parameter estimation. Generally not called directly by
users.

## Usage

``` r
mixed_subjects_loss(pars, q_observed, q_predicted, q_llm, lambda = 0)
```

## Arguments

- pars:

  A vector of item parameters in slope-intercept form. Passed with
  discriminations first, location second, in item order.

- q_observed:

  A quadrature object constructed from running mixed_subjects_quadrature
  on LLM-generated item responses

- q_predicted:

  A quadrature object constructed from running mixed_subjects_quadrature
  on predicted item responses for the human calibration sample

- lambda:

  The power tuning parameter. Takes values between 0 and 1. When lambda
  is zero, disregards predicted data completely

## Value

Loss for the mixed subjects estimator at the values of `pars`
