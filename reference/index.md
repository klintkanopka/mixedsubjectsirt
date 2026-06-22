# Package index

## Fitting mixed-subjects models

Estimators for mixed-subjects IRT calibration.
[`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
is the recommended marginal-likelihood estimator;
[`fit_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
and
[`fit_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
fit human-only baselines using a `mirt` backend

- [`fit_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects.md)
  : Fit a mixed-subjects 2PL calibration
- [`fit_mixed_subjects_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_1pl.md)
  : Fit a mixed-subjects 1PL calibration (frozen expected-count)
- [`fit_mixed_subjects_from_quadrature()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_from_quadrature.md)
  : Fit from precomputed quadrature summaries
- [`fit_mixed_subjects_iterative()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_iterative.md)
  : Fit a mixed-subjects 2PL calibration with iterative EM
- [`fit_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md)
  : Fit a mixed-subjects 2PL calibration via marginal maximum likelihood
- [`fit_mixed_subjects_mml_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml_1pl.md)
  : Fit a mixed-subjects 1PL calibration via marginal maximum likelihood
- [`fit_mixed_subjects_split()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_split.md)
  : Fit a split-sample mixed-subjects 2PL calibration
- [`fit_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_2pl.md)
  : Fit a unidimensional 2PL IRT model
- [`fit_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/fit_1pl.md)
  : Fit a 1PL (one-parameter logistic) model

## Choosing lambda

Power-tuning the mixed-subjects correction.
[`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
is the recommended practical criterion (downstream ability-score risk);
[`tune_lambda_ppi_score()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
is a theoretical PPI++ diagnostic;
[`diagnose_lambda_grid()`](https://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
is a sensitivity check

- [`tune_lambda_ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
  : Tune lambda by downstream ability-score risk
- [`tune_lambda_ability_risk_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_1pl.md)
  : Tune lambda by downstream ability-score risk for a 1PL model
- [`tune_lambda_ability_risk_crossfit()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_crossfit.md)
  : Cross-fit ability-score-risk lambda tuning
- [`tune_lambda_ability_risk_item()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk_item.md)
  : Per-item ability-risk lambda tuning via coordinate descent
- [`tune_lambda_ppi_score()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
  : Plug-in PPI++ optimal tuning parameter
- [`tune_lambda_ppi_score_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score_1pl.md)
  : Plug-in PPI++ optimal tuning parameter for a 1PL model
- [`tune_lambda_ppi_score_item()`](https://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score_item.md)
  : Per-item PPI++ optimal tuning parameters
- [`diagnose_lambda_grid()`](https://klintkanopka.com/mixedsubjectsirt/reference/diagnose_lambda_grid.md)
  : Diagnose lambda values over a grid

## Covariance and uncertainty

Sandwich covariance for fitted models, used through the
[`vcov()`](https://rdrr.io/r/stats/vcov.html) S3 method

- [`vcov_mixed_subjects()`](https://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects.md)
  : Sandwich covariance for a mixed-subjects fit
- [`vcov_mixed_subjects_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_1pl.md)
  : Sandwich covariance for a 1PL mixed-subjects fit
- [`vcov_mixed_subjects_mml()`](https://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md)
  : Marginal-MML sandwich covariance for a mixed-subjects fit

## Ability scoring and risk

Estimating abilities and propagating item-parameter uncertainty into
scores

- [`score_theta()`](https://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md)
  : Estimate ability scores from a 2PL calibration
- [`ability_gradient()`](https://klintkanopka.com/mixedsubjectsirt/reference/ability_gradient.md)
  : Gradient of ML ability scores with respect to item parameters
- [`ability_gradient_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/ability_gradient_1pl.md)
  : Gradient of ML ability scores w.r.t. 1PL item parameters
- [`ability_risk()`](https://klintkanopka.com/mixedsubjectsirt/reference/ability_risk.md)
  : Propagated ability risk from item-parameter uncertainty
- [`ability_risk_1pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/ability_risk_1pl.md)
  : Propagated ability risk for a 1PL fit

## Scale linking

Using established linking procedures to place LLM item parameters onto
the human scale

- [`link_item_parameters()`](https://klintkanopka.com/mixedsubjectsirt/reference/link_item_parameters.md)
  : Link item parameters onto a target scale

## Simulation and low-level building blocks

Data simulation, Gauss-Hermite quadrature, posterior weights, and
expected-count summaries

- [`simulate_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/simulate_2pl.md)
  : Simulate 2PL item responses
- [`make_quadrature()`](https://klintkanopka.com/mixedsubjectsirt/reference/make_quadrature.md)
  : Create a standard-normal Gauss-Hermite quadrature grid
- [`posterior_weights_2pl()`](https://klintkanopka.com/mixedsubjectsirt/reference/posterior_weights_2pl.md)
  : Compute posterior quadrature weights for a 2PL model
- [`summarize_expected_counts()`](https://klintkanopka.com/mixedsubjectsirt/reference/summarize_expected_counts.md)
  : Summarize response data as expected quadrature counts
- [`mixed_subjects_quadrature()`](https://klintkanopka.com/mixedsubjectsirt/reference/mixed_subjects_quadrature.md)
  : Convert responses to quadrature form
- [`mixed_subjects_loss()`](https://klintkanopka.com/mixedsubjectsirt/reference/mixed_subjects_loss.md)
  : Mixed-subjects objective function
