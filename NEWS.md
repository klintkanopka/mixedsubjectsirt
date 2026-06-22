# mixedsubjectsirt 1.0.0

Initial CRAN release.

* Mixed-subjects 2PL/1PL IRT calibration that augments human responses with
  LLM-generated responses through a PPI++ marginal-MML estimator
  (`fit_mixed_subjects_mml()` and relatives). The estimator is anchored to the
  human data and is asymptotically unbiased for the human item parameters at any
  tuning weight.
* Power tuning by **ability-score risk** (`tune_lambda_ability_risk()`), which
  selects the tuning weight by direct 1-D optimization of propagated
  ability-recovery risk (pass `method = "grid"` to scan a grid instead). Also
  included: a theoretical PPI++ score diagnostic (`tune_lambda_ppi_score()`),
  cross-fitted tuning (`tune_lambda_ability_risk_crossfit()`, the recommended
  workflow for reported analyses), and experimental per-item tuning
  (`tune_lambda_ability_risk_item()`). All non-experimental tuners use the
  marginal-MML estimator by default; the frozen expected-count estimator remains
  available via `fit_fn` but is discouraged.
* Louis-corrected marginal sandwich covariance through the `vcov()` S3 method
  (`vcov_mixed_subjects_mml()`), with ability scoring and item-parameter
  uncertainty propagation (`score_theta()`, `ability_risk()`).
* Vignettes covering the recommended workflow, lambda tuning, the 1PL model,
  per-item tuning, scale linking, and a simulation-validation study; an
  `R-CMD-check` GitHub Actions workflow.
* Currently `predicted` and `generated` data must be **binary 0/1 responses** in
  the high-level fitting and PPI-score functions; the low-level quadrature
  utilities accept fractional input.