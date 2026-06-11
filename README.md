
# mixedsubjectsirt

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/klintkanopka/mixedsubjectsirt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/klintkanopka/mixedsubjectsirt/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`mixedsubjectsirt` is a package that facilitates augmenting human pilot data with LLM-generated item responses in psychometric calibration studies. We do this by implementing the Mixed-Subject Design[^1][^2] for latent variable measurement models. This package ports the Prediction Powered Inference (PPI)[^3] and PPI++[^4] paradigms to EM-based estimation procedures that don't have the clear independent and dependent variables usually thought of in these PPI-based workflows. The goal is item-parameter estimates that retain the human-data target while using synthetic responses only when they appear informative. This works because the estimator is anchored to the human responses and the LLM contribution is down-weighted when it does not help.

The strength of this method is that it tunes the contribution of the LLM-generated responses based on how informative they are. This is done through a procedure called _power tuning_ (derived from PPI++) with one key deviation: Instead of selecting a tuning parameter to minimize the standard errors of the estimated model parameters, we minimize _ability risk_, a quantification of the expected measurement error in downstream ability estimation, integrated over the assumed ability distribution. This allows our method to target parts of the scale where reductions in item parameter uncertainty are the most valuable, increasing operational measurement precision. Additionally, this approach guards against poor-quality synthetic data: ability-risk tuning can shrink the tuning parameter λ toward zero when synthetic responses do not improve downstream scoring precision, so estimation leans on the human responses where the LLM is uninformative. This means that whenever users are able to produce better quality predictions (through the use of using auxiliary data, better prompting, stronger models, or other new and unforeseen advances in LLMs or response prediction), the utility of this method increases in kind.

Implemented here are methods for standard dichotomous 2PL and 1PL IRT models. There are multiple options for estimation, with the recommended approach being Marginal Maximum Likelihood-based EM cross-fit to split samples.[^5] Other options include approximations based upon quadrature-based expected count regressions and iterated expected counts. This package is under active development, with experimental features such as per-item power tuning available for users to try.

## What should I use?

| Goal | Recommended function |
| --- | --- |
| Complete calibration workflow, including cross-fit λ tuning | `tune_lambda_ability_risk_crossfit()` |
| Complete workflow without cross-fit λ tuning | `tune_lambda_ability_risk()` |
| Fitting models with user-specified λ value | `fit_mixed_subjects_mml()` |
| Experimental item-specific λ tuning | `tune_lambda_ability_risk_item()` |


See the [Mixed-Subjects Workflow](https://klintkanopka.com/mixedsubjectsirt/articles/mixed-subjects-workflow.html) vignette for the recommended end-to-end pipeline.

## Installation

Interested users can install the development version using:

``` r
devtools::install_github('klintkanopka/mixedsubjectsirt')
```

[^1]: [Broska, D., Howes, M., & van Loon, A. (2025). The mixed subjects design: Treating large language models as potentially informative observations. _Sociological Methods & Research_, 54(3), 1074-1109.](https://journals.sagepub.com/doi/abs/10.1177/00491241251326865)
[^2]: [Van Loon, A., & Kanopka, K. (2026). Using large language models as a source of human behavioral data in social science experiments. _SocArXiv preprint_.](https://osf.io/preprints/socarxiv/y74mu)
[^3]: [Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., & Zrnic, T. (2023). Prediction-powered inference. _Science, 382_(6671), 669-674.](https://www.science.org/doi/abs/10.1126/science.adi6000)
[^4]: [Angelopoulos, A. N., Duchi, J. C., & Zrnic, T. (2023). Ppi++: Efficient prediction-powered inference. _arXiv preprint arXiv:2311.01453_.](https://arxiv.org/abs/2311.01453)
[^5]: [Mani, P., Xu, P., Lipton, Z. C., & Oberst, M. (2025). No free lunch: Non-asymptotic analysis of prediction-powered inference. _arXiv preprint arXiv:2505.20178_.](https://arxiv.org/abs/2505.20178)

