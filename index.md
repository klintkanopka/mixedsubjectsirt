# mixedsubjectsirt

`mixedsubjectsirt` is a package that facilitates augmenting human pilot
data with LLM-generated item responses in psychometric calibration
studies. We do this by implementing the Mixed-Subject Design[^1][^2] for
latent variable measurement models. This package ports the Prediction
Powered Inference (PPI)[^3] and PPI++[^4] paradigms to EM-based
estimation procedures that don’t have the clear independent and
dependent variables usually thought of in these PPI-based workflows.
**The result is unbiased item parameter estimates with improved
precision compared to the human sample alone.**

The strength of this method is that it tunes the contribution of the
LLM-generated responses based on how informative they are. This is done
through a procedure called *power tuning* (derived from PPI++) with one
key deviation: Instead of selecting a tuning parameter to minimize the
standard errors of the estimated model parameters, we minimize *ability
risk*, a quantification of the expected measurement error in downstream
ability estimation, integrated over the assumed ability distribution.
This allows our method to target parts of the scale where reductions in
item parameter uncertainty are the most valuable, increasing operational
measurement precision. Additionally, this approach also protects users
from poor quality synthetic data by **completely disregarding
LLM-generated responses when they are poorly aligned or unhelpful**. In
these cases, item parameter estimation is automatically carried out
using only human responses. This means that whenever users are able to
produce better quality predictions (through the use of using auxiliary
data, better prompting, stronger models, or other new and unforeseen
advances in LLMs or response prediction), the utility of this method
increases in kind.

Implemented here are methods for standard dichotomous 2PL and 1PL IRT
models. There are multiple options for estimation, with the recommended
approach being Marginal Maximum Likelihood-based EM cross-fit to split
samples.[^5] Other options include approximations based upon
quadrature-based expected count regressions and iterated expected
counts. This package is under active development, with experimental
features such as per-item power tuning available for users to try.

## Installation

Interested users can install the development version using:

``` r

devtools::install_github('klintkanopka/mixedsubjectsirt')
```

[^1]: [Broska, D., Howes, M., & van Loon, A. (2025). The mixed subjects
    design: Treating large language models as potentially informative
    observations. *Sociological Methods & Research*, 54(3),
    1074-1109.](https://journals.sagepub.com/doi/abs/10.1177/00491241251326865)

[^2]: [Van Loon, A., & Kanopka, K. (2026). Using large language models
    as a source of human behavioral data in social science experiments.
    *SocArXiv preprint*.](https://osf.io/preprints/socarxiv/y74mu)

[^3]: [Angelopoulos, A. N., Bates, S., Fannjiang, C., Jordan, M. I., &
    Zrnic, T. (2023). Prediction-powered inference. *Science,
    382*(6671),
    669-674.](https://www.science.org/doi/abs/10.1126/science.adi6000)

[^4]: [Angelopoulos, A. N., Duchi, J. C., & Zrnic, T. (2023). Ppi++:
    Efficient prediction-powered inference. *arXiv preprint
    arXiv:2311.01453*.](https://arxiv.org/abs/2311.01453)

[^5]: [Mani, P., Xu, P., Lipton, Z. C., & Oberst, M. (2025). No free
    lunch: Non-asymptotic analysis of prediction-powered inference.
    *arXiv preprint
    arXiv:2505.20178*.](https://arxiv.org/abs/2505.20178)
