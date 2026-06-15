# Understanding Ability-Risk Tuning

## Why this vignette exists

This vignette is meant to alleviate some common misunderstandings of
mixed-subjects IRT calibration. While the core idea of augmenting human
data with LLM-generated data and estimating a $`\lambda`$ value to tune
how much you want the LLM-generated to contribute is straightforward,
the nuances of the ability-risk objective are not. In this vignette, we
will go into the data requirements for mixed-subject calibrations and
the ability-risk objective itself to help set users up for success. If
there is one major takeaway from this document, it is that the tuning
parameter $`\lambda`$ is not asking whether LLM-generated responses look
human in the aggregate. It is, instead, asking whether an LLM-based
response-generation procedure can predict the **row-level** human
response structure enough to reduce downstream ability-estimation error.
The most important object is therefore the **paired prediction matrix**,
$`P`$.

## Key Intuition

### The three response matrices

Mixed-subjects IRT requires three item response matrices. Let $`J`$ be
the number of items, $`n`$ be the number of observed human respondents,
and $`N`$ be the number of additional generated respondents.

#### 1. Observed human responses: $`O`$

This is the real human pilot calibration response matrix, with structure

``` math
O \in \{0,1\}^{n \times J}.
```

Each row is the full observed response string from one human respondent.
Each column corresponds to an item. Each entry is an observed,
dichotomously scored, response.

#### 2. Paired LLM-predicted human responses: $`P`$

This is the LLM-predicted response matrix for the **same human
respondents** in $`O`$, with structure

``` math
P \in \{0,1\}^{n \times J}.
```

Pay special attention to the requirement of “same human respondents.”
Row $`i`$ in $`P`$ must be predicted responses for the respondent in row
$`i`$ in $`O`$. Column $`j`$ in $`P`$ must correspond to column $`j`$ in
$`O`$. This is the diagnostic matrix. It tells the method whether the
LLM response-generation procedure is informative about the human
response process and defines the magnitude and confidence in the
mixed-subjects correction.

Importantly, this means that when generating $`P`$, users need to
transmit both the content of the items being responded to and some sort
of information about the rows of $`O`$ to the LLM generating the
predicted responses. This can be in the form of narrative or covariate
information about the respondents, held out item responses, or something
else.

#### 3. Additional LLM-generated responses: $`G`$

This is the larger synthetic or LLM-generated response matrix, with
structure

``` math
G \in \{0,1\}^{N \times J}.
```

Typically, $`N \gg n`$ to maximize the potential improvement in
post-calibration ability estimation precision. The rows in $`G`$ are
*not* paired with rows in $`O`$, but instead additional generated
respondents. However, the crucial requirement is that $`G`$ is meant to
be sampled from the same distribution as $`P`$, meaning it should be
created using a procedure that is as close as possible to the procedure
used to create $`P`$ and some amount of information about the ability
distribution in $`P`$.

A useful way to remember the design is:

| Matrix | Shape | Rows | Purpose |
|----|----|----|----|
| $`O`$ | $`n \times J`$ | observed human respondents | anchor the calibration target to the human population |
| $`P`$ | $`n \times J`$ | LLM predictions for those same human respondents | estimate how the LLM procedure agrees with or deviates from human responses |
| $`G`$ | $`N \times J`$ | additional LLM-generated rows | provide extra precision derived from synthetic information after correction |

Another helpful rule of thumb is that if $`P`$ is not row-aligned with
$`O`$, you should expect $`\lambda \to 0`$.

### The mixed-subjects IRT objective

The recommended scalar-$`\lambda`$ workflow in `mixedsubjectsirt` uses a
marginal maximum likelihood (MML) objective of the form

``` math
\hat\gamma_\lambda =\arg\min_\gamma\big[L_O^{\mathrm{marg}}(\gamma)+\lambda\big(L_G^{\mathrm{marg}}(\gamma)-L_P^{\mathrm{marg}}(\gamma)\big)\big].
```

Here $`\gamma = \{a_1, \ldots, a_J, d_1, \ldots, d_J \}`$ is the vector
of item parameters. For a 2PL model, we write item $`j`$’s response
probability as

``` math
p_j(\theta;\gamma_j) =\mathrm{logit}^{-1}(d_j + a_j \theta),
```

where $`d_j`$ is the intercept and $`a_j`$ is the discrimination.

The three pieces of the objective each have different jobs in the
mixed-subjects loss:
``` math
L_O^{\mathrm{marg}}(\gamma)
```
is the usual human-data marginal IRT likelihood evaluated at $`\gamma`$.
This term anchors the calibration to humans.
``` math
L_G^{\mathrm{marg}}(\gamma)
```
is the likelihood contribution from the additional generated response
rows.
``` math
L_P^{\mathrm{marg}}(\gamma)
```
is the correction term. It estimates what the same LLM-generation
procedure says about the humans whose actual responses are observed.

The term $`L_G^{\mathrm{marg}} - L_P^{\mathrm{marg}}`$ is the
prediction-powered correction. The generated matrix $`G`$ adds
information, while the paired prediction matrix $`P`$ allows us to
estimate the LLM procedure’s bias, noise, and covariance with the human
response process.

### What lambda is learning

The tuning parameter $`\lambda`$ weights how much the LLM-generated
responses are allowed to contribute to parameter estimation. When
$`\lambda = 0`$,

``` math
L_O^{\mathrm{marg}}(\gamma)
+
\lambda
\{L_G^{\mathrm{marg}}(\gamma)-L_P^{\mathrm{marg}}(\gamma)\}
=
L_O^{\mathrm{marg}}(\gamma),
```

so the method falls back to human-only calibration. When
$`\lambda > 0`$, the method borrows information from $`G`$, corrected by
$`P`$.

This means $`\lambda=0`$ is not a failure of the software, instead it is
a protective outcome that occurs when the LLM responses do not reduce
expected ability-estimation error. A high $`\lambda`$ requires more than
plausible synthetic response rows. It requires the paired LLM
predictions in $`P`$ to track the human responses in $`O`$ in a way that
is consistent across respondents.

### Ability-risk tuning

There are two related (but distinct) tuning ideas in the package. The
original PPI++ paper minimizes the standard errors of the estimated
parameters (the trace of the covariance matrix). The function
[`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
implements this.

The function
[`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md)
asks a more practical psychometric question: Which value of $`\lambda`$
minimizes *expected downstream ability-estimation error*?

The approximate target is
``` math
\widehat R(\lambda)
=
\frac{1}{M}
\sum_{m=1}^{M}
g_m^\top
\widehat\Sigma_{\gamma,\lambda}
g_m,
```
where:

- $`y_m`$ is a target response pattern;
- $`\hat\theta(y_m;\hat\gamma_\lambda)`$ is the ability estimate
  produced from that response pattern;
- $`g_m = \nabla_\gamma \hat\theta(y_m;\hat\gamma_\lambda)`$ is gradient
  of the ability estimate, which captures the sensitivity of the ability
  estimate to item-parameter error;
- $`\widehat\Sigma_{\gamma,\lambda}`$ is the full covariance matrix of
  the item-parameter estimates under the mixed-subjects fit.

$`g_m^\top \widehat\Sigma_{\gamma,\lambda} g_m`$ thus propagates
item-parameter uncertainty and covariance structure into
ability-estimation uncertainty.

This matters because ability-risk tuning is not the same as minimizing
average item-parameter standard errors. The off-diagonal elements of the
item-parameter covariance matrix now matter, describing how uncertainty
in one item parameter moves with uncertainty in another. Some covariance
patterns may cancel out for ability scoring; others may distort the
scale in high-information regions. Ability-risk tuning weights these
covariance patterns by their downstream impact on $`\hat\theta`$,
averaged over an expected ability distribution.

### Why row alignment matters

The easiest way to understand $`\lambda`$ is to compare three cases.

#### Case A: perfect paired prediction

Suppose

``` math
P = O.
```

Here, the paired LLM prediction is exactly equal to the human response
matrix. From the perspective of maximizing the contribution of $`G`$ to
estimation, this is the best possible version of $`P`$. In this case,
$`\lambda`$ should be large, though not necessarily $`\lambda = 1`$. The
finite-$`N`$ benchmark is

``` math
\lambda_{\max}
=
\frac{1}{1+n/N}
=
\frac{N}{n+N}.
```

If $`n=400`$ and $`N=1200`$, we see

``` math
\lambda_{\max} = \frac{1200}{400+1200}=0.75.
```

So even when $`P=O`$, if $`N \not\gg n`$, $`\lambda < 1`$ should be
expected.

#### Case B: row-shuffled perfect predictions

Now suppose

``` math
P = \mathrm{shuffle\_rows}(O).
```

The marginal item means are identical. The total-score distribution is
identical. The item difficulty information is identical. The only thing
that’s changed is row $`i`$ in $`P`$ is no longer a prediction for row
$`i`$ in $`O`$. This should produce $`\lambda \approx 0`$, because the
row-aligned covariance structure has been destroyed. To observe this for
both $`\lambda`$ tuning objectives, you can run:

``` r

predicted_shuffled <- observed[sample(nrow(observed)), ]

lambda_shuffled_ppi <- tune_lambda_ppi_score(
  observed    = observed,
  predicted   = predicted_shuffled,
  item_pars   = human_fit$pars,
  n_generated = nrow(generated)
)

lambda_shuffled_ppi$lambda

lambda_shuffled_ability_risk <- tune_lambda_ability_risk(
  observed    = observed,
  predicted   = predicted_shuffled,
  item_pars   = human_fit$pars,
  n_generated = nrow(generated)
)

lambda_shuffled_ability_risk$lambda
```

#### Case C: same DGP, fresh Bernoulli draw

Suppose both $`O`$ and $`P`$ are generated from the same IRT model:

``` math
O_{ij} \mid \theta_i \sim \mathrm{Bernoulli}\{p_j(\theta_i)\},
```

``` math
P_{ij} \mid \theta_i \sim \mathrm{Bernoulli}\{p_j(\theta_i)\}.
```

This is “same DGP,” but $`P`$ is still a fresh stochastic response. It
is not the same as $`O`$. The two matrices share person ability and item
parameters, but they do not share response noise.

For a single item,

``` math
\operatorname{Cov}(O_{ij},P_{ij})
=
\operatorname{Var}_\theta[p_j(\theta)].
```

But

``` math
\operatorname{Var}(P_{ij})
=
\operatorname{Var}_\theta[p_j(\theta)]
+
\mathbb{E}_\theta[p_j(\theta)\{1-p_j(\theta)\}].
```

The Bernoulli noise in $`P`$ dilutes the control-variate signal. As a
result, a fresh same-DGP draw may produce a modest $`\lambda`$, not a
large one. This distinction is important: Merely producing the same item
parameters does not imply strong paired prediction.

### What kind of LLM data produces higher lambda?

A useful $`P`$ matrix has to predict row-level response structure. A
good $`P`$ should have these properties:

1.  **Same rows as $`O`$**: row $`i`$ in $`P`$ predicts row $`i`$ in
    $`O`$.
2.  **Same item columns as $`O`$**.
3.  **Target response not leaked**: when predicting $`P_{ij}`$, the
    prompt must not include $`O_{ij}`$.
4.  **Construct-relevant respondent information**: covariates or context
    should help infer the respondent’s likely response.
5.  **Within-person response structure**: the LLM should be able to
    infer something about the respondent and their knowledge or ability
    level from other responses or covariates.
6.  **Same procedure for $`P`$ and $`G`$**: the generated matrix should
    be produced by an analogous procedure to the paired predicted
    matrix.

#### One approach to row alignment: leave-one-item-out prediction

When you have a human response matrix $`O`$, an approach to build $`P`$
is leave-one-item-out response prediction. For each respondent $`i`$ and
item $`j`$:

1.  Mask response $`O_{ij}`$.
2.  Give the LLM the text and responses for the other items
    $`O_{i,-j}`$.
3.  Give the LLM the text for item $`j`$.
4.  Ask it to predict the response to item $`j`$.
5.  Store the result as $`P_{ij}`$.

#### Another approach: covariate-based prediction

If you have construct-relevant covariates, you can build $`P`$ without
using other item responses or augment the LOO-response prediction
approach outlined above.

Example covariates include:

- grade level;
- age;
- prior achievement;
- language background;
- prior placement scores;
- response-time or engagement indicators;
- classroom, school, or instructional context;
- demographic variables, where appropriate and ethically justified.

#### Something that probably won’t work: item-text-only generation

Item-text-only generation usually predicts column properties (like item
parameters and relative item spacings), not row-aligned responses.
Row-aligned responses are important because they allow the method to
link the underlying ability distributions between human respondents and
LLM-generated respondents. This is why an item-text-only approach may
produce synthetic data that looks plausible in aggregate, but still
produces $`\lambda=0`$.

### How to generate $`G`$

The generated matrix $`G`$ should be produced using a procedure that
mirrors the procedure used to produce $`P`$.

If $`P`$ is generated with leave-one-item-out prompts, then $`G`$ should
also be generated with leave-one-item-out or masked-item prompts.

One possible procedure:

1.  Sample or resample a respondent profile.
2.  Create a partial response context.
3.  Mask one target item.
4.  Ask the LLM to predict the masked response.
5.  Repeat until a full generated response row is built.

The most important rule is that $`G`$ and $`P`$ should be generated by
the same prediction mechanism. If $`P`$ is generated with row-aligned
covariates and response history, but $`G`$ is generated by asking the
LLM to invent full response strings from scratch, then
$`L_G^{\mathrm{marg}} - L_P^{\mathrm{marg}}`$ may not be zero in
expectation, producing asymptotically biased parameter estimates.

### Summary

Mixed-subjects IRT is useful when the LLM response-generation procedure
captures respondent-level response structure. The generated matrix $`G`$
matters, but the paired matrix $`P`$ is what lets the method learn
whether $`G`$ should be trusted. The key understanding is that $`G`$
supplies extra synthetic rows, but $`P`$ tells us how much to trust
them.

If $`P`$ is row-aligned and predictive of $`O`$, $`\lambda`$ can be
positive and the generated data can improve calibration. If $`P`$ only
reproduces marginal item difficulty, or if its rows are not aligned with
$`O`$, then $`\lambda`$ should shrink toward zero. This is a *feature*,
not a bug.

## Technical Explanation

The [Choosing
Lambda](http://klintkanopka.com/mixedsubjectsirt/articles/lambda-tuning.md)
vignette explains which tuning function to use and when. This section
derives the mathematics those functions implement: how item-parameter
uncertainty is propagated into ability scores, and why minimizing the
resulting ability risk is a fundamentally different objective from the
original PPI++ trace criterion.

Throughout, let $`\gamma = \{a_1, \dots, a_J, d_1, \dots, d_J\}`$
collect the $`2J`$ item parameters of a $`J`$-item 2PL model, ordered as
all discriminations followed by all intercepts. This is the ordering
convention used by package functions like `fit$par`,
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`ability_gradient()`](http://klintkanopka.com/mixedsubjectsirt/reference/ability_gradient.md).
The item response function is again

``` math
p_{j}(\theta;\gamma_j) \;=\; \Pr(Y_j = 1 \mid \theta)
            \;=\; \operatorname{logit}^{-1}\!\big(d_j + a_j \theta\big).
```

### Overview: four objects, one objective

Ability-risk tuning chains four quantities together:

1.  The mixed-subjects estimator $`\hat\gamma(\lambda)`$, a function of
    the tuning parameter $`\lambda`$.
2.  Its sandwich covariance
    $`\Sigma_\gamma(\lambda) = \operatorname{Cov}(\hat\gamma)`$.
3.  The ability estimate $`\hat\theta_i(\gamma)`$ for a response pattern
    $`y_i`$, together with its gradient
    $`g_i = \partial \hat\theta_i / \partial \gamma`$.
4.  The propagated risk $`g_i' \Sigma_\gamma(\lambda)\, g_i`$, averaged
    over a target population.

Tuning chooses $`\lambda`$ to minimize that average. The sections below
build up each link in turn.

### 1. The estimator and its estimating equation

The mixed-subjects estimator minimizes a PPI++-style combined objective
over human (observed), paired-predicted, and generated responses,

``` math
L_\lambda(\gamma)
  \;=\; L_O(\gamma)
      \;+\; \lambda\,\big[L_G(\gamma)
                         - L_P(\gamma)\big],
```

where each $`L`$ is a marginal (or Bock–Aitkin expected-count) negative
log-likelihood. The estimator $`\hat\gamma(\lambda)`$ solves the
estimating equation

``` math
\Psi_\lambda(\gamma)
  \;=\; \psi_O(\gamma)
      + \lambda\,\big[\psi_G(\gamma) - \psi_P(\gamma)\big]
  \;=\; 0,
\qquad \psi = \nabla_\gamma L.
```

Setting $`\lambda = 0`$ recovers the human-only calibration;
$`\lambda > 0`$ borrows strength from the LLM responses while the
$`-\psi_P`$ term de-biases that contribution at the population level
(the PPI++ correction). The per-person score contributions to $`\psi`$
are, for item $`j`$,

``` math
s_{ij}^{a} = (y_{ij} - \bar p_{ij})\,\bar\theta_i, \qquad
s_{ij}^{d} = (y_{ij} - \bar p_{ij}),
```

evaluated under the posterior over $`\theta`$.

### 2. The sandwich covariance of $`\hat\gamma`$

$`\hat\gamma(\lambda)`$ has asymptotic covariance of the form

``` math
\Sigma_\gamma(\lambda)
  \;=\; A_\lambda^{-1}\, B_\lambda\, A_\lambda^{-1},
```

with bread $`A_\lambda = \mathbb{E}[\nabla_\gamma \Psi_\lambda]`$ and
meat $`B_\lambda = \operatorname{Cov}(\Psi_\lambda)`$.

**Bread.** Normally we would combine the three Hessians block by block,

``` math
A_\lambda \;=\; H_O
              + \lambda\,\big(H_G - H_P \big),
```

where each $`H`$ is the appropriate Hessian. Since he MML estimator
marginalizes over $`\theta`$, its bread must use Louis’s (1982)
observed-information identity
$`A^{\text{marg}} = H^{\text{comp}} - I^{\text{miss}}`$, subtracting the
missing information `louis_missing_info()` from the complete-data
Hessian. Using the complete-data Hessian alone would overstate
efficiency. This is what the package computes by default.

**Meat.** The meat is the covariance of the labeled correction plus the
independent generated contribution,

``` math
B_\lambda \;=\; \frac{1}{n}\operatorname{Cov}\!\big(S_{\text{obs}} - \lambda S_{\text{pred}}\big)
              \;+\; \frac{\lambda^2}{N}\operatorname{Cov}\!\big(S_{\text{gen}}\big),
```

with $`n`$ labeled and $`N`$ generated subjects. The
[`vcov()`](https://rdrr.io/r/stats/vcov.html) S3 method dispatches
automatically.

### 3. Ability scoring and the implicit gradient

Given item parameters $`\gamma`$, the bounded maximum-likelihood ability
estimate for response pattern $`y_i`$ solves the scoring equation

``` math
S(\theta;\gamma, y_i)
  \;=\; \sum_{j} a_j\,\big(y_{ij} - p_j(\theta)\big) \;=\; 0,
```

which
[`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md)
finds by 1-D optimization on the interval `bounds`. The risk machinery
needs the sensitivity of that solution $`\hat\theta_i`$ to the item
parameters, $`g_i = \partial \hat\theta_i / \partial \gamma`$. Because
$`\hat\theta_i`$ is defined *implicitly* by
$`S(\hat\theta_i; \gamma) = 0`$, the implicit function theorem gives

``` math
\frac{\partial \hat\theta_i}{\partial \gamma_k}
  \;=\; -\,\Big(\frac{\partial S}{\partial \theta}\Big)^{-1}
          \frac{\partial S}{\partial \gamma_k}.
```

The denominator is the (negative) test information at $`\hat\theta_i`$,

``` math
\frac{\partial S}{\partial \theta}
  \;=\; -\sum_j a_j^2\, p_j(1 - p_j),
```

and the numerators, for the discrimination and intercept of item $`j`$,
are

``` math
\frac{\partial S}{\partial a_j} = (y_{ij} - p_j) - a_j\, p_j(1 - p_j)\,\hat\theta_i,
\qquad
\frac{\partial S}{\partial d_j} = -\,a_j\, p_j(1 - p_j).
```

**Where the gradient is undefined.** The implicit-function argument
requires an *interior* optimum. At a boundary estimate (all-correct or
all-incorrect patterns push $`\hat\theta_i`$ to a bound), the score
equation does not hold and the gradient is theoretically undefined;
[`ability_gradient()`](http://klintkanopka.com/mixedsubjectsirt/reference/ability_gradient.md)
returns `NA` for those rows, and they drop out of the risk average via
`na.rm = TRUE`. Rows with vanishing test information
$`(|\partial S/\partial\theta| < \varepsilon)`$ are treated the same
way.

### 4. Delta-method propagation and the risk

With $`g_i`$ in hand, the delta method propagates item-parameter
uncertainty into the score:

``` math
\rho_i(\lambda) = \operatorname{Var}\big(\hat\theta_i\big)
  \;\approx\; g_i'\, \Sigma_\gamma(\lambda)\, g_i.
```

Averaging over a target population of $`M`$ response patterns gives the
scalar objective that the tuners minimize,

``` math
R(\lambda)
  \;=\; \frac{1}{M}\sum_{i=1}^{M} \rho_i(\lambda)
  \;=\; \mathbb{E}_{\text{target}}\!\big[\, g'\,\Sigma_\gamma(\lambda)\, g \,\big].
```

The expectation is over the *target* population’s ability distribution,
which is why `target_resp` matters: tuning for the observed calibration
sample (`target_resp = observed`) and tuning for a separate operational
scoring population generally give different $`\lambda`$.

### 5. Why this differs from the PPI++ trace objective

The original PPI++ tuning rule minimizes the trace of the item-parameter
covariance,

``` math
\lambda^{\star}_{\text{PPI}} = \arg\min_\lambda \operatorname{Tr}\big[\Sigma_\gamma(\lambda)\big],
```

which
[`tune_lambda_ppi_score()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ppi_score.md)
evaluates in closed form. Writing $`\Sigma_\gamma = (\sigma_{kl})`$, the
two objectives expand as

``` math
\operatorname{Tr}(\Sigma_\gamma) = \sum_k \sigma_{kk},
\qquad
g'\Sigma_\gamma g = \sum_{k} g_k^2\,\sigma_{kk}
                   + \sum_{k \neq l} g_k g_l\,\sigma_{kl}.
```

The trace sees only the diagonal variances and weights every parameter
equally. The ability risk weights each variance by $`g_k^2`$ (how much
that particular parameter moves the score) and uses the off-diagonal
covariances $`\sigma_{kl}`$, which encode the scale/identification
structure of the 2PL. Errors in $`a_j`$ and $`d_j`$ that are correlated
in a direction that leaves $`\hat\theta`$ unchanged are penalized by the
trace but (correctly) ignored by the ability risk. The two criteria
therefore select different $`\lambda`$ in general; use the ability risk
for operational scoring and the trace as a theoretical diagnostic.

### Summary

| Symbol | Meaning | Computed by |
|----|----|----|
| $`\hat\gamma(\lambda)`$ | Mixed-subjects item parameters | [`fit_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/fit_mixed_subjects_mml.md) |
| $`\Sigma_\gamma(\lambda)`$ | Sandwich covariance of $`\hat\gamma`$ | [`vcov()`](https://rdrr.io/r/stats/vcov.html) → [`vcov_mixed_subjects_mml()`](http://klintkanopka.com/mixedsubjectsirt/reference/vcov_mixed_subjects_mml.md) |
| $`\hat\theta_i`$ | Bounded ML ability score | [`score_theta()`](http://klintkanopka.com/mixedsubjectsirt/reference/score_theta.md) |
| $`g_i = \partial\hat\theta_i/\partial\gamma`$ | Implicit ability gradient | [`ability_gradient()`](http://klintkanopka.com/mixedsubjectsirt/reference/ability_gradient.md) |
| $`g_i'\Sigma_\gamma g_i`$ | Propagated score variance | [`ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/ability_risk.md) |
| $`R(\lambda) = \mathbb{E}[g'\Sigma_\gamma g]`$ | Ability-risk objective | [`tune_lambda_ability_risk()`](http://klintkanopka.com/mixedsubjectsirt/reference/tune_lambda_ability_risk.md) |

For the practical workflow built on these pieces — cross-fitting, target
populations, and the choice between estimators — see the [Choosing
Lambda](http://klintkanopka.com/mixedsubjectsirt/articles/lambda-tuning.md)
vignette.
