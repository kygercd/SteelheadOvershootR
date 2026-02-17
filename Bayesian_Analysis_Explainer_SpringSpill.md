# Steelhead Overshoot Bayesian Analysis

## A Guide to Interpreting Bayesian Results with Frequentist Analogies

**Upper Columbia River Steelhead — Wells Dam Overshoot Study**

---

## Table of Contents

1. Introduction and Study Overview
2. The Model: Bayesian Hierarchical Logistic Regression
3. Priors — Starting Assumptions
4. Posterior Distributions — The Main Results
5. Posterior Mean (β) and Standard Deviation
6. 95% Highest Density Interval (HDI)
7. Probability of Direction — P(β > 0) or P(β < 0)
8. Odds Ratios from Posterior
9. Random Effects (Year Intercepts) and σ_alpha
10. MCMC Diagnostics: R-hat, ESS, Trace Plots
11. Posterior Predictive Checks
12. Prior Sensitivity Analysis
13. LOO Cross-Validation (LOO-CV)
14. Summary Comparison Table
15. How to Read the Results CSV Files

---

## 1. Introduction and Study Overview

This document explains the Bayesian statistical analysis used in the Steelhead Overshoot study at Wells Dam on the Upper Columbia River. The analysis investigates whether dam spill conditions and harvest exposure influence whether overshooting steelhead return downstream to their natal tributaries or remain at Wells Dam.

Because Bayesian statistics uses different terminology and produces different output than the frequentist methods more commonly taught in introductory statistics courses, each section below explains a component of the Bayesian results and then describes how it maps to familiar frequentist concepts such as p-values, confidence intervals, and maximum-likelihood estimates.

**Study at a Glance**

- 1,355 fish across 20 spill years (2005–2024)
- 388 fish (28.6%) returned downstream; 967 fish (71.4%) remained at Wells Dam
- 1,149 fish (84.8%) classified as harvest-exposed

---

## 2. The Model: Bayesian Hierarchical Logistic Regression

### Model Specification

The core model is a hierarchical (mixed-effects) logistic regression. For fish *i* observed in year *j*:

```
yᵢⱼ ~ Bernoulli(pᵢⱼ)
logit(pᵢⱼ) = αⱼ + β_spill × Spillⱼ + β_harvest × Harvestᵢ
αⱼ ~ Normal(μ_α, σ_α)     [year-specific random intercepts]
```

Where y = 1 if the fish returned downstream, y = 0 if it remained at Wells Dam. The αⱼ terms are year-level random intercepts that account for unmeasured annual variation in baseline return rates.

### Frequentist Analogy

This model is directly analogous to a frequentist generalized linear mixed model (GLMM) with a logit link function, as would be fitted using `glmer()` in R's lme4 package or PROC GLIMMIX in SAS. The fixed effects (spill, harvest) serve the same role, and the year-level random intercepts serve the same role as random effects in a mixed model. The key difference is in how the parameters are estimated: frequentist GLMMs use maximum likelihood or restricted maximum likelihood (REML), whereas the Bayesian approach uses Markov Chain Monte Carlo (MCMC) sampling to characterize the full posterior distribution of each parameter.

---

## 3. Priors — Starting Assumptions

### What Are Priors?

In Bayesian analysis, a prior distribution encodes what we believe about a parameter before seeing the data. The prior is then updated by the data (via the likelihood) to produce the posterior distribution. When the sample size is large, the data typically overwhelm the prior, and results will be similar regardless of prior choice. This analysis uses "weakly informative" priors — broad enough not to strongly constrain the results, but informative enough to regularize the model and help MCMC sampling.

### Priors Used in This Analysis

| Parameter | Prior | Interpretation |
|-----------|-------|----------------|
| μ_α (population intercept) | Normal(0, 2) | Allows baseline return rates from ~12% to ~88% |
| σ_α (year SD) | HalfNormal(0, 1) | Allows moderate year-to-year variation |
| β_spill | Normal(0, 1) | Centered on zero, allows moderate effects |
| β_harvest | Normal(0, 1) | Centered on zero, allows moderate effects |

### Frequentist Analogy

Frequentist statistics does not use priors. The closest analogy is regularization (penalized regression). A Normal(0, 1) prior on a coefficient is mathematically equivalent to an L2 (ridge) penalty in penalized logistic regression. In practice, with 1,355 observations, these weakly informative priors have negligible influence on the results — as confirmed by the prior sensitivity analysis, which showed that varying priors from narrow (Normal(0, 0.5)) to wide (Normal(0, 5)) changed the harvest coefficient by less than 0.3 units.

---

## 4. Posterior Distributions — The Main Results

### What Is a Posterior Distribution?

The posterior distribution is the central output of a Bayesian analysis. It represents the updated probability distribution for each parameter after combining the prior with the observed data via Bayes' theorem:

```
Posterior ∝ Likelihood × Prior
```

Rather than producing a single point estimate and a p-value, the Bayesian approach produces a full probability distribution for each parameter. This distribution tells you: given the data and model, how probable are different values of this parameter? The posterior is summarized by its mean (or median), its spread (standard deviation), and a credible interval.

### Frequentist Analogy

In frequentist statistics, the analogy to the posterior distribution is the sampling distribution of the estimator. However, the interpretation is fundamentally different. A frequentist sampling distribution describes how estimates would vary across hypothetical repeated experiments. A Bayesian posterior directly describes the probability of different parameter values given this data. In practice, when priors are weak (as in this analysis), the posterior mean will be very close to the maximum-likelihood estimate (MLE), and the posterior standard deviation will be very close to the frequentist standard error.

---

## 5. Posterior Mean (β) and Standard Deviation

### In This Analysis

Each model reports `Beta_Spill_Mean`, `Beta_Spill_SD`, `Beta_Harvest_Mean`, and `Beta_Harvest_SD`. These are the mean and standard deviation of the posterior distribution for each regression coefficient, computed from the 12,000 MCMC samples (4 chains × 3,000 post-warmup iterations).

### Example from the Spring TotalSpill Model

| Parameter | Posterior Mean | Posterior SD |
|-----------|----------------|--------------|
| β_spill (Spring TotalSpill) | 0.03 | 0.13 |
| β_harvest | −2.17 | 0.21 |

### Interpretation

The posterior mean of β_harvest = −2.17 means that, after accounting for spill and year-to-year variation, harvest-exposed fish have an estimated 2.17-unit decrease in log-odds of returning downstream compared to non-harvest-exposed fish. The posterior SD of 0.21 describes the remaining uncertainty about this estimate.

The posterior mean of β_spill = 0.03 for spring total spill is essentially zero, indicating no detectable relationship between spring spill volume and downstream return probability.

### Frequentist Analogy

The posterior mean is analogous to the maximum-likelihood estimate (MLE) or the coefficient estimate from `glmer()`. The posterior SD is analogous to the standard error (SE) of that coefficient. In frequentist output, these would appear as the "Estimate" and "Std. Error" columns in a coefficient table. With weakly informative priors and a large sample, the posterior mean ≈ MLE and the posterior SD ≈ SE. The difference is in interpretation: the frequentist SE describes sampling variability, while the Bayesian posterior SD describes remaining uncertainty about the parameter given the observed data.

---

## 6. 95% Highest Density Interval (HDI)

### What Is an HDI?

The 95% Highest Density Interval (HDI) is the narrowest interval that contains 95% of the posterior distribution's probability mass. Every point inside the HDI has higher posterior density than any point outside it. It is the Bayesian analogue of a confidence interval.

### Example from the Results

| Parameter | 95% HDI |
|-----------|---------|
| β_harvest | [−2.59, −1.77] |
| β_spill (Spring TotalSpill) | [−0.22, 0.28] |

The harvest HDI of [−2.59, −1.77] is entirely below zero, indicating strong evidence that harvest exposure reduces the probability of returning downstream. In contrast, the spring spill HDI of [−0.22, 0.28] straddles zero, indicating that the data are consistent with both positive and negative spill effects (i.e., no clear effect).

### Frequentist Analogy

The HDI is directly analogous to the 95% confidence interval (CI) from frequentist analysis, and when priors are weak the numerical values will be nearly identical. However, the interpretation differs:

- **Frequentist 95% CI**: "If we repeated the experiment many times, 95% of the calculated intervals would contain the true parameter value."
- **Bayesian 95% HDI**: "There is a 95% probability that the true parameter value lies within this interval, given the data."

The Bayesian HDI allows the intuitive interpretation that most people incorrectly apply to frequentist confidence intervals: "I am 95% sure the true value is in this range." A related Bayesian interval is the 95% credible interval (or equal-tailed interval, ETI), which clips the lowest and highest 2.5% of the posterior. The HDI is preferred here because it is the shortest interval containing 95% of the mass, which is more informative for skewed distributions.

---

## 7. Probability of Direction — P(β > 0) or P(β < 0)

### What Is It?

The probability of direction is the proportion of the posterior distribution that falls above (or below) zero. It answers: "Given the data and model, what is the probability that this effect is positive (or negative)?" This is reported in the results as `P_Spill_Positive` and `P_Harvest_Positive`.

### Example Values

| Parameter | Probability | Interpretation |
|-----------|-------------|----------------|
| P(β_harvest > 0) | 0.000 | 100% probability that harvest reduces return |
| P(β_spill > 0) for Spring TotalSpill | 0.609 | Weak, inconclusive evidence (~61% positive) |

### Frequentist Analogy

The probability of direction is loosely analogous to a one-sided p-value, but with a critical distinction in interpretation:

- **Frequentist p-value**: "The probability of observing data this extreme if the null hypothesis (β = 0) were true."
- **Bayesian P(β > 0)**: "The probability that the true effect is positive, given the data."

The Bayesian probability of direction provides the intuitive answer most researchers want: How likely is it that this effect is real and in this direction? A frequentist p-value answers a different question: How surprising is the data if there were no effect at all?

A rough (but imperfect) mapping: P(β > 0) ≈ 0.975 corresponds roughly to a two-sided p-value of ~0.05; P(β > 0) ≈ 0.995 corresponds roughly to p ≈ 0.01. However, this mapping is only approximate and depends on the prior.

---

## 8. Odds Ratios from the Posterior

### Calculation

Because the model uses a logit link, each coefficient β represents a change in log-odds. Exponentiating the posterior mean gives the odds ratio (OR):

```
OR = exp(β)
```

### Example

| Parameter | β | Odds Ratio | Interpretation |
|-----------|---|------------|----------------|
| Harvest exposure | −2.17 | 0.11 | 89% reduction in odds |
| Spring TotalSpill (1 SD increase) | 0.03 | 1.03 | 3% increase in odds (negligible) |

The harvest OR of 0.11 means that harvest-exposed fish have only 11% the odds of returning downstream compared to non-exposed fish — an 89% reduction in odds. Similarly, the entire 95% HDI can be exponentiated to obtain a 95% credible interval for the odds ratio: exp([−2.59, −1.77]) = [0.075, 0.170].

The spring spill OR of 1.03 indicates essentially no relationship — a one standard deviation increase in spring total spill volume is associated with only a 3% increase in odds of returning downstream, which is indistinguishable from no effect.

### Frequentist Analogy

Odds ratios are computed identically in frequentist logistic regression. The interpretation is the same: they describe the multiplicative change in odds for a one-unit increase in the predictor. The 95% HDI for the OR is analogous to the 95% confidence interval for the OR from a frequentist model. In both frameworks, an OR whose interval excludes 1.0 indicates a "significant" effect; an interval that includes 1.0 indicates no clear effect.

---

## 9. Random Effects (Year Intercepts) and σ_α

### What Are They?

The model includes year-specific random intercepts (αⱼ) drawn from a normal distribution: αⱼ ~ Normal(μ_α, σ_α). Each year gets its own baseline log-odds of returning downstream, and these baselines are partially pooled toward the grand mean. The parameter σ_α (reported as `Sigma_Alpha` in the results CSV) quantifies how much baseline return rates vary from year to year.

### Values from This Analysis

Across the six models, σ_α ranged from 0.52 to 0.69 on the log-odds scale. For the Spring TotalSpill model, σ_α = 0.69. This indicates substantial year-to-year variation in baseline return rates beyond what is explained by spill and harvest. On the probability scale, this translates to annual baseline return rates varying by roughly ±12–17 percentage points around the grand mean.

### Partial Pooling

The Bayesian hierarchical model implements "partial pooling," where year estimates are shrunk toward the grand mean. Years with few observations are shrunk more heavily (borrowed strength from other years), while years with many observations are shrunk less. This is a major advantage over treating year as a fixed effect, which would overfit years with small samples.

### Frequentist Analogy

The year random intercepts are directly analogous to the random effects in a frequentist GLMM (e.g., the `(1 | Year)` term in `glmer()`). The σ_α corresponds to the estimated random-effect standard deviation. Frequentist mixed models also perform partial pooling via best linear unbiased prediction (BLUP). The Bayesian approach differs in that it provides a full posterior distribution for σ_α and each αⱼ, rather than a single point estimate with an approximate standard error. This is particularly useful when the number of groups (years) is moderate (20 years here) — where frequentist estimates of variance components can be unreliable or boundary-constrained.

---

## 10. MCMC Diagnostics: R-hat, ESS, and Trace Plots

Unlike frequentist methods that find parameter estimates by optimization (maximizing the likelihood), Bayesian methods use Markov Chain Monte Carlo (MCMC) sampling to draw samples from the posterior distribution. Because MCMC is an iterative algorithm, we need diagnostics to confirm that the chains have converged and that the samples are reliable.

### R-hat (Potential Scale Reduction Factor)

R-hat compares the variance within each MCMC chain to the variance between chains. If all chains have converged to the same distribution, R-hat ≈ 1.0. Values above 1.01 suggest the chains have not converged and the results may not be trustworthy.

**This analysis**: All R-hat values were below 1.01 (range: 0.9998–1.0005 for fixed effects), indicating excellent convergence.

### Effective Sample Size (ESS)

Because MCMC samples are autocorrelated (each sample depends on the previous one), the effective sample size is smaller than the nominal number of samples. ESS estimates how many independent samples the chain is worth. Guidelines suggest ESS > 400 for reliable posterior summaries and ESS > 1,000 for stable tail estimates.

**This analysis**: Bulk ESS exceeded 3,000 for all fixed effects parameters, well above recommended minimums.

### Trace Plots

Trace plots show the sampled values of each parameter across iterations for all chains. Well-mixed chains look like "fuzzy caterpillars" — stationary, well-overlapping, and without trends or long excursions. The trace plots for this analysis show excellent mixing across all four chains.

### Divergent Transitions

Divergent transitions are numerical problems during sampling that indicate the sampler had difficulty exploring the posterior. They can signal model misspecification or difficult geometry. Zero divergent transitions were observed in this analysis.

### Frequentist Analogy

Frequentist methods generally do not require convergence diagnostics because they use closed-form solutions or deterministic optimization algorithms. The closest frequentist analogy is checking that the optimizer in `glmer()` or `optim()` has converged (e.g., no convergence warnings, the Hessian is positive definite). MCMC diagnostics serve the same purpose — confirming that the estimation procedure worked correctly — but are more involved because the algorithm is stochastic. If a frequentist optimization fails to converge, you cannot trust the estimates; similarly, if MCMC diagnostics show problems (high R-hat, low ESS, divergences), the Bayesian estimates should not be trusted.

---

## 11. Posterior Predictive Checks (PPC)

### What Are They?

Posterior predictive checks assess model adequacy by generating new datasets from the fitted model and comparing them to the observed data. If the model fits well, simulated datasets should look similar to the real data. This analysis performed three types of PPCs:

1. **Overall proportion check**: Does the model reproduce the observed overall downstream return rate (28.6%)?
2. **Year-level calibration**: For each year, does the observed return rate fall within the model's 95% posterior predictive interval?
3. **Harvest-group calibration**: Does the model reproduce the return rates separately for harvest-exposed and non-exposed fish?

### Results

All PPCs indicated adequate model fit. Year-level calibration showed years falling within the 95% posterior predictive intervals. Bayesian p-values near 0.5 confirmed that the model generates data consistent with observations.

### Frequentist Analogy

Posterior predictive checks are analogous to frequentist goodness-of-fit tests, such as the Hosmer-Lemeshow test for logistic regression or residual diagnostics (deviance residuals, Pearson residuals). The key difference is that PPCs account for parameter uncertainty by averaging over the posterior, rather than conditioning on a single point estimate. In frequentist terms, a PPC is like running a chi-squared goodness-of-fit test at every plausible parameter value and summarizing the results. This makes PPCs more robust than frequentist GOF tests, which can be overly sensitive or insensitive depending on sample size.

---

## 12. Prior Sensitivity Analysis

### What Is It?

A prior sensitivity analysis refits the model under several different prior specifications to check whether the results are driven by the choice of priors or by the data. If results are robust, conclusions should not change substantially across reasonable priors.

### Prior Specifications Tested

| Specification | β priors | σ_α prior |
|--------------|----------|-----------|
| Narrow | Normal(0, 0.5) | HalfNormal(0, 0.5) |
| Default | Normal(0, 1) | HalfNormal(0, 1) |
| Wide | Normal(0, 2) | HalfNormal(0, 2) |
| Very Wide | Normal(0, 5) | HalfNormal(0, 5) |

### Results

The harvest coefficient varied by less than 0.3 units across all four prior specifications, and the spill coefficient varied by less than 0.2 units. All specifications led to the same substantive conclusions. This confirms that the results are data-driven, not prior-driven.

### Frequentist Analogy

Frequentist analysis has no direct analogue to prior sensitivity analysis because it does not use priors. The closest analogy is running sensitivity analyses by varying model specifications (e.g., different link functions, different covariance structures, or different variable transformations) to see if conclusions are robust. Another partial analogy is varying the regularization penalty in penalized regression (e.g., lasso or ridge) and checking whether significant predictors remain significant across penalty values.

---

## 13. LOO Cross-Validation (LOO-CV)

### What Is It?

Leave-One-Out Cross-Validation (LOO-CV) estimates a model's out-of-sample predictive performance by approximating what would happen if each observation were left out of the fitting data and then predicted. In Bayesian analysis, this is efficiently approximated using Pareto Smoothed Importance Sampling (PSIS-LOO), which avoids actually refitting the model n times.

### Pareto k Diagnostics

The Pareto k statistic for each observation indicates how influential that observation is. Values above 0.7 suggest that the LOO estimate for that observation is unreliable (the observation is too influential). In this analysis, no observations had Pareto k > 0.7, indicating that the LOO estimates are reliable and no single observation is unduly influential.

### Frequentist Analogy

LOO-CV is directly analogous to the frequentist concept of cross-validation and, more specifically, to information criteria used for model comparison:

| Bayesian | Frequentist Analogue |
|----------|---------------------|
| LOO-CV (ELPD) | AIC, BIC, or k-fold CV |
| Pareto k diagnostic | Cook's distance (influential observations) |

LOO-CV is generally considered superior to AIC/BIC because it accounts for the full posterior uncertainty rather than relying on asymptotic approximations. The Pareto k diagnostic serves the same purpose as Cook's distance — identifying influential observations — but within the Bayesian framework.

---

## 14. Summary Comparison Table

The following table provides a comprehensive mapping between Bayesian results components and their frequentist analogues.

| Bayesian Concept | Frequentist Analogue | Notes |
|-----------------|---------------------|-------|
| Posterior mean | MLE / coefficient estimate | Nearly identical with weak priors |
| Posterior SD | Standard error | Nearly identical with weak priors |
| 95% HDI | 95% confidence interval | Different interpretation |
| P(β > 0) | One-sided p-value (loosely) | Different interpretation |
| Prior | Regularization penalty | Mathematical equivalence |
| Random effects | BLUP / random effects | Full posterior vs. point estimate |
| σ_α | Random effect SD | Full posterior vs. point estimate |
| R-hat | Convergence check | MCMC-specific |
| ESS | N/A | MCMC-specific |
| PPC | Hosmer-Lemeshow, residuals | Accounts for parameter uncertainty |
| LOO-CV | AIC, BIC, cross-validation | Accounts for parameter uncertainty |

---

## 15. How to Read the Results CSV Files

### Main Results: `full_hierarchical_results_with_harvest_R.csv`

This file contains one row per model (6 models total). The columns are:

| Column | Description |
|--------|-------------|
| Window | Seasonal window (Fall or Spring) |
| Variable | Spill metric (SpillHours, SpillDays, TotalSpill) |
| Beta_Spill_Mean | Posterior mean of spill coefficient |
| Beta_Spill_SD | Posterior SD of spill coefficient |
| Spill_HDI_Lower | Lower bound of 95% HDI for spill |
| Spill_HDI_Upper | Upper bound of 95% HDI for spill |
| P_Spill_Positive | P(β_spill > 0) |
| Beta_Harvest_Mean | Posterior mean of harvest coefficient |
| Beta_Harvest_SD | Posterior SD of harvest coefficient |
| Harvest_HDI_Lower | Lower bound of 95% HDI for harvest |
| Harvest_HDI_Upper | Upper bound of 95% HDI for harvest |
| P_Harvest_Positive | P(β_harvest > 0) |
| Sigma_Alpha | Year random effect SD |

### Quick Interpretation Guide

- **If HDI excludes zero**: strong evidence of an effect (analogous to p < 0.05).
- **If P_Harvest_Positive = 0.000**: essentially certain that harvest reduces return probability.
- **If P_Spill_Positive is near 0.5**: no directional evidence for a spill effect.
- **To get an odds ratio**: compute exp(Beta_Mean). For harvest: exp(−2.17) = 0.11.
- **To get an OR confidence interval**: compute exp(HDI_Lower) and exp(HDI_Upper).
- **Sigma_Alpha > 0** indicates meaningful year-to-year variation beyond what the fixed effects explain.

### Example: Reading Spring TotalSpill Results

From the CSV:
- β_spill = 0.03, HDI = [−0.22, 0.28], P(β > 0) = 0.61 → **No spill effect**
- β_harvest = −2.17, HDI = [−2.59, −1.77], P(β > 0) = 0.00 → **Strong harvest effect**
- Harvest OR = exp(−2.17) = 0.11 → **89% reduction in odds**

---

*This document was prepared as a companion guide to the Steelhead Overshoot Bayesian Analysis at Wells Dam.*

*Analysis conducted February 2026*

*Statistical software: R 4.3.2, brms 2.21 (Stan backend)*
