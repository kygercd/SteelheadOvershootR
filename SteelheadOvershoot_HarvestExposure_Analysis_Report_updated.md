# Effects of Seasonal Spill and Harvest Exposure on Steelhead Overshoot Fate at Wells Dam

## Bayesian Hierarchical Analysis Using R/brms: 2005–2024

---

## Abstract

Steelhead (*Oncorhynchus mykiss*) that "overshoot" their natal tributaries and ascend past Wells Dam on the Upper Columbia River face uncertain fates, with potential consequences for population productivity and ESA recovery goals. This study employed Bayesian hierarchical logistic regression implemented in R using the brms package to evaluate whether dam spill conditions and individual fish harvest exposure influence the probability that overshooting steelhead return downstream versus remain at Wells Dam. Analysis of 1,321 PIT-tagged steelhead across 20 years (2005–2024) revealed that harvest exposure—defined as hatchery fish detected during years when harvest was open—was the dominant predictor of fate, with harvest-exposed fish showing approximately 90% lower odds of returning downstream (β = −2.28; 95% HDI: −2.72, −1.87). Neither fall nor spring spill significantly influenced downstream return probability after accounting for harvest exposure. These findings have important implications for understanding the mechanisms underlying steelhead overshoot mortality and informing recovery strategies for this ESA-threatened species.

---

## 1. Introduction

### 1.1 Steelhead Overshoot Phenomenon

Anadromous salmonids exhibit remarkable natal homing fidelity, returning to their birth streams to spawn after years at sea. However, a phenomenon known as "overshoot" occurs when adult fish migrate upstream past their natal tributary confluence before subsequently attempting to relocate olfactory cues and fall back downstream (Keefer and Caudill 2014). While overshooting is believed to be relatively rare in free-flowing rivers, it has become pervasive throughout the Columbia River basin, with some steelhead populations exhibiting overshoot rates exceeding 50% (Richins and Skalski 2018; Murdoch et al. 2022).

The Columbia River basin presents unique challenges for returning adult salmonids. Eight to nine mainstem dams impede migration between the ocean and spawning grounds, and adult fish must navigate multiple fish ladders while experiencing altered thermal and hydrological conditions. Richins and Skalski (2018) documented that steelhead overshoot rates were strongly correlated with water temperatures near natal streams—the probability of direct natal stream entry declined from over 90% to less than 25% as water temperatures increased from 10°C to 20°C. This finding supports the hypothesis that overshooting steelhead may be seeking coldwater thermal refugia upstream (WDFW 2024).

### 1.2 Conservation Significance

Upper Columbia River steelhead are listed as threatened under the Endangered Species Act, having been originally listed as endangered in 1997 and reclassified to threatened in 2006 (NOAA Fisheries 2024). The Upper Columbia Spring-run Chinook Salmon and Steelhead Recovery Plan, adopted in 2007, identifies 306 actions to improve habitat conditions, ensure safe fish passage, and allow sustainable harvest opportunities.

The fate of overshoot steelhead has direct implications for population recovery. Murdoch et al. (2022) estimated that approximately 40% of wild steelhead passing Priest Rapids Dam are overshoot fish, and critically, 40% of these overshoots fail to return downstream. Only 3% are subsequently detected on spawning grounds in the Upper Columbia. This substantial mortality represents a significant drain on population productivity. Wells Dam had the second-highest overshoot rate in the upper Columbia River, with 20% of wild steelhead passing the dam being overshoot fish.

### 1.3 Dam Passage and Spill Operations

Dam spill—the release of water over dam spillways rather than through turbines—affects multiple aspects of salmonid migration. While the effects of spill on juvenile downstream passage have been extensively studied (Muir et al. 2001), the influence of spill on adult upstream migrants and their subsequent downstream fallback behavior is less well understood. Spill operations alter flow patterns, acoustic environments, and potentially the attractiveness of fishway entrances (Keefer et al. 2008). Whether spill conditions influence the probability that overshoot steelhead successfully navigate back downstream to their natal tributaries has not been systematically evaluated.

The location of fish ladders at dams is also important—adult steelhead using ladders on the opposite shore from their natal streams overshot their home waters 18% more often than those using same-shore ladders (Richins and Skalski 2018). This suggests that physical dam configuration and hydraulic conditions play a role in overshoot behavior.

### 1.4 Harvest Exposure

Fisheries for hatchery steelhead in the Upper Columbia have been subject to varying regulations over the study period. Harvest was open from 2000–2015 and again in 2024–2025, while emergency harvest closures were implemented from 2016–2023 due to low wild steelhead returns. Mark-selective fisheries, where adipose fin-clipped hatchery fish can be retained while unmarked wild fish must be released, are the primary management tool (WDFW 2024). ESA impacts to wild steelhead in non-treaty fisheries are limited to 2% for most stocks.

Hatchery fish may also differ behaviorally from wild fish in ways that affect overshoot fate. Richins and Skalski (2018) found that segregated hatchery stocks had 15% higher overshoot rates than wild or integrated stocks, and outplanted juveniles showed 65% higher overshoot rates than wild fish. These behavioral differences, combined with direct harvest mortality during open seasons, could influence the probability that hatchery fish are detected returning downstream.

### 1.5 Study Objectives

This study tested two primary hypotheses:

1. **Spill Hypothesis**: Spill conditions at Wells Dam during the fall arrival period (August–November) or the subsequent spring holding period (January–March) influence the probability that overshoot steelhead return downstream.

2. **Harvest Exposure Hypothesis**: Hatchery steelhead detected during years when harvest was open are less likely to be detected returning downstream than other fish, either due to direct harvest mortality or behavioral differences.

We employed Bayesian hierarchical logistic regression to account for both individual-level (harvest exposure) and year-level (spill) predictors while appropriately modeling inter-annual variation in baseline return rates.

---

## 2. Methods

### 2.1 Study Area

The study focused on steelhead from the Wenatchee and Entiat River basins that were detected at Wells Dam on the Upper Columbia River. Wells Dam (river kilometer 830) is located upstream of these natal tributaries, and fish detected there have "overshot" their spawning destination. Fish fates were classified based on their last PIT-tag detection: at Wells Dam (failure to return downstream), downstream of Wells (successful fallback), or above Wells (continued upstream migration to the Okanogan or Methow basins).

### 2.2 Data Sources

**Fish Detection Data**: PIT-tag detection records for steelhead that overshot their natal tributaries and were detected at Wells Dam (2001–2025) were obtained from the Columbia Basin PIT Tag Information System (PTAGIS). Detection sites were classified by geographic zone relative to Wells Dam using an updated site classification that correctly assigned NAU (Upper Nason Creek) to the downstream zone. A small number of tags identified as potentially compromised or censored were excluded from the analysis prior to model fitting.

**Spill Data**: Hourly spill discharge measurements (KCFS) at Wells Dam from 2000–2024 were obtained from Columbia River DART (Data Access in Real Time).

**Rear Type and Harvest Status**: Fish rearing origin (hatchery or wild) was obtained from PTAGIS mark records. Annual harvest status for the study area was compiled from Washington Department of Fish and Wildlife regulations, with harvest open 2000–2015 and 2024–2025, and closed 2016–2023.

### 2.3 Outcome Variable and Fate Classification

Fish were classified into three fates based on their last detection location: (1) **Downstream**—last detected below Wells Dam, indicating successful fallback toward natal tributaries; (2) **At Wells**—last detected at Wells Dam, suggesting failure to return downstream; (3) **Above**—last detected upstream of Wells Dam in non-natal basins. For the primary analysis, we modeled the binary outcome of Downstream (1) versus At Wells (0), excluding fish with Above fates.

### 2.4 Predictor Variables

**Spill Year Assignment**: Fish were assigned to a "spill year" using a June 1 cutoff: fish first detected at Wells Dam from June of year Y through May of year Y+1 were assigned to SpillYear Y+1. Most overshooting steelhead arrive at Wells in late summer through fall, so this grouping pairs each fish with the fall spill conditions present during their arrival and spring spill conditions during the subsequent holding period.

**Seasonal Spill Metrics**: For each spill year, three metrics were calculated for two seasonal windows:
- **Fall spill** (August 1–November 30): Total spill hours, number of spill days, and cumulative spill volume (KCFS-hours)
- **Spring spill** (January 1–March 31): Same three metrics

**Harvest Exposure**: A binary individual-level covariate was created:
- **Harvest Exposure = 1**: Hatchery fish (Rear Type Code = "H") detected at Wells Dam in a year when harvest was open (2000–2015, 2024–2025)
- **Harvest Exposure = 0**: Wild fish, or hatchery fish detected in closed harvest years (2016–2023)

### 2.5 Statistical Analysis

#### 2.5.1 Bayesian Hierarchical Logistic Regression

We fitted hierarchical Bayesian logistic regression models using the **brms** package (Bürkner 2017, 2018) in R version 4.3, which provides an interface to Stan for full Bayesian inference via Hamiltonian Monte Carlo. The brms package is widely used in ecological and fisheries research for fitting complex hierarchical models (Ogle 2024; Thorson and Kristensen 2024).

The model structure was:

$$
\begin{aligned}
y_{ij} &\sim \text{Bernoulli}(p_{ij}) \\
\text{logit}(p_{ij}) &= \alpha_j + \beta_{\text{spill}} \cdot \text{Spill}_j + \beta_{\text{harvest}} \cdot \text{Harvest}_i \\
\alpha_j &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha)
\end{aligned}
$$

Where:
- $y_{ij}$ = 1 if fish $i$ in year $j$ returned downstream, 0 otherwise
- $\alpha_j$ = year-specific random intercept (partial pooling across years)
- $\text{Spill}_j$ = standardized spill metric for year $j$
- $\text{Harvest}_i$ = harvest exposure indicator for fish $i$

#### 2.5.2 Prior Specifications

Weakly informative priors were specified to regularize estimates while allowing the data to dominate inference:

- $\mu_\alpha \sim \text{Normal}(0, 2)$ — population mean intercept
- $\sigma_\alpha \sim \text{HalfNormal}(0, 1)$ — year-level standard deviation
- $\beta_{\text{spill}} \sim \text{Normal}(0, 1)$ — spill effect
- $\beta_{\text{harvest}} \sim \text{Normal}(0, 1)$ — harvest exposure effect

These priors are consistent with recommendations for logistic regression (Gelman et al. 2008) and ensure stable estimation with limited sample sizes in some years.

#### 2.5.3 MCMC Sampling

Models were fitted using 4 chains with 4,500 iterations each (1,500 warmup + 3,000 sampling), yielding 12,000 posterior samples per parameter. The target acceptance rate was set to 0.95 (`adapt_delta = 0.95`) to minimize divergent transitions. Convergence was assessed using the Gelman-Rubin diagnostic ($\hat{R}$), with values <1.01 indicating convergence, and effective sample size (ESS), with values >1,000 considered adequate.

#### 2.5.4 Models Fitted

Six hierarchical models were fitted, crossing three spill metrics with two seasonal windows:
1. Fall SpillHours + Harvest Exposure
2. Fall SpillDays + Harvest Exposure
3. Fall TotalSpill + Harvest Exposure
4. Spring SpillHours + Harvest Exposure
5. Spring SpillDays + Harvest Exposure
6. Spring TotalSpill + Harvest Exposure

#### 2.5.5 Inference

Posterior distributions were summarized using means, standard deviations, and 95% highest density intervals (HDI). The probability of direction, P(β > 0) or P(β < 0), was calculated to assess evidence for positive or negative effects. Effect sizes were also expressed as odds ratios by exponentiating the harvest coefficient: OR = exp(β_harvest).

### 2.6 Software and Reproducibility

All analyses were conducted in R version 4.3.2 (R Core Team 2023) using the following packages:
- **brms** version 2.21 (Bürkner 2017) for Bayesian model fitting
- **tidyverse** (Wickham et al. 2019) for data manipulation
- **bayesplot** (Gabry and Mahr 2022) for MCMC diagnostics
- **posterior** (Bürkner et al. 2022) for posterior summarization
- **patchwork** (Pedersen 2024) for figure composition

Model objects and analysis data are saved in `hierarchical_models_with_harvest_R.RData` for reproducibility. The complete analysis script is available as `run_full_hierarchical_analysis_with_harvest.R`.

---

## 3. Results

### 3.1 Sample Summary

The analysis included 1,321 fish across 20 spill years (2005–2024). Of these, 381 fish (28.8%) were last detected downstream of Wells Dam and 940 fish (71.2%) were last detected at Wells Dam. A total of 1,126 fish (85.2%) were classified as harvest-exposed (hatchery fish in open harvest years).

Sample sizes varied substantially across years, from a single fish in 2022, 2023, and 2024 to 286 fish in 2010 (Table 1). The proportion returning downstream ranged from 4.0% (2005) to 100% (2022, 2023, 2024), with higher return rates generally observed during the harvest closure period (2016–2023).

**Table 1. Fish fates, harvest exposure, and spill metrics by year**

| Year | Downstream | At Wells | Total | % Downstream | Harvest Exposed | Fall Spill (KCFS-hr) | Spring Spill (KCFS-hr) |
|------|------------|----------|-------|--------------|-----------------|---------------------|----------------------|
| 2005 | 1 | 24 | 25 | 4.0% | 25 | 5,546 | 564 |
| 2006 | 4 | 53 | 57 | 7.0% | 57 | 6,330 | 1,611 |
| 2007 | 8 | 99 | 107 | 7.5% | 107 | 4,682 | 6,310 |
| 2008 | 11 | 80 | 91 | 12.1% | 88 | 4,960 | 153 |
| 2009 | 22 | 108 | 130 | 16.9% | 127 | 5,951 | 256 |
| 2010 | 80 | 206 | 286 | 28.0% | 257 | 3,866 | 1 |
| 2011 | 54 | 137 | 191 | 28.3% | 178 | 4,091 | 1,892 |
| 2012 | 46 | 84 | 130 | 35.4% | 114 | 6,876 | 2,877 |
| 2013 | 36 | 78 | 114 | 31.6% | 104 | 9,667 | 1,117 |
| 2014 | 19 | 16 | 35 | 54.3% | 21 | 4,207 | 2,903 |
| 2015 | 22 | 15 | 37 | 59.5% | 13 | 4,187 | 4,606 |
| 2016 | 35 | 27 | 62 | 56.5% | 35 | 4,098 | 165 |
| 2017 | 12 | 4 | 16 | 75.0% | 0 | 3,623 | 17,709 |
| 2018 | 8 | 3 | 11 | 72.7% | 0 | 4,904 | 14,037 |
| 2019 | 4 | 3 | 7 | 57.1% | 0 | 8,117 | 847 |
| 2020 | 8 | 2 | 10 | 80.0% | 0 | 5,134 | 3,547 |
| 2021 | 8 | 1 | 9 | 88.9% | 0 | 8,756 | 2,215 |
| 2022 | 1 | 0 | 1 | 100.0% | 0 | 3,912 | 7,001 |
| 2023 | 1 | 0 | 1 | 100.0% | 0 | 4,276 | 840 |
| 2024 | 1 | 0 | 1 | 100.0% | 0 | 724 | 306 |
| **Total** | **381** | **940** | **1,321** | **28.8%** | **1,126** | — | — |

*Note: Harvest was closed 2016–2023.*

### 3.2 Model Convergence

All six models achieved excellent convergence. The Gelman-Rubin diagnostic ($\hat{R}$) was ≤1.01 for all parameters, and bulk effective sample sizes exceeded 3,000 for all fixed effects. Visual inspection of trace plots confirmed adequate mixing across chains (Figure 1).

### 3.3 Harvest Exposure Effect

The harvest exposure effect was the dominant predictor of downstream return probability, with highly consistent results across all six models (Table 2). Harvest-exposed fish had substantially lower odds of returning downstream:

- **Fall TotalSpill Model**: β_harvest = −2.29 (95% HDI: −2.71, −1.88), P(β < 0) = 1.000, Odds Ratio = 0.10
- **Spring TotalSpill Model**: β_harvest = −2.28 (95% HDI: −2.72, −1.87), P(β < 0) = 1.000, Odds Ratio = 0.10

The harvest exposure coefficient was negative and excluded zero in all models, with posterior probability of a negative effect equal to 1.000. The odds ratio of approximately 0.10 indicates that harvest-exposed fish had about 90% lower odds of returning downstream compared to fish without harvest exposure.

**Table 2. Posterior summaries for all hierarchical models**

| Window | Spill Metric | β_Spill (Mean) | Spill 95% HDI | P(β_Spill > 0) | β_Harvest (Mean) | Harvest 95% HDI | Harvest OR | σ_year |
|--------|--------------|----------------|---------------|----------------|------------------|-----------------|------------|--------|
| Fall | SpillHours | −0.30 | [−0.54, −0.07] | 0.005 | −2.18 | [−2.58, −1.77] | 0.11 | 0.49 |
| Fall | SpillDays | −0.16 | [−0.41, 0.09] | 0.093 | −2.33 | [−2.74, −1.91] | 0.10 | 0.55 |
| Fall | TotalSpill | −0.03 | [−0.39, 0.33] | 0.446 | −2.29 | [−2.71, −1.88] | 0.10 | 0.66 |
| Spring | SpillHours | 0.09 | [−0.16, 0.35] | 0.752 | −2.26 | [−2.68, −1.86] | 0.10 | 0.63 |
| Spring | SpillDays | 0.09 | [−0.20, 0.38] | 0.742 | −2.27 | [−2.68, −1.85] | 0.10 | 0.65 |
| Spring | TotalSpill | 0.03 | [−0.22, 0.28] | 0.589 | −2.28 | [−2.72, −1.87] | 0.10 | 0.65 |

*Note: β coefficients on logit scale. Spill standardized (mean 0, SD 1). OR = Odds Ratio = exp(β_harvest). σ_year = standard deviation of year random intercepts.*

### 3.4 Spill Effects

Neither fall nor spring spill showed a consistent significant effect on downstream return probability after accounting for harvest exposure (Table 2):

- **Fall Spill**: The Fall SpillHours model showed a marginally negative effect (β = −0.30, 95% HDI excludes zero), but this was not replicated by Fall SpillDays or Fall TotalSpill, whose HDIs included zero.

- **Spring Spill**: All three spring spill metrics showed estimated effects near zero with 95% HDIs that included zero. P(β > 0) ranged from 0.59 to 0.75, indicating weak and inconsistent evidence for any spill effect.

The spill effects were an order of magnitude smaller than the harvest exposure effect in all models.

### 3.5 Year Random Effects

Year-level random intercepts captured substantial inter-annual variation in baseline downstream return probability not explained by the fixed effects (Figure 2). The estimated standard deviation of year effects ranged from σ = 0.49 to σ = 0.66 across models. Years during the harvest closure period (2017–2023) generally showed higher baseline return rates (positive random effects), while early study years (2005–2009) showed lower baseline rates (negative random effects).

### 3.6 Model Comparison

All models yielded nearly identical estimates for the harvest exposure effect, demonstrating robustness to the choice of spill metric. The primary difference among models was the estimated spill coefficient, which varied depending on the specific metric used. The TotalSpill metric showed the weakest spill effect (closest to zero), while SpillHours showed somewhat stronger effects in opposite directions for fall (negative) and spring (positive).

---

## 4. Discussion

### 4.1 Harvest Exposure as the Dominant Predictor

The most striking finding of this analysis is the overwhelming dominance of harvest exposure in predicting steelhead overshoot fate. Hatchery fish detected during open harvest years had approximately 90% lower odds of returning downstream compared to other fish. This effect was:

1. **Large in magnitude**: The odds ratio of 0.10 represents nearly an order of magnitude difference in downstream return probability
2. **Highly precise**: The 95% HDI excluded zero with P(β < 0) = 1.000
3. **Robust across models**: All six models yielded virtually identical harvest effect estimates regardless of which spill metric was included

These findings are consistent with Murdoch et al. (2022), who estimated that 40% of overshoot steelhead at Priest Rapids Dam fail to return downstream. Our results suggest that a substantial portion of this non-return may be attributable to harvest exposure.

### 4.2 Potential Mechanisms for the Harvest Effect

The observed harvest exposure effect could reflect several non-mutually exclusive mechanisms:

**Direct Harvest Mortality**: During open harvest years, hatchery steelhead are legal to harvest. Anglers may encounter overshooting fish as they hold at or migrate through the study area. Mark-selective fisheries rely on adipose fin clips to identify hatchery fish, and WDFW estimates that 50–65% of hatchery fish in the Upper Columbia are adipose-clipped and thus harvestable (WDFW 2024).

**Behavioral Differences**: Hatchery and wild steelhead may differ in behaviors that influence detection probability or downstream migration success. Richins and Skalski (2018) found that hatchery stocks had 15–65% higher overshoot rates than wild fish, suggesting fundamental behavioral differences in navigation and homing.

**Temporal Confounding**: The harvest closure period (2016–2023) coincided with changes in environmental conditions, fish abundance, and potentially other management actions. While the harvest exposure covariate attempts to isolate the harvest effect, some confounding with period effects may remain.

Distinguishing among these mechanisms would require additional data on actual harvest encounters, fish movements between detection sites, and individual fate assignments. Radio or acoustic telemetry studies could provide finer-resolution movement data to assess harvest interactions.

### 4.3 Spill Effects are Negligible After Accounting for Harvest

In contrast to harvest exposure, spill effects were small and inconsistent across models and seasonal windows. The one potentially significant finding—a negative effect of Fall SpillHours—was not replicated by other fall spill metrics and may represent a spurious correlation or confounding with other factors that vary with SpillHours specifically.

This finding has important management implications. It suggests that spill operations during the fall or spring periods do not substantially influence whether overshoot steelhead return downstream. Previous hypotheses that spill creates cues or conditions facilitating fallback are not supported by this analysis. However, our analysis cannot rule out spill effects on other aspects of overshoot biology, such as migration timing, route selection, or survival during downstream passage.

### 4.4 Year-to-Year Variation

The substantial year random effects (σ ≈ 0.5–0.7) indicate that factors beyond spill and harvest exposure influence annual variation in downstream return rates. Potential sources of this variation include:

- Water temperature conditions affecting fish physiology and behavior
- Flow conditions affecting dam passage
- Fish condition and health varying by cohort
- Changes in detection infrastructure affecting apparent fates
- Stochastic variation due to small sample sizes in some years

The hierarchical model structure appropriately accounts for this variation through partial pooling, preventing any single unusual year from unduly influencing the overall inference.

### 4.5 Comparison with Previous Research

Our findings complement the growing body of research on Columbia Basin steelhead overshoot. Murdoch et al. (2022) documented the abundance and migration success of overshoot steelhead but did not evaluate harvest exposure as a predictor. Richins and Skalski (2018) characterized overshoot patterns and identified temperature as a key driver but focused on overshoot occurrence rather than subsequent fate.

This study adds to the literature by demonstrating that post-overshoot fate is strongly associated with harvest exposure status. The use of Bayesian hierarchical methods, as recommended for conservation applications (Thorson and Kristensen 2024; Scheuerell et al. 2015), provides principled uncertainty quantification and appropriate handling of the complex data structure.

### 4.6 Limitations

Several limitations should be considered when interpreting these results:

1. **Detection-based inference**: Fish fates are inferred from last detection locations, not from direct observation. Fish last detected at Wells Dam may have died there, emigrated undetected, or simply passed detection sites without being recorded.

2. **Incomplete covariate information**: Not all fish had rear type data available, potentially introducing bias if missing data is non-random.

3. **Sample size variation**: Later years (2017–2024) had substantially smaller sample sizes, increasing uncertainty for those periods.

4. **Correlation between harvest status and time**: The harvest closure period (2016–2023) represents a contiguous block of years, making it impossible to fully separate harvest effects from other temporal trends.

### 4.7 Management Implications

These findings suggest several considerations for steelhead management:

1. **Harvest exposure is associated with reduced downstream return**: Whether through direct mortality or other mechanisms, harvest-open conditions are associated with substantially lower detection of downstream fallback. Managers should consider this when evaluating harvest effects on overshoot steelhead.

2. **Spill operations may not be effective for facilitating fallback**: The lack of consistent spill effects suggests that modifying spill operations is unlikely to improve downstream return rates for overshoot fish.

3. **Year effects warrant further investigation**: The substantial unexplained year-to-year variation suggests that additional factors influence overshoot fate. Temperature, flow, and fish condition merit further study.

---

## 5. Conclusions

1. **Harvest exposure is the primary predictor of steelhead overshoot fate**, explaining the vast majority of predictable variation in downstream return probability. Harvest-exposed fish (hatchery fish in open harvest years) had 90% lower odds of returning downstream compared to other fish.

2. **Neither fall nor spring spill has a detectable effect on downstream return** after accounting for harvest exposure. Previous hypotheses about spill facilitating fallback are not supported by this analysis.

3. **Results are robust across different spill metrics**, with the harvest effect virtually identical regardless of which spill variable was included in the model.

4. **Substantial year-to-year variation remains unexplained**, suggesting that additional factors beyond spill and harvest influence overshoot fate. Future research should investigate temperature, flow, and fish condition.

5. **Bayesian hierarchical modeling provides appropriate inference** for this data structure, accounting for inter-annual variation while estimating population-level effects of management-relevant predictors.

---

## References

Bürkner, P.-C. 2017. brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software* 80(1):1–28.

Bürkner, P.-C. 2018. Advanced Bayesian multilevel modeling with the R package brms. *The R Journal* 10(1):395–411.

Gelman, A., A. Jakulin, M. G. Pittau, and Y.-S. Su. 2008. A weakly informative default prior distribution for logistic and other regression models. *Annals of Applied Statistics* 2(4):1360–1383.

Keefer, M. L., and C. C. Caudill. 2014. Homing and straying by anadromous salmonids: a review of mechanisms and rates. *Reviews in Fish Biology and Fisheries* 24:333–368.

Keefer, M. L., C. A. Peery, T. C. Bjornn, M. A. Jepson, and L. C. Stuehrenberg. 2008. Hydrosystem, dam, and reservoir passage rates of adult Chinook salmon and steelhead in the Columbia and Snake rivers. *Transactions of the American Fisheries Society* 137:1185–1207.

Muir, W. D., S. G. Smith, J. G. Williams, and B. P. Sandford. 2001. Survival of juvenile salmonids passing through bypass systems, turbines, and spillways with and without flow deflectors at Snake River dams. *North American Journal of Fisheries Management* 21:135–146.

Murdoch, A. R., C. O. Ochieng, and C. M. Moffett. 2022. Abundance and migration success of overshoot steelhead in the Upper Columbia River. *North American Journal of Fisheries Management* 42:1066–1079.

NOAA Fisheries. 2024. Upper Columbia River steelhead. https://www.fisheries.noaa.gov/west-coast/endangered-species-conservation/upper-columbia-river-steelhead

Ogle, D. H. 2024. Bayesian LVB I - brms. fishR. https://fishr-core-team.github.io/fishR/blog/posts/2024-2-5_LVB_brms/

R Core Team. 2023. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

Richins, S. M., and J. R. Skalski. 2018. Steelhead overshoot and fallback rates in the Columbia–Snake River basin and the influence of hatchery and hydrosystem operations. *North American Journal of Fisheries Management* 38:1122–1137.

Scheuerell, M. D., C. P. Ruff, J. H. Anderson, and E. M. Beamer. 2015. Analyzing large-scale conservation interventions with Bayesian hierarchical models: a case study of supplementing threatened Pacific salmon. *Ecology and Evolution* 5:2115–2125.

Thorson, J. T., and K. Kristensen. 2024. Hierarchical Bayesian models in fisheries stock assessment. *ICES Journal of Marine Science* 81:fsad199.

WDFW (Washington Department of Fish and Wildlife). 2024. Steelhead management in Columbia and Snake River basins. https://wdfw.wa.gov/fishing/management/columbia-river

---

## Appendix: Supplementary Information

### A.1 Software Versions

- R version 4.3.2 (2023-10-31)
- brms version 2.21.0
- Stan version 2.32.2
- tidyverse version 2.0.0
- bayesplot version 1.10.0
- posterior version 1.5.0

### A.2 Data Files

| File | Description |
|------|-------------|
| `all_overshoot_fish_2001_2025.csv` | Primary fish detection records |
| `OvershootRearType.csv` | Rear type classifications |
| `Harvest_Open_Closed.csv` | Annual harvest status |
| `merged_spill_2000_2024.csv` | Hourly spill data |
| `site_classifications.csv` | Detection site zone assignments |

### A.3 Output Files

| File | Description |
|------|-------------|
| `full_hierarchical_results_with_harvest_R.csv` | Parameter estimates for all 6 models |
| `full_year_summary_with_harvest_R.csv` | Year-level data summary |
| `analysis_fish_with_harvest_R.csv` | Final analysis dataset |
| `hierarchical_models_with_harvest_R.RData` | Saved R model objects |
| `full_data_overview_with_harvest_R.png` | Data summary figure |
| `full_spill_effects_with_harvest_R.png` | Spill effect posteriors |
| `full_harvest_effects_R.png` | Harvest effect posteriors |
| `full_year_effects_with_harvest_R.png` | Year random effects |
| `full_trace_plots_R.png` | MCMC diagnostics |

---

*Analysis conducted February 2026; updated March 2026*

*Statistical software: R 4.3.2, brms 2.21 (Stan backend)*

*Data sources: PTAGIS, Columbia River DART, WDFW*
