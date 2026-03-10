# Multinomial Analysis of Wild Steelhead Overshoot Fate at Wells Dam: Effects of Spill and Harvest Season

## Bayesian Hierarchical Multinomial Logistic Regression Using R/brms: 2007–2023

---

## Abstract

Wild steelhead (*Oncorhynchus mykiss*) that overshoot their natal tributaries and are detected at Wells Dam on the Upper Columbia River face three potential fates: fallback downstream toward natal spawning grounds, continued upstream movement into non-natal basins (Okanogan or Methow), or remaining at Wells Dam without subsequent detection. This study employed Bayesian hierarchical multinomial logistic regression implemented in R using the brms package to evaluate whether dam spill conditions and year-level harvest season status influence the fate distribution of 208 wild overshoot steelhead observed across 17 years (2007–2023). The multinomial model estimated separate log-odds contrasts for each outcome (downstream fallback and above-Wells migration) relative to the reference category (remaining at Wells Dam). Neither fall spill nor spring spill showed a detectable effect on either fate contrast. Similarly, whether harvest was open or closed in a given year did not significantly predict the fate of wild fish—a result consistent with the selective nature of mark-selective fisheries, which target adipose-clipped hatchery fish rather than wild fish. The majority of wild overshooting steelhead (68.3%) were last detected downstream of Wells Dam, indicating high fallback success for this group. Year-level random effects captured residual inter-annual variation. These findings contrast with a parallel analysis of hatchery fish, where harvest exposure was the dominant predictor of fate, and reinforce the biological distinction between wild and hatchery overshoot fish.

---

## 1. Introduction

### 1.1 Steelhead Overshoot Phenomenon

Anadromous salmonids exhibit strong natal homing fidelity, returning to their birth streams after years at sea. However, a proportion of fish "overshoot" their natal tributary confluence during upstream migration, ascending beyond their spawning destination before attempting to fall back (Keefer and Caudill 2014). Overshoot has become pervasive in the Columbia River Basin, with some steelhead populations exhibiting overshoot rates exceeding 50% (Richins and Skalski 2018; Murdoch et al. 2022). Critically, not all overshoot fish successfully fall back to their natal stream—a portion fail to return downstream, representing a meaningful source of mortality for threatened populations.

Upper Columbia River steelhead are listed as threatened under the Endangered Species Act. Wells Dam (river kilometer 830) sits upstream of the Wenatchee and Entiat River confluences, meaning any steelhead detected there has overshot its spawning destination. Murdoch et al. (2022) estimated that Wells Dam had the second-highest overshoot rate in the Upper Columbia, with approximately 20% of wild steelhead passing the dam being overshoot fish. Among those overshoot fish, a proportion failed to be detected returning downstream.

### 1.2 Multinomial Fate Structure

Unlike analyses that treat fallback as a binary outcome (downstream vs. remaining at dam), overshoot steelhead at Wells Dam exhibit a three-way fate structure:

1. **Downstream fallback (Below)**: Last detected below Wells Dam, indicating successful fallback toward natal tributaries—the "success" outcome from a recovery perspective.
2. **Above-Wells migration (Above)**: Last detected upstream of Wells Dam in non-natal basins (Okanogan or Methow rivers). These fish continued upstream migration rather than falling back, potentially entering non-natal spawning areas.
3. **At Wells Dam (Wells Dam)**: Last detected at Wells Dam without subsequent observation downstream or further upstream. These fish may have died at or near the dam, emigrated undetected, or remained in the tailwater.

A multinomial model is necessary to simultaneously characterize all three fate probabilities and properly account for the compositional nature of the outcome.

### 1.3 Wild Fish–Specific Considerations

This analysis focuses exclusively on wild steelhead (Rear Type Code = "W"), distinct from a parallel analysis of hatchery fish. This distinction is ecologically and statistically important for several reasons:

**Harvest exposure**: Upper Columbia mark-selective fisheries target adipose fin-clipped hatchery fish. Wild fish (unmarked) must be released under most circumstances. Consequently, "harvest open" years are expected to have little or no direct harvest effect on wild steelhead fates, in contrast to hatchery fish for which open harvest years dramatically reduce downstream detection probability.

**Behavioral differences**: Wild steelhead exhibit stronger homing fidelity and lower baseline overshoot rates than hatchery-origin fish (Richins and Skalski 2018). Their navigation behavior and post-overshoot fate may differ systematically from hatchery fish.

**Recovery implications**: Wild fish survival and successful spawning are paramount to ESA recovery goals. Understanding the drivers of wild overshoot fate—independently of hatchery fish dynamics—is directly relevant to recovery planning.

### 1.4 Study Objectives

This study addressed two primary questions for wild overshoot steelhead at Wells Dam:

1. **Spill Hypothesis**: Do dam spill conditions during the fall arrival period (August–November) or the subsequent spring holding period (January–March) influence the multinomial fate distribution of wild overshoot steelhead?

2. **Harvest Season Hypothesis**: Does year-level harvest season status (open vs. closed) predict wild fish fate, or are wild fish fate distributions independent of whether harvest is occurring in that year?

---

## 2. Methods

### 2.1 Study Area and Species

The study focused on wild steelhead from the Wenatchee and Entiat River basins detected at Wells Dam (river kilometer 830) on the Upper Columbia River. Fish detected at Wells Dam have passed their natal tributary confluences and are classified as overshoot fish. Fate was determined from PIT-tag detection records at monitoring sites distributed above and below Wells Dam.

### 2.2 Data Sources

**Fish Detection Data**: PIT-tag detection records for wild steelhead (Rear Type Code = "W") detected at Wells Dam from 2001–2025 were obtained from the Columbia Basin PIT Tag Information System (PTAGIS). Rearing origin was determined from PTAGIS mark records. Site classifications assigned each detection location to one of three geographic zones relative to Wells Dam: Below (downstream), Above (upstream, non-natal basins), or Wells Dam.

**Spill Data**: Hourly spill discharge measurements (KCFS) at Wells Dam from 2000–2024 were obtained from Columbia River DART (Data Access in Real Time). Seasonal spill metrics were calculated for the fall arrival window (August 1–November 30) and the spring holding period (January 1–March 31).

**Harvest Status**: Annual harvest status was compiled from Washington Department of Fish and Wildlife regulations. Harvest was open from 2000–2015 and again in 2024–2025; emergency closures were in effect from 2016–2023. For wild fish, this variable captures year-level harvest season context rather than individual harvest exposure.

### 2.3 Fish Fate Classification

Fish were assigned to one of three fate categories based on their last PIT-tag detection location:

- **Below**: Last detected at a site downstream of Wells Dam (fallback toward natal tributaries)
- **Above**: Last detected at a site upstream of Wells Dam in non-natal basins (Okanogan or Methow)
- **Wells Dam**: Last detected at Wells Dam facilities without evidence of subsequent directional movement

Fish with unknown or unclassifiable final detection locations were excluded. The analysis included only fish with complete covariate data (spill metrics and harvest status).

### 2.4 Predictor Variables

**Spill Year Assignment**: Fish were assigned to a spill year using a June 1 cutoff: fish first detected at Wells Dam between June of year Y and May of year Y+1 were assigned to SpillYear Y+1. This pairs each fish's arrival period with the fall spill conditions present during overshoot and the subsequent spring spill conditions during holding.

**Seasonal Spill Metric**: Total cumulative spill volume (KCFS-hours; "TotalSpill") was calculated for each seasonal window:
- Fall spill: August 1–November 30
- Spring spill: January 1–March 31

Spill metrics were standardized (mean 0, SD 1) prior to model fitting to facilitate prior specification and coefficient comparison.

**Harvest Open**: A binary year-level covariate indicating whether directed harvest on hatchery steelhead was open (1) or closed (0) in a given spill year. For wild fish, this variable does not represent individual harvest exposure but captures the broader management context of each year.

### 2.5 Statistical Analysis

#### 2.5.1 Bayesian Hierarchical Multinomial Logistic Regression

We fitted Bayesian hierarchical multinomial logistic regression models using the **brms** package (Bürkner 2017, 2018) in R version 4.3, with Stan as the computational backend. The multinomial model estimates simultaneous log-odds contrasts for K−1 outcome categories relative to a reference category. Wells Dam was designated the reference category (representing the "failure" outcome), yielding two sets of contrasts:

$$
\begin{aligned}
\mathbf{y}_{ij} &\sim \text{Categorical}(\mathbf{p}_{ij}) \\
\log\left(\frac{p_{ij}^{[\text{Below}]}}{p_{ij}^{[\text{Wells}]}}\right) &= \alpha_j^{[\text{Below}]} + \beta_{\text{spill}}^{[\text{Below}]} \cdot \text{Spill}_j + \beta_{\text{harvest}}^{[\text{Below}]} \cdot \text{HarvestOpen}_j \\
\log\left(\frac{p_{ij}^{[\text{Above}]}}{p_{ij}^{[\text{Wells}]}}\right) &= \alpha_j^{[\text{Above}]} + \beta_{\text{spill}}^{[\text{Above}]} \cdot \text{Spill}_j + \beta_{\text{harvest}}^{[\text{Above}]} \cdot \text{HarvestOpen}_j \\
\alpha_j^{[k]} &\sim \text{Normal}(\mu_\alpha^{[k]}, \sigma_\alpha^{[k]})
\end{aligned}
$$

Where:
- $p_{ij}^{[k]}$ = probability of outcome category $k$ for fish $i$ in year $j$
- $\alpha_j^{[k]}$ = year-specific random intercept for outcome $k$ (partial pooling across years)
- $\text{Spill}_j$ = standardized seasonal spill metric for year $j$
- $\text{HarvestOpen}_j$ = binary harvest status for year $j$

#### 2.5.2 Prior Specifications

Weakly informative priors were specified:

- $\mu_\alpha^{[k]} \sim \text{Normal}(0, 2)$ — population mean intercept for each outcome contrast
- $\sigma_\alpha^{[k]} \sim \text{HalfNormal}(0, 1)$ — year-level standard deviation for each outcome contrast
- $\beta_{\text{spill}}^{[k]} \sim \text{Normal}(0, 1)$ — spill effect for each outcome contrast
- $\beta_{\text{harvest}}^{[k]} \sim \text{Normal}(0, 1)$ — harvest season effect for each outcome contrast

Prior sensitivity was assessed by refitting the preferred model under three alternative prior specifications: wider priors (σ = 2), narrower priors (σ = 0.5), and skeptical (tight) priors.

#### 2.5.3 MCMC Sampling and Convergence

Models were fitted using 4 chains with sufficient warmup and sampling iterations to yield large effective sample sizes. Convergence was assessed using the Gelman-Rubin diagnostic (R̂), with values ≤1.01 indicating convergence. Effective sample sizes (ESS) exceeding 1,000 (and ideally >10,000) for all parameters were targeted. Visual inspection of trace plots, density plots, rank plots, and MCMC pair plots was used to confirm mixing and absence of pathological sampling behavior.

#### 2.5.4 Models Fitted

Two primary models were fitted, one for each seasonal spill window:

1. **Fall TotalSpill Model**: log-odds of Below vs. Wells and Above vs. Wells modeled as functions of fall cumulative spill volume and harvest season status, with year random intercepts
2. **Spring TotalSpill Model**: Same structure using spring cumulative spill volume

#### 2.5.5 Model Diagnostics

Posterior predictive checks (PPC) were conducted to assess calibration: the proportion of years for which the observed outcome count fell within the 95% posterior predictive interval was calculated for each outcome category. Leave-one-out cross-validation (LOO-CV) was used to assess predictive performance, with Pareto-k diagnostic values flagging potentially influential observations.

#### 2.5.6 Inference

Posterior distributions were summarized using means, standard deviations, and 95% highest density intervals (HDI). The probability of direction, P(β > 0), was calculated for each parameter to assess evidence for positive effects. Coefficients near zero with HDIs spanning zero and P(β > 0) near 0.5 indicate absence of evidence for an effect.

### 2.6 Software and Reproducibility

All analyses were conducted in R version 4.3.2 (R Core Team 2023) using:

- **brms** version 2.21 (Bürkner 2017) for Bayesian multinomial model fitting
- **tidyverse** (Wickham et al. 2019) for data manipulation
- **bayesplot** (Gabry and Mahr 2022) for MCMC and PPC diagnostics
- **posterior** (Bürkner et al. 2022) for posterior summarization
- **tidybayes** (Kay 2023) for tidy posterior extraction
- **patchwork** (Pedersen 2024) for figure composition

Model objects and analysis data are saved in `wild_fish_multinomial_models.RData` for full reproducibility. The complete analysis script is `multinomial_wild_fish_analysis.R` and the diagnostic script is `model_diagnostics_and_sensitivity.R`.

---

## 3. Results

### 3.1 Sample Summary

The analysis included **208 wild steelhead** across **17 spill years (2007–2023)**. Fate distribution across the full study period was:

- **Below (downstream fallback)**: 142 fish (68.3%)
- **Above (non-natal upstream)**: 19 fish (9.1%)
- **Wells Dam (no fallback detected)**: 47 fish (22.6%)

Wild fish demonstrated substantially higher downstream return rates than the hatchery-dominated full dataset (68.3% vs. 28.8%), consistent with better homing behavior and absence of direct harvest.

Annual sample sizes ranged from 1 fish (2007, 2023) to 33 fish (2016). Harvest was open during 2007–2016 (163 fish, 78.4%) and closed during 2017–2023 (45 fish, 21.6%).

**Table 1. Wild steelhead fate counts by year, harvest season status, and seasonal spill volume**

| Year | N | Below | Above | Wells | % Below | Harvest Open | Fall TotalSpill (KCFS-hr) | Spring TotalSpill (KCFS-hr) |
|------|---|-------|-------|-------|---------|--------------|--------------------------|----------------------------|
| 2007 | 1 | 0 | 0 | 1 | 0.0% | Yes | 4,682 | 6,310 |
| 2008 | 3 | 2 | 0 | 1 | 66.7% | Yes | 4,960 | 153 |
| 2009 | 3 | 2 | 0 | 1 | 66.7% | Yes | 5,951 | 256 |
| 2010 | 31 | 23 | 4 | 4 | 74.2% | Yes | 3,866 | 1 |
| 2011 | 19 | 11 | 3 | 5 | 57.9% | Yes | 4,091 | 1,892 |
| 2012 | 19 | 12 | 1 | 6 | 63.2% | Yes | 6,876 | 2,877 |
| 2013 | 11 | 6 | 1 | 4 | 54.5% | Yes | 9,667 | 1,117 |
| 2014 | 16 | 15 | 0 | 1 | 93.8% | Yes | 4,207 | 2,903 |
| 2015 | 27 | 18 | 1 | 8 | 66.7% | Yes | 4,187 | 4,606 |
| 2016 | 33 | 23 | 4 | 6 | 69.7% | Yes | 4,098 | 165 |
| 2017 | 16 | 10 | 3 | 3 | 62.5% | No | 3,623 | 17,709 |
| 2018 | 4 | 3 | 0 | 1 | 75.0% | No | 4,904 | 14,037 |
| 2019 | 4 | 2 | 0 | 2 | 50.0% | No | 8,117 | 847 |
| 2020 | 8 | 6 | 0 | 2 | 75.0% | No | 5,134 | 3,547 |
| 2021 | 10 | 8 | 1 | 1 | 80.0% | No | 8,756 | 2,215 |
| 2022 | 2 | 1 | 0 | 1 | 50.0% | No | 3,912 | 7,001 |
| 2023 | 1 | 0 | 1 | 0 | 0.0% | No | 4,276 | 840 |
| **Total** | **208** | **142** | **19** | **47** | **68.3%** | — | — | — |

*Note: Harvest closed 2017–2023. % Below = proportion of fish last detected downstream.*

### 3.2 Model Convergence and Diagnostics

Both models achieved excellent convergence and sampling performance:

- **R̂ (Gelman-Rubin)**: All fixed-effect parameters in the range [0.9998, 1.0005], well below the ≤1.01 threshold
- **Effective Sample Size**: All fixed-effect ESS in the range [11,914, 14,842], far exceeding the 1,000 threshold and indicating highly efficient sampling
- **Divergent transitions**: 0 in both models, indicating no pathological geometry in the posterior
- **LOO Pareto-k diagnostics**: 0 of 208 observations with Pareto-k > 0.7, indicating no unduly influential data points and reliable LOO estimates

**Posterior predictive calibration** was excellent: both the downstream (Below) and tributary (Above) contrast models correctly bracketed the observed annual counts within 95% predictive intervals for all 17 years (17/17 for both outcome contrasts). Visual diagnostics including trace plots, density plots, rank bar plots, and pairwise parameter scatter plots confirmed adequate mixing and absence of multicollinearity or identification problems.

### 3.3 Spill Effects

Neither fall nor spring total spill showed a detectable effect on wild steelhead fate in either outcome contrast (Table 2). All spill coefficients had 95% HDIs that included zero, and the probability of direction was near 0.5 for most estimates, indicating near-complete uncertainty about the sign of any effect.

**Fall TotalSpill model:**

- **Downstream vs. Wells**: β_spill = −0.17 (95% HDI: −0.54, 0.20), P(β > 0) = 0.17
- **Above vs. Wells**: β_spill = −0.38 (95% HDI: −1.02, 0.24), P(β > 0) = 0.11

**Spring TotalSpill model:**

- **Downstream vs. Wells**: β_spill = −0.04 (95% HDI: −0.49, 0.42), P(β > 0) = 0.42
- **Above vs. Wells**: β_spill = −0.01 (95% HDI: −0.68, 0.65), P(β > 0) = 0.49

The marginally negative point estimates for fall spill on both outcomes are small in magnitude and have credible intervals that span zero by a wide margin. There is no consistent evidence that higher fall or spring spill volumes shift wild steelhead fate probabilities in either direction relative to remaining at Wells Dam.

### 3.4 Harvest Season Effects

Year-level harvest season status (open vs. closed) did not significantly predict wild steelhead fate in either outcome contrast (Table 2):

**Fall TotalSpill model:**

- **Downstream vs. Wells**: β_harvest = −0.04 (95% HDI: −0.87, 0.79), P(β > 0) = 0.47
- **Above vs. Wells**: β_harvest = −0.25 (95% HDI: −1.46, 0.85), P(β > 0) = 0.33

**Spring TotalSpill model:**

- **Downstream vs. Wells**: β_harvest = −0.05 (95% HDI: −1.05, 0.90), P(β > 0) = 0.46
- **Above vs. Wells**: β_harvest = −0.15 (95% HDI: −1.47, 1.21), P(β > 0) = 0.41

All harvest season coefficients had wide HDIs spanning zero and probability-of-direction values near 0.5, indicating negligible information about whether harvest season status influences wild fish fate. This result is consistent with biological expectations: wild steelhead are protected under mark-selective regulations and are rarely harvested in directed fisheries.

**Table 2. Posterior summaries for multinomial model parameters (wild fish)**

| Window | Outcome Contrast | β_Spill (Mean) | Spill 95% HDI | P(β_Spill > 0) | β_Harvest (Mean) | Harvest 95% HDI | P(β_Harvest > 0) |
|--------|-----------------|----------------|---------------|----------------|-----------------|-----------------|-----------------|
| Fall | Downstream vs. Wells | −0.17 | [−0.54, 0.20] | 0.17 | −0.04 | [−0.87, 0.79] | 0.47 |
| Fall | Above vs. Wells | −0.38 | [−1.02, 0.24] | 0.11 | −0.25 | [−1.46, 0.85] | 0.33 |
| Spring | Downstream vs. Wells | −0.04 | [−0.49, 0.42] | 0.42 | −0.05 | [−1.05, 0.90] | 0.46 |
| Spring | Above vs. Wells | −0.01 | [−0.68, 0.65] | 0.49 | −0.15 | [−1.47, 1.21] | 0.41 |

*Note: β coefficients on log-odds scale. Spill standardized (mean 0, SD 1). Reference category: Wells Dam. Harvest = year-level binary indicator of open harvest season.*

### 3.5 Prior Sensitivity Analysis

Prior sensitivity was assessed by refitting the fall TotalSpill model under four prior specifications (Table 3). Results were consistent across all specifications:

- **Downstream contrast, β_spill**: Range of means from −0.09 to −0.18 across all priors; all HDIs included zero
- **Downstream contrast, β_harvest**: Range of means from −0.007 to −0.076; all HDIs included zero
- **Above contrast, β_spill**: Range of means from −0.15 to −0.43; all HDIs included zero
- **Above contrast, β_harvest**: Range of means from −0.04 to −0.33; all HDIs included zero

The maximum range across priors was 0.299 for the tributary harvest effect—modest variation that does not alter qualitative conclusions. Results are robust to reasonable choices of prior width.

**Table 3. Prior sensitivity analysis: downstream contrast β estimates**

| Prior Specification | β_Spill (Downstream) | β_Spill 95% HDI | β_Harvest (Downstream) | β_Harvest 95% HDI |
|--------------------|---------------------|-----------------|----------------------|------------------|
| Default (σ = 1) | −0.168 | [−0.525, 0.188] | −0.038 | [−0.818, 0.789] |
| Wider (σ = 2) | −0.178 | [−0.541, 0.218] | −0.076 | [−0.991, 0.856] |
| Narrower (σ = 0.5) | −0.136 | [−0.466, 0.191] | −0.007 | [−0.626, 0.624] |
| Skeptical (tight) | −0.091 | [−0.370, 0.191] | 0.006 | [−0.412, 0.423] |

### 3.6 Year Random Effects

Year-level random intercepts accounted for residual inter-annual variation in fate probabilities not attributable to the fixed-effect predictors. The year random effects confirmed that fate distributions varied across years in ways not captured by spill or harvest season. Sample size limitations in many years (N = 1–4 fish in some years) necessitated the hierarchical partial-pooling structure to obtain stable estimates.

---

## 4. Discussion

### 4.1 High Downstream Return Rate for Wild Fish

The most notable descriptive finding is that 68.3% of wild overshoot steelhead at Wells Dam were last detected downstream—a substantially higher rate than the 28.8% observed in the full (predominantly hatchery) dataset analyzed in the companion hatchery report. This difference is consistent with two non-exclusive explanations: (1) wild fish are not subject to direct harvest mortality during open seasons, and (2) wild fish may exhibit stronger homing behavior and better post-overshoot navigation than hatchery fish.

Richins and Skalski (2018) documented that wild and integrated-stock fish had 15% lower overshoot rates than segregated hatchery stocks, and outplanted juvenile-origin fish had 65% higher overshoot rates than wild fish. If similar behavioral differences extend to post-overshoot navigation, we would expect wild fish to return downstream at higher rates—consistent with what we observe.

### 4.2 Absence of Harvest Season Effect on Wild Fish

The finding that year-level harvest season status does not predict wild steelhead fate is both unsurprising and informative. Wild steelhead (unmarked fish) are protected under mark-selective harvest regulations; directed retention of wild fish is prohibited or extremely limited. Therefore, whether harvest is "open" or "closed" in a given year has little direct bearing on individual wild fish mortality from fishing.

This result stands in sharp contrast to the companion hatchery analysis, where harvest-exposed hatchery fish had 90% lower odds of returning downstream (β_harvest = −2.28, OR ≈ 0.10, P(β < 0) = 1.000). The divergence between wild and hatchery fish responses to harvest season status provides strong evidence that the harvest effect observed for hatchery fish reflects actual harvest-related processes (direct mortality or detection changes due to harvest), rather than confounding with other year-level factors that co-vary with harvest season.

### 4.3 Absence of Spill Effects

As in the companion hatchery analysis, neither fall nor spring spill predicted wild fish fate. The fall spill coefficients showed slightly negative point estimates for both outcome contrasts (β ≈ −0.17 for downstream, β ≈ −0.38 for Above), but these estimates were imprecise and widely credible intervals spanned zero.

The biological mechanism by which spill might facilitate fallback—altered hydraulic conditions, cue dispersion, or fishway entrance attractiveness—does not appear to operate strongly enough to be detectable in these data. Alternatively, if spill does influence fallback, the effect may be too small relative to inter-annual variability in other unmeasured factors to detect with the available sample size of 208 fish distributed across 17 years.

It is important to note that the power to detect spill effects is limited by the relatively small wild fish sample size and the cross-sectional, year-level nature of the spill predictor. Spill effects on other aspects of overshoot behavior (e.g., the probability of overshoot occurrence, or the timing of fallback) are not evaluated here.

### 4.4 Above-Wells Migration: A Distinct Fate

The 9.1% of wild fish last detected in above-Wells basins (Okanogan or Methow) represent a distinct biological outcome not captured in binary analyses. These fish continued migrating upstream past Wells Dam into non-natal basins rather than falling back. Whether these fish eventually strayed to spawn in non-natal tributaries, perished without returning, or were simply lost from the PIT-tag detection network is unknown from these data.

The Above contrast showed slightly larger (though still non-significant) negative point estimates for both spill and harvest season effects than the Downstream contrast. This could suggest that higher spill conditions might weakly reduce above-Wells migration relative to remaining at Wells, but the evidence is too weak to draw firm conclusions.

### 4.5 Comparison with Full-Dataset Analysis

The multinomial wild-fish model differs from the binary hatchery-focused hierarchical model in several structural ways:

| Feature | Wild Fish Multinomial | Full Dataset Binary |
|---------|----------------------|---------------------|
| Fish origin | Wild only | All (predominantly hatchery) |
| Outcome | 3-category multinomial | Binary (Downstream vs. Wells) |
| Key predictor | Year-level harvest season | Individual harvest exposure |
| N fish | 208 | 1,321 |
| N years | 17 | 20 |
| Harvest effect | None detected | Very large (OR ≈ 0.10) |
| Downstream rate | 68.3% | 28.8% |

The complementary analyses underscore the importance of analyzing wild and hatchery fish separately, and of using appropriate outcome models (multinomial vs. binary) to capture the full fate distribution.

### 4.6 Limitations

Several limitations should be considered:

1. **Small sample size**: 208 wild fish across 17 years yields limited statistical power, particularly for year-level predictors like spill and harvest season, which have only 17 independent observations at the year level.

2. **Detection-based inference**: Fates are inferred from last detection locations, not from direct observation. Fish last detected at Wells Dam may have died, emigrated undetected, or moved through unmonitored sites.

3. **Year-level confounding**: Harvest season status and year co-vary by definition; the harvest closure period (2016–2023) is a contiguous block. Separating harvest season effects from other temporal trends is inherently difficult.

4. **Spill metric limitations**: Cumulative spill volume may not capture the specific hydraulic conditions experienced by individual fish during their residence at or passage through the dam.

5. **Incomplete above-Wells fate resolution**: Fish last detected in non-natal basins (Okanogan, Methow) have indeterminate final outcomes—they may have perished, strayed to spawn, or subsequently returned without detection.

### 4.7 Management Implications

The findings for wild fish have several management-relevant implications:

1. **Wild fish fallback appears robust to harvest season context**: The absence of a harvest season effect on wild fish fate suggests that current mark-selective regulations are effectively protecting wild steelhead from direct harvest-related fate consequences.

2. **Spill management is unlikely to improve wild overshoot fallback rates**: Consistent with the hatchery analysis, spill operations do not appear to influence downstream return probability for wild fish.

3. **Wild overshoot fish are a relatively resilient sub-group**: The 68.3% downstream return rate suggests that most wild fish that overshoot Wells Dam eventually return. Management focus should remain on reducing overshoot occurrence in the first place (e.g., addressing thermal conditions) and ensuring that the 22.6% that remain at Wells Dam without subsequent detection receive continued attention.

4. **Monitoring of above-Wells wild fish**: The 9.1% of wild fish last detected in non-natal above-Wells basins warrants monitoring. If these fish stray to spawn in non-natal watersheds, it could have genetic and demographic consequences for Upper Columbia steelhead recovery.

---

## 5. Conclusions

1. **Wild overshoot steelhead at Wells Dam fall back downstream at substantially higher rates (68.3%) than hatchery-dominated samples (28.8%)**, consistent with better homing behavior and absence of direct harvest effects.

2. **Neither fall nor spring spill conditions had a detectable effect on wild steelhead fate**, for either the downstream fallback or above-Wells migration contrasts.

3. **Year-level harvest season status did not predict wild steelhead fate**, confirming that mark-selective regulations protect wild fish from harvest-season-associated fate changes. This result provides indirect support for the interpretation that the large harvest effect seen in hatchery fish reflects actual harvest mortality.

4. **Model convergence and diagnostics were excellent**: R̂ ≤ 1.001, ESS > 11,000, zero divergent transitions, and perfect posterior predictive calibration (17/17 years) for both outcome contrasts.

5. **Results were robust to prior specification**, with maximum parameter ranges of 0.08–0.30 across four prior sensitivity specifications, none altering qualitative conclusions.

6. **The multinomial modeling framework** appropriately accounts for the three-way fate structure of wild overshoot fish, providing a more complete characterization of outcomes than binary approaches and enabling simultaneous inference on both downstream fallback and above-Wells migration.

---

## References

Bürkner, P.-C. 2017. brms: An R package for Bayesian multilevel models using Stan. *Journal of Statistical Software* 80(1):1–28.

Bürkner, P.-C. 2018. Advanced Bayesian multilevel modeling with the R package brms. *The R Journal* 10(1):395–411.

Gelman, A., A. Jakulin, M. G. Pittau, and Y.-S. Su. 2008. A weakly informative default prior distribution for logistic and other regression models. *Annals of Applied Statistics* 2(4):1360–1383.

Kay, M. 2023. tidybayes: Tidy Data and Geoms for Bayesian Models. R package version 3.0.6. https://doi.org/10.5281/zenodo.1308151

Keefer, M. L., and C. C. Caudill. 2014. Homing and straying by anadromous salmonids: a review of mechanisms and rates. *Reviews in Fish Biology and Fisheries* 24:333–368.

Murdoch, A. R., C. O. Ochieng, and C. M. Moffett. 2022. Abundance and migration success of overshoot steelhead in the Upper Columbia River. *North American Journal of Fisheries Management* 42:1066–1079.

NOAA Fisheries. 2024. Upper Columbia River steelhead. https://www.fisheries.noaa.gov/west-coast/endangered-species-conservation/upper-columbia-river-steelhead

R Core Team. 2023. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.

Richins, S. M., and J. R. Skalski. 2018. Steelhead overshoot and fallback rates in the Columbia–Snake River basin and the influence of hatchery and hydrosystem operations. *North American Journal of Fisheries Management* 38:1122–1137.

Scheuerell, M. D., C. P. Ruff, J. H. Anderson, and E. M. Beamer. 2015. Analyzing large-scale conservation interventions with Bayesian hierarchical models: a case study of supplementing threatened Pacific salmon. *Ecology and Evolution* 5:2115–2125.

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
- tidybayes version 3.0.6

### A.2 Input Data Files

| File | Description |
|------|-------------|
| `all_overshoot_fish_2001_2025.csv` | Primary fish detection records |
| `OvershootRearType.csv` | Rear type classifications (W = Wild) |
| `Harvest_Open_Closed.csv` | Annual harvest status |
| `merged_spill_2000_2024.csv` | Hourly spill data |
| `site_classifications.csv` | Detection site zone assignments |

### A.3 Output Files (Wild/ folder)

| File | Description |
|------|-------------|
| `wild_fish_analysis_data.csv` | Final analysis dataset (208 wild fish) |
| `wild_fish_year_summary.csv` | Year-level fate counts and spill metrics |
| `wild_fish_multinomial_results.csv` | Posterior summaries for all model parameters |
| `wild_fish_diagnostics_summary.csv` | Convergence diagnostics and calibration results |
| `wild_fish_prior_sensitivity_results.csv` | Prior sensitivity parameter estimates |
| `wild_fish_multinomial_models.RData` | Saved brms model objects |
| `wild_fish_sensitivity_models.RData` | Prior sensitivity brms model objects |
| `wild_fish_data_overview.png` | Data summary figure |
| `wild_fish_model_results.png` | Posterior coefficient plot |
| `wild_fish_trace_plots.png` | MCMC trace plots |
| `wild_fish_ppc_diagnostics.png` | Posterior predictive checks |
| `wild_fish_prior_sensitivity.png` | Prior sensitivity comparison |
| `wild_fish_mcmc_trace.png` | Extended MCMC trace diagnostics |
| `wild_fish_mcmc_density.png` | MCMC density plots |
| `wild_fish_mcmc_ranks.png` | MCMC rank bar plots |
| `wild_fish_mcmc_pairs.png` | MCMC pairwise scatter plots |

---

*Analysis conducted February–March 2026*

*Statistical software: R 4.3.2, brms 2.21 (Stan backend)*

*Data sources: PTAGIS, Columbia River DART, WDFW*
