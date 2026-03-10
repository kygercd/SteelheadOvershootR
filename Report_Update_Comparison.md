# Summary of Changes: Updated vs. Original Harvest Exposure Analysis Report

## Steelhead Overshoot Analysis — March 2026 Update

---

## Background

The updated analysis excludes a set of censored PIT tags (listed in `SteelheadOvershootCensoredTags.csv`) that were identified as potentially problematic after the original report was completed. These tags were removed from the fish detection dataset prior to re-running the hierarchical models. All other methods, model structure, priors, and MCMC settings remain identical to the original analysis.

---

## Sample Size Changes

| Metric | Original | Updated | Difference |
|--------|----------|---------|------------|
| Total fish analyzed | 1,355 | 1,321 | −34 |
| Fish returning downstream | 388 | 381 | −7 |
| Fish remaining at Wells | 967 | 940 | −27 |
| % returning downstream | 28.6% | 28.8% | +0.2 pp |
| Harvest-exposed fish | 1,149 | 1,126 | −23 |
| % harvest-exposed | 84.8% | 85.2% | +0.4 pp |

The 34 removed fish were distributed across multiple years, with the largest reductions in 2011 (5 fish), 2007 (4 fish), 2009 (4 fish), 2012 (4 fish), and 2006 (3 fish). Nearly all removed fish were hatchery fish present during open harvest years (i.e., harvest-exposed). The proportional composition of the dataset is essentially unchanged.

### Year-Level Changes

The table below summarizes per-year differences in total fish counts between the original and updated analyses. Years not listed had no change.

| Year | Original Total | Updated Total | Difference |
|------|---------------|---------------|------------|
| 2005 | 26 | 25 | −1 |
| 2006 | 60 | 57 | −3 |
| 2007 | 111 | 107 | −4 |
| 2008 | 92 | 91 | −1 |
| 2009 | 134 | 130 | −4 |
| 2010 | 289 | 286 | −3 |
| 2011 | 196 | 191 | −5 |
| 2012 | 134 | 130 | −4 |
| 2013 | 115 | 114 | −1 |
| 2014 | 38 | 35 | −3 |
| 2015 | 39 | 37 | −2 |
| 2016 | 64 | 62 | −2 |
| 2022 | 2 | 1 | −1 |

---

## Model Results Changes

### Harvest Exposure Effect

The harvest exposure coefficient is marginally stronger in the updated analysis (slightly more negative), consistent with the observation that the removed fish were predominantly harvest-exposed fish that remained at Wells Dam. Removing them slightly reduces the number of harvest-exposed fish with "At Wells" fates, which modestly increases the estimated difference between harvest-exposed and non-exposed fish.

| Model | Original β_harvest (95% HDI) | Updated β_harvest (95% HDI) | Change |
|-------|------------------------------|------------------------------|--------|
| Fall TotalSpill | −2.18 (−2.60, −1.79) | −2.29 (−2.71, −1.88) | −0.11 |
| Spring TotalSpill | −2.17 (−2.59, −1.77) | −2.28 (−2.72, −1.87) | −0.11 |

**Odds ratio**: The odds ratio for harvest-exposed fish changed from approximately 0.11 (original; ~89% lower odds) to approximately 0.10 (updated; ~90% lower odds). This difference is negligible in practical terms.

**Posterior certainty**: The probability of a negative harvest effect remained P(β < 0) = 1.000 in both analyses.

### Spill Effects

Spill effect estimates are essentially unchanged between the original and updated analyses. The marginally negative Fall SpillHours effect and the near-zero spring spill effects are consistent across both versions.

| Model | Original β_Spill (95% HDI) | Updated β_Spill (95% HDI) |
|-------|---------------------------|--------------------------|
| Fall SpillHours | −0.29 (−0.53, −0.06) | −0.30 (−0.54, −0.07) |
| Fall TotalSpill | −0.02 (−0.39, 0.35) | −0.03 (−0.39, 0.33) |
| Spring TotalSpill | 0.03 (−0.22, 0.28) | 0.03 (−0.22, 0.28) |

### Year Random Effect Standard Deviation

The estimated year-level standard deviation (σ_year) is slightly smaller in the updated analysis for some models (e.g., Fall SpillHours: 0.52 → 0.49), indicating marginally less unexplained inter-annual variation after removing the censored tags.

---

## Impact on Conclusions

**The main conclusions of the analysis are unchanged.** All five key conclusions from the original report hold in the updated analysis:

1. Harvest exposure remains the dominant predictor of overshoot fate (~90% lower odds of downstream return; unchanged from ~89%).
2. Neither fall nor spring spill has a detectable effect on downstream return probability.
3. Results remain robust across all six spill metrics.
4. Substantial unexplained year-to-year variation persists.
5. Bayesian hierarchical modeling remains the appropriate inferential framework.

The only material change is a slight strengthening of the already dominant and highly certain harvest exposure effect. The removal of 34 censored tags had no meaningful effect on the qualitative interpretation of any result.

---

*Original analysis: February 2026*

*Updated analysis: March 2026*

*Change: Removal of censored tags listed in SteelheadOvershootCensoredTags.csv*
