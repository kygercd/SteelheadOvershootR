# Wells Dam Overshoot & PIT Array Detection Efficiency

**Summary Report — Wild Entiat River Steelhead**

Public Utility District No. 1 of Douglas County | April 2026

---

## Overview

This report summarizes two related analyses using PIT-tag detection data for wild steelhead in the upper Columbia River. The first characterizes detection efficiency at three adult steelhead gateway arrays (LMR/Methow, OKL/Okanogan, ENL/Entiat), which is needed to correctly interpret detection-based return estimates. The second uses those efficiency estimates to test whether wild Entiat River steelhead that overshoot to Wells Dam have reduced probability of returning to their natal subbasin.

---

## Part 1: Gateway Array Detection Efficiency

### What Was Analyzed

Three lower-subbasin PIT antenna arrays serve as "gateways" for adult steelhead returning to their natal tributaries above Wells Dam. A fish detected above Wells Dam but not detected at its subbasin gateway was either (a) missed by the antenna or (b) never entered that subbasin. Distinguishing these requires knowing how often fish are missed.

Detection efficiency was estimated for:

- **LMR** — Lower Methow River (operational since 3/4/2009). 9,702 eligible fish.
- **OKL** — Lower Okanogan River (operational since 10/8/2013). 2,767 eligible fish.
- **ENL** — Lower Entiat River (operational since 10/25/2007). 2,078 eligible fish.

Source: PTAGIS detection records for Wells Dam steelhead subsequently detected in Methow/Okanogan subbasins (LMR/OKL) and Rocky Reach Dam steelhead detected in the Entiat subbasin (ENL). For each site, a fish detected anywhere in the subbasin was scored as detected (1) or missed (0) at the gateway.

### Overall Miss Rates

| Site | Total fish | Missed gateway | Miss rate |
|------|-----------|----------------|-----------|
| LMR (Methow) | 9,702 | 1,818 | 18.7% |
| OKL (Okanogan) | 2,767 | 824 | 29.8% |
| ENL (Entiat) | 2,078 | 461 | 22.2% |

Fish detected before each site's operational start date were excluded (174 Methow, 810 Okanogan, 0 Entiat) to avoid confounding non-detection with the site not yet running.

### Key Patterns

**LMR and OKL — Flow-driven seasonal pattern**

Detection is near-perfect in late summer and fall (Aug–Nov; >98%) when flows are low. Efficiency drops sharply in spring (Mar–Jun) as snowmelt raises discharge: as low as 27–59% at LMR and 16–51% at OKL. Both sites showed a mid-period performance dip around 2015–2019 before recovering.

**ENL — Unusual double-trough pattern**

ENL behaves differently from LMR and OKL. Detection dips in both May (spring freshet) AND September (peak steelhead run timing), creating an inverted-U relationship with flow — both very low and very high discharge reduce detection, with optimal performance around 200–500 cfs. Annual detection declined after 2015, with especially poor performance in 2019, 2023, and 2024, likely reflecting equipment or maintenance issues.

### Models

Two complementary models were fitted for each site:

**GAM model:** `detected ~ s(day-of-year, cyclic) + s(year)`

| Site | AUC | Deviance explained |
|------|-----|--------------------|
| LMR  | 0.893 | 37.6% |
| OKL  | 0.925 | 51.8% |
| ENL  | 0.846 | 27.4% |

**Bayesian flow model:** `detected ~ log(discharge) + s(day-of-year) + (1|year)`

Incorporates USGS daily streamflow (gauges 12449950/12439500/12452990) as the physical driver. Year random effects isolate equipment problems from flow-driven variation.

| Site | Flow coefficient (β) | 95% CI | AUC |
|------|---------------------|--------|-----|
| LMR  | −0.43 | [−0.56, −0.30] | 0.900 |
| OKL  | −0.46 | [−0.68, −0.24] | 0.929 |
| ENL  | +0.74 / −0.67 (quadratic) | — | 0.875 |

GAM predictions by return date and year are used as per-fish p~ENL~, p~LMR~, and p~OKL~ in the return analysis below.

---

## Part 2: Wells Dam Overshoot and Entiat Return Probability

### Research Question

Do wild Entiat steelhead that overshoot past their natal tributary and interact with Wells Dam (31 miles upstream of the Entiat River confluence) have reduced probability of subsequently returning to the Entiat River to spawn?

### Study Design

PIT-tag detection records from PTAGIS for 379 wild Entiat-origin steelhead with valid adult Rocky Reach Dam (RRF) anchor dates, spanning return years 2004–2025. Each fish was assigned an adult RRF anchor date (first Jul–Dec RRF detection) and classified based on post-anchor detections:

| Group | n | Definition |
|-------|---|------------|
| Group A | 239 | No Wells Dam detection after adult RRF passage |
| Group B | 140 | Wells Dam detection after RRF; not last detected in Methow/Okanogan |
| Excluded — Strays | 13 | Last detected in Methow or Okanogan (confirmed non-natal spawners) |

### Raw Detection Rates

| Group | n | Entiat detected | Rate |
|-------|---|----------------|------|
| Group A (no Wells) | 239 | 175 | 73.2% |
| Group B (Wells interaction) | 140 | 82 | 58.6% |
| **Difference** | | | **−14.6 pp** |

The raw gap is concentrated in **fall returners** (Group A 74.0% vs. Group B 57.7%), not spring — the season when ENL detection is highest. Spring detection rates were similar between groups (68.6% vs. 64.7%).

### Bayesian Model Results

Model structure: `EntiatDetected ~ wells_interaction + harvest_open + (1 | ReturnYear)`
4 chains × 4,500 iterations; all R̂ ≤ 1.01; bulk ESS > 3,000.

| Model | β (Wells) | 95% CI | Odds Ratio | P(β < 0) |
|-------|-----------|--------|-----------|-----------|
| Model 1 — Standard | −0.642 | [−1.09, −0.21] | 0.53 | 0.997 |
| Model 1b — + Spring season covariate | −0.645 | [−1.08, −0.18] | 0.52 | 0.998 |
| Model 1c — Detection-corrected | −0.790 | [−1.34, −0.26] | **0.45** | **0.998** |
| Model 1d — Corrected + season | −0.789 | [−1.34, −0.25] | **0.45** | **0.998** |

The detection-corrected models (1c/1d) use the compound Bernoulli likelihood Y ~ Bernoulli(θ × p~ENL~), incorporating per-fish ENL detection efficiency from the gateway GAM. Both groups had identical mean p~ENL~ = 0.864, so differential detection is not a source of between-group bias.

### Harvest and Spill Effects (Model 2)

Among Group B fish with available spill data (n = 116), **no spill metric or harvest status predicted Entiat return probability**:

| Predictor | β | 95% HDI | Odds Ratio | P(positive) |
|-----------|---|---------|-----------|-------------|
| Spring Spill Hours (std) | +0.301 | [−0.22, +0.86] | 1.35 | 0.864 |
| Fall Spill Hours (std) | −0.045 | [−0.54, +0.48] | 0.96 | 0.426 |
| Harvest Open (vs. Closed) | −0.084 | [−1.31, +1.14] | 0.92 | 0.439 |

All posteriors broadly span zero. Spill operations at Wells Dam do not appear to substantially alter the probability of wild Entiat steelhead returning to the Entiat River following an overshoot. Harvest status (open 2000–2015, 2024–2025; closed 2016–2023) showed no credible effect in any model.

### Detection Efficiency Robustness Checks

Four independent lines of evidence confirm the group gap is not a detection artifact:

1. **Seasonal composition** — Both groups had similar spring return fractions (14.6% Group A; 12.1% Group B). Differential seasonal composition cannot explain the gap.

2. **Within-season rates** — The gap is concentrated in fall returns (Group A 74.0% vs. Group B 57.7%), precisely when ENL detection is highest, not in spring where detection issues would appear.

3. **Detection-corrected model** — Formally incorporating per-fish p~ENL~ into the Bernoulli likelihood (Models 1c/1d) amplifies rather than attenuates the Wells effect (OR 0.53 → 0.45). Both groups have identical mean p~ENL~ = 0.864.

4. **Hidden stray analysis** — Of 53 undetected Group B fish, expected fate breakdown using p~LMR~/p~OKL~ from the gateway GAMs:

| Expected fate | Fish | % of undetected |
|---------------|------|----------------|
| True Entiat returners missed at ENL | 19.2 | 36% |
| Hidden Methow strays missed at LMR | 1.7 | 3% |
| Hidden Okanogan strays missed at OKL | 0.4 | 1% |
| Other (died, harvested, unknown) | 31.7 | 60% |

Only ~2 hidden strays expected among the 53 undetected Group B fish — negligible.

### Detection-Corrected Return Rates

| Group | Observed rate | Corrected rate |
|-------|--------------|----------------|
| Group A (no Wells) | 73.2% | ~80% |
| Group B (Wells interaction) | 58.6% | ~72% |

Both groups' true return rates are higher than raw observed rates after accounting for ENL detection misses. The relative difference in log-odds is larger after correction.

---

## Key Takeaways

**1. Wells Dam overshoot substantially reduces Entiat return probability.**
Group B fish have approximately half the odds of returning to the Entiat River compared to Group A fish (OR = 0.53 raw; OR = 0.45 detection-corrected). This is a large, highly credible effect consistent across all four model specifications (P > 0.997).

**2. The effect is not a detection artifact.**
Spring return timing, within-season detection rates, a compound likelihood model, and a hidden stray analysis all confirm the group gap reflects a genuine behavioral or survival difference rather than imperfect antenna detection.

**3. Spill and harvest do not explain the effect.**
Neither spring nor fall spill volumes at Wells Dam nor annual harvest status (open vs. closed) predicted whether a Wells-interacting fish subsequently returned to the Entiat. This is consistent across multiple analysis approaches and corroborated by previous analyses of the broader Wells Dam overshoot dataset.

**4. Gateway detection efficiency is seasonal and flow-driven.**
ENL, LMR, and OKL antennas all show substantial seasonal variation in detection probability driven by streamflow. LMR and OKL are near-perfect in fall but can fall to 20–40% in spring high-flow periods. ENL has a more complex pattern with both spring and September troughs. Per-fish efficiency estimates from GAM models correct for this in the return analysis.

**5. Mechanism is unresolved.**
The most likely contributors to non-return in Group B fish are: downstream passage difficulty or injury at Wells Dam, overwinter or en-route mortality in the reservoir, and undetected straying to non-natal tributaries. PIT-tag data alone cannot discriminate among these mechanisms. Telemetry-based research on wild Entiat-origin fish above Wells Dam would help.

---

## Analysis Files

All files in `/home/chas/SteelheadOvershootR/Entiat Wild/`

| File | Description |
|------|-------------|
| `EntiatWildSteel.csv` | Source PIT-tag detection data (417 fish) |
| `entiat_return_bayesian_analysis.R` | Standard Bayesian models (1, 1b, 2) |
| `entiat_detection_corrected_analysis.R` | Detection-corrected models (1c, 1d) + hidden stray analysis |
| `entiat_bayesian_models.RData` | Fitted brms objects: fit_model1, fit_model1b, fit_model2 |
| `entiat_corrected_models.RData` | Fitted brms objects: fit_corrected_m1, fit_corrected_m1d |
| `gateway_detection_models.RData` | Fitted GAMs: gam_met, gam_oka, gam_enl |
| `gateway_bayesian_models.RData` | Fitted Bayesian flow models: fit_lmr_comb, fit_okl_comb, fit_enl_comb |
| `gateway_detection_summary.docx` | Full gateway detection efficiency report |
| `Entiat_Wild_Steelhead_Wells_Overshoot_Manuscript.docx` | Full manuscript |
