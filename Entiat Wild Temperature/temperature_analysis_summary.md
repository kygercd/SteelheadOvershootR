# Temperature Analysis Summary: Rocky Reach Forebay and Entiat Wild Steelhead

**Analysis:** Rocky Reach Forebay (RRH) water temperature and its relationship to (1) Wells Dam overshoot probability and (2) Entiat wild steelhead return probability (criterion C2).

**Data:** 20 return years (2006–2025), n = 404 fish. Temperature data from USACE CROHMS via DART WQM; gap years 2006–2011 imputed from calibrated surrogate stations (WEL, PRD, CHJ). Models are Bayesian logistic regressions (brms, Bernoulli logit) with random year intercepts and weakly informative priors. All continuous predictors standardized (mean = 0, SD = 1).

---

## Temperature Data Coverage and Trends

All 20 return years have complete Aug–Sep coverage. Aug–Nov metrics (≥85 days required) are available for 19 of 20 years; 2007 falls just above the threshold at 87 days with November data absent. Years 2006–2011 rely partially on surrogate stations calibrated against RRH over complete-coverage years (2012–2016; R² = 0.94–0.98).

Rocky Reach forebay temperatures show a warming trend over the study period. The 2025 return year was the warmest on record, with a mean Aug–Sep temperature of 20.2°C and 73 days above 18°C during Aug–Nov.

---

## Section 7 — Annual Temperature and Overshoot Probability

None of the five annual temperature metrics credibly predicted Wells Dam overshoot probability. All 95% credible intervals spanned zero.

| Metric | β | OR | 95% CI (β) | P(β > 0) |
|--------|---:|----|:----------:|:--------:|
| Annual max daily mean | 0.20 | 1.22 | −0.04, 0.43 | 0.95 |
| Mean Aug–Sep | 0.11 | 1.12 | −0.16, 0.37 | 0.81 |
| Mean Aug–Nov | 0.04 | 1.04 | −0.25, 0.33 | 0.62 |
| Cumulative degree-days Aug–Nov | 0.00 | 1.00 | −0.27, 0.27 | 0.52 |
| Days >18°C Aug–Nov | 0.03 | 1.03 | −0.26, 0.30 | 0.61 |

Annual maximum temperature had the strongest posterior lean (95% of draws positive) but the credible interval still included zero. Harvest status (open/closed) was not a credible predictor of overshoot in any model.

---

## Section 8 — Annual Temperature and Return Probability (C2)

No annual temperature metric credibly predicted return probability once wells group and harvest status were controlled. Across all five models the wells group effect was dominant and nearly identical.

| Metric | β (temp) | P(β < 0) | Wells group β | Wells group OR |
|--------|--------:|:--------:|-------------:|:--------------:|
| Annual max daily mean | −0.10 | 0.68 | −0.82 | 0.44 |
| Mean Aug–Sep | −0.06 | 0.59 | −0.83 | 0.44 |
| Mean Aug–Nov | −0.13 | 0.71 | −0.83 | 0.44 |
| Cumulative degree-days Aug–Nov | +0.15 | 0.21 | −0.83 | 0.44 |
| Days >18°C Aug–Nov | −0.25 | 0.88 | −0.83 | 0.44 |

The days >18°C Aug–Nov metric showed the strongest signal (88% posterior probability of a negative effect on return), but its credible interval spanned zero (β = −0.70 to 0.19). The wells group effect on return probability was robust across all temperature specifications (OR ≈ 0.44, equivalent to Group B fish having roughly 56% lower odds of meeting criterion C2 compared to Group A fish).

---

## Section 9 — Individual Arrival Temperature and Overshoot Probability

Individual fish arrival temperature at Rocky Reach was strongly and positively associated with Wells Dam overshoot probability. This is the most consistent finding of the analysis.

| Metric | β | OR | 95% CI (β) | P(β > 0) |
|--------|---:|----|:----------:|:--------:|
| Arrival temp (day of arrival) | 0.52 | 1.68 | 0.27, 0.79 | 1.00 |
| Arrival temp (±15-day mean) | 0.66 | 1.94 | 0.37, 0.99 | 1.00 |

Fish arriving at Rocky Reach during warmer thermal conditions were substantially more likely to continue upstream past Wells Dam. The 30-day window metric (±15 days around the anchor date) produced a slightly larger effect estimate than the single-day arrival temperature (OR = 1.94 vs. 1.68), suggesting the broader thermal context at arrival is at least as informative as point temperature.

**Arrival temperature did not predict return probability (C2).** After controlling for wells group, the effect of arrival temperature on return was negligible:

| Metric | β | OR | 95% CI (β) | P(β < 0) |
|--------|---:|----|:----------:|:--------:|
| Arrival temp (day of arrival) | −0.11 | 0.89 | −0.38, 0.13 | 0.81 |
| Arrival temp (±15-day mean) | −0.05 | 0.95 | −0.32, 0.21 | 0.63 |

The wells group coefficient remained stable at OR ≈ 0.45 in both models.

---

## Section 10 — Temperature × Wells Group Interaction on Return Probability

No temperature metric showed a credible interaction with wells group on return probability (C2 or C3). All interaction terms were small and 95% credible intervals spanned zero.

| Interaction term | β | OR | 95% CI (β) | P(β < 0) |
|-----------------|---:|----|:----------:|:--------:|
| Mean Aug–Sep × wells group | −0.05 | 0.96 | −0.49, 0.40 | 0.58 |
| Mean Aug–Nov × wells group | −0.14 | 0.87 | −0.55, 0.28 | 0.74 |
| Arrival temp × wells group (C2) | −0.14 | 0.87 | −0.66, 0.38 | 0.70 |
| Arrival temp × wells group (C3) | −0.28 | 0.75 | −0.96, 0.35 | 0.81 |

Interaction terms were consistently negative in direction — suggesting that warmer temperatures may reduce return probability somewhat more for Group B (overshooting) fish — but none cleared a credibility threshold. The wells group main effect remained negative and credible across all interaction models.

---

## Summary

| Question | Finding |
|----------|---------|
| Does annual forebay temperature predict overshoot? | No — no credible effect across all five annual metrics |
| Does annual forebay temperature predict return probability? | No — days >18°C Aug–Nov shows the strongest signal (88% posterior) but CI spans zero |
| Does individual arrival temperature predict overshoot? | **Yes — strongly.** OR = 1.68–1.94; essentially all posterior mass positive |
| Does individual arrival temperature predict return probability? | No — effect disappears after controlling for wells group |
| Does temperature moderate the wells group effect on return? | No — interaction terms are consistently near zero |
| Is the wells group effect robust to temperature? | **Yes** — OR ≈ 0.44 across all 18 models |

The principal finding is that warm arrival temperatures at Rocky Reach Forebay trigger overshoot of Wells Dam, but temperature does not independently determine whether an overshooting fish ultimately returns to the Entiat. The overshoot event itself — not the thermal conditions that precipitate it — is the primary driver of reduced return probability.
