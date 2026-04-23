# Entiat Wild Steelhead Overshoot Analysis — Results Summary

**Date:** April 2026  
**Analyst:** Chas Kyger

---

## Overview

This analysis tests whether wild Entiat steelhead that interact with Wells Dam (overshoot past Rocky Reach Dam and are subsequently detected at Wells Dam) have a reduced probability of returning to the Entiat River subbasin to spawn. Data span return years 2006–2025 using PIT tag detections from PTAGIS.

---

## Study Design

**Anchor event:** First Rocky Reach adult fishway (RRF) detection per fish in July–December (fall) or January–May (spring). Two fish with confirmed adult returns that were missed by RRF detection (detected at Priest Rapids and Rock Island adult ladders) were anchored to their Priest Rapids detection date.

**Group assignment:**

- **Group A** (n = 223): No Wells Dam detection following RRF anchor
- **Group B** (n = 181): Wells Dam detection following RRF anchor; last detection not in Methow or Okanogan subbasin

Fish whose last post-anchor detection was in the Methow or Okanogan subbasin were excluded as strays.

**Return seasons:** 395 fall (July–December RRF), 9 spring (January–May RRF). All Group B fish are fall returners; all 9 spring fish are Group A.

**Return years:** 2006–2025 (20 years)

**Harvest status:** Tributary harvest was open 2006–2015 and 2024–2025; closed 2016–2023.

---

## Outcome Definitions

### Original outcome: EntiatDetected

Any detection at any Entiat River subbasin antenna following the RRF anchor, with no time restriction.

| Group | n | Detected | Rate |
|-------|---|----------|------|
| A (no Wells) | 223 | 200 | 89.7% |
| B (Wells) | 181 | 123 | 68.0% |

### Strict outcome: EntiatSpawned

A more conservative definition requiring evidence of reaching spawning habitat within 18 months of the RRF anchor:

- **Fall returners:** Detection at ENL (Lower Entiat River, rkm 1–2) AND at least one upstream Entiat site, both within 18 months of anchor
- **Spring returners:** Detection at ENL alone within 18 months of anchor

ENL-only fall fish are excluded because ENL is at the mouth of the Entiat with no spawning habitat below; brief or incomplete ENL detections in fall fish likely represent fish that entered but did not reach spawning grounds.

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A (no Wells) | 223 | 134 | 60.1% |
| B (Wells) | 181 | 86 | 47.5% |

---

## Data Quality

### Tag reuse

Release-date checks (age at anchor = days from release to RRF / 365.25) identified three fish with anomalous detection histories. All three were reviewed against their complete PTAGIS detection records and retained with corrections:

- **3D9.1C2C451B01** — Valid 2009 fall return (age at anchor 1.35 yr). A single 2017 ENS detection was confirmed as a data error and excluded from outcome scoring; EntiatDetected = 0.
- **3D9.1C2CDD70AD** — Tagged as a juvenile smolt (age at anchor 0.27 yr from tagging to first RRF detection). Confirmed 2011 adult return via BON → MCN → PRA → RIA → ENL 2011-09-26; upstream Entiat detections spring 2012. Re-anchored to PRA 2011-09-17; EntiatSpawned = 1.
- **3D9.1C2D45A06D** — Tagged at ENL as juvenile; RRF detection in 2010 was juvenile upstream passage. Confirmed 2012 adult return via BON → MCN → PRA → RIA; MAD/TLT/ENA detections spring 2013; kelt detected at RRJ April 2013. Re-anchored to PRA 2012-08-29; EntiatSpawned = 1.

---

## Statistical Models

All models are Bayesian hierarchical logistic regressions fit with brms. All include a random intercept for return year (20 levels). MCMC: 4 chains × 3,000 post-warmup iterations (iter = 4,500, warmup = 1,500), adapt_delta = 0.95. Priors: Normal(0, 1.5) on all fixed effects; Exponential(1) on year SD.

The detection-corrected models use a custom compound Bernoulli likelihood that accounts for imperfect ENL detection efficiency. For each fish, the probability of observing a spawning detection is modeled as π_obs = π_true × p_ENL, where p_ENL is the posterior mean predicted ENL detection efficiency from a separate Bayesian flow model (quadratic log-discharge + year random effect). Detection efficiency values were unavailable for 9 fish in the strict outcome models; those fish are excluded from detection-corrected model fits (n = 395).

---

## Results

### Original Outcome (EntiatDetected, n = 404)

#### Model 1: Wells interaction + harvest + year RE

| Parameter | β | OR | 95% CI | P(β < 0) |
|-----------|---|----|--------|----------|
| Wells (Group B) | −1.348 | 0.260 | [−1.881, −0.834] | 1.000 |

**Model 1b** (wells + year RE, no harvest covariate): β = −1.354, OR = 0.258, CI [−1.886, −0.829], P(β < 0) = 1.000. Results are nearly identical, confirming the Wells effect is not confounded by harvest timing.

#### Model 1c: Detection-corrected

| Parameter | β | OR | 95% CI | P(β < 0) |
|-----------|---|----|--------|----------|
| Wells (Group B) | −1.836 | 0.160 | [−2.737, −1.067] | 1.000 |

Correcting for imperfect ENL detection efficiency strengthens the estimated effect, as expected.

#### Model 1e: Season × Wells interaction (standard)

| Season | β | OR | 95% CI | P(β < 0) |
|--------|---|----|--------|----------|
| Fall | −1.402 | 0.246 | [−1.952, −0.885] | 1.000 |
| Spring | −1.407 | 0.245 | [−3.448, +0.651] | 0.907 |
| Interaction (spring − fall) | −0.006 | — | [−1.979, +1.975] | — |

The season interaction is negligible and the CI spans zero completely, though the spring estimate is imprecise (n = 9 spring fish, all Group A; spring Wells effect is unidentifiable from data).

#### Model 1f: Season × Wells interaction (detection-corrected)

| Season | β | OR | 95% CI | P(β < 0) |
|--------|---|----|--------|----------|
| Fall | −1.831 | 0.160 | [−2.708, −1.076] | 1.000 |
| Spring | −1.835 | 0.160 | [−3.990, +0.278] | 0.956 |

---

### Model 2: Spill and Harvest Effects on Group B (n = 160)

Tests whether Wells Dam spill conditions during the overshoot period predict subsequent Entiat return probability among Group B fish. Spring spill = January–March spill hours at Wells; fall spill = August–November spill hours at Wells (prior year, assigned to the following return year).

| Parameter | β | OR | 95% CI | P(β > 0) |
|-----------|---|----|--------|----------|
| Spring spill hours (std) | +0.063 | 1.066 | [−0.401, +0.554] | 0.601 |
| Fall spill hours (std) | −0.045 | 0.956 | [−0.505, +0.421] | 0.413 |
| Harvest open | +0.122 | 1.130 | [−0.817, +1.053] | — |

No credible effect of spring or fall spill hours on Entiat return probability among Group B fish. Posteriors span zero broadly for both spill predictors.

---

### Strict Outcome (EntiatSpawned, n = 404)

#### Model S1: Standard (wells + harvest + year RE)

| Parameter | β | OR | 95% CI | P(β < 0) |
|-----------|---|----|--------|----------|
| Wells (Group B) | −0.593 | 0.553 | [−1.016, −0.182] | 0.998 |

#### Model S2: Detection-corrected (n = 395)

| Parameter | β | OR | 95% CI | P(β < 0) |
|-----------|---|----|--------|----------|
| Wells (Group B) | −0.723 | 0.485 | [−1.251, −0.212] | 0.997 |

#### Model S3: Season × Wells interaction

| Season | β | OR | 95% CI | P(β < 0) |
|--------|---|----|--------|----------|
| Fall | −0.599 | 0.550 | [−1.015, −0.187] | 0.997 |
| Spring | −0.601 | 0.548 | [−3.518, +2.316] | 0.661 |
| Interaction | −0.002 | — | [−2.885, +2.852] | — |

No season interaction on the strict outcome. The spring estimate is unidentifiable (no Group B spring fish).

---

## Summary Across Outcomes

| Model | Outcome | β (Wells) | OR | 95% CI | P(β < 0) |
|-------|---------|-----------|-----|--------|----------|
| M1 standard | ENL detected | −1.348 | 0.260 | [−1.881, −0.834] | 1.000 |
| M1c det-corrected | ENL detected | −1.836 | 0.160 | [−2.737, −1.067] | 1.000 |
| S1 standard | Strict spawned | −0.593 | 0.553 | [−1.016, −0.182] | 0.998 |
| S2 det-corrected | Strict spawned | −0.723 | 0.485 | [−1.251, −0.212] | 0.997 |

The Wells Dam interaction effect is highly credible under all model specifications. The effect is stronger on the original ENL-detection outcome (OR ≈ 0.16–0.26) than on the strict spawning outcome (OR ≈ 0.49–0.55). The true effect on spawning success likely falls between these bounds, depending on what fraction of ENL-only fall fish (detected at the river mouth but not upstream) successfully spawned in the unmonitored lower Entiat reach (rkm 1–~10, below the ENM antenna).

---

## Post-Spawn Kelt Detection

Among the 220 fish scored EntiatSpawned = 1 in the strict outcome (groups A and B combined), a substantial fraction were subsequently detected at Rocky Reach Dam juvenile bypass (RRJ/BCC) within weeks of their last Entiat detection. This is consistent with normal post-spawn kelt outmigration and provides independent confirmation of spawning.

---

## Notes on Detection Efficiency

ENL detection efficiency varies substantially with river discharge. The Bayesian flow model (quadratic log-discharge + year random effect) predicts mean efficiency across return years and typical fall flows of approximately 60–80%. Years with high discharge show lower efficiency, which is the primary motivation for the detection-corrected models. The detection-corrected results (M1c, S2) are considered the primary analysis.
