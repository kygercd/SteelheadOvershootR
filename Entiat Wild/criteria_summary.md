# Entiat Wild Steelhead Return Analysis: Outcome Criteria Summary

## Overview

This document summarizes the four outcome criteria used to define successful return/spawning for wild Entiat River steelhead, along with the Bayesian model results for the effect of Wells Dam interaction on each outcome. All analyses use 404 fish (Group A: n=223, no Wells Dam detection post-anchor; Group B: n=181, Wells Dam detection post-anchor) with fall Rocky Reach fishway (RRF) detection as the anchor event. Models include a random year intercept and harvest status as a covariate. Each criterion is modeled under two likelihood specifications: standard Bernoulli and a detection-corrected compound likelihood that accounts for imperfect detection at the Entiat River mouth array (ENL).

---

## Study Groups

- **Group A** (n=223): Fish with no Wells Dam detection after the fall RRF anchor. These fish are assumed to have turned back before reaching Wells.
- **Group B** (n=181): Fish with a Wells Dam detection after the fall RRF anchor (excluding confirmed strays to the Methow or Okanogan subbasins). These fish overshot past Wells Dam before returning.

---

## Criterion 1 — Broad Detection

**Definition:** Any detection at any Entiat River subbasin site after the fall RRF anchor date, with no time window restriction.

**Rationale:** The broadest possible measure of Entiat return. A fish is considered to have returned if it was detected anywhere in the subbasin at any point after anchoring.

**Return rates:**

| Group | n | Returned | Rate |
|-------|---|----------|------|
| A | 223 | 200 | 89.7% |
| B | 181 | 123 | 68.0% |

**Model results — Wells Dam effect (β < 0 = reduced return probability):**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −1.348 | 0.26 | [−1.881, −0.834] | 1.000 |
| Detection-corrected | −1.836 | 0.16 | [−2.737, −1.067] | 1.000 |

Wells Dam interaction is strongly associated with reduced broad return probability. The detection-corrected model accounts for the possibility that some Group A fish returned but were not detected at ENL, which would bias the observed Group A rate downward; correcting for this makes the Wells effect appear larger.

---

## Criterion 2 — Original Strict (All-Fall Rule)

**Definition:** All fish are treated as fall returners regardless of their actual Entiat entry date. A fish is considered to have spawned if it was detected at ENL **and** at least one upstream Entiat site (ENM or higher) within 18 months of the fall RRF anchor.

**Rationale:** The strict criterion requires evidence of upstream passage, not just arrival at the river mouth. The 18-month window captures one complete spawning season after the anchor event. Because all fish have fall RRF anchors, applying the fall criterion uniformly (ENL + upstream required) is conservative and consistent. However, it misclassifies fish that actually entered the Entiat in spring — those fish typically spawn closer to the mouth and may not be detected at upstream fixed arrays, even when present.

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 136 | 61.0% |
| B | 181 | 91 | 50.3% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.517 | 0.60 | [−0.937, −0.095] | 0.992 |
| Detection-corrected | −0.671 | 0.51 | [−1.191, −0.148] | 0.995 |

The Wells effect is credible and robust. The detection-corrected estimate is modestly larger, reflecting upward correction for Group A fish that returned in spring (when ENL detection efficiency is lower due to higher flows) but were not observed.

---

## Criterion 3 — Reclassified Strict (Actual Entry Season)

**Definition:** Each fish is reclassified by its actual Entiat entry season — the calendar month of its first post-anchor Entiat site detection — rather than the RRF anchor season (which is always fall by design). Spring-entry fish (first Entiat detection January–June) are scored under a spring criterion: detection at ENL within 18 months is sufficient for success. Fall-entry Group A fish retain the strict fall criterion (ENL + upstream within 18 months). Fall-entry Group B fish are credited as successful if they have any Entiat detection at all, since navigating back past Wells Dam to reach the Entiat already demonstrates homing success; ENL detection alone is considered sufficient evidence.

For spring-entry fish, the ENL detection probability used in the detection-corrected model is recomputed at the actual spring entry date's flow (rather than the fall anchor date's flow), since spring flows at ENL are typically higher and detection efficiency lower.

**Rationale:** The all-fall rule (Criterion 2) underestimates spawning success for the roughly half of fish that overwinter between Rocky Reach and the Entiat mouth and then enter the system in spring. Those fish spawn in the lower Entiat reach and are correctly scored as successful when detected at ENL, without requiring upstream array detections they may never encounter. Reclassifying by actual entry season produces a more biologically accurate outcome.

**Key reclassification results (Group A):**

- Of 223 Group A fish, 113 (50.7%) entered the Entiat in spring; 87 (39.0%) entered in fall; 23 (10.3%) were never detected in the Entiat.
- 16 spring-entry Group A fish gained a spawned=1 coding under C3 that they lacked under C2 (previously required upstream detection they did not have).
- 2 fall-entry Group A fish that were incorrectly scored 0 under C2 due to a date-window calculation error were corrected to 1.

**Key reclassification results (Group B):**

- 4 Group B fall-entry fish (2 ENL-only, 2 with upstream detections but no ENL) gained spawned=1 coding; under Criterion 2 these fish lacked ENL + upstream in window, but under the Group B rule any Entiat detection constitutes success.
- 5 additional Group B spring-entry fish that were incorrectly scored 0 under C2 due to a date-window calculation error were corrected to 1.

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 150 | 67.3% |
| B | 181 | 108 | 59.7% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.372 | 0.69 | [−0.801, 0.051] | 0.958 |
| Detection-corrected | −0.518 | 0.60 | [−1.095, 0.055] | 0.964 |

The Wells effect is attenuated relative to Criterion 2, and the 95% credible interval includes zero under both model specifications. The reclassification raises Group B's success rate more than Group A's — the Group B fall-entry fix and the correction of misclassified spring-entry Group B fish together add proportionally more to Group B. The effect remains directionally consistent (P(β < 0) = 0.958–0.964), but is less strongly supported at this outcome definition.

---

## Criterion 4 — ENL Directionality Sensitivity (Upward)

**Definition:** Builds on Criterion 3. Additionally credits Group A fall-entry fish that have ENL-only detections (no upstream Entiat site in window) but whose ENL antenna records show upstream-directed movement. Two antenna patterns are credited: (1) *Entering* — detection sequence progressed from the lower/downstream ENL antennas to the upper/upstream antennas, indicating inriver passage; (2) *Upper/upstream-only with multiple observations* — fish detected exclusively on upstream antennas on more than one occasion, suggesting they had already passed through. Single-observation upper-antenna detections are not credited due to ambiguity. All other fish are unchanged from Criterion 3.

**Rationale:** ENL has two antenna tiers oriented to the direction of fish movement. A Group A fall-entry fish detected only at ENL (scored 0 under C3 because no upstream array confirmation) but whose antenna sequence shows it moving inriver may well have spawned in the unmonitored lower Entiat reach (approximately rkm 1–10, between ENL and the next fixed array ENM). Crediting these fish as an upward sensitivity tests whether the C3 result is conservative with respect to this unmonitored reach.

**Fish credited (0 → 1 under C4):**

Five Group A fall-entry fish were recoded from spawned=0 to spawned=1:

| Tag Code | Year | Antenna Pattern |
|----------|------|-----------------|
| 3D9.1C2C517ED3 | 2010 | Entering (Lower→Upper) |
| 3D9.1C2C51DAD4 | 2010 | Entering (Lower→Upper) |
| 3D9.1C2D774C1E | 2014 | Entering (Lower→Upper) |
| 3D9.1C2D8D4137 | 2013 | Entering (Lower→Upper) |
| 3D9.1C2DDAD709 | 2014 | Upper/Upstream only (multi-obs) |

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 155 | 69.5% |
| B | 181 | 108 | 59.7% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.476 | 0.62 | [−0.898, −0.058] | 0.987 |
| Detection-corrected | −0.691 | 0.50 | [−1.279, −0.113] | 0.991 |

Crediting the entering fish increases Group A's rate (69.5%) without changing Group B (59.7%), widening the gap and strengthening the Wells effect relative to C3. The 95% credible intervals exclude zero under both specifications, and P(β < 0) returns to 0.987–0.991. The ENL sensitivity is therefore a *conservative* check — if anything, applying directional evidence makes the Wells Dam effect appear larger, not smaller.

---

## Cross-Criteria Comparison

**Return rates by criterion:**

| Criterion | Group A (n=223) | Group B (n=181) | Gap (A−B) |
|-----------|-----------------|-----------------|-----------|
| C1 — Broad | 200 (89.7%) | 123 (68.0%) | 21.7 pp |
| C2 — Original strict | 136 (61.0%) | 91 (50.3%) | 10.7 pp |
| C3 — Reclassified | 150 (67.3%) | 108 (59.7%) | 7.6 pp |
| C4 — ENL sensitivity | 155 (69.5%) | 108 (59.7%) | 9.8 pp |

**Wells Dam effect (β) across criteria and model types:**

| Criterion | Standard β | OR | CI | P(β<0) | Det-corr β | OR | CI | P(β<0) |
|-----------|-----------|----|----|--------|-----------|----|----|--------|
| C1 — Broad | −1.348 | 0.26 | [−1.881, −0.834] | 1.000 | −1.836 | 0.16 | [−2.737, −1.067] | 1.000 |
| C2 — Original strict | −0.517 | 0.60 | [−0.937, −0.095] | 0.992 | −0.671 | 0.51 | [−1.191, −0.148] | 0.995 |
| C3 — Reclassified | −0.372 | 0.69 | [−0.801, 0.051] | 0.958 | −0.518 | 0.60 | [−1.095, 0.055] | 0.964 |
| C4 — ENL sensitivity | −0.476 | 0.62 | [−0.898, −0.058] | 0.987 | −0.691 | 0.50 | [−1.279, −0.113] | 0.991 |

**Key patterns:**

- The negative Wells Dam effect is present and directionally consistent across all four criteria and both model types.
- C1 (broad) shows the strongest effect because the raw detection gap between groups is largest when no movement threshold is required.
- C3 (reclassified) shows the weakest effect, with the 95% credible interval slightly straddling zero. Reclassification raises Group B's rate more than Group A's, narrowing the gap to 7.6 percentage points. This is the most biologically defensible outcome definition, but also the most conservative estimate of the Wells effect.
- C4 (ENL sensitivity) recovers most of the C2 effect magnitude and the credible interval returns to excluding zero. Because the directionality credit raises Group A but not Group B, the gap widens back toward C2 levels.
- Detection correction consistently strengthens the estimated Wells effect, reflecting that imperfect ENL detection (especially during high-flow spring conditions) more frequently causes Group A spring-entry fish to appear unsuccessful when they are not.
- The consistency of direction and the P(β < 0) values of 0.958–1.000 across all criteria provide robust evidence that Wells Dam interaction is associated with reduced Entiat return probability, regardless of which outcome definition is applied.

---

## Notes on Detection Correction

ENL detection efficiency was modeled as a function of daily discharge using a Bayesian GAM fit to paired mark-recapture data at the ENL array. Predicted detection probabilities (p_ENL) range from approximately 0.55 at peak spring flows to near 1.0 during low summer flows. In the detection-corrected models, the observed binary outcome is modeled with a compound likelihood: a fish scored as spawned=0 could either have truly not spawned or could have spawned but been missed at ENL. For Criteria 3 and 4, p_ENL is recomputed at the actual spring entry date's discharge for reclassified spring-entry fish (rather than the fall anchor date's flow), which appropriately accounts for the lower detection efficiency they would have faced on entry.
