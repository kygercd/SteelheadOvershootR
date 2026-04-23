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
| Detection-corrected | −2.046 | 0.13 | [−3.262, −0.923] | 0.999 |

Wells Dam interaction is strongly associated with reduced broad return probability. The detection-corrected model accounts for the possibility that some Group A fish returned but were not detected at ENL, which would bias the observed Group A rate downward; correcting for this makes the Wells effect appear larger.

---

## Criterion 2 — Original Strict (Confirmed Upstream Passage)

**Definition:** A fish is considered to have spawned if it was detected at any Entiat subbasin site AND at least one upstream Entiat site (ENM or higher) within 18 months of the fall RRF anchor. Fish confirmed above ENM (with or without an ENL gateway detection) are credited as successful, as confirmed upstream passage constitutes strong evidence of spawning.

**Rationale:** The strict criterion requires evidence of upstream passage, not just arrival at the river mouth or subbasin. The 18-month window captures one complete spawning season after the anchor event. This is more conservative than C1 (broad detection) because ENL-only detections in fall fish — where ENL is at the river mouth with no spawning habitat below — are not sufficient for success unless upstream confirmation follows.

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 168 | 75.3% |
| B | 181 | 108 | 59.7% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.801 | 0.45 | [−1.246, −0.352] | 1.000 |
| Detection-corrected | −1.115 | 0.33 | [−1.782, −0.488] | 1.000 |

The Wells effect is highly credible under both model specifications. The detection-corrected estimate is larger, reflecting upward correction for Group A fish that returned but were missed at ENL — particularly spring-entry fish encountering higher flows.

---

## Criterion 3 — Reclassified Strict (Actual Entry Season)

**Definition:** Each fish is reclassified by its actual Entiat entry season — the calendar month of its first post-anchor Entiat site detection — rather than the RRF anchor season (which is always fall by design). Spring-entry fish (first Entiat detection January–June) are scored under a spring criterion: detection at ENL within 18 months is sufficient for success. Fall-entry Group A fish retain the strict criterion (confirmed upstream Entiat passage within 18 months). Fall-entry Group B fish are credited as successful if they have any Entiat detection at all, since navigating back past Wells Dam to reach the Entiat already demonstrates homing success.

For spring-entry fish, the ENL detection probability used in the detection-corrected model is recomputed at the actual spring entry date's flow (rather than the fall anchor date's flow), since spring flows at ENL are typically higher and detection efficiency lower.

**Rationale:** Some fish overwinter between Rocky Reach and the Entiat mouth and enter the system in spring. Those fish spawn in the lower Entiat reach and may not be detected at upstream fixed arrays, even when present. Reclassifying by actual entry season and applying the appropriate seasonal criterion produces a more biologically accurate outcome. However, this criterion requires ENL detection specifically for spring-entry fish, which means fish confirmed upstream without an ENL gateway detection are scored more conservatively than under C2.

**Key reclassification results (Group A):**

- Of 223 Group A fish, 113 (50.7%) entered the Entiat in spring; 87 (39.0%) entered in fall; 23 (10.3%) were never detected in the Entiat.
- Spring-entry Group A fish with ENL detections are credited under the spring criterion (ENL alone is sufficient), gaining credit that fall-entry criteria would withhold.
- Spring-entry Group A fish with only upstream detections (no ENL gateway detection) are scored 0 under C3 (ENL required for spring criterion), whereas under C2 confirmed upstream passage was sufficient.

**Key reclassification results (Group B):**

- Fall-entry Group B fish are credited for any Entiat detection under C3, relaxing the upstream confirmation requirement.
- Spring-entry Group B fish with ENL detections are credited under the spring criterion.
- Group B fish confirmed upstream in the Entiat but lacking an ENL detection are scored 0 under C3's spring criterion if spring-entry, producing a small set of reclassification losses that offset gains elsewhere. Net Group B count is unchanged from C2 (n=108).

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 153 | 68.6% |
| B | 181 | 108 | 59.7% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.444 | 0.64 | [−0.864, −0.015] | 0.979 |
| Detection-corrected | −0.629 | 0.53 | [−1.235, −0.048] | 0.983 |

The Wells effect is attenuated relative to C2 and the credible intervals barely exclude zero under both specifications. P(β < 0) remains high (0.979–0.983), indicating the effect is directionally consistent though weaker than under the broader criteria. The attenuation reflects that reclassification by actual entry season credits some spring-entry Group A fish with ENL-only detections that were missed under the fall-entry rule (raising Group A), while simultaneously removing credit from Group A fish confirmed upstream without an ENL gateway detection that were counted under C2.

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
| A | 223 | 158 | 70.9% |
| B | 181 | 108 | 59.7% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.548 | 0.58 | [−0.976, −0.126] | 0.994 |
| Detection-corrected | −0.772 | 0.46 | [−1.390, −0.183] | 0.995 |

Crediting the entering fish increases Group A's rate (70.9%) without changing Group B (59.7%), widening the gap and strengthening the Wells effect relative to C3. The 95% credible intervals exclude zero under both specifications, and P(β < 0) rises to 0.994–0.995. The ENL sensitivity is therefore a *conservative* check — if anything, applying directional evidence makes the Wells Dam effect appear larger, not smaller.

---

## Cross-Criteria Comparison

**Return rates by criterion:**

| Criterion | Group A (n=223) | Group B (n=181) | Gap (A−B) |
|-----------|-----------------|-----------------|-----------|
| C1 — Broad | 200 (89.7%) | 123 (68.0%) | 21.7 pp |
| C2 — Strict (upstream confirmed) | 168 (75.3%) | 108 (59.7%) | 15.6 pp |
| C3 — Reclassified | 153 (68.6%) | 108 (59.7%) | 8.9 pp |
| C4 — ENL sensitivity | 158 (70.9%) | 108 (59.7%) | 11.2 pp |

**Wells Dam effect (β) across criteria and model types:**

| Criterion | Standard β | OR | CI | P(β<0) | Det-corr β | OR | CI | P(β<0) |
|-----------|-----------|----|----|--------|-----------|----|----|--------|
| C1 — Broad | −1.348 | 0.26 | [−1.881, −0.834] | 1.000 | −2.046 | 0.13 | [−3.262, −0.923] | 0.999 |
| C2 — Strict upstream | −0.801 | 0.45 | [−1.246, −0.352] | 1.000 | −1.115 | 0.33 | [−1.782, −0.488] | 1.000 |
| C3 — Reclassified | −0.444 | 0.64 | [−0.864, −0.015] | 0.979 | −0.629 | 0.53 | [−1.235, −0.048] | 0.983 |
| C4 — ENL sensitivity | −0.548 | 0.58 | [−0.976, −0.126] | 0.994 | −0.772 | 0.46 | [−1.390, −0.183] | 0.995 |

**Key patterns:**

- The negative Wells Dam effect is present and directionally consistent across all four criteria and both model types.
- C1 (broad) shows the strongest effect because the raw detection gap between groups is largest when no movement threshold is required.
- C3 (reclassified) shows the weakest effect — the 95% credible intervals barely exclude zero under both model specifications. Reclassifying by actual entry season narrows the gap to 8.9 percentage points; the spring criterion credits some Group A fish with ENL-only detections but simultaneously removes credit for Group A fish confirmed upstream without an ENL gateway detection (previously counted under C2), producing a net reduction in Group A relative to C2.
- C4 (ENL sensitivity) recovers part of the C2–C3 gap: crediting fall-entry Group A fish with upstream-directed ENL antenna patterns raises Group A to 70.9% (11.2 pp gap), and the credible interval excludes zero under both specifications (P(β < 0) = 0.994–0.995).
- Detection correction consistently strengthens the estimated Wells effect across all criteria, reflecting that imperfect ENL detection — especially during high-flow spring conditions — more frequently causes Group A spring-entry fish to appear unsuccessful when they are not.
- The consistency of direction and the P(β < 0) values of 0.979–1.000 across all criteria provide robust evidence that Wells Dam interaction is associated with reduced Entiat return probability, regardless of which outcome definition is applied.

---

## Notes on Detection Correction

ENL detection efficiency was modeled as a function of daily discharge using a Bayesian GAM fit to paired mark-recapture data at the ENL array. Predicted detection probabilities (p_ENL) range from approximately 0.55 at peak spring flows to near 1.0 during low summer flows. In the detection-corrected models, the observed binary outcome is modeled with a compound likelihood: a fish scored as spawned=0 could either have truly not spawned or could have spawned but been missed at ENL. For Criteria 3 and 4, p_ENL is recomputed at the actual spring entry date's discharge for reclassified spring-entry fish (rather than the fall anchor date's flow), which appropriately accounts for the lower detection efficiency they would have faced on entry.
