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
| Detection-corrected | −2.035 | 0.13 | [−3.259, −0.929] | 0.999 |

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
| Detection-corrected | −1.120 | 0.33 | [−1.767, −0.491] | 1.000 |

The Wells effect is highly credible under both model specifications. The detection-corrected estimate is larger, reflecting upward correction for fish that returned but were missed at ENL — particularly spring-entry fish encountering higher flows. Detection correction is applied to both Group A and Group B fish based on each fish's actual ENL entry date (or seasonal mean flow when no ENL detection is available). Fish detected at ENL in both fall and spring — indicating two separate entry attempts — receive a compound detection probability: p_detect = 1 − (1 − p_fall)(1 − p_spring), which is always greater than either single-crossing probability.

---

## Criterion 3 — Reclassified Strict (Actual Entry Season)

**Definition:** Each fish is reclassified by its actual Entiat entry season — the calendar month of any post-anchor Entiat detection — rather than the RRF anchor season (which is always fall by design). A fish with any spring detection (January–June) at ENL or any upstream Entiat site is classified as a spring migrant. Only fish with exclusively fall Entiat detections and no subsequent spring activity are treated as fall fish.

Outcome scoring:

- **Spring-entry Group A fish**: detection at ENL in spring (months 1–6) OR at any upstream Entiat site within 18 months = success. ENL-only fall detections with no spring return or upstream confirmation = failure.
- **Fall-entry Group A fish**: confirmed upstream Entiat passage within 18 months = success (same as C2).
- **Group B fish** (spring or fall entry): any Entiat detection within 18 months = success. These fish navigated back past Wells Dam to reach the Entiat subbasin, which already demonstrates homing success regardless of how far upstream they traveled.

For spring-entry fish, the ENL detection probability used in the detection-corrected model is recomputed at the actual spring entry date's flow (rather than the fall anchor date's flow), since spring flows at ENL are typically higher and detection efficiency lower.

**Rationale:** Some fish overwinter between Rocky Reach and the Entiat mouth and enter the system in spring. Those fish spawn in the lower Entiat reach and may not be detected at upstream fixed arrays, even when present. Reclassifying by actual entry season and applying the appropriate seasonal criterion produces a more biologically accurate outcome. For Group B, requiring upstream confirmation is unnecessarily strict — a fish that overshot to Wells Dam and then successfully homed back to the Entiat has already demonstrated the key navigation behavior of interest.

**Key reclassification results (Group A):**

- Of 223 Group A fish, 113 (50.7%) entered the Entiat in spring; 87 (39.0%) entered in fall; 23 (10.3%) were never detected in the Entiat.
- Spring-entry Group A fish with spring ENL detections are credited under the spring criterion (spring ENL alone is sufficient), gaining credit that C2's fall-entry strict rule would withhold.
- Fall-entry Group A fish retain the C2 criterion (upstream confirmation required within 18 months).

**Key reclassification results (Group B):**

- Of 181 Group B fish, 104 (57.5%) entered the Entiat in spring; 18 (9.9%) entered in fall; 59 (32.6%) were never detected in the Entiat.
- All Group B fish with any Entiat detection are credited as successful under C3.
- Group B fish never detected in the Entiat (n=59) are scored 0; some of these may represent hidden strays to the Methow or Okanogan that were missed at those gateway arrays.

**Return rates:**

| Group | n | Spawned | Rate |
|-------|---|---------|------|
| A | 223 | 184 | 82.5% |
| B | 181 | 123 | 68.0% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.832 | 0.44 | [−1.305, −0.364] | 1.000 |
| Detection-corrected | −1.238 | 0.29 | [−2.031, −0.515] | 1.000 |

The Wells effect is strongly supported under both model specifications, with 95% credible intervals solidly excluding zero and P(β < 0) = 1.000. The reclassification by actual entry season raises Group A's rate (82.5%) relative to C2 (75.3%) because spring-entry Group A fish with spring ENL detections now receive credit. Group B's rate also rises relative to C2 (68.0% vs. 59.7%) because Group B fish with ENL-only detections — who navigated back past Wells Dam to reach the Entiat — are now credited. Despite both groups gaining, the gap narrows modestly (14.5 pp vs. 15.6 pp under C2), and the model-estimated Wells effect remains large and highly credible.

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
| A | 223 | 189 | 84.8% |
| B | 181 | 123 | 68.0% |

**Model results — Wells Dam effect:**

| Model | β | OR | 95% CI | P(β < 0) |
|-------|---|----|--------|----------|
| Standard | −0.548 | 0.58 | [−0.976, −0.126] | 0.994 |
| Detection-corrected | −0.772 | 0.46 | [−1.390, −0.183] | 0.995 |

Adding the five directionally confirmed fish raises Group A from 82.5% to 84.8%, while Group B remains at 68.0%. The raw detection gap widens to 16.8 percentage points. Both credible intervals exclude zero and P(β < 0) remains high (0.994–0.995), confirming that the directional ENL evidence does not undermine the C3 conclusion. The ENL directionality sensitivity therefore constitutes a conservative check — the Wells Dam effect is present regardless of whether these five fish are credited.

---

## Cross-Criteria Comparison

**Return rates by criterion:**

| Criterion | Group A (n=223) | Group B (n=181) | Gap (A−B) |
|-----------|-----------------|-----------------|-----------|
| C1 — Broad | 200 (89.7%) | 123 (68.0%) | 21.7 pp |
| C2 — Strict (upstream confirmed) | 168 (75.3%) | 108 (59.7%) | 15.6 pp |
| C3 — Reclassified | 184 (82.5%) | 123 (68.0%) | 14.5 pp |
| C4 — ENL sensitivity | 189 (84.8%) | 123 (68.0%) | 16.8 pp |

**Wells Dam effect (β) across criteria and model types:**

| Criterion | Standard β | OR | CI | P(β<0) | Det-corr β | OR | CI | P(β<0) |
|-----------|-----------|----|----|--------|-----------|----|----|--------|
| C1 — Broad | −1.348 | 0.26 | [−1.881, −0.834] | 1.000 | −2.035 | 0.13 | [−3.259, −0.929] | 0.999 |
| C2 — Strict upstream | −0.801 | 0.45 | [−1.246, −0.352] | 1.000 | −1.120 | 0.33 | [−1.767, −0.491] | 1.000 |
| C3 — Reclassified | −0.832 | 0.44 | [−1.305, −0.364] | 1.000 | −1.238 | 0.29 | [−2.031, −0.515] | 1.000 |
| C4 — ENL sensitivity | −0.548 | 0.58 | [−0.976, −0.126] | 0.994 | −0.772 | 0.46 | [−1.390, −0.183] | 0.995 |

**Key patterns:**

- The negative Wells Dam effect is present and directionally consistent across all four criteria and both model types.
- C1 (broad) shows the strongest effect because the raw detection gap between groups is largest when no movement threshold is required.
- C2 and C3 produce similar model-estimated effects (β ≈ −0.80 to −0.83 standard; β ≈ −1.12 to −1.24 detection-corrected). The reclassification raises both Group A and Group B rates — Group A gains spring-entry fish credited by the spring ENL criterion, and Group B gains fish previously penalized under C2's upstream-confirmation requirement. The net gap narrows slightly (14.5 pp vs. 15.6 pp), but the model estimates remain very similar because the hierarchical model accounts for within-year structure.
- C4 (ENL sensitivity) credits five additional Group A fish based on upstream-directed antenna patterns, raising the raw detection gap to 16.8 pp. The model-estimated effect (β = −0.548 standard; β = −0.772 detection-corrected) remains firmly negative with credible intervals excluding zero, confirming that the C3 result does not depend on conservative treatment of the unmonitored lower Entiat reach.
- Detection correction consistently strengthens the estimated Wells effect across all criteria, reflecting that imperfect ENL detection — especially during high-flow spring conditions — more frequently causes fish to appear unsuccessful when they are not. For fish with two ENL crossing opportunities (fall entry followed by spring re-entry), a compound detection probability is used: p_detect = 1 − (1 − p_fall)(1 − p_spring).
- The consistency of direction and the P(β < 0) values of 0.994–1.000 across all criteria provide robust evidence that Wells Dam interaction is associated with reduced Entiat return probability, regardless of which outcome definition is applied.

---

## Notes on Detection Correction

ENL detection efficiency was modeled as a function of daily discharge using a Bayesian GAM fit to paired mark-recapture data at the ENL array. Predicted detection probabilities (p_ENL) range from approximately 0.55 at peak spring flows to near 1.0 during low summer flows. Detection correction is applied to both Group A and Group B fish.

In the detection-corrected models, the observed binary outcome is modeled with a compound likelihood: a fish scored as spawned=0 could either have truly not spawned or could have spawned but been missed at ENL. The per-fish p_ENL is computed at the flow conditions prevailing on the fish's actual ENL entry date. For fish with no ENL detection (upstream passage only), p_ENL is set to the year-season mean for the fish's actual entry season.

For fish detected at ENL in both fall and spring — indicating two separate inriver entry attempts — a compound detection probability is applied: p_detect = 1 − (1 − p_fall)(1 − p_spring), where p_fall and p_spring are the predicted detection efficiencies at the respective crossing dates. This compound probability is always greater than either single-crossing probability, appropriately reflecting the reduced chance of missing a fish on either of two detection opportunities.

For Criteria 3 and 4, p_ENL is recomputed at the actual spring entry date's discharge for reclassified spring-entry fish (rather than the fall anchor date's flow), which appropriately accounts for the lower detection efficiency they would have faced on entry.
