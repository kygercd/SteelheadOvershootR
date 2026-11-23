## Robustness to Sparse Recent-Year Data (2022--2024)

The three most recent spill years (2022, 2023, and 2024) each contributed
only a single fish, all of which returned downstream (100% apparent return
rate). Because these years sit at the extremes of 0% or 100% by
construction, they could in principle exert disproportionate leverage on
the model, particularly on the year-level variance and on any spill
coefficient tied to a year-specific spill value. The hierarchical model
structure is designed to guard against this: each year's intercept
(α_j) is drawn from a shared Normal(μ_α, σ_α) distribution, so a
single-fish year is partially pooled toward the overall mean rather than
being treated as direct evidence of a 0% or 100% true return rate. The
wide posterior uncertainty contributed by these years is absorbed into
σ_year rather than distorting the fixed-effect estimates.

To confirm this empirically, we refit all six hierarchical models after
excluding the 2022--2024 spill years (n = 1,318 fish across 17 years),
using identical priors and MCMC settings. The harvest exposure effect was
essentially unchanged across all six models (β_harvest ranged from −2.18
to −2.30 versus −2.18 to −2.33 in the full dataset; OR ≈ 0.10--0.11 in
both cases), and remained strongly credible (95% HDI excluding zero in
all models). Spill effects were similarly stable: the Fall SpillHours
model remained the only spill metric with a 95% HDI excluding zero
(β = −0.27 [−0.52, −0.04] versus −0.30 [−0.54, −0.07] with the full
dataset), and all other spill metrics continued to span zero with
estimates differing from the full-data results by less than 0.05 on the
logit scale. These results indicate that the three sparse, perfectly
separated recent years do not meaningfully bias either the harvest or
spill effect estimates, and that the partial-pooling structure of the
hierarchical model appropriately limits their influence on the overall
inference.
