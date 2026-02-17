# Steelhead Overshoot Hierarchical Analysis with Harvest (R Version)

## Overview

This analysis examines whether steelhead that "overshoot" their natal tributaries (Wenatchee and Entiat rivers) and reach Wells Dam on the Upper Columbia River are more likely to return downstream versus remain at Wells Dam. The model tests the effects of dam spill conditions and harvest exposure on downstream return probability.

## Analysis Script

### `run_full_hierarchical_analysis_with_harvest.R`

A Bayesian hierarchical logistic regression implemented in R using `brms` (Stan backend)

**Model structure:**

```
y_ij ~ Bernoulli(p_ij)
logit(p_ij) = alpha_j + beta_spill * Spill_j + beta_harvest * Harvest_i

alpha_j ~ Normal(mu_alpha, sigma_alpha)   [year random intercepts]
mu_alpha ~ Normal(0, 2)
sigma_alpha ~ HalfNormal(0, 1)
beta_spill ~ Normal(0, 1)                 [population-level spill effect]
beta_harvest ~ Normal(0, 1)               [individual-level harvest effect]
```

Where:
- `i` = individual fish
- `j` = spill year
- `y_ij` = 1 if fish returned downstream, 0 if last detected at Wells Dam
- `alpha_j` = year-specific baseline (partial pooling via random intercepts)
- `Spill_j` = standardized spill metric for that year
- `Harvest_i` = 1 if hatchery fish in a harvest-open year, 0 otherwise

**What the script does:**

1. Loads and merges raw data files (fish detections, rear type, harvest status, spill, site classifications)
2. Reclassifies fish fates using corrected site zone assignments (notably NAU = downstream)
3. Creates the binary harvest exposure covariate (hatchery fish in open-harvest years)
4. Assigns fish to a "spill year" using a June 1 cutoff: fish first detected at Wells Dam from June of year Y through May of year Y+1 are assigned to SpillYear Y+1. Most overshooting steelhead arrive at Wells in late summer through fall, so this grouping ensures each fish is paired with the fall spill conditions present when they arrived and the spring spill conditions of the following year while they may still be holding at the dam. For example, a fish arriving in October 2017 is assigned SpillYear 2018, and is matched with fall spill from Aug-Nov 2017 (the fall it arrived) and spring spill from Jan-Mar 2018 (the subsequent spring).
5. Calculates spill metrics for two seasonal windows:
   - Fall: August 1 - November 30 of year Y, assigned to SpillYear Y+1
   - Spring: January 1 - March 31 of year Y, assigned to SpillYear Y
6. Fits 6 hierarchical Bayesian models (3 spill metrics x 2 seasons):
   - SpillHours, SpillDays, TotalSpill for each window
7. Reports posterior summaries: means, 95% HDI, P(beta > 0), odds ratios
8. Generates diagnostic and visualization figures (PNG)
9. Saves results to CSV and model objects to RData

**MCMC settings:** 4 chains, 4500 iterations (1500 warmup + 3000 sampling), adapt_delta = 0.95.

**Required R packages:** `brms`, `tidyverse`, `lubridate`, `bayesplot`, `posterior`, `patchwork`, `tidybayes` (optionally `bayestestR` for HDI computation). The script will attempt to install missing packages automatically.

## Source Data Files

### `all_overshoot_fish_2001_2025.csv`

Primary detection records for 9,688 steelhead that overshot their natal tributaries and were detected at Wells Dam (2001-2025).

| Column | Description |
|--------|-------------|
| TagCode | PIT tag identifier |
| Subbasin | Natal subbasin (Wenatchee, Upper Columbia-Entiat) |
| ReleaseSite | Tagging/release location |
| FirstWellsDate | Date/time of first Wells Dam detection |
| LastDetectionDate | Date/time of last detection anywhere |
| LastDetectionSite | Site code of final detection |
| Fate | Original fate classification (Downstream, At Wells, Above) |
| Basin | Simplified basin (Wenatchee, Entiat) |
| SpillYear | Assigned spill year (June-May cycle) |

### `OvershootRearType.csv`

Rearing origin classification for each tagged fish.

| Column | Description |
|--------|-------------|
| Tag Code | PIT tag identifier |
| Rear Type Code | W = Wild/Natural, H = Hatchery |
| Rear Type Name | Full rearing type description |
| Mark Site Name | Location where fish was marked |
| Release Site Name | Location where fish was released |

### `Harvest_Open_Closed.csv`

Annual harvest regulation status for hatchery steelhead in the study area.

| Column | Description |
|--------|-------------|
| Year | Calendar year |
| Harvest_Status | "Open" (2000-2015, 2024-2025) or "Closed" (2016-2023) |

The harvest closure period (2016-2023) eliminated directed harvest of hatchery steelhead in the Upper Columbia tributaries.

### `merged_spill_2000_2024.csv`

Hourly spill discharge records at Wells Dam (2000-2024).

| Column | Description |
|--------|-------------|
| DateTime | Hourly timestamp |
| SPILL (KCFS) | Spill discharge in thousand cubic feet per second |

Used to calculate three spill metrics per seasonal window: total hours of spill, number of days with any spill, and cumulative spill volume.

### `site_classifications.csv`

Lookup table mapping PIT tag detection site codes to geographic zones relative to Wells Dam.

| Column | Description |
|--------|-------------|
| Site | Detection site code |
| Zone | Above, Below, or Wells Dam |
| Fish_Count | Number of fish last detected at this site |

This file contains corrected classifications, notably assigning NAU (Upper Nason Creek) to the "Below" zone, which corrects the original fate assignment for fish detected there.

## Output Files (generated by the script)

- `full_hierarchical_results_with_harvest_R.csv` - Model parameter estimates for all 6 models
- `full_year_summary_with_harvest_R.csv` - Year-level data summary
- `analysis_fish_with_harvest_R.csv` - Final analysis dataset
- `hierarchical_models_with_harvest_R.RData` - Saved R model objects
- `full_data_overview_with_harvest_R.png` - Data summary panels
- `full_spill_effects_with_harvest_R.png` - Spill effect posterior distributions
- `full_harvest_effects_R.png` - Harvest effect posterior distributions
- `full_year_effects_with_harvest_R.png` - Year random effect estimates
- `full_trace_plots_R.png` - MCMC trace plot diagnostics

## Running the Analysis

```r
source("run_full_hierarchical_analysis_with_harvest.R")
```

Or open in RStudio and run interactively. Stan compilation occurs on the first model fit and is cached for subsequent runs.
