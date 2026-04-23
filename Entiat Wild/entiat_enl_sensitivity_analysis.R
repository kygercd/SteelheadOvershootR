# ENL Directionality Sensitivity Analysis
#
# Builds on the reclassified strict analysis (entiat_reclassified_strict_analysis.R),
# which scores all fish by their actual Entiat entry season:
#   - Spring-entry fish (Groups A and B): ENL alone within 18 months = spawned
#   - Fall-entry fish: ENL + ≥1 upstream site within 18 months = spawned
#
# After that reclassification, some Group A fall-entry fish still have ENL-only
# detections (scored 0) but their antenna records show they passed through ENL
# in the upstream direction. These fish may have spawned in the unmonitored
# lower Entiat reach (rkm 1–~10, between ENL and ENM).
#
# This sensitivity credits those fish as spawned and compares against the
# reclassified strict base. Group B is not modified — any Entiat detection
# post-Wells is already treated as a successful return.
#
# Sensitivity outcome (EntiatSpawned_s):
#   Group A fall-entry ENL-only fish with "Entering (Lower→Upper)" or confirmed
#   "Upper/Upstream only" antenna pattern → recoded from 0 to 1 in fish_r.
#   All other fish: unchanged from EntiatSpawned_r.
#
# Models:
#   S1s: EntiatSpawned_s ~ wells + harvest + (1|year)    [standard]
#   S2s: detection-corrected compound likelihood          [det-corr]
#
# Comparison: original strict → reclassified strict → ENL-sensitivity

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(lubridate)
  library(tidybayes)
  library(purrr)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

MCMC_ITER   <- 4500
MCMC_WARMUP <- 1500
MCMC_CHAINS <- 4
MCMC_CORES  <- 4
ADAPT_DELTA <- 0.95

# ── Load reclassified strict base ──────────────────────────────────────────────
load("entiat_reclassified_strict_models.RData")   # fish_r, fit_s1r, fit_s2r
load("entiat_cleaned_strict_models.RData")        # fit_s1c, fit_s2c for comparison

# ── Load directionality classifications ────────────────────────────────────────
dir_df <- read_csv("~/ptagis_data/enl_directionality_ga_enl_only_29fish.csv",
                   show_col_types = FALSE) |>
  select(TagCode, direction_class, n_obs)

# Credit fish that passed through ENL in the upstream direction
entering_classes <- c(
  "Entering (Lower→Upper)",
  "Upper/Upstream only"     # confirmed multi-obs only; single-obs excluded
)

dir_df <- dir_df |>
  mutate(enl_entering = direction_class %in% entering_classes)

cat("Directionality summary for Group A fall-entry ENL-only fish:\n")
print(count(dir_df, direction_class, enl_entering))

# ── Identify eligible fish: Group A, fall-entry, ENL-only, EntiatSpawned_r == 0
eligible <- fish_r |>
  filter(group == "A", spring_return_r == 0, EntiatSpawned_r == 0) |>
  pull(TagCode)
cat(sprintf("\nEligible Group A fall-entry ENL-only fish (EntiatSpawned_r=0): %d\n",
            length(eligible)))

dir_eligible <- dir_df |> filter(TagCode %in% eligible)
cat(sprintf("Of these, with entering directionality (to credit): %d\n",
            sum(dir_eligible$enl_entering, na.rm = TRUE)))

# ── Build sensitivity dataset ──────────────────────────────────────────────────
fish_s <- fish_r |>
  left_join(dir_df |> select(TagCode, enl_entering), by = "TagCode") |>
  mutate(
    enl_entering    = replace_na(enl_entering, FALSE),
    EntiatSpawned_s = case_when(
      group == "A" & spring_return_r == 0 & EntiatSpawned_r == 0 & enl_entering ~ 1L,
      TRUE ~ EntiatSpawned_r
    )
  )

cat("\nOutcome change from directionality sensitivity:\n")
fish_s |>
  count(group, EntiatSpawned_r, EntiatSpawned_s) |>
  mutate(change = case_when(
    EntiatSpawned_r == EntiatSpawned_s ~ "unchanged",
    EntiatSpawned_r == 0 & EntiatSpawned_s == 1 ~ "0 → 1 (credited)",
    TRUE ~ "other"
  )) |>
  print()

cat("\nGroup × season × outcome (sensitivity):\n")
fish_s |>
  mutate(SeasonR = if_else(spring_return_r == 1, "Spring", "Fall")) |>
  group_by(group, SeasonR) |>
  summarise(n = n(),
            spawned_reclassified = sum(EntiatSpawned_r),
            spawned_sens         = sum(EntiatSpawned_s),
            pct_reclassified     = round(mean(EntiatSpawned_r) * 100, 1),
            pct_sens             = round(mean(EntiatSpawned_s) * 100, 1),
            .groups = "drop") |>
  print()

# ── brms custom family (detection-corrected) ───────────────────────────────────
det_corrected_family <- custom_family(
  "det_corrected",
  dpars = "mu", links = "logit", type = "int",
  vars  = "vreal1[n]"
)
stan_functions <- "
  real det_corrected_lpmf(int y, real mu, real vreal1) {
    real pi_obs = mu * vreal1;
    pi_obs = fmax(1e-10, fmin(1 - 1e-10, pi_obs));
    return bernoulli_lpmf(y | pi_obs);
  }
"
sv <- stanvar(scode = stan_functions, block = "functions")

model_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

# ── Model S1s: standard ────────────────────────────────────────────────────────
cat("\nMODEL S1s: EntiatSpawned_s ~ wells + harvest + (1|year)\n")
fit_s1s <- brm(
  EntiatSpawned_s ~ wells_interaction + harvest_open + (1 | ReturnYearFactor),
  data    = fish_s,
  family  = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500,
  file    = "fit_s1s"
)

post_s1s <- as_draws_df(fit_s1s)
b_s1s    <- post_s1s$b_wells_interaction
ci_s1s   <- quantile(b_s1s, c(0.025, 0.975))
cat(sprintf("  β=%.3f  OR=%.2f  95%%CI=[%.3f, %.3f]  P(β<0)=%.3f\n",
            mean(b_s1s), exp(mean(b_s1s)), ci_s1s[1], ci_s1s[2], mean(b_s1s < 0)))

# ── Model S2s: detection-corrected ────────────────────────────────────────────
cat("\nMODEL S2s: detection-corrected\n")
fit_s2s <- brm(
  bf(EntiatSpawned_s | vreal(p_enl_r_safe) ~
       wells_interaction + harvest_open + (1 | ReturnYearFactor)),
  data      = fish_s |> filter(!is.na(p_enl_r_safe)),
  family    = det_corrected_family,
  stanvars  = sv,
  prior     = model_priors,
  iter      = MCMC_ITER, warmup = MCMC_WARMUP,
  chains    = MCMC_CHAINS, cores  = MCMC_CORES,
  seed      = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent    = 2, refresh = 500,
  file      = "fit_s2s"
)

post_s2s <- as_draws_df(fit_s2s)
b_s2s    <- post_s2s$b_wells_interaction
ci_s2s   <- quantile(b_s2s, c(0.025, 0.975))
cat(sprintf("  β=%.3f  OR=%.2f  95%%CI=[%.3f, %.3f]  P(β<0)=%.3f\n",
            mean(b_s2s), exp(mean(b_s2s)), ci_s2s[1], ci_s2s[2], mean(b_s2s < 0)))

# ── Comparison table ───────────────────────────────────────────────────────────
make_row <- function(label, b_vec, n) {
  ci <- quantile(b_vec, c(0.025, 0.975))
  tibble(Model = label, n = n,
         Beta  = round(mean(b_vec), 3),
         OR    = round(exp(mean(b_vec)), 2),
         CI_lo = round(ci[1], 3),
         CI_hi = round(ci[2], 3),
         P_neg = round(mean(b_vec < 0), 3))
}

comparison_s <- bind_rows(
  make_row("Strict (original) — standard",          as_draws_df(fit_s1c)$b_wells_interaction,   nrow(fit_s1c$data)),
  make_row("Strict (original) — det-corr",           as_draws_df(fit_s2c)$b_wells_interaction,   nrow(fit_s2c$data)),
  make_row("Strict (reclassified) — standard",       as_draws_df(fit_s1r)$b_wells_interaction,   nrow(fit_s1r$data)),
  make_row("Strict (reclassified) — det-corr",       as_draws_df(fit_s2r)$b_wells_interaction,   nrow(fit_s2r$data)),
  make_row("Strict (ENL-sensitivity) — standard",    b_s1s, nrow(fit_s1s$data)),
  make_row("Strict (ENL-sensitivity) — det-corr",    b_s2s, nrow(fit_s2s$data))
)

cat("\nCOMPARISON: Original → Reclassified → ENL-sensitivity\n")
print(comparison_s, n = Inf)

# ── Save ──────────────────────────────────────────────────────────────────────
save(fit_s1s, fit_s2s, fish_s, comparison_s,
     file = "entiat_enl_sensitivity_models.RData")
cat("\nSaved: entiat_enl_sensitivity_models.RData\n")

write_csv(comparison_s, "entiat_enl_sensitivity_comparison.csv")
cat("Saved: entiat_enl_sensitivity_comparison.csv\n")
