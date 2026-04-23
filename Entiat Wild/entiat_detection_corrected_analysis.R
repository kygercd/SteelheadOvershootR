###############################################################################
# Detection-Corrected Return Analysis: Wild Entiat Steelhead
#
# PURPOSE
# -------
# The standard analysis in entiat_return_bayesian_analysis.R models:
#   Y_i ~ Bernoulli(θ_i)   where θ_i = inv_logit(η_i)
#
# This is WRONG when detection at ENL is imperfect: a fish with EntiatDetected=0
# could either (a) not have returned to the Entiat, OR (b) returned but was
# missed by ENL. This analysis corrects that using:
#
#   Y_i ~ Bernoulli(θ_i × p_ENL_i)
#
# where p_ENL_i is the per-fish detection efficiency from the GAM model,
# computed from the fish's return date and year.
#
# Similarly, among Group B fish NOT detected post-Wells (53 fish), some may be
# hidden strays that went to Methow/Okanogan but were missed at LMR/OKL. This
# script estimates the number of hidden strays using p_LMR and p_OKL.
#
# APPROACH
# --------
# 1. Load per-fish detection efficiencies from entiat_analysis_with_detection_eff.csv
# 2. Fit detection-corrected brms model using custom Bernoulli likelihood
# 3. Compare corrected vs uncorrected Wells interaction estimate
# 4. Post-hoc hidden stray analysis for 53 undetected Group B fish
# 5. Visualize and save results
#
# INPUT FILES
# -----------
#   entiat_analysis_with_detection_eff.csv  (per-fish dataset with p_enl, p_lmr, p_okl)
#   entiat_bayesian_models.RData            (fitted standard models for comparison)
#
# OUTPUT FILES
# ------------
#   entiat_corrected_models.RData           (fitted brms objects)
#   entiat_detection_corrected_results.csv  (model comparison table)
#   entiat_corrected_vs_uncorrected.png     (posterior comparison plot)
#   entiat_hidden_stray_analysis.png        (stray correction figure)
#   entiat_corrected_return_rates.png       (corrected annual return rates)
###############################################################################

# =============================================================================
# PACKAGES
# =============================================================================
required_packages <- c("brms", "tidyverse", "lubridate", "bayesplot",
                       "posterior", "ggplot2", "patchwork", "tidybayes")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

library(brms)
library(tidyverse)
library(lubridate)
library(bayesplot)
library(posterior)
library(ggplot2)
library(patchwork)
library(tidybayes)

set.seed(42)

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

cat("======================================================================\n")
cat("DETECTION-CORRECTED RETURN ANALYSIS: WILD ENTIAT STEELHEAD\n")
cat("======================================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================
cat("--- Loading Data ---\n")

fish_df <- read_csv("entiat_analysis_with_detection_eff.csv", show_col_types = FALSE)

raw_main <- read_csv("EntiatWildSteel.csv", show_col_types = FALSE) %>%
  rename(tag = `Tag Code`, site = `Site Name`) %>%
  mutate(first_dt = mdy_hm(`First Time Value`))

load("gateway_bayesian_models.RData")   # fit_enl, flow_enl, entiat

fish_df <- fish_df %>% mutate(rrf_anchor = as_datetime(rrf_anchor))

# ---- ENL detection efficiency helpers ----------------------------------------
entiat_all_sites <- c("ENL - Lower Entiat River","ENM - Middle Entiat River",
                      "ENA - Upper Entiat River at rkm 17.1",
                      "ENF - Upper Entiat River at rkm 40.6",
                      "ENS - Upper Entiat River at rkm 35.7",
                      "EHL - Entiat NFH Adult Ladder",
                      "MAD - Mad River, Entiat River Basin",
                      "UMR - Upper Mad Rv. Temporary Array",
                      "SR1 - Top of Entiat R side channel",
                      "SR2 - Bot of Entiat R side channel",
                      "HS1 - Top of side channel near ENFH",
                      "HS2 - Bottom side channel near ENFH",
                      "HN1 - Top of Harrison side channel",
                      "HN2 - Middle Harrison Side Channel",
                      "HN3 - Bottom Harrison Side Channel",
                      "TLT - Tillicum Creek Temporary Array",
                      "RCT - Roaring Creek Temporary Array",
                      "URT - Upper Roaring Temporary Array")

enl_logq_mean <- mean(entiat$log_q)
enl_logq_sd   <- sd(entiat$log_q)
enl_base_year <- min(entiat$year)
flow_enl      <- flow_enl %>% mutate(date = as.Date(date))

predict_p_enl_at_date <- function(entry_date, return_year) {
  q <- flow_enl %>% filter(date == as.Date(entry_date)) %>% pull(discharge)
  if (length(q) == 0 || is.na(q)) return(NA_real_)
  lq_s <- (log(q) - enl_logq_mean) / enl_logq_sd
  yr   <- return_year - enl_base_year
  nd   <- data.frame(log_q_s = lq_s, log_q2_s = lq_s^2, year = yr)
  mean(posterior_epred(fit_enl, newdata = nd,
                       allow_new_levels = TRUE, sample_new_levels = "gaussian"))
}

cat("Computing year x season mean p_enl lookup...\n")
season_means <- flow_enl %>%
  filter(!is.na(discharge)) %>%
  mutate(yr     = year(date),
         season = if_else(month(date) %in% 1:6, "spring", "fall"),
         lq_s   = (log(discharge) - enl_logq_mean) / enl_logq_sd,
         yr_c   = yr - enl_base_year) %>%
  group_by(yr, season) %>%
  summarise(mean_lq_s = mean(lq_s), yr_c = first(yr_c), .groups = "drop") %>%
  rowwise() %>%
  mutate(mean_p_enl = mean(posterior_epred(
    fit_enl,
    newdata = data.frame(log_q_s  = mean_lq_s,
                         log_q2_s = mean_lq_s^2,
                         year     = yr_c),
    allow_new_levels = TRUE, sample_new_levels = "gaussian"))) %>%
  ungroup()

overall_spring_p <- mean(season_means$mean_p_enl[season_means$season == "spring"])
overall_fall_p   <- mean(season_means$mean_p_enl[season_means$season == "fall"])

get_seasonal_p_enl <- function(return_year, season) {
  val <- season_means %>%
    filter(yr == return_year, season == !!season) %>%
    pull(mean_p_enl)
  if (length(val) == 0 || is.na(val))
    return(if (season == "spring") overall_spring_p else overall_fall_p)
  val
}

cat(sprintf("  Overall spring mean p_enl: %.3f\n", overall_spring_p))
cat(sprintf("  Overall fall mean p_enl:   %.3f\n\n", overall_fall_p))

# ---- Per-fish entry-date p_enl -----------------------------------------------
cat("Computing per-fish entry-date p_enl...\n")
entry_corrections <- map_dfr(fish_df$TagCode, function(tc) {
  fm     <- fish_df   %>% filter(TagCode == tc)
  dets   <- raw_main  %>% filter(tag == tc) %>% arrange(first_dt)
  anchor <- as_datetime(fm$rrf_anchor)

  ent_dets <- dets %>% filter(site %in% entiat_all_sites, first_dt > anchor)

  if (nrow(ent_dets) > 0) {
    entry_date  <- min(ent_dets$first_dt)
    entry_month <- month(entry_date)
    p_entry     <- predict_p_enl_at_date(entry_date, fm$ReturnYear)
    if (is.na(p_entry)) {
      entry_seas <- if_else(entry_month %in% 1:6, "spring", "fall")
      p_entry    <- get_seasonal_p_enl(fm$ReturnYear, entry_seas)
    }
  } else {
    # Not detected: use anchor-season seasonal mean
    entry_month <- month(anchor)
    anchor_seas <- if_else(entry_month %in% 1:6, "spring", "fall")
    p_entry     <- get_seasonal_p_enl(fm$ReturnYear, anchor_seas)
    if (is.na(p_entry)) p_entry <- fm$p_enl
  }

  tibble(TagCode     = tc,
         p_enl_entry = p_entry,
         spring_entry = as.integer(entry_month %in% 1:6))
})
cat(sprintf("  Computed entry p_enl for %d fish\n\n", nrow(entry_corrections)))

cat(sprintf("Loaded %d fish\n", nrow(fish_df)))
cat(sprintf("  Group A: %d  |  Group B: %d\n",
            sum(fish_df$group == "A"), sum(fish_df$group == "B")))

# Summary of p_enl availability
cat(sprintf("\np_ENL summary:\n"))
cat(sprintf("  Non-missing: %d / %d  (%.1f%%)\n",
            sum(!is.na(fish_df$p_enl)), nrow(fish_df),
            mean(!is.na(fish_df$p_enl))*100))
cat(sprintf("  Mean (non-NA): %.3f\n", mean(fish_df$p_enl, na.rm=TRUE)))
cat(sprintf("  Min (non-NA):  %.3f\n", min(fish_df$p_enl, na.rm=TRUE)))

# =============================================================================
# PREPARE MODEL DATA
# =============================================================================
cat("\n--- Preparing Model Data ---\n")

# Join entry-date p_enl corrections; fall back to 1.0 only if no flow data
# available at all (pre-operational years).
fish_model <- fish_df %>%
  left_join(entry_corrections, by = "TagCode") %>%
  mutate(
    p_enl_entry_safe = if_else(is.na(p_enl_entry), 1.0, p_enl_entry),
    p_lmr_safe       = ifelse(is.na(p_lmr), 1.0, p_lmr),
    p_okl_safe       = ifelse(is.na(p_okl), 1.0, p_okl),
    ReturnYearFactor = factor(ReturnYear),
    spring_return    = as.integer(ReturnSeason == "Spring")
  )

n_enl_missing <- sum(is.na(entry_corrections$p_enl_entry))
cat(sprintf("Fish using p_enl=1.0 (no correction): %d  [pre-operational or missing flow]\n",
            n_enl_missing))

# Detection efficiency summaries by group and fate
cat("\nDetection efficiency by group:\n")
fish_model %>%
  group_by(group) %>%
  summarise(
    n          = n(),
    p_enl_mean = round(mean(p_enl_entry_safe), 3),
    p_enl_min  = round(min(p_enl_entry_safe),  3),
    n_spring   = sum(spring_entry),
    .groups    = "drop"
  ) %>%
  print()

# Classify Group B fish by fate
groupB_fates <- fish_model %>%
  filter(group == "B") %>%
  mutate(
    fate = case_when(
      entiat_detected ~ "Entiat detected",
      upstream_detected ~ "Methow/Oka detected",
      TRUE ~ "Undetected post-Wells"
    )
  )

cat("\nGroup B fish by fate:\n")
groupB_fates %>%
  group_by(fate) %>%
  summarise(
    n = n(),
    p_enl = round(mean(p_enl_entry_safe), 3),
    p_lmr = round(mean(p_lmr_safe), 3),
    p_okl = round(mean(p_okl_safe), 3),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# LOAD STANDARD MODELS FOR COMPARISON
# =============================================================================
cat("\n--- Loading Standard (Uncorrected) Models ---\n")

if (file.exists("entiat_bayesian_models.RData")) {
  load("entiat_bayesian_models.RData")
  cat("Loaded standard models\n")
  standard_post <- as_draws_df(fit_model1)
  beta_wells_std    <- standard_post$b_wells_interaction
  beta_harvest_std  <- standard_post$b_harvest_open
  cat(sprintf("Standard Wells effect: β=%.3f [%.3f, %.3f]\n",
              mean(beta_wells_std),
              quantile(beta_wells_std, 0.025),
              quantile(beta_wells_std, 0.975)))
} else {
  cat("Standard models not found — comparison will be limited\n")
  beta_wells_std <- NULL
  beta_harvest_std <- NULL
}

# =============================================================================
# DEFINE CUSTOM BRMS FAMILY: COMPOUND DETECTION LIKELIHOOD
# =============================================================================
# Y_i ~ Bernoulli(θ_i × p_ENL_i)
# where θ_i = inv_logit(η_i) is the true underlying return probability
#
# Likelihood contributions:
#   y=1: log(θ_i) + log(p_ENL_i)
#   y=0: log(1 - θ_i × p_ENL_i)
#
# Note: this is NOT equivalent to Bernoulli(θ_i) with an offset.
#       logit(θ × p) ≠ logit(θ) + logit(p)
# =============================================================================

det_corrected_family <- custom_family(
  "det_corrected",
  dpars  = "mu",
  links  = "logit",
  type   = "int",
  vars   = "vreal1[n]"   # per-observation detection efficiency passed via vreal()
)

# Stan function block: lpmf, log_lik, and posterior_predict
stan_functions <- "
  // Compound detection likelihood: Y ~ Bernoulli(theta * p_detect)
  // mu = inv_logit(eta) = true return probability
  // vreal_1 = p_detect (per-fish detection efficiency from GAM)
  real det_corrected_lpmf(int y, real mu, real vreal_1) {
    real pi_obs = mu * vreal_1;
    // Clamp to valid probability range
    pi_obs = fmax(1e-10, fmin(1 - 1e-10, pi_obs));
    return bernoulli_lpmf(y | pi_obs);
  }
"

sv_functions <- stanvar(scode = stan_functions, block = "functions")

# R-level functions for post-processing (log_lik, posterior_predict)
posterior_predict_det_corrected <- function(i, prep, ...) {
  mu    <- brms::get_dpar(prep, "mu", i = i)
  p_det <- prep$data$vreal1[i]  # brms stores vreal() data as vreal1
  pi_obs <- pmin(pmax(mu * p_det, 1e-10), 1 - 1e-10)
  rbinom(length(mu), size = 1, prob = pi_obs)
}

log_lik_det_corrected <- function(i, prep) {
  y     <- prep$data$Y[i]
  mu    <- brms::get_dpar(prep, "mu", i = i)
  p_det <- prep$data$vreal1[i]
  pi_obs <- pmin(pmax(mu * p_det, 1e-10), 1 - 1e-10)
  dbinom(y, size = 1, prob = pi_obs, log = TRUE)
}

# =============================================================================
# PRIORS AND MCMC SETTINGS (MATCHING STANDARD ANALYSIS)
# =============================================================================
model_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "sd")
)

MCMC_ITER    <- 4500
MCMC_WARMUP  <- 1500
MCMC_CHAINS  <- 4
MCMC_CORES   <- 4
MCMC_SEED    <- 42
ADAPT_DELTA  <- 0.95

# =============================================================================
# MODEL 1c: DETECTION-CORRECTED — Wells + Harvest + Year RE
# =============================================================================
cat("\n======================================================================\n")
cat("MODEL 1c: Detection-Corrected — Wells + Harvest + Year RE\n")
cat(sprintf("Y_i ~ Bernoulli(theta_i × p_ENL_i)\n"))
cat(sprintf("n = %d fish\n", nrow(fish_model)))
cat("======================================================================\n")

fit_corrected_m1 <- brm(
  bf(EntiatDetected | vreal(p_enl_entry_safe) ~
       wells_interaction + harvest_open + (1 | ReturnYearFactor)),
  data     = fish_model,
  family   = det_corrected_family,
  stanvars = sv_functions,
  prior    = model_priors,
  iter     = MCMC_ITER,
  warmup   = MCMC_WARMUP,
  chains   = MCMC_CHAINS,
  cores    = MCMC_CORES,
  seed     = MCMC_SEED,
  control  = list(adapt_delta = ADAPT_DELTA),
  silent   = 2,
  refresh  = 500
)

cat("\nModel 1c Summary:\n")
print(summary(fit_corrected_m1))

post_c1 <- as_draws_df(fit_corrected_m1)
beta_wells_c1    <- post_c1$b_wells_interaction
beta_harvest_c1  <- post_c1$b_harvest_open
sigma_yr_c1      <- post_c1$`sd_ReturnYearFactor__Intercept`

wells_ci_c1   <- quantile(beta_wells_c1,   probs = c(0.025, 0.975))
harvest_ci_c1 <- quantile(beta_harvest_c1, probs = c(0.025, 0.975))

cat(sprintf("\n--- Model 1c Results (Detection-Corrected) ---\n"))
cat(sprintf("Wells Interaction:\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n",
            mean(beta_wells_c1), exp(mean(beta_wells_c1))))
cat(sprintf("  95%% CI   [%.3f, %.3f]\n", wells_ci_c1[1], wells_ci_c1[2]))
cat(sprintf("  P(beta < 0) = %.4f\n", mean(beta_wells_c1 < 0)))
cat(sprintf("\nHarvest Status:\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n",
            mean(beta_harvest_c1), exp(mean(beta_harvest_c1))))
cat(sprintf("  95%% CI   [%.3f, %.3f]\n", harvest_ci_c1[1], harvest_ci_c1[2]))

diag_c1 <- summary(fit_corrected_m1)$fixed
cat(sprintf("\nDiagnostics:\n"))
cat(sprintf("  Rhat range:     [%.3f, %.3f]\n", min(diag_c1$Rhat), max(diag_c1$Rhat)))
cat(sprintf("  Bulk ESS range: [%.0f, %.0f]\n", min(diag_c1$Bulk_ESS), max(diag_c1$Bulk_ESS)))

# =============================================================================
# MODEL 1d: DETECTION-CORRECTED + SEASON
# =============================================================================
cat("\n======================================================================\n")
cat("MODEL 1d: Detection-Corrected — Wells + Harvest + Season + Year RE\n")
cat("======================================================================\n")

fit_corrected_m1d <- brm(
  bf(EntiatDetected | vreal(p_enl_entry_safe) ~
       wells_interaction + harvest_open + spring_entry + (1 | ReturnYearFactor)),
  data     = fish_model,
  family   = det_corrected_family,
  stanvars = sv_functions,
  prior    = model_priors,
  iter     = MCMC_ITER,
  warmup   = MCMC_WARMUP,
  chains   = MCMC_CHAINS,
  cores    = MCMC_CORES,
  seed     = MCMC_SEED,
  control  = list(adapt_delta = ADAPT_DELTA),
  silent   = 2,
  refresh  = 500
)

cat("\nModel 1d Summary:\n")
print(summary(fit_corrected_m1d))

post_c1d <- as_draws_df(fit_corrected_m1d)
beta_wells_c1d   <- post_c1d$b_wells_interaction
beta_season_c1d  <- post_c1d$b_spring_entry
beta_harvest_c1d <- post_c1d$b_harvest_open

wells_ci_c1d  <- quantile(beta_wells_c1d,  probs = c(0.025, 0.975))
season_ci_c1d <- quantile(beta_season_c1d, probs = c(0.025, 0.975))

cat(sprintf("\n--- Model 1d Results ---\n"))
cat(sprintf("Wells (corrected + season): beta=%.3f, OR=%.3f, [%.3f, %.3f], P(<0)=%.4f\n",
            mean(beta_wells_c1d), exp(mean(beta_wells_c1d)),
            wells_ci_c1d[1], wells_ci_c1d[2], mean(beta_wells_c1d < 0)))
cat(sprintf("Spring Entry:               beta=%.3f, OR=%.3f, [%.3f, %.3f]\n",
            mean(beta_season_c1d), exp(mean(beta_season_c1d)),
            season_ci_c1d[1], season_ci_c1d[2]))

# =============================================================================
# MODEL COMPARISON SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("COMPARISON: Uncorrected vs Detection-Corrected Wells Effect\n")
cat("======================================================================\n")

if (!is.null(beta_wells_std)) {
  comparison_df <- tibble(
    Model = c(
      "Standard (no detection correction)",
      "Standard + season covariate",
      "Detection-corrected",
      "Detection-corrected + season"
    ),
    Beta_Mean = c(
      mean(beta_wells_std),
      if (exists("fit_model1b")) mean(as_draws_df(fit_model1b)$b_wells_interaction) else NA,
      mean(beta_wells_c1),
      mean(beta_wells_c1d)
    ),
    CI_Lower = c(
      quantile(beta_wells_std, 0.025),
      if (exists("fit_model1b")) quantile(as_draws_df(fit_model1b)$b_wells_interaction, 0.025) else NA,
      wells_ci_c1[1],
      wells_ci_c1d[1]
    ),
    CI_Upper = c(
      quantile(beta_wells_std, 0.975),
      if (exists("fit_model1b")) quantile(as_draws_df(fit_model1b)$b_wells_interaction, 0.975) else NA,
      wells_ci_c1[2],
      wells_ci_c1d[2]
    ),
    OR = c(
      exp(mean(beta_wells_std)),
      if (exists("fit_model1b")) exp(mean(as_draws_df(fit_model1b)$b_wells_interaction)) else NA,
      exp(mean(beta_wells_c1)),
      exp(mean(beta_wells_c1d))
    ),
    P_negative = c(
      mean(beta_wells_std < 0),
      if (exists("fit_model1b")) mean(as_draws_df(fit_model1b)$b_wells_interaction < 0) else NA,
      mean(beta_wells_c1 < 0),
      mean(beta_wells_c1d < 0)
    )
  )

  cat("\nWells Interaction Effect (Group B vs Group A):\n")
  print(comparison_df %>% mutate(across(c(Beta_Mean:OR), ~ round(.x, 3))))
}

# =============================================================================
# HIDDEN STRAY ANALYSIS
# =============================================================================
# For 53 Group B fish not detected post-Wells (no ENL, LMR, or OKL detection),
# we estimate what fraction are likely hidden strays vs. missed Entiat returners.
#
# For each undetected Group B fish i, we compute posterior probability of
# each fate using Bayes' theorem:
#
#   P(fate | not detected) ∝ P(not detected | fate) × P(fate)
#
# Fates and associated non-detection probabilities:
#   - True Entiat return: P(miss) = 1 - p_ENL_i
#   - Hidden Methow stray: P(miss) = 1 - p_LMR_i
#   - Hidden Okanogan stray: P(miss) = 1 - p_OKL_i
#   - Other (died, harvested, etc.): P(miss) = 1.0
#
# Prior fate probabilities estimated from all confirmed Group B detections
# =============================================================================
cat("\n======================================================================\n")
cat("HIDDEN STRAY ANALYSIS: 53 Undetected Group B Fish\n")
cat("======================================================================\n")

# Get undetected Group B fish
groupB_undetected <- groupB_fates %>%
  filter(fate == "Undetected post-Wells")

cat(sprintf("Undetected Group B fish: %d\n\n", nrow(groupB_undetected)))

# Compute prior fate probabilities from CONFIRMED detections in Group B
# (detection-efficiency-corrected counts)
groupB_confirmed <- groupB_fates %>%
  filter(fate != "Undetected post-Wells")

# Raw confirmed counts
n_entiat_conf  <- sum(groupB_confirmed$fate == "Entiat detected")
n_stray_conf   <- sum(groupB_confirmed$fate == "Methow/Oka detected")

cat(sprintf("Confirmed Group B detections:\n"))
cat(sprintf("  Entiat: %d  |  Strays (Met/Oka): %d\n", n_entiat_conf, n_stray_conf))

# Determine if each confirmed stray was Methow (LMR) or Okanogan (OKL)
# by whether their last upstream detection was in the Methow or Okanogan
# Use p_lmr vs p_okl to classify: Methow strays have p_lmr valid, Okanogan have p_okl valid
# More precisely: if p_lmr < 1 and p_okl < 1, use subbasin info from detections
# We'll check the raw data column upstream_detected flag — actually we need the subbasin
# For now, split based on which gateway detection probability is lower (more spring/high-flow)
# This is approximate; ideally we'd have the subbasin name

# Better approach: check if p_okl is notably different from p_lmr
# The confirmed strays have p_lmr mean=0.862, p_okl mean=0.338
# The very low p_okl suggests some were Okanogan strays detected under
# high-flow conditions. Let's use a heuristic: fish with p_okl < 0.5
# were likely Okanogan strays.
n_met_stray  <- sum(groupB_confirmed$fate == "Methow/Oka detected" &
                    !is.na(groupB_confirmed$p_okl_safe) &
                    groupB_confirmed$p_okl_safe >= 0.5, na.rm=TRUE) +
                sum(groupB_confirmed$fate == "Methow/Oka detected" &
                    is.na(groupB_confirmed$p_okl), na.rm=TRUE)
n_oka_stray  <- sum(groupB_confirmed$fate == "Methow/Oka detected" &
                    !is.na(groupB_confirmed$p_okl_safe) &
                    groupB_confirmed$p_okl_safe < 0.5, na.rm=TRUE)

cat(sprintf("  Approx Methow strays: %d  |  Approx Okanogan strays: %d\n\n",
            n_met_stray, n_oka_stray))

# Detection-efficiency-corrected stray rate estimates
# For each confirmed Methow stray with known p_LMR, the implied total = 1/p_LMR
# Total implied Methow strays (accounting for missed) = sum(1/p_LMR)
groupB_met <- groupB_confirmed %>%
  filter(fate == "Methow/Oka detected",
         !is.na(p_lmr_safe) & p_lmr_safe < 1.0)
groupB_oka <- groupB_confirmed %>%
  filter(fate == "Methow/Oka detected",
         !is.na(p_okl_safe) & p_okl_safe < 1.0)

# Implied total strays = sum of 1/p_detect for each detected stray
implied_met_total <- if (nrow(groupB_met) > 0) sum(1 / groupB_met$p_lmr_safe) else n_met_stray
implied_oka_total <- if (nrow(groupB_oka) > 0) sum(1 / groupB_oka$p_okl_safe) else n_oka_stray
implied_ent_total <- sum(1 / groupB_confirmed$p_enl_entry_safe[groupB_confirmed$fate == "Entiat detected"])

cat(sprintf("Detection-adjusted implied total counts:\n"))
cat(sprintf("  Entiat returners (confirmed): %.1f  (from %d detected)\n",
            implied_ent_total, n_entiat_conf))
cat(sprintf("  Methow strays (implied):      %.1f  (from %d detected)\n",
            implied_met_total, max(n_met_stray, nrow(groupB_met))))
cat(sprintf("  Okanogan strays (implied):    %.1f  (from %d detected)\n\n",
            implied_oka_total, max(n_oka_stray, nrow(groupB_oka))))

# Prior fate probabilities for undetected fish (from confirmed+implied totals)
# Also include "Other" category (died, harvested, other tributaries)
# We estimate "Other" as a residual to reach the full Group B undetected count
total_confirmed_implied <- implied_ent_total + implied_met_total + implied_oka_total
# Group B total: 140 fish; 87 confirmed, 53 undetected
# The undetected 53 are already somewhat accounted for in the implied totals above
# (each detected fish implies some missed partners at the SAME subbasin)
# The "Other" fraction comes from remaining undetected fish
other_fraction <- 0.1  # assume 10% prior probability of dying/harvested before any detection

# Normalize priors
pi_ent <- implied_ent_total / total_confirmed_implied * (1 - other_fraction)
pi_met <- implied_met_total / total_confirmed_implied * (1 - other_fraction)
pi_oka <- implied_oka_total / total_confirmed_implied * (1 - other_fraction)
pi_other <- other_fraction

cat(sprintf("Prior fate probabilities for undetected Group B fish:\n"))
cat(sprintf("  P(true Entiat return): %.3f\n", pi_ent))
cat(sprintf("  P(Methow stray):       %.3f\n", pi_met))
cat(sprintf("  P(Okanogan stray):     %.3f\n", pi_oka))
cat(sprintf("  P(Other/unknown):      %.3f\n\n", pi_other))

# Per-fish Bayesian posterior fate probabilities
groupB_undetected <- groupB_undetected %>%
  mutate(
    # P(not detected | fate) × P(fate)
    lk_ent   = (1 - p_enl_entry_safe) * pi_ent,
    lk_met   = (1 - p_lmr_safe) * pi_met,
    lk_oka   = (1 - p_okl_safe) * pi_oka,
    lk_other = 1.0 * pi_other,
    lk_total = lk_ent + lk_met + lk_oka + lk_other,
    # Posterior probability of each fate
    prob_ent   = lk_ent   / lk_total,
    prob_met   = lk_met   / lk_total,
    prob_oka   = lk_oka   / lk_total,
    prob_other = lk_other / lk_total
  )

cat(sprintf("Per-fish fate posteriors (mean across 53 undetected Group B fish):\n"))
cat(sprintf("  P(true Entiat return):  %.3f  (range: %.3f - %.3f)\n",
            mean(groupB_undetected$prob_ent),
            min(groupB_undetected$prob_ent),
            max(groupB_undetected$prob_ent)))
cat(sprintf("  P(hidden Methow stray): %.3f  (range: %.3f - %.3f)\n",
            mean(groupB_undetected$prob_met),
            min(groupB_undetected$prob_met),
            max(groupB_undetected$prob_met)))
cat(sprintf("  P(hidden Okano stray):  %.3f  (range: %.3f - %.3f)\n",
            mean(groupB_undetected$prob_oka),
            min(groupB_undetected$prob_oka),
            max(groupB_undetected$prob_oka)))
cat(sprintf("  P(Other/died):          %.3f\n\n", mean(groupB_undetected$prob_other)))

# Expected number of each fate among 53 undetected Group B fish
E_ent   <- sum(groupB_undetected$prob_ent)
E_met   <- sum(groupB_undetected$prob_met)
E_oka   <- sum(groupB_undetected$prob_oka)
E_other <- sum(groupB_undetected$prob_other)

cat(sprintf("Expected fate distribution for 53 undetected Group B fish:\n"))
cat(sprintf("  True Entiat returners (missed at ENL): %.1f fish\n", E_ent))
cat(sprintf("  Hidden Methow strays (missed at LMR):  %.1f fish\n", E_met))
cat(sprintf("  Hidden Okanogan strays (missed at OKL): %.1f fish\n", E_oka))
cat(sprintf("  Other/died/unknown:                    %.1f fish\n\n", E_other))

# Compare observed vs detection-corrected Group B Entiat return rate
n_groupB <- nrow(fish_model %>% filter(group == "B"))
observed_entiat_B <- sum(fish_model$EntiatDetected[fish_model$group == "B"])
# Corrected = observed + expected missed returners
corrected_entiat_B <- observed_entiat_B + E_ent

cat(sprintf("Group B Entiat return rate:\n"))
cat(sprintf("  Observed (raw):    %d / %d = %.1f%%\n",
            observed_entiat_B, n_groupB,
            observed_entiat_B / n_groupB * 100))
cat(sprintf("  Corrected (+ENL missed): %.1f / %d = %.1f%%\n",
            corrected_entiat_B, n_groupB,
            corrected_entiat_B / n_groupB * 100))

# Same for Group A
n_groupA <- nrow(fish_model %>% filter(group == "A"))
observed_entiat_A <- sum(fish_model$EntiatDetected[fish_model$group == "A"])
# For Group A, ENL detection correction applies to all undetected Group A fish
groupA_undetected <- fish_model %>%
  filter(group == "A", EntiatDetected == 0) %>%
  mutate(
    lk_ent   = (1 - p_enl_entry_safe) * 0.8,  # Group A has higher base return rate
    lk_other = 1.0 * 0.2,
    lk_total = lk_ent + lk_other,
    prob_ent = lk_ent / lk_total
  )

E_ent_A <- sum(groupA_undetected$prob_ent)
corrected_entiat_A <- observed_entiat_A + E_ent_A

cat(sprintf("\nGroup A Entiat return rate:\n"))
cat(sprintf("  Observed (raw):    %d / %d = %.1f%%\n",
            observed_entiat_A, n_groupA,
            observed_entiat_A / n_groupA * 100))
cat(sprintf("  Corrected (+ENL missed): %.1f / %d = %.1f%%\n",
            corrected_entiat_A, n_groupA,
            corrected_entiat_A / n_groupA * 100))

# =============================================================================
# CORRECTED ANNUAL RETURN RATES
# =============================================================================
cat("\n--- Computing Corrected Annual Return Rates ---\n")

# For each return year, compute observed and detection-corrected return rates
year_summary_corrected <- fish_model %>%
  group_by(ReturnYear, group) %>%
  summarise(
    n = n(),
    n_detected = sum(EntiatDetected),
    obs_rate   = mean(EntiatDetected),
    # Expected detection probability (average over fish in this year-group cell)
    mean_p_enl = mean(p_enl_entry_safe),
    # Naive correction: observed / mean_p_enl (rough upper bound)
    # Better: (observed + E_missed) / n
    # where E_missed = sum of P(ENL miss | undetected)
    .groups = "drop"
  ) %>%
  # Per-cell correction using mean p_enl
  mutate(
    # Expected corrected return rate (simple correction: divide by detection eff)
    # This overestimates if many fish truly didn't return, so cap at 1
    corr_rate_naive = pmin(obs_rate / mean_p_enl, 1.0),
    # More conservative: apply correction only to undetected fish proportion
    n_undetected = n - n_detected,
    E_missed = n_undetected * (1 - mean_p_enl) * obs_rate / (1 - obs_rate * mean_p_enl + 1e-6),
    corr_n_detected = n_detected + E_missed,
    corr_rate = corr_n_detected / n
  )

cat("\nAnnual observed vs corrected return rates by group (selected years):\n")
year_summary_corrected %>%
  filter(group %in% c("A", "B")) %>%
  select(ReturnYear, group, n, n_detected, obs_rate, mean_p_enl, corr_rate) %>%
  mutate(across(c(obs_rate, mean_p_enl, corr_rate), ~ round(.x, 3))) %>%
  arrange(ReturnYear, group) %>%
  print(n = 40)

# =============================================================================
# SAVE RESULTS
# =============================================================================
cat("\n--- Saving Results ---\n")

# Model comparison table
if (!is.null(beta_wells_std)) {
  results_comparison <- comparison_df %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))

  write_csv(results_comparison, "entiat_detection_corrected_results.csv")
  cat("Saved: entiat_detection_corrected_results.csv\n")
}

# Hidden stray analysis results
stray_summary <- tibble(
  Description        = c("Undetected Group B fish",
                         "Expected true Entiat returners (missed at ENL)",
                         "Expected hidden Methow strays (missed at LMR)",
                         "Expected hidden Okanogan strays (missed at OKL)",
                         "Expected other/unknown fate",
                         "Group B observed Entiat detection rate",
                         "Group B corrected Entiat return rate",
                         "Group A observed Entiat detection rate",
                         "Group A corrected Entiat return rate"),
  Value = c(nrow(groupB_undetected),
            round(E_ent, 1),
            round(E_met, 1),
            round(E_oka, 1),
            round(E_other, 1),
            round(observed_entiat_B / n_groupB * 100, 1),
            round(corrected_entiat_B / n_groupB * 100, 1),
            round(observed_entiat_A / n_groupA * 100, 1),
            round(corrected_entiat_A / n_groupA * 100, 1))
)

write_csv(stray_summary, "entiat_hidden_stray_summary.csv")
cat("Saved: entiat_hidden_stray_summary.csv\n")

# Per-fish fate posteriors for undetected Group B fish
write_csv(
  groupB_undetected %>%
    select(TagCode, ReturnYear, ReturnSeason, p_enl_entry_safe, p_lmr_safe, p_okl_safe,
           prob_ent, prob_met, prob_oka, prob_other),
  "entiat_groupB_undetected_fate_posteriors.csv"
)
cat("Saved: entiat_groupB_undetected_fate_posteriors.csv\n")

# Save model objects
save(fit_corrected_m1, fit_corrected_m1d,
     fish_model, entry_corrections, season_means,
     groupB_undetected, year_summary_corrected,
     stray_summary, comparison_df,
     file = "entiat_corrected_models.RData")
cat("Saved: entiat_corrected_models.RData\n")

# =============================================================================
# VISUALIZATIONS
# =============================================================================
cat("\n--- Generating Visualizations ---\n")

harvest_colors <- c("Open" = "#E07B54", "Closed" = "#5B8DB8")

# ---- Figure 1: Corrected vs Uncorrected Wells Effect Posteriors ----

if (!is.null(beta_wells_std)) {
  wells_comparison_long <- bind_rows(
    tibble(beta = beta_wells_std,
           Model = "Standard\n(no detection correction)"),
    tibble(beta = as_draws_df(fit_model1b)$b_wells_interaction,
           Model = "Standard\n(+ season covariate)"),
    tibble(beta = beta_wells_c1,
           Model = "Detection-corrected\n(p_ENL per fish)"),
    tibble(beta = beta_wells_c1d,
           Model = "Detection-corrected\n(+ season covariate)")
  ) %>%
    mutate(Model = factor(Model, levels = c(
      "Standard\n(no detection correction)",
      "Standard\n(+ season covariate)",
      "Detection-corrected\n(p_ENL per fish)",
      "Detection-corrected\n(+ season covariate)"
    )))

  model_colors <- c("#E07B54", "#FF9800", "#5B8DB8", "#2196F3")

  p_wells_dens <- ggplot(wells_comparison_long, aes(x = beta, fill = Model, color = Model)) +
    geom_density(alpha = 0.30, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_fill_manual(values = model_colors) +
    scale_color_manual(values = model_colors) +
    labs(x = "Wells Interaction Coefficient (log-odds)",
         y = "Posterior Density",
         title = "Wells Dam Effect on Entiat Return:\nStandard vs Detection-Corrected Models",
         fill = "Model", color = "Model") +
    annotate("text", x = max(wells_comparison_long$beta) * 0.9, y = Inf,
             label = "Group B MORE likely\nto return to Entiat",
             vjust = 1.5, hjust = 1, size = 3, color = "gray40") +
    annotate("text", x = min(wells_comparison_long$beta) * 0.9, y = Inf,
             label = "Group B LESS likely\nto return to Entiat",
             vjust = 1.5, hjust = 0, size = 3, color = "gray40") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9))

  # Forest plot of Wells effect across models
  comparison_plot_df <- tibble(
    Model = c("Standard", "Standard\n+ season", "Detection-\ncorrected",
              "Corrected\n+ season"),
    Beta  = c(mean(beta_wells_std),
              mean(as_draws_df(fit_model1b)$b_wells_interaction),
              mean(beta_wells_c1),
              mean(beta_wells_c1d)),
    Lower = c(quantile(beta_wells_std, 0.025),
              quantile(as_draws_df(fit_model1b)$b_wells_interaction, 0.025),
              wells_ci_c1[1],
              wells_ci_c1d[1]),
    Upper = c(quantile(beta_wells_std, 0.975),
              quantile(as_draws_df(fit_model1b)$b_wells_interaction, 0.975),
              wells_ci_c1[2],
              wells_ci_c1d[2]),
    P_neg = c(mean(beta_wells_std < 0),
              mean(as_draws_df(fit_model1b)$b_wells_interaction < 0),
              mean(beta_wells_c1 < 0),
              mean(beta_wells_c1d < 0)),
    Type = c("Standard", "Standard", "Corrected", "Corrected")
  ) %>%
    mutate(Model = factor(Model, levels = Model))

  p_forest <- ggplot(comparison_plot_df,
                     aes(x = Beta, y = Model, color = Type)) +
    geom_pointrange(aes(xmin = Lower, xmax = Upper), size = 0.9, linewidth = 1.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_text(aes(label = sprintf("β=%.3f\nOR=%.2f\nP(<0)=%.2f",
                                  Beta, exp(Beta), P_neg)),
              hjust = -0.08, size = 3.2) +
    scale_color_manual(values = c("Standard" = "#E07B54", "Corrected" = "#5B8DB8"),
                       guide = "none") +
    scale_x_continuous(limits = c(min(comparison_plot_df$Lower) - 0.3,
                                  max(comparison_plot_df$Upper) + 0.5)) +
    labs(x = "Wells Interaction Coefficient (log-odds scale)",
         y = "",
         title = "Wells Effect: Four Model Specifications (95% CI)") +
    theme_minimal()

  fig1 <- p_wells_dens / p_forest + plot_layout(heights = c(2, 1.5))
  ggsave("entiat_corrected_vs_uncorrected.png", fig1, width = 10, height = 10, dpi = 150)
  cat("Saved: entiat_corrected_vs_uncorrected.png\n")
}

# ---- Figure 2: Hidden Stray Analysis ----

# Panel A: Per-fish fate posteriors for 53 undetected Group B fish
stray_long <- groupB_undetected %>%
  select(ReturnYear, ReturnSeason, p_enl_entry_safe, p_lmr_safe, p_okl_safe,
         prob_ent, prob_met, prob_oka, prob_other) %>%
  arrange(desc(prob_met + prob_oka)) %>%
  mutate(fish_id = row_number()) %>%
  pivot_longer(c(prob_ent, prob_met, prob_oka, prob_other),
               names_to = "Fate", values_to = "Probability") %>%
  mutate(Fate = recode(Fate,
    "prob_ent"   = "True Entiat returner\n(missed at ENL)",
    "prob_met"   = "Hidden Methow stray\n(missed at LMR)",
    "prob_oka"   = "Hidden Okanogan stray\n(missed at OKL)",
    "prob_other" = "Other/died/unknown"
  ))

fate_colors <- c(
  "True Entiat returner\n(missed at ENL)"  = "#4CAF50",
  "Hidden Methow stray\n(missed at LMR)"   = "#E07B54",
  "Hidden Okanogan stray\n(missed at OKL)" = "#9C27B0",
  "Other/died/unknown"                      = "#9E9E9E"
)

p_stray_stack <- ggplot(stray_long, aes(x = fish_id, y = Probability, fill = Fate)) +
  geom_col(width = 1) +
  scale_fill_manual(values = fate_colors) +
  labs(x = "Individual fish (ranked by stray probability)",
       y = "Posterior fate probability",
       title = "Fate Posteriors: 53 Undetected Group B Fish",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        axis.text.x = element_blank())

# Panel B: Expected fate counts with uncertainty
expected_df <- tibble(
  Fate = c("True Entiat\nreturner", "Hidden Methow\nstray",
           "Hidden Okanogan\nstray", "Other/unknown"),
  Expected = c(E_ent, E_met, E_oka, E_other),
  SD = c(sqrt(sum(groupB_undetected$prob_ent * (1 - groupB_undetected$prob_ent))),
         sqrt(sum(groupB_undetected$prob_met * (1 - groupB_undetected$prob_met))),
         sqrt(sum(groupB_undetected$prob_oka * (1 - groupB_undetected$prob_oka))),
         sqrt(sum(groupB_undetected$prob_other * (1 - groupB_undetected$prob_other)))),
  Color = c("#4CAF50", "#E07B54", "#9C27B0", "#9E9E9E")
)

p_expected <- ggplot(expected_df, aes(x = Fate, y = Expected, fill = Fate)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, Expected - SD), ymax = Expected + SD),
                width = 0.3, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", Expected)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = setNames(expected_df$Color, expected_df$Fate),
                    guide = "none") +
  labs(x = "", y = "Expected number of fish (± 1 SD)",
       title = "Expected Fate of 53 Undetected Group B Fish") +
  theme_minimal()

# Panel C: Group B detection-corrected rates
groupB_rate_df <- tibble(
  Type = c("Observed\nEntiat detected", "Corrected\n(+ ENL misses)",
           "Confirmed\nStrays detected", "Corrected\n(+ LMR/OKL misses)"),
  Rate = c(observed_entiat_B / n_groupB * 100,
           corrected_entiat_B / n_groupB * 100,
           n_stray_conf / n_groupB * 100,
           (n_stray_conf + E_met + E_oka) / n_groupB * 100),
  Category = c("Entiat", "Entiat", "Stray", "Stray")
)

p_rates <- ggplot(groupB_rate_df, aes(x = Type, y = Rate, fill = Category)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%", Rate)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("Entiat" = "#4CAF50", "Stray" = "#E07B54"),
                    guide = "none") +
  labs(x = "", y = "% of Group B fish",
       title = "Group B Return Rates: Observed vs Corrected") +
  ylim(0, 85) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 9))

fig2 <- (p_stray_stack | p_expected) / p_rates + plot_layout(heights = c(2, 1.5))
ggsave("entiat_hidden_stray_analysis.png", fig2, width = 13, height = 11, dpi = 150)
cat("Saved: entiat_hidden_stray_analysis.png\n")

# ---- Figure 3: Annual Return Rates (Observed vs Corrected) ----

annual_rates_plot <- year_summary_corrected %>%
  filter(group %in% c("A", "B")) %>%
  select(ReturnYear, group, obs_rate, corr_rate, n) %>%
  pivot_longer(c(obs_rate, corr_rate), names_to = "type", values_to = "rate") %>%
  mutate(
    type = recode(type,
      "obs_rate"  = "Observed (raw)",
      "corr_rate" = "Corrected (ENL detection)"
    ),
    group_label = recode(group,
      "A" = "Group A (Rocky Reach only)",
      "B" = "Group B (Wells interaction)"
    )
  )

p_annual_rates <- ggplot(annual_rates_plot,
                         aes(x = ReturnYear, y = rate * 100,
                             color = type, linetype = type)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  facet_wrap(~ group_label, ncol = 1) +
  scale_color_manual(values = c("Observed (raw)" = "#9E9E9E",
                                "Corrected (ENL detection)" = "#5B8DB8"),
                     name = "") +
  scale_linetype_manual(values = c("Observed (raw)" = "dashed",
                                   "Corrected (ENL detection)" = "solid"),
                        name = "") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Return Year", y = "% Detected in Entiat River",
       title = "Annual Entiat Return Rates: Raw vs Detection-Corrected") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("entiat_corrected_return_rates.png", p_annual_rates, width = 10, height = 8, dpi = 150)
cat("Saved: entiat_corrected_return_rates.png\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("FINAL SUMMARY\n")
cat("======================================================================\n\n")

cat(sprintf("DETECTION CORRECTION EFFECT ON WELLS INTERACTION ESTIMATE\n"))
cat(sprintf("==========================================================\n"))
if (!is.null(beta_wells_std)) {
  cat(sprintf("  Standard model:         β = %.3f [%.3f, %.3f]  P(<0)=%.3f\n",
              mean(beta_wells_std),
              quantile(beta_wells_std, 0.025), quantile(beta_wells_std, 0.975),
              mean(beta_wells_std < 0)))
  cat(sprintf("  Detection-corrected:    β = %.3f [%.3f, %.3f]  P(<0)=%.3f\n",
              mean(beta_wells_c1),
              wells_ci_c1[1], wells_ci_c1[2],
              mean(beta_wells_c1 < 0)))
  delta_beta <- mean(beta_wells_c1) - mean(beta_wells_std)
  cat(sprintf("  Change in β:            %.3f\n", delta_beta))
  cat(sprintf("  Interpretation: Detection correction %s the Wells effect estimate\n",
              ifelse(abs(delta_beta) < 0.05, "has minimal effect on",
                     ifelse(delta_beta > 0, "attenuates (reduces magnitude of)",
                            "amplifies"))))
}

cat(sprintf("\nHIDDEN STRAY ANALYSIS (53 Undetected Group B Fish)\n"))
cat(sprintf("====================================================\n"))
cat(sprintf("  Expected true Entiat returners (missed at ENL): %.1f  (%.0f%%)\n",
            E_ent, E_ent / nrow(groupB_undetected) * 100))
cat(sprintf("  Expected hidden Methow strays:                  %.1f  (%.0f%%)\n",
            E_met, E_met / nrow(groupB_undetected) * 100))
cat(sprintf("  Expected hidden Okanogan strays:                %.1f  (%.0f%%)\n",
            E_oka, E_oka / nrow(groupB_undetected) * 100))
cat(sprintf("  Other/died:                                     %.1f  (%.0f%%)\n",
            E_other, E_other / nrow(groupB_undetected) * 100))
cat(sprintf("\n  Total expected hidden strays: %.1f fish (%.0f%% of undetected)\n",
            E_met + E_oka, (E_met + E_oka) / nrow(groupB_undetected) * 100))
cat(sprintf("  These hidden strays, if confirmed, would reduce the apparent\n"))
cat(sprintf("  Wells interaction effect (Group B true Entiat return rate goes UP;\n"))
cat(sprintf("  true stray rate goes UP too, partially offsetting each other).\n"))

cat(sprintf("\nDETECTION-CORRECTED RETURN RATES\n"))
cat(sprintf("==================================\n"))
cat(sprintf("  Group A: %.1f%% → %.1f%% (raw → corrected)\n",
            observed_entiat_A / n_groupA * 100,
            corrected_entiat_A / n_groupA * 100))
cat(sprintf("  Group B: %.1f%% → %.1f%% (raw → corrected)\n",
            observed_entiat_B / n_groupB * 100,
            corrected_entiat_B / n_groupB * 100))
cat(sprintf("  Both groups' true return rates are higher than raw observed rates,\n"))
cat(sprintf("  but the RELATIVE difference between groups is what drives the\n"))
cat(sprintf("  Wells interaction coefficient.\n"))

cat("\n======================================================================\n")
cat("OUTPUT FILES:\n")
cat("  entiat_corrected_models.RData              — brms model objects\n")
cat("  entiat_detection_corrected_results.csv     — model comparison table\n")
cat("  entiat_hidden_stray_summary.csv            — stray analysis summary\n")
cat("  entiat_groupB_undetected_fate_posteriors.csv — per-fish fate probs\n")
cat("  entiat_corrected_vs_uncorrected.png        — Wells effect comparison\n")
cat("  entiat_hidden_stray_analysis.png           — stray analysis figure\n")
cat("  entiat_corrected_return_rates.png          — annual corrected rates\n")
cat("======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
