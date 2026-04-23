###############################################################################
# Reclassified strict-outcome analysis
#
# Extends the cleaned strict analysis by reclassifying ALL fish based on
# their ACTUAL Entiat entry season rather than their RRF anchor season.
# Both Group A and Group B fish have fall RRF anchors, but many first enter
# the Entiat in spring after overwintering. Those fish are reclassified as
# spring returners and scored with the spring criterion (ENL within 18 months
# is sufficient), and their p_ENL is recomputed at the actual spring entry
# date's flow rather than the fall anchor date's flow (spring flows are
# typically higher → lower detection efficiency).
#
# Fall-entry fish of both groups retain the strict fall criterion
# (ENL + ≥1 upstream site within 18 months).
#
# Models fit:
#   S1r: EntiatSpawned_r ~ wells + harvest + (1|year)         [standard]
#   S2r: detection-corrected compound likelihood               [det-corr]
#   S3r: wells * spring_return_r interaction                  [season]
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(brms)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

MCMC_ITER   <- 4500
MCMC_WARMUP <- 1500
MCMC_CHAINS <- 4
MCMC_CORES  <- 4
ADAPT_DELTA <- 0.95

# ---- Load data ---------------------------------------------------------------
raw_main <- read_csv("EntiatWildSteel.csv", show_col_types = FALSE) %>%
  rename(tag = `Tag Code`, site = `Site Name`) %>%
  mutate(first_dt     = mdy_hm(`First Time Value`),
         release_date = mdy(`Release Date MMDDYYYY`))

load("entiat_cleaned_strict_models.RData")   # fish_clean, bad_tags
load("gateway_bayesian_models.RData")        # fit_enl, flow_enl, entiat

entiat_sites <- c("ENL - Lower Entiat River","ENM - Middle Entiat River",
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
entiat_upstream <- setdiff(entiat_sites, "ENL - Lower Entiat River")

# ---- ENL detection efficiency scaling parameters ----------------------------
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

# =============================================================================
# RECLASSIFY ALL FISH BY ACTUAL ENTIAT ENTRY SEASON
# =============================================================================
cat("Reclassifying all fish by actual Entiat entry season...\n")

reclassify_results <- map_dfr(fish_clean$TagCode, function(tc) {
  fm      <- fish_clean %>% filter(TagCode == tc)
  dets    <- raw_main  %>% filter(tag == tc) %>% arrange(first_dt)
  anchor  <- as_datetime(fm$rrf_anchor)
  win_end <- anchor %m+% months(18)

  # Find first post-anchor Entiat detection (any Entiat site)
  first_ent <- dets %>%
    filter(site %in% entiat_sites, first_dt > anchor) %>%
    slice_min(first_dt, n = 1)

  if (nrow(first_ent) == 0) {
    # Never detected in Entiat — fall, outcome 0
    return(tibble(TagCode = tc,
                  spring_return_r   = 0L,
                  EntiatSpawned_r   = 0L,
                  p_enl_r           = fm$p_enl_corrected,
                  entiat_entry_date = as.Date(NA),
                  reclassified      = FALSE))
  }

  entry_mon       <- month(first_ent$first_dt)
  is_spring_entry <- entry_mon %in% 1:6

  if (!is_spring_entry) {
    # Fall entry
    # Group B: any Entiat detection = successful return — these fish already
    #   navigated back past Wells Dam, so ENL detection alone is sufficient.
    # Group A: strict fall criterion (ENL + ≥1 upstream site required)
    spawned_fall <- if (fm$group == "B") 1L else fm$EntiatSpawned
    return(tibble(TagCode = tc,
                  spring_return_r   = 0L,
                  EntiatSpawned_r   = spawned_fall,
                  p_enl_r           = fm$p_enl_corrected,
                  entiat_entry_date = as.Date(first_ent$first_dt),
                  reclassified      = FALSE))
  }

  # Spring entry: apply spring criterion (ENL within window is sufficient)
  enl_in    <- dets %>% filter(site == "ENL - Lower Entiat River",
                                first_dt > anchor, first_dt <= win_end)
  spawned_r <- as.integer(nrow(enl_in) > 0)

  # Recompute p_enl at actual spring entry date
  entry_date <- as.Date(first_ent$first_dt)
  p_enl_new  <- predict_p_enl_at_date(entry_date, fm$ReturnYear)
  if (is.na(p_enl_new)) p_enl_new <- fm$p_enl_corrected   # fallback to entry-date corrected value

  tibble(TagCode = tc,
         spring_return_r   = 1L,
         EntiatSpawned_r   = spawned_r,
         p_enl_r           = p_enl_new,
         entiat_entry_date = entry_date,
         reclassified      = TRUE)
})

cat(sprintf("Fish reclassified as spring entry: %d total (Group A: %d, Group B: %d)\n",
            sum(reclassify_results$reclassified),
            sum(reclassify_results$reclassified &
                  fish_clean$group[match(reclassify_results$TagCode, fish_clean$TagCode)] == "A"),
            sum(reclassify_results$reclassified &
                  fish_clean$group[match(reclassify_results$TagCode, fish_clean$TagCode)] == "B")))

# Summarise outcome change by group
cat("\nOutcome change from reclassification:\n")
fish_clean %>%
  left_join(reclassify_results %>% select(TagCode, EntiatSpawned_r, reclassified),
            by = "TagCode") %>%
  count(group, EntiatSpawned, EntiatSpawned_r) %>%
  mutate(change = case_when(
    EntiatSpawned == EntiatSpawned_r ~ "unchanged",
    EntiatSpawned == 0 & EntiatSpawned_r == 1 ~ "0 -> 1 (gained)",
    TRUE ~ "1 -> 0 (lost)")) %>%
  print()

# Merge reclassified columns into fish_clean
fish_r <- fish_clean %>%
  left_join(reclassify_results %>%
              select(TagCode, spring_return_r, EntiatSpawned_r, p_enl_r,
                     reclassified),
            by = "TagCode") %>%
  mutate(
    p_enl_r_safe      = pmax(0.01, pmin(0.99, p_enl_r)),
    ReturnYearFactor  = factor(ReturnYear)
  )

# ---- Descriptive summary -----------------------------------------------------
cat("\n======================================================================\n")
cat("RECLASSIFIED STRICT OUTCOME\n")
cat(sprintf("n = %d fish\n", nrow(fish_r)))
cat("======================================================================\n\n")

cat("Group × reclassified season breakdown:\n")
fish_r %>%
  mutate(SeasonR = if_else(spring_return_r == 1, "Spring", "Fall")) %>%
  group_by(group, SeasonR) %>%
  summarise(n = n(), spawned = sum(EntiatSpawned_r),
            pct = round(mean(EntiatSpawned_r) * 100, 1), .groups = "drop") %>%
  print()

cat("\nOverall by group:\n")
fish_r %>%
  group_by(group) %>%
  summarise(n = n(), spawned = sum(EntiatSpawned_r),
            pct = round(mean(EntiatSpawned_r) * 100, 1), .groups = "drop") %>%
  print()

cat("\np_enl comparison for reclassified Group B fish:\n")
fish_r %>% filter(group == "B", reclassified) %>%
  summarise(mean_p_enl_orig  = round(mean(p_enl,   na.rm = TRUE), 3),
            mean_p_enl_spring = round(mean(p_enl_r, na.rm = TRUE), 3)) %>%
  print()

# ---- MCMC priors -------------------------------------------------------------
model_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

# ---- Model S1r: Standard -----------------------------------------------------
cat("\n======================================================================\n")
cat("MODEL S1r: EntiatSpawned_r ~ wells + harvest + (1|year)\n")
cat("======================================================================\n")

fit_s1r <- brm(
  EntiatSpawned_r ~ wells_interaction + harvest_open + (1 | ReturnYearFactor),
  data    = fish_r, family = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s1r <- as_draws_df(fit_s1r)
b_s1r    <- post_s1r$b_wells_interaction
ci_s1r   <- quantile(b_s1r, c(0.025, 0.975))
cat(sprintf("\nWells: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_s1r), exp(mean(b_s1r)), ci_s1r[1], ci_s1r[2], mean(b_s1r < 0)))

# ---- Model S2r: Detection-corrected ------------------------------------------
cat("\n======================================================================\n")
cat("MODEL S2r: detection-corrected, p_enl at actual entry date\n")
cat("======================================================================\n")

det_corrected_family <- custom_family(
  "det_corrected", dpars = "mu", links = "logit", type = "int",
  vars = "vreal1[n]")
stan_functions <- "
  real det_corrected_lpmf(int y, real mu, real vreal1) {
    real pi_obs = mu * vreal1;
    pi_obs = fmax(1e-10, fmin(1 - 1e-10, pi_obs));
    return bernoulli_lpmf(y | pi_obs);
  }
"
sv <- stanvar(scode = stan_functions, block = "functions")

fit_s2r <- brm(
  bf(EntiatSpawned_r | vreal(p_enl_r_safe) ~
       wells_interaction + harvest_open + (1 | ReturnYearFactor)),
  data    = fish_r %>% filter(!is.na(p_enl_r_safe)),
  family  = det_corrected_family, stanvars = sv,
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s2r <- as_draws_df(fit_s2r)
b_s2r    <- post_s2r$b_wells_interaction
ci_s2r   <- quantile(b_s2r, c(0.025, 0.975))
cat(sprintf("\nWells: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f  n=%d\n",
            mean(b_s2r), exp(mean(b_s2r)), ci_s2r[1], ci_s2r[2], mean(b_s2r < 0),
            nrow(fit_s2r$data)))

# ---- Model S3r: Season interaction -------------------------------------------
cat("\n======================================================================\n")
cat("MODEL S3r: wells * spring_return_r interaction\n")
cat("======================================================================\n")

fit_s3r <- brm(
  EntiatSpawned_r ~ wells_interaction * spring_return_r + harvest_open
                  + (1 | ReturnYearFactor),
  data    = fish_r, family = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s3r    <- as_draws_df(fit_s3r)
b_fall_s3r  <- post_s3r$b_wells_interaction
b_inter_s3r <- post_s3r$`b_wells_interaction:spring_return_r`
b_spr_s3r   <- b_fall_s3r + b_inter_s3r
ci_fall_r   <- quantile(b_fall_s3r, c(0.025, 0.975))
ci_spr_r    <- quantile(b_spr_s3r,  c(0.025, 0.975))
ci_inter_r  <- quantile(b_inter_s3r, c(0.025, 0.975))

cat(sprintf("\nWells FALL:   beta=%.3f, OR=%.3f, CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_fall_s3r),  exp(mean(b_fall_s3r)),  ci_fall_r[1],  ci_fall_r[2],
            mean(b_fall_s3r < 0)))
cat(sprintf("Wells SPRING: beta=%.3f, OR=%.3f, CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_spr_s3r),   exp(mean(b_spr_s3r)),   ci_spr_r[1],   ci_spr_r[2],
            mean(b_spr_s3r < 0)))
cat(sprintf("Interaction:  beta=%.3f, CI [%.3f, %.3f], P(>0)=%.4f\n",
            mean(b_inter_s3r), ci_inter_r[1], ci_inter_r[2],
            mean(b_inter_s3r > 0)))

# ---- Comparison table --------------------------------------------------------
cat("\n======================================================================\n")
cat("COMPARISON: Original strict vs. Reclassified strict\n")
cat("======================================================================\n")

load("entiat_cleaned_strict_models.RData")

make_row <- function(label, b_vec, n = NULL) {
  ci <- quantile(b_vec, c(0.025, 0.975))
  tibble(Model = label, n = n, Beta = mean(b_vec), OR = exp(mean(b_vec)),
         CI_low = ci[1], CI_high = ci[2], P_neg = mean(b_vec < 0))
}

comparison_r <- bind_rows(
  make_row("Strict (original) — standard",    as_draws_df(fit_s1c)$b_wells_interaction, nrow(fit_s1c$data)),
  make_row("Strict (original) — det-corr",    as_draws_df(fit_s2c)$b_wells_interaction, nrow(fit_s2c$data)),
  make_row("Strict (reclassified) — standard", b_s1r, nrow(fit_s1r$data)),
  make_row("Strict (reclassified) — det-corr", b_s2r, nrow(fit_s2r$data))
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

print(comparison_r)

# ---- Save --------------------------------------------------------------------
save(fit_s1r, fit_s2r, fit_s3r, fish_r, reclassify_results, comparison_r,
     file = "entiat_reclassified_strict_models.RData")
cat("\nSaved: entiat_reclassified_strict_models.RData\n")

write_csv(comparison_r, "entiat_reclassified_strict_comparison.csv")
cat("Saved: entiat_reclassified_strict_comparison.csv\n")

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
