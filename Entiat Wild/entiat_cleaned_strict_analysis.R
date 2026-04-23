###############################################################################
# Cleaned strict-outcome analysis
# Removes tag-reuse / biologically implausible fish, then fits:
#   S1c: EntiatSpawned ~ wells + harvest + (1|year)          [standard]
#   S2c: detection-corrected compound likelihood              [det-corr]
#   S3c: wells * actual_entry_season interaction              [season]
#
# Strict outcome definition:
#   Any detection at ≥1 upstream Entiat site within the 18-month anchor
#   window constitutes success. ENL detection alone (no upstream) = 0.
#   Upstream detection implies passage through ENL regardless of whether
#   the fish was detected there.
#
# Detection correction:
#   Fish with ENL detection: p_enl computed at the actual ENL entry date
#     and return year (not the RRF anchor date).
#   Fish with upstream-only detection (no ENL): p_enl is the mean for the
#     entry season (spring = Jan–Jun, fall = Jul–Dec) of that return year.
#   Both Group A and Group B use actual Entiat entry season.
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
         last_dt      = mdy_hm(`Last Time Value`),
         release_date = mdy(`Release Date MMDDYYYY`))

fish_model <- read_csv("entiat_analysis_with_detection_eff.csv",
                       show_col_types = FALSE) %>%
  mutate(rrf_anchor = as_datetime(rrf_anchor))

load("gateway_bayesian_models.RData")   # fit_enl, flow_enl, entiat

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

# ---- ENL detection efficiency helpers ----------------------------------------
enl_logq_mean <- mean(entiat$log_q)
enl_logq_sd   <- sd(entiat$log_q)
enl_base_year <- min(entiat$year)
flow_enl      <- flow_enl %>% mutate(date = as.Date(date))

# p_enl at a specific date and return year
predict_p_enl_at_date <- function(entry_date, return_year) {
  q <- flow_enl %>% filter(date == as.Date(entry_date)) %>% pull(discharge)
  if (length(q) == 0 || is.na(q)) return(NA_real_)
  lq_s <- (log(q) - enl_logq_mean) / enl_logq_sd
  yr   <- return_year - enl_base_year
  nd   <- data.frame(log_q_s = lq_s, log_q2_s = lq_s^2, year = yr)
  mean(posterior_epred(fit_enl, newdata = nd,
                       allow_new_levels = TRUE, sample_new_levels = "gaussian"))
}

# Pre-compute year x season mean p_enl (for upstream-only fish)
# Predicts at the mean log-discharge for each year x season combination
cat("Computing year x season mean p_enl lookup...\n")
season_means <- flow_enl %>%
  filter(!is.na(discharge)) %>%
  mutate(yr       = year(date),
         season   = if_else(month(date) %in% 1:6, "spring", "fall"),
         lq_s     = (log(discharge) - enl_logq_mean) / enl_logq_sd,
         yr_c     = yr - enl_base_year) %>%
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

# Fallback: overall season mean across all years
overall_spring_p <- mean(season_means$mean_p_enl[season_means$season == "spring"])
overall_fall_p   <- mean(season_means$mean_p_enl[season_means$season == "fall"])

get_seasonal_p_enl <- function(return_year, season) {
  val <- season_means %>%
    filter(yr == return_year, season == !!season) %>%
    pull(mean_p_enl)
  if (length(val) == 0 || is.na(val)) {
    return(if (season == "spring") overall_spring_p else overall_fall_p)
  }
  val
}

cat(sprintf("  Overall spring mean p_enl: %.3f\n", overall_spring_p))
cat(sprintf("  Overall fall mean p_enl:   %.3f\n", overall_fall_p))

# ---- Identify bad fish (tag reuse / implausible) ----------------------------
# These three fish are manually verified as legitimate adult returns:
#   3D9.1C2C451B01 — valid 2009 return; 2017 ENS detection is a data error
#   3D9.1C2CDD70AD — RRF-missed 2011 adult return; manually re-anchored below
#   3D9.1C2D45A06D — RRF-missed 2012 adult return; manually re-anchored below
verified_fish <- c("3D9.1C2C451B01", "3D9.1C2CDD70AD", "3D9.1C2D45A06D")

bad_tags <- map_dfr(fish_model$TagCode, function(tc) {
  if (tc %in% verified_fish) return(NULL)
  fm       <- fish_model %>% filter(TagCode == tc)
  dets     <- raw_main  %>% filter(tag == tc)
  anchor   <- as_datetime(fm$rrf_anchor)
  win_end  <- anchor %m+% months(18)

  ent_in  <- dets %>% filter(site %in% entiat_sites,
                              first_dt > anchor, first_dt <= win_end)
  ent_out <- dets %>% filter(site %in% entiat_sites, first_dt > win_end)

  if (nrow(ent_in) > 0 || nrow(ent_out) == 0) return(NULL)

  rel_date        <- dets %>% pull(release_date) %>% first()
  late_ent_dt     <- min(ent_out$first_dt, na.rm = TRUE)
  age_at_anchor   <- as.numeric(difftime(anchor, as_datetime(rel_date),
                                          units = "days")) / 365.25
  age_at_late_ent <- as.numeric(difftime(late_ent_dt, as_datetime(rel_date),
                                          units = "days")) / 365.25

  is_bad <- age_at_anchor < 1 | age_at_late_ent > 9
  if (!is_bad) return(NULL)
  tibble(TagCode = tc, age_at_anchor = round(age_at_anchor, 2),
         age_at_late_ent = round(age_at_late_ent, 2))
})

cat(sprintf("Bad fish identified: %d (will be excluded)\n", nrow(bad_tags)))
cat(sprintf("Remaining: %d of %d\n\n",
            nrow(fish_model) - nrow(bad_tags), nrow(fish_model)))

# ---- Build cleaned dataset with strict outcome ------------------------------
fish_clean <- map_dfr(
  fish_model %>% filter(!TagCode %in% bad_tags$TagCode) %>% pull(TagCode),
  function(tc) {
    fm      <- fish_model %>% filter(TagCode == tc)
    dets    <- raw_main   %>% filter(tag == tc) %>% arrange(first_dt)
    anchor  <- as_datetime(fm$rrf_anchor)
    win_end <- anchor %m+% months(18)

    enl_in <- dets %>% filter(site == "ENL - Lower Entiat River",
                               first_dt > anchor, first_dt <= win_end)
    up_in  <- dets %>% filter(site %in% entiat_upstream,
                               first_dt > anchor, first_dt <= win_end)

    has_enl      <- nrow(enl_in) > 0
    has_upstream <- nrow(up_in)  > 0

    # Outcome: upstream detection = success; ENL-only = failure
    spawned <- as.integer(has_upstream)

    # p_enl correction: use actual ENL entry date, or seasonal mean if no ENL
    if (has_enl) {
      enl_date    <- min(enl_in$first_dt)
      p_enl_corr  <- predict_p_enl_at_date(enl_date, fm$ReturnYear)
      if (is.na(p_enl_corr)) p_enl_corr <- fm$p_enl
      entry_month <- month(enl_date)
    } else if (has_upstream) {
      up_date     <- min(up_in$first_dt)
      entry_month <- month(up_date)
      entry_seas  <- if_else(entry_month %in% 1:6, "spring", "fall")
      p_enl_corr  <- get_seasonal_p_enl(fm$ReturnYear, entry_seas)
    } else {
      # Never detected in Entiat — fallback (spawned=0, p_enl irrelevant)
      entry_month <- month(anchor)
      p_enl_corr  <- fm$p_enl
    }

    fm %>% mutate(
      EntiatSpawned  = spawned,
      p_enl_corrected = p_enl_corr,
      spring_entry    = as.integer(entry_month %in% 1:6),
      p_enl_safe      = pmax(0.01, pmin(0.99, p_enl_corr))
    )
  }) %>%
  mutate(ReturnYearFactor = factor(ReturnYear))

# ---- Manual patch: RRF-missed adults ----------------------------------------
# Two fish had confirmed adult returns (BON→MCN→PRA→RIA→Entiat) but were not
# detected at Rocky Reach adult fishway. Re-anchor to PRA detection.
# 3D9.1C2CDD70AD: 2011-09-17 PRA; ENL 2011-09-26 + upstream spring 2012
# 3D9.1C2D45A06D: 2012-08-29 PRA; MAD/TLT/ENA spring 2013 + RRJ kelt

# p_enl corrections for patched fish
p_enl_70AD <- predict_p_enl_at_date("2011-09-26", 2011)   # actual ENL date
if (is.na(p_enl_70AD)) p_enl_70AD <- get_seasonal_p_enl(2011, "fall")
p_enl_A06D <- get_seasonal_p_enl(2013, "spring")           # upstream-only, spring

fish_clean <- fish_clean %>%
  mutate(
    rrf_anchor     = case_when(
      TagCode == "3D9.1C2CDD70AD" ~ as_datetime("2011-09-17 16:23:00"),
      TagCode == "3D9.1C2D45A06D" ~ as_datetime("2012-08-29 07:54:00"),
      TRUE ~ rrf_anchor),
    ReturnYear     = case_when(
      TagCode == "3D9.1C2CDD70AD" ~ 2011,
      TagCode == "3D9.1C2D45A06D" ~ 2012,
      TRUE ~ ReturnYear),
    ReturnSeason   = case_when(
      TagCode %in% c("3D9.1C2CDD70AD", "3D9.1C2D45A06D") ~ "Fall",
      TRUE ~ ReturnSeason),
    EntiatSpawned  = case_when(
      TagCode %in% c("3D9.1C2CDD70AD", "3D9.1C2D45A06D") ~ 1L,
      TRUE ~ EntiatSpawned),
    p_enl_corrected = case_when(
      TagCode == "3D9.1C2CDD70AD" ~ p_enl_70AD,
      TagCode == "3D9.1C2D45A06D" ~ p_enl_A06D,
      TRUE ~ p_enl_corrected),
    spring_entry   = case_when(
      TagCode == "3D9.1C2CDD70AD" ~ 0L,   # fall ENL entry (Sep)
      TagCode == "3D9.1C2D45A06D" ~ 1L,   # spring upstream entry (Apr)
      TRUE ~ spring_entry),
    p_enl_safe     = pmax(0.01, pmin(0.99, p_enl_corrected))
  ) %>%
  mutate(ReturnYearFactor = factor(ReturnYear),
         spring_return    = spring_entry)   # align spring_return with actual entry

cat("Manual patch applied: 2 fish re-anchored to PRA with EntiatSpawned=1\n")

# ---- Descriptive summary ----------------------------------------------------
cat("======================================================================\n")
cat("CLEANED + STRICT OUTCOME ANALYSIS\n")
cat(sprintf("n = %d fish (%d excluded as tag-reuse/implausible)\n",
            nrow(fish_clean), nrow(bad_tags)))
cat("======================================================================\n\n")

cat("Group x actual entry season breakdown:\n")
fish_clean %>%
  mutate(EntrySeasonLabel = if_else(spring_entry == 1, "Spring", "Fall")) %>%
  group_by(group, EntrySeasonLabel) %>%
  summarise(n = n(), spawned = sum(EntiatSpawned),
            pct = round(mean(EntiatSpawned)*100, 1), .groups = "drop") %>%
  print()

cat("\nOverall by group:\n")
fish_clean %>%
  group_by(group) %>%
  summarise(n = n(), spawned = sum(EntiatSpawned),
            pct = round(mean(EntiatSpawned)*100, 1), .groups = "drop") %>%
  print()

cat("\nMean p_enl_corrected by group and entry season:\n")
fish_clean %>%
  mutate(EntrySeasonLabel = if_else(spring_entry == 1, "Spring", "Fall")) %>%
  group_by(group, EntrySeasonLabel) %>%
  summarise(mean_p = round(mean(p_enl_corrected, na.rm=TRUE), 3),
            .groups = "drop") %>%
  print()

# ---- MCMC priors -------------------------------------------------------------
model_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

# ---- Model S1c: Standard -----------------------------------------------------
cat("\n======================================================================\n")
cat("MODEL S1c: EntiatSpawned ~ wells + harvest + (1|year)\n")
cat("======================================================================\n")

fit_s1c <- brm(
  EntiatSpawned ~ wells_interaction + harvest_open + (1 | ReturnYearFactor),
  data    = fish_clean, family = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s1c <- as_draws_df(fit_s1c)
b_s1c    <- post_s1c$b_wells_interaction
ci_s1c   <- quantile(b_s1c, c(0.025, 0.975))
cat(sprintf("\nWells: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_s1c), exp(mean(b_s1c)), ci_s1c[1], ci_s1c[2], mean(b_s1c<0)))

# ---- Model S2c: Detection-corrected ------------------------------------------
cat("\n======================================================================\n")
cat("MODEL S2c: detection-corrected, p_enl at actual entry date/season\n")
cat("======================================================================\n")

det_corrected_family <- custom_family(
  "det_corrected", dpars="mu", links="logit", type="int", vars="vreal1[n]")
stan_functions <- "
  real det_corrected_lpmf(int y, real mu, real vreal1) {
    real pi_obs = mu * vreal1;
    pi_obs = fmax(1e-10, fmin(1 - 1e-10, pi_obs));
    return bernoulli_lpmf(y | pi_obs);
  }
"
sv <- stanvar(scode = stan_functions, block = "functions")

fit_s2c <- brm(
  bf(EntiatSpawned | vreal(p_enl_safe) ~
       wells_interaction + harvest_open + (1 | ReturnYearFactor)),
  data    = fish_clean %>% filter(!is.na(p_enl_safe)),
  family  = det_corrected_family, stanvars = sv,
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s2c <- as_draws_df(fit_s2c)
b_s2c    <- post_s2c$b_wells_interaction
ci_s2c   <- quantile(b_s2c, c(0.025, 0.975))
cat(sprintf("\nWells: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f  n=%d\n",
            mean(b_s2c), exp(mean(b_s2c)), ci_s2c[1], ci_s2c[2], mean(b_s2c<0),
            nrow(fit_s2c$data)))

# ---- Model S3c: Actual entry season interaction ------------------------------
cat("\n======================================================================\n")
cat("MODEL S3c: wells * spring_entry interaction (actual Entiat entry season)\n")
cat("======================================================================\n")

fit_s3c <- brm(
  EntiatSpawned ~ wells_interaction * spring_entry + harvest_open
                + (1 | ReturnYearFactor),
  data    = fish_clean, family = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER, warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores  = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

post_s3c    <- as_draws_df(fit_s3c)
b_fall_s3c  <- post_s3c$b_wells_interaction
b_inter_s3c <- post_s3c$`b_wells_interaction:spring_entry`
b_spr_s3c   <- b_fall_s3c + b_inter_s3c
ci_fall   <- quantile(b_fall_s3c,  c(0.025, 0.975))
ci_spr    <- quantile(b_spr_s3c,   c(0.025, 0.975))
ci_inter  <- quantile(b_inter_s3c, c(0.025, 0.975))

cat(sprintf("\nWells FALL:   beta=%.3f, OR=%.3f, CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_fall_s3c),  exp(mean(b_fall_s3c)),  ci_fall[1],  ci_fall[2],
            mean(b_fall_s3c < 0)))
cat(sprintf("Wells SPRING: beta=%.3f, OR=%.3f, CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_spr_s3c),   exp(mean(b_spr_s3c)),   ci_spr[1],   ci_spr[2],
            mean(b_spr_s3c < 0)))
cat(sprintf("Interaction:  beta=%.3f, CI [%.3f, %.3f], P(>0)=%.4f\n",
            mean(b_inter_s3c), ci_inter[1], ci_inter[2],
            mean(b_inter_s3c > 0)))

# ---- Comparison table --------------------------------------------------------
cat("\n======================================================================\n")
cat("COMPARISON TABLE\n")
cat("======================================================================\n")

make_row <- function(label, b_vec, n = NULL) {
  ci <- quantile(b_vec, c(0.025, 0.975))
  tibble(Model = label, n = n, Beta = mean(b_vec), OR = exp(mean(b_vec)),
         CI_low = ci[1], CI_high = ci[2], P_neg = mean(b_vec < 0))
}

comparison <- bind_rows(
  make_row("Strict (original) — standard",  b_s1c, nrow(fit_s1c$data)),
  make_row("Strict (original) — det-corr",  b_s2c, nrow(fit_s2c$data))
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))

print(comparison)

# ---- Save --------------------------------------------------------------------
save(fit_s1c, fit_s2c, fit_s3c, fish_clean, bad_tags, comparison, season_means,
     file = "entiat_cleaned_strict_models.RData")
cat("\nSaved: entiat_cleaned_strict_models.RData\n")
write_csv(comparison, "entiat_cleaned_strict_comparison.csv")
cat("Saved: entiat_cleaned_strict_comparison.csv\n")

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
