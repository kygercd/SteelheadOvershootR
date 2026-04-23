###############################################################################
# Strict outcome analysis: EntiatSpawned
#   Fall returners: ENL + upstream Entiat, both in 18-month window
#   Spring returners: ENL alone in 18-month window
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(brms)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

MCMC_ITER    <- 4500
MCMC_WARMUP  <- 1500
MCMC_CHAINS  <- 4
MCMC_CORES   <- 4
ADAPT_DELTA  <- 0.95

raw <- read_csv("EntiatWildSteel.csv", show_col_types = FALSE) %>%
  rename(tag = `Tag Code`, site = `Site Name`) %>%
  mutate(first_dt = mdy_hm(`First Time Value`), last_dt = mdy_hm(`Last Time Value`))

fish_model <- read_csv("entiat_analysis_with_detection_eff.csv", show_col_types = FALSE) %>%
  mutate(rrf_anchor = as_datetime(rrf_anchor))

entiat_upstream <- c("ENM - Middle Entiat River",
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

# Build new outcome
fish <- map_dfr(fish_model$TagCode, function(tc) {
  fm      <- fish_model %>% filter(TagCode == tc)
  dets    <- raw %>% filter(tag == tc) %>% arrange(first_dt)
  anchor  <- fm$rrf_anchor
  win_end <- anchor + months(18)

  enl_in <- dets %>% filter(site == "ENL - Lower Entiat River",
                              first_dt > anchor, first_dt <= win_end)
  up_in  <- dets %>% filter(site %in% entiat_upstream,
                              first_dt > anchor, first_dt <= win_end)

  has_enl      <- nrow(enl_in) > 0
  has_upstream <- nrow(up_in)  > 0
  is_spring    <- fm$ReturnSeason == "Spring"

  fm %>% mutate(
    EntiatSpawned = as.integer(if_else(is_spring, has_enl, has_enl & has_upstream)),
    has_enl_window      = has_enl,
    has_upstream_window = has_upstream
  )
}) %>%
  mutate(
    ReturnYearFactor = factor(ReturnYear),
    spring_return    = as.integer(ReturnSeason == "Spring"),
    p_enl_safe       = pmax(0.01, pmin(0.99, p_enl))
  )

# ---- Descriptive summary ----
cat("======================================================================\n")
cat("STRICT OUTCOME: EntiatSpawned\n")
cat("Fall = ENL + upstream (in-window); Spring = ENL only (in-window)\n")
cat("======================================================================\n\n")

cat("Overall:\n")
fish %>% group_by(group) %>%
  summarise(n = n(), orig = sum(EntiatDetected), orig_pct = mean(EntiatDetected)*100,
            spawned = sum(EntiatSpawned), spawned_pct = mean(EntiatSpawned)*100,
            dropped = orig - spawned, .groups = "drop") %>%
  print()

cat("\nBy group × season:\n")
fish %>% group_by(group, ReturnSeason) %>%
  summarise(n = n(), spawned = sum(EntiatSpawned), pct = mean(EntiatSpawned)*100,
            .groups = "drop") %>%
  print()

cat("\nRaw log-odds differences:\n")
for (seas in c("Fall", "Spring")) {
  la <- fish %>% filter(group=="A", ReturnSeason==seas) %>% pull(EntiatSpawned)
  lb <- fish %>% filter(group=="B", ReturnSeason==seas) %>% pull(EntiatSpawned)
  diff_lo <- qlogis(mean(lb)) - qlogis(mean(la))
  cat(sprintf("  %s: A=%.1f%% (n=%d), B=%.1f%% (n=%d), diff=%.1f pp, log-odds diff=%.3f\n",
              seas, mean(la)*100, length(la), mean(lb)*100, length(lb),
              (mean(lb)-mean(la))*100, diff_lo))
}

# ---- Model priors ----
model_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

# ---- Model S1: Standard — new outcome ----
cat("\n======================================================================\n")
cat("MODEL S1: Standard Bayesian — EntiatSpawned ~ wells + harvest + (1|year)\n")
cat("======================================================================\n")

fit_s1 <- brm(
  EntiatSpawned ~ wells_interaction + harvest_open + (1 | ReturnYearFactor),
  data     = fish,
  family   = bernoulli(link = "logit"),
  prior    = model_priors,
  iter     = MCMC_ITER, warmup = MCMC_WARMUP,
  chains   = MCMC_CHAINS, cores  = MCMC_CORES,
  seed     = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent   = 2, refresh = 500
)

print(summary(fit_s1))
post_s1 <- as_draws_df(fit_s1)
b_wells_s1 <- post_s1$b_wells_interaction
ci_s1 <- quantile(b_wells_s1, c(0.025, 0.975))
cat(sprintf("\nWells effect: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_wells_s1), exp(mean(b_wells_s1)), ci_s1[1], ci_s1[2], mean(b_wells_s1<0)))

# ---- Model S2: Detection-corrected — new outcome ----
cat("\n======================================================================\n")
cat("MODEL S2: Detection-corrected — EntiatSpawned | vreal(p_enl)\n")
cat("======================================================================\n")

det_corrected_family <- custom_family(
  "det_corrected",
  dpars  = "mu",
  links  = "logit",
  type   = "int",
  vars   = "vreal1[n]"
)
stan_functions <- "
  real det_corrected_lpmf(int y, real mu, real vreal1) {
    real pi_obs = mu * vreal1;
    pi_obs = fmax(1e-10, fmin(1 - 1e-10, pi_obs));
    return bernoulli_lpmf(y | pi_obs);
  }
"
sv <- stanvar(scode = stan_functions, block = "functions")

fit_s2 <- brm(
  bf(EntiatSpawned | vreal(p_enl_safe) ~
       wells_interaction + harvest_open + (1 | ReturnYearFactor)),
  data     = fish,
  family   = det_corrected_family,
  stanvars = sv,
  prior    = model_priors,
  iter     = MCMC_ITER, warmup = MCMC_WARMUP,
  chains   = MCMC_CHAINS, cores  = MCMC_CORES,
  seed     = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent   = 2, refresh = 500
)

print(summary(fit_s2))
post_s2 <- as_draws_df(fit_s2)
b_wells_s2 <- post_s2$b_wells_interaction
ci_s2 <- quantile(b_wells_s2, c(0.025, 0.975))
cat(sprintf("\nWells effect: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_wells_s2), exp(mean(b_wells_s2)), ci_s2[1], ci_s2[2], mean(b_wells_s2<0)))

# ---- Model S3: Season interaction — new outcome ----
cat("\n======================================================================\n")
cat("MODEL S3: Season interaction — EntiatSpawned ~ wells * season + harvest\n")
cat("======================================================================\n")

fit_s3 <- brm(
  EntiatSpawned ~ wells_interaction * spring_return + harvest_open + (1 | ReturnYearFactor),
  data     = fish,
  family   = bernoulli(link = "logit"),
  prior    = model_priors,
  iter     = MCMC_ITER, warmup = MCMC_WARMUP,
  chains   = MCMC_CHAINS, cores  = MCMC_CORES,
  seed     = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent   = 2, refresh = 500
)

print(summary(fit_s3))
post_s3 <- as_draws_df(fit_s3)
b_fall_s3   <- post_s3$b_wells_interaction
b_inter_s3  <- post_s3$`b_wells_interaction:spring_return`
b_spring_s3 <- b_fall_s3 + b_inter_s3
ci_fall_s3   <- quantile(b_fall_s3,   c(0.025, 0.975))
ci_spring_s3 <- quantile(b_spring_s3, c(0.025, 0.975))
ci_inter_s3  <- quantile(b_inter_s3,  c(0.025, 0.975))
cat(sprintf("\nWells FALL:   beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_fall_s3),   exp(mean(b_fall_s3)),   ci_fall_s3[1],   ci_fall_s3[2],   mean(b_fall_s3<0)))
cat(sprintf("Wells SPRING: beta=%.3f, OR=%.3f, 95%% CI [%.3f, %.3f], P(<0)=%.4f\n",
            mean(b_spring_s3), exp(mean(b_spring_s3)), ci_spring_s3[1], ci_spring_s3[2], mean(b_spring_s3<0)))
cat(sprintf("Interaction:  beta=%.3f, 95%% CI [%.3f, %.3f], P(>0)=%.4f\n",
            mean(b_inter_s3), ci_inter_s3[1], ci_inter_s3[2], mean(b_inter_s3>0)))

# ---- Comparison table ----
cat("\n======================================================================\n")
cat("MODEL COMPARISON: Original EntiatDetected vs Strict EntiatSpawned\n")
cat("======================================================================\n")

load("entiat_corrected_models.RData")  # fit_corrected_m1 (original)
post_orig <- as_draws_df(fit_corrected_m1)
b_orig    <- post_orig$b_wells_interaction
ci_orig   <- quantile(b_orig, c(0.025, 0.975))

comparison <- tribble(
  ~Outcome,         ~Model,                    ~Beta,            ~OR,                   ~CI_low,       ~CI_high,       ~P_neg,
  "Original (ENL)", "Standard (Model 1)",      -0.642,           0.526,                 -1.09,          -0.21,          0.997,
  "Original (ENL)", "Det-corrected (Model 1c)", mean(b_orig),    exp(mean(b_orig)),      ci_orig[1],    ci_orig[2],     mean(b_orig<0),
  "Strict (spawned)","Standard (S1)",           mean(b_wells_s1),exp(mean(b_wells_s1)), ci_s1[1],      ci_s1[2],       mean(b_wells_s1<0),
  "Strict (spawned)","Det-corrected (S2)",      mean(b_wells_s2),exp(mean(b_wells_s2)), ci_s2[1],      ci_s2[2],       mean(b_wells_s2<0)
)

print(comparison %>% mutate(across(where(is.numeric), ~ round(.x, 3))))

# ---- Save ----
save(fit_s1, fit_s2, fit_s3, fish,
     file = "entiat_strict_outcome_models.RData")
cat("\nSaved: entiat_strict_outcome_models.RData\n")

write_csv(comparison, "entiat_strict_outcome_comparison.csv")
cat("Saved: entiat_strict_outcome_comparison.csv\n")

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
