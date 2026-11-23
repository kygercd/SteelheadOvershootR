###############################################################################
# Temperature Analysis: Rocky Reach Forebay Water Temperature and
#   (1) Wells Dam Overshoot Probability
#   (2) Entiat Wild Steelhead Return Probability (criteria C2 and C3)
#
# Temperature source: USACE CROHMS, Rocky Reach Forebay (RRH)
#   https://public.crohms.org/ftppub/water_quality/tdg/RRH{MonAbbYear}_132.csv
#   Hourly data, aggregated to daily mean here.
#
# Annual temperature metrics (per return year):
#   temp_max_ann       : annual max daily mean temp (°C), calendar year
#   temp_mean_augsep   : mean daily mean temp, Aug–Sep  [available all years]
#   temp_mean_augnov   : mean daily mean temp, Aug–Nov  [NA if <90 days coverage]
#   temp_dd_augnov     : cumulative thermal units (sum daily mean), Aug–Nov  [NA if <90 days]
#   temp_days18_augnov : days with daily mean > 18 °C, Aug–Nov  [NA if <90 days]
#
# Data coverage note: CROHMS data for RRH is complete (122/122 days) from 2012 onward
#   but has gaps in 2006–2011 (34–90 days). Aug–Nov metrics are set to NA for years
#   with < 90 days coverage. Aug–Sep metrics are available for all years.
#   MIN_DAYS_AUGNOV = 90 (threshold below)
#
# Individual fish metrics (matched to rrf_anchor date):
#   arrival_temp       : daily mean temp on anchor date
#   arrival_temp_30d   : mean daily temp ±15 days centred on anchor date
#
# Models (brms, Bernoulli logit, random year intercept):
#   Annual metrics: 4 overshoot models + 4 return (C2) models
#   Individual temp: 2 overshoot models + 2 return (C2) models
#   Interaction: temp × wells group on return (best annual metric, arrival temp)
#   All continuous predictors standardised (mean=0, sd=1).
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(readr)
  library(httr)
  library(brms)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild Temperature")

MCMC_ITER       <- 4500
MCMC_WARMUP     <- 1500
MCMC_CHAINS     <- 4
MCMC_CORES      <- 4
ADAPT_DELTA     <- 0.95
MIN_DAYS_AUGNOV <- 85    # minimum Aug-Nov days to compute aug-nov metrics
                         # 2007 has 87 days (Nov missing from all sources)

# ============================================================
# SECTION 1: LOAD FISH DATA
# ============================================================

fish_model <- read_csv("../Entiat Wild/entiat_analysis_with_detection_eff.csv",
                       show_col_types = FALSE) %>%
  mutate(rrf_anchor = as_datetime(rrf_anchor))

outcomes <- read_csv("../Entiat Wild/entiat_outcome_criteria_summary.csv",
                     show_col_types = FALSE)

fish <- fish_model %>%
  left_join(
    outcomes %>% select(TagCode, criterion_2_orig_strict, criterion_3_reclassified),
    by = "TagCode"
  ) %>%
  mutate(
    wells_group   = as.integer(wells_detected),   # 1 = Group B (went to Wells)
    anchor_date   = as.Date(rrf_anchor),
    anchor_year   = ReturnYear                     # pre-existing return year
  )

cat("Fish loaded: n =", nrow(fish), "\n")
cat("Return years:", min(fish$anchor_year), "–", max(fish$anchor_year), "\n")

# ============================================================
# SECTION 2: DOWNLOAD AND GAP-FILL ROCKY REACH FOREBAY TEMPERATURE
#
# Primary source: DART WQM daily data for RRH (Rocky Reach Forebay)
#   USACE CROHMS sensor via CBR DART mg.php API
#   URL: https://www.cbr.washington.edu/dart/cs/php/rpt/mg.php
#
# Data coverage: RRH WQM is complete (122/122 Aug-Nov days) from 2012 onward
#   but has sensor gaps in 2006-2011 (34-90 days).
#
# Gap-fill surrogates (calibrated against RRH using 2012-2016 Aug-Nov overlap):
#   WEL — Wells Forebay (~30km upstream)      R²=0.973  RMSE=0.44°C
#   PRD — Priest Rapids Forebay (downstream)  R²=0.976  RMSE=0.42°C
#   CHJ — Chief Joseph Forebay (~100km up)    R²=0.941  RMSE=0.52°C
#
# Imputation priority: RRH → WEL → PRD → CHJ → NA
# After imputation all years reach ≥85 days; 2007 = 87 days (Nov unobserved)
# ============================================================

DART_BASE  <- "https://www.cbr.washington.edu/dart/cs/php/rpt/mg.php"
TEMP_CACHE <- "rrh_temperature_daily.csv"
CAL_YEARS  <- 2012:2016   # complete-coverage years used for calibration

return_years <- sort(unique(fish$anchor_year))

# Fetch one full year of daily WQM temperature from DART for a given location
fetch_dart_year <- function(loc, year) {
  url <- paste0(
    DART_BASE,
    "?sc=1&mgconfig=river&outputFormat=csv",
    "&year%5B%5D=", year,
    "&loc%5B%5D=", loc,
    "&data%5B%5D=Temp+%28WQM%29",
    "&startdate=1%2F1&enddate=12%2F31",
    "&avgyear=0&consolidate=1&grid=1"
  )
  resp <- tryCatch(httr::GET(url, httr::timeout(30)), error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)
  txt <- strsplit(httr::content(resp, "text", encoding = "UTF-8"), "\n")[[1]]
  if (length(txt) < 2) return(NULL)
  dat <- tryCatch(
    read.csv(text = paste(txt, collapse = "\n"), stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(dat) || ncol(dat) < 2) return(NULL)
  names(dat)[1:2] <- c("mmdd", "temp_raw")
  dat %>%
    mutate(
      date = suppressWarnings(as.Date(paste0(year, "/", mmdd), "%Y/%m/%d")),
      temp = suppressWarnings(as.numeric(temp_raw))
    ) %>%
    filter(!is.na(date), !is.na(temp), temp >= 0, temp <= 30) %>%
    select(date, temp)
}

# Fetch all years for one station
fetch_station_all <- function(loc, years) {
  cat("  Downloading", loc, "for", length(years), "years ...")
  d <- purrr::map_dfr(years, function(yr) {
    df <- fetch_dart_year(loc, yr)
    if (is.null(df)) tibble(date = as.Date(character()), temp = numeric()) else df
  })
  cat(nrow(d), "days\n")
  d
}

fetch_temperature_data <- function(cache_path = TEMP_CACHE, force = FALSE) {
  if (!force && file.exists(cache_path)) {
    cat("Loading cached temperature data from", cache_path, "\n")
    return(read_csv(cache_path, show_col_types = FALSE) %>%
             mutate(date = as.Date(date)))
  }

  cat("\nDownloading temperature data from DART WQM...\n")

  # Primary: RRH all return years
  rrh_raw <- fetch_station_all("RRH", return_years)

  # Surrogates: calibration years (2012-2016) + gap years (2006-2011)
  sur_years <- sort(union(CAL_YEARS, 2006:2011))
  wel_raw   <- fetch_station_all("WEL", sur_years)
  prd_raw   <- fetch_station_all("PRD", sur_years)
  chj_raw   <- fetch_station_all("CHJ", sur_years)

  # ---- Fit calibration models (Aug-Nov, complete years only) ----
  cat("\nFitting calibration models (", paste(range(CAL_YEARS), collapse="-"), "Aug-Nov):\n")
  rrh_cal <- rrh_raw %>%
    filter(year(date) %in% CAL_YEARS, month(date) %in% 8:11) %>%
    rename(rrh = temp)

  fit_cal_model <- function(sur_df, sur_name) {
    d <- sur_df %>%
      filter(year(date) %in% CAL_YEARS, month(date) %in% 8:11) %>%
      rename(sur = temp) %>%
      inner_join(rrh_cal, by = "date")
    m <- lm(rrh ~ sur, data = d)
    cat(sprintf("  %s: n=%d  RRH = %.4f + %.4f*%s  R²=%.4f  RMSE=%.3f°C\n",
                sur_name, nrow(d), coef(m)[1], coef(m)[2], sur_name,
                summary(m)$r.squared, sqrt(mean(resid(m)^2))))
    m
  }
  m_wel <- fit_cal_model(wel_raw, "WEL")
  m_prd <- fit_cal_model(prd_raw, "PRD")
  m_chj <- fit_cal_model(chj_raw, "CHJ")

  # ---- Build imputed daily series ----
  # Full date grid from June of first year to December of last year
  date_grid <- tibble(date = seq(
    as.Date(paste0(min(return_years), "-06-01")),
    as.Date(paste0(max(return_years), "-12-31")),
    by = "day"
  ))

  daily <- date_grid %>%
    left_join(rrh_raw %>% rename(rrh = temp), by = "date") %>%
    left_join(wel_raw %>% rename(wel = temp), by = "date") %>%
    left_join(prd_raw %>% rename(prd = temp), by = "date") %>%
    left_join(chj_raw %>% rename(chj = temp), by = "date") %>%
    mutate(
      temp_c = rrh,
      source = if_else(!is.na(rrh), "RRH", NA_character_),
      # Priority 2: WEL (closest upstream, R²=0.973)
      temp_c = if_else(is.na(temp_c) & !is.na(wel),
                       coef(m_wel)[1] + coef(m_wel)[2] * wel, temp_c),
      source = if_else(is.na(source) & !is.na(wel), "WEL_cal", source),
      # Priority 3: PRD (downstream, R²=0.976)
      temp_c = if_else(is.na(temp_c) & !is.na(prd),
                       coef(m_prd)[1] + coef(m_prd)[2] * prd, temp_c),
      source = if_else(is.na(source) & !is.na(prd), "PRD_cal", source),
      # Priority 4: CHJ (far upstream, R²=0.941)
      temp_c = if_else(is.na(temp_c) & !is.na(chj),
                       coef(m_chj)[1] + coef(m_chj)[2] * chj, temp_c),
      source = if_else(is.na(source) & !is.na(chj), "CHJ_cal", source),
      source = replace_na(source, "missing")
    ) %>%
    filter(!is.na(temp_c), temp_c >= 0, temp_c <= 30) %>%
    select(date, temp_c, source)

  cat("\nFinal Aug-Nov daily temperature coverage by year and source:\n")
  daily %>%
    filter(month(date) %in% 8:11, year(date) %in% return_years) %>%
    mutate(yr = year(date)) %>%
    group_by(yr, source) %>% summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = source, values_from = n, values_fill = 0L) %>%
    mutate(total = rowSums(across(-yr))) %>%
    print(n = Inf)

  write_csv(daily, cache_path)
  daily
}

daily_temp <- fetch_temperature_data()

# ============================================================
# SECTION 3: ANNUAL TEMPERATURE METRICS
# ============================================================

annual_metrics <- purrr::map_dfr(return_years, function(yr) {
  yr_data <- daily_temp %>% filter(year(date) == yr)
  augsep  <- yr_data %>% filter(month(date) %in% 8:9)
  augnov  <- yr_data %>% filter(month(date) %in% 8:11)
  n_augnov     <- nrow(augnov)
  n_rrh_direct <- sum(augnov$source == "RRH",     na.rm = TRUE)
  n_imputed    <- sum(augnov$source != "RRH",     na.rm = TRUE)

  tibble(
    ReturnYear           = yr,
    n_days_augnov        = n_augnov,
    n_days_rrh_direct    = n_rrh_direct,
    n_days_imputed       = n_imputed,
    temp_max_ann         = if (nrow(yr_data) > 0)
                             max(yr_data$temp_c, na.rm = TRUE) else NA_real_,
    # Aug-Sep: available for all years (no coverage filter)
    temp_mean_augsep     = if (nrow(augsep) >= 20)
                             mean(augsep$temp_c, na.rm = TRUE) else NA_real_,
    # Aug-Nov: only computed when coverage is sufficient
    temp_mean_augnov     = if (n_augnov >= MIN_DAYS_AUGNOV)
                             mean(augnov$temp_c, na.rm = TRUE) else NA_real_,
    temp_dd_augnov       = if (n_augnov >= MIN_DAYS_AUGNOV)
                             sum(augnov$temp_c, na.rm = TRUE) else NA_real_,
    temp_days18_augnov   = if (n_augnov >= MIN_DAYS_AUGNOV)
                             sum(augnov$temp_c > 18, na.rm = TRUE) else NA_real_
  )
})

cat("\nAnnual temperature metrics:\n")
print(annual_metrics, n = Inf)

# Coverage summary
n_full   <- sum(annual_metrics$n_days_augnov >= 122, na.rm = TRUE)
n_thresh <- sum(annual_metrics$n_days_augnov >= MIN_DAYS_AUGNOV, na.rm = TRUE)
cat(sprintf("\nData coverage: %d/%d years have full Aug-Nov coverage (122 days)\n",
            n_full, nrow(annual_metrics)))
cat(sprintf("             %d/%d years meet threshold (%d days) for Aug-Nov metrics\n",
            n_thresh, nrow(annual_metrics), MIN_DAYS_AUGNOV))
low_cov <- annual_metrics %>% filter(n_days_augnov < MIN_DAYS_AUGNOV)
if (nrow(low_cov) > 0) {
  cat("Years with insufficient coverage (Aug-Nov metrics set to NA):",
      paste(low_cov$ReturnYear, collapse = ", "), "\n")
}

write_csv(annual_metrics, "rrh_annual_temp_metrics.csv")

# ============================================================
# SECTION 4: INDIVIDUAL FISH TEMPERATURE METRICS
# ============================================================

fish_base <- fish %>%
  select(TagCode, ReturnYear, anchor_date, wells_detected, wells_group,
         harvest_open, criterion_2_orig_strict, criterion_3_reclassified)

# arrival_temp: vectorized join on anchor date
arrival_direct <- daily_temp %>% select(date, arrival_temp = temp_c)

# arrival_temp_30d: build lookup over unique dates to avoid redundant scans
unique_dates <- sort(unique(fish_base$anchor_date))
rolling_lookup <- purrr::map_dfr(unique_dates, function(d) {
  r <- daily_temp %>%
    filter(date >= d - 15L, date <= d + 15L) %>%
    pull(temp_c)
  tibble(anchor_date = d,
         arrival_temp_30d = if (length(r) < 5) NA_real_ else mean(r, na.rm = TRUE))
})

fish_temp <- fish_base %>%
  left_join(arrival_direct,  by = c("anchor_date" = "date")) %>%
  left_join(rolling_lookup,  by = "anchor_date") %>%
  left_join(annual_metrics,  by = "ReturnYear")

cat("\nIndividual fish temperature metrics (first 6 rows):\n")
print(head(fish_temp %>%
  select(TagCode, ReturnYear, anchor_date, arrival_temp, arrival_temp_30d,
         temp_mean_augnov, temp_days18_augnov)))

cat("\nMissing arrival_temp:", sum(is.na(fish_temp$arrival_temp)), "fish\n")
cat("Missing arrival_temp_30d:", sum(is.na(fish_temp$arrival_temp_30d)), "fish\n")

write_csv(fish_temp, "rrh_fish_temp_metrics.csv")

# ============================================================
# SECTION 5: STANDARDISE PREDICTORS
# ============================================================

# Standardise all continuous temperature predictors
std_metrics <- c("temp_max_ann", "temp_mean_augsep", "temp_mean_augnov",
                 "temp_dd_augnov", "temp_days18_augnov",
                 "arrival_temp", "arrival_temp_30d")

means_sds <- fish_temp %>%
  select(all_of(std_metrics)) %>%
  summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE),
                                      sd   = ~sd(.x,   na.rm = TRUE))))

standardise <- function(x, col) {
  m <- means_sds[[paste0(col, "_mean")]]
  s <- means_sds[[paste0(col, "_sd")]]
  (x - m) / s
}

for (col in std_metrics) {
  fish_temp[[paste0(col, "_std")]] <- standardise(fish_temp[[col]], col)
}

cat("\nStandardisation parameters:\n")
for (col in std_metrics) {
  cat(sprintf("  %-25s mean = %6.2f  sd = %5.2f\n",
              col,
              means_sds[[paste0(col, "_mean")]],
              means_sds[[paste0(col, "_sd")]]))
}

# ============================================================
# SECTION 6: PRIOR AND MODEL-FITTING HELPERS
# ============================================================

# Weakly informative priors consistent with the main analysis
model_priors <- c(
  prior(normal(0, 2.5), class = b),
  prior(normal(0, 2.5), class = Intercept),
  prior(normal(0, 1),   class = sd)
)

fit_temp_model <- function(formula, data, label) {
  cat("\n--- Fitting:", label, "(n =", nrow(data), ") ---\n")
  dir.create("model_cache", showWarnings = FALSE)
  brm(
    formula   = formula,
    data      = data,
    family    = bernoulli(link = "logit"),
    prior     = model_priors,
    iter      = MCMC_ITER,
    warmup    = MCMC_WARMUP,
    chains    = MCMC_CHAINS,
    cores     = MCMC_CORES,
    control   = list(adapt_delta = ADAPT_DELTA),
    seed      = 42,
    silent    = 2,
    file      = file.path("model_cache", label)
  )
}

extract_results <- function(fit, label) {
  draws <- as_draws_df(fit)
  # Find all b_ parameters
  b_params <- names(draws)[startsWith(names(draws), "b_")]
  purrr::map_dfr(b_params, function(p) {
    vals   <- draws[[p]]
    param  <- sub("^b_", "", p)
    tibble(
      model      = label,
      parameter  = param,
      mean_beta  = mean(vals),
      OR         = exp(mean(vals)),
      ci_lo_beta = quantile(vals, 0.025),
      ci_hi_beta = quantile(vals, 0.975),
      ci_lo_OR   = exp(quantile(vals, 0.025)),
      ci_hi_OR   = exp(quantile(vals, 0.975)),
      p_positive = mean(vals > 0),
      p_negative = mean(vals < 0)
    )
  })
}

# ============================================================
# SECTION 7: ANNUAL TEMPERATURE MODELS — OVERSHOOT (wells_detected)
# ============================================================

cat("\n============================================================\n")
cat("SECTION 7: Annual temp × overshoot models\n")
cat("============================================================\n")

# 5 models, one per annual metric; both groups included
# temp_mean_augsep uses all 20 years; aug-nov metrics limited to years with ≥90 days coverage
overshoot_annual_models <- list(
  max_ann     = "temp_max_ann_std",
  mean_augsep = "temp_mean_augsep_std",
  mean_augnov = "temp_mean_augnov_std",
  dd_augnov   = "temp_dd_augnov_std",
  days18      = "temp_days18_augnov_std"
)

overshoot_annual_fits   <- list()
overshoot_annual_results <- tibble()

for (nm in names(overshoot_annual_models)) {
  pred  <- overshoot_annual_models[[nm]]
  label <- paste0("overshoot_annual_", nm)
  fmla  <- as.formula(
    paste("wells_group ~", pred, "+ harvest_open + (1 | ReturnYear)")
  )
  dat   <- fish_temp %>% filter(!is.na(.data[[pred]]))
  fit   <- fit_temp_model(fmla, dat, label)
  overshoot_annual_fits[[nm]]   <- fit
  overshoot_annual_results       <- bind_rows(overshoot_annual_results,
                                               extract_results(fit, label))
}

cat("\nAnnual temperature × overshoot results:\n")
print(overshoot_annual_results %>% arrange(model, parameter))

# ============================================================
# SECTION 8: ANNUAL TEMPERATURE MODELS — RETURN (C2)
# ============================================================

cat("\n============================================================\n")
cat("SECTION 8: Annual temp × return (C2) models\n")
cat("============================================================\n")

# C2 outcome with Wells group + annual temp + harvest
return_annual_fits    <- list()
return_annual_results  <- tibble()

for (nm in names(overshoot_annual_models)) {
  pred  <- overshoot_annual_models[[nm]]
  label <- paste0("return_c2_annual_", nm)
  fmla  <- as.formula(
    paste("criterion_2_orig_strict ~", pred,
          "+ wells_group + harvest_open + (1 | ReturnYear)")
  )
  dat   <- fish_temp %>% filter(!is.na(.data[[pred]]))
  fit   <- fit_temp_model(fmla, dat, label)
  return_annual_fits[[nm]]    <- fit
  return_annual_results        <- bind_rows(return_annual_results,
                                             extract_results(fit, label))
}

cat("\nAnnual temperature × return (C2) results:\n")
print(return_annual_results %>% arrange(model, parameter))

# ============================================================
# SECTION 9: INDIVIDUAL ARRIVAL TEMPERATURE MODELS
# ============================================================

cat("\n============================================================\n")
cat("SECTION 9: Individual arrival temp models\n")
cat("============================================================\n")

ind_metrics <- list(
  arrival      = "arrival_temp_std",
  arrival_30d  = "arrival_temp_30d_std"
)

# — Overshoot —
overshoot_ind_fits    <- list()
overshoot_ind_results  <- tibble()

for (nm in names(ind_metrics)) {
  pred  <- ind_metrics[[nm]]
  label <- paste0("overshoot_ind_", nm)
  fmla  <- as.formula(
    paste("wells_group ~", pred, "+ harvest_open + (1 | ReturnYear)")
  )
  dat   <- fish_temp %>% filter(!is.na(.data[[pred]]))
  fit   <- fit_temp_model(fmla, dat, label)
  overshoot_ind_fits[[nm]]    <- fit
  overshoot_ind_results        <- bind_rows(overshoot_ind_results,
                                             extract_results(fit, label))
}

# — Return (C2) —
return_ind_fits    <- list()
return_ind_results  <- tibble()

for (nm in names(ind_metrics)) {
  pred  <- ind_metrics[[nm]]
  label <- paste0("return_c2_ind_", nm)
  fmla  <- as.formula(
    paste("criterion_2_orig_strict ~", pred,
          "+ wells_group + harvest_open + (1 | ReturnYear)")
  )
  dat   <- fish_temp %>% filter(!is.na(.data[[pred]]))
  fit   <- fit_temp_model(fmla, dat, label)
  return_ind_fits[[nm]]    <- fit
  return_ind_results        <- bind_rows(return_ind_results,
                                          extract_results(fit, label))
}

cat("\nIndividual temp × overshoot results:\n")
print(overshoot_ind_results %>% arrange(model, parameter))
cat("\nIndividual temp × return (C2) results:\n")
print(return_ind_results %>% arrange(model, parameter))

# ============================================================
# SECTION 10: INTERACTION MODELS (temp × wells_group)
# ============================================================

cat("\n============================================================\n")
cat("SECTION 10: Temperature × Wells group interaction models\n")
cat("============================================================\n")

# Use the annual metric with the strongest return association (from Section 8)
# and also arrival_temp as individual metric.

interaction_fits    <- list()
interaction_results  <- tibble()

# Annual interaction: temp_mean_augsep (all years) and temp_mean_augnov × wells_group
for (pred_nm in c("mean_augsep", "mean_augnov", "arrival")) {
  pred  <- switch(pred_nm,
    mean_augsep = "temp_mean_augsep_std",
    mean_augnov = "temp_mean_augnov_std",
    arrival     = "arrival_temp_std"
  )
  label <- paste0("return_c2_interact_", pred_nm)
  fmla  <- as.formula(
    paste("criterion_2_orig_strict ~", pred,
          "* wells_group + harvest_open + (1 | ReturnYear)")
  )
  dat <- fish_temp %>% filter(!is.na(.data[[pred]]))
  fit <- fit_temp_model(fmla, dat, label)
  interaction_fits[[pred_nm]]    <- fit
  interaction_results             <- bind_rows(interaction_results,
                                                extract_results(fit, label))
}

# Also run for C3 (reclassified) with arrival temp
fit_c3_interact <- fit_temp_model(
  formula = criterion_3_reclassified ~ arrival_temp_std * wells_group +
    harvest_open + (1 | ReturnYear),
  data  = fish_temp %>% filter(!is.na(arrival_temp_std)),
  label = "return_c3_interact_arrival"
)
interaction_results <- bind_rows(
  interaction_results,
  extract_results(fit_c3_interact, "return_c3_interact_arrival")
)
interaction_fits[["c3_arrival"]] <- fit_c3_interact

cat("\nInteraction model results:\n")
print(interaction_results %>% arrange(model, parameter))

# ============================================================
# SECTION 11: SAVE ALL MODELS
# ============================================================

save(
  overshoot_annual_fits,
  return_annual_fits,
  overshoot_ind_fits,
  return_ind_fits,
  interaction_fits,
  file = "entiat_temperature_models.RData"
)
cat("\nAll models saved to entiat_temperature_models.RData\n")

# ============================================================
# SECTION 12: COMPILE AND SAVE RESULTS
# ============================================================

all_results <- bind_rows(
  overshoot_annual_results %>% mutate(outcome = "overshoot", pred_type = "annual"),
  return_annual_results    %>% mutate(outcome = "return_c2", pred_type = "annual"),
  overshoot_ind_results    %>% mutate(outcome = "overshoot", pred_type = "individual"),
  return_ind_results       %>% mutate(outcome = "return_c2", pred_type = "individual"),
  interaction_results      %>% mutate(outcome = "interaction", pred_type = "mixed")
) %>%
  select(outcome, pred_type, model, parameter, mean_beta, OR,
         ci_lo_beta, ci_hi_beta, ci_lo_OR, ci_hi_OR, p_positive, p_negative) %>%
  mutate(across(c(mean_beta, OR, ci_lo_beta, ci_hi_beta,
                  ci_lo_OR, ci_hi_OR), ~round(.x, 4)),
         across(c(p_positive, p_negative), ~round(.x, 4)))

write_csv(all_results, "entiat_temperature_results.csv")

cat("\n=== FINAL RESULTS SUMMARY ===\n")
cat("\n--- Temperature coefficients (temp predictors only) ---\n")
temp_params <- c("temp_max_ann_std", "temp_mean_augsep_std",
                 "temp_mean_augnov_std", "temp_dd_augnov_std",
                 "temp_days18_augnov_std", "arrival_temp_std", "arrival_temp_30d_std",
                 "temp_mean_augsep_std:wells_group",
                 "temp_mean_augnov_std:wells_group", "arrival_temp_std:wells_group")
all_results %>%
  filter(parameter %in% temp_params) %>%
  select(outcome, model, parameter, mean_beta, OR, ci_lo_OR, ci_hi_OR,
         p_positive, p_negative) %>%
  print(n = Inf)

# ============================================================
# SECTION 13: PLOTS
# ============================================================

# ---- Plot 1: Annual temperature metrics over time ----------------------------
annual_plot_df <- annual_metrics %>%
  pivot_longer(
    cols      = c(temp_max_ann, temp_mean_augsep, temp_mean_augnov, temp_days18_augnov),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(metric = recode(metric,
    temp_max_ann       = "Annual max (°C)",
    temp_mean_augsep   = "Mean Aug–Sep (°C, all years)",
    temp_mean_augnov   = "Mean Aug–Nov (°C, 2012+)",
    temp_days18_augnov = "Days >18°C Aug–Nov (2012+)"
  ))

p1 <- ggplot(annual_plot_df, aes(x = ReturnYear, y = value)) +
  geom_line(colour = "#E07B54", linewidth = 0.8) +
  geom_point(colour = "#E07B54", size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = seq(2006, 2025, by = 2)) +
  labs(
    title = "Rocky Reach Forebay Temperature by Return Year",
    x     = "Return year",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )
ggsave("fig_temp_annual_metrics.png", p1, width = 7, height = 8, dpi = 180)

# ---- Plot 2: Annual overshoot rate vs. Aug–Nov mean temp ---------------------
year_summary <- fish_temp %>%
  group_by(ReturnYear) %>%
  summarise(
    n              = n(),
    overshoot_pct  = mean(wells_group) * 100,
    return_pct     = mean(criterion_2_orig_strict, na.rm = TRUE) * 100,
    temp_augsep    = first(temp_mean_augsep),
    temp_augnov    = first(temp_mean_augnov),
    days18         = first(temp_days18_augnov),
    .groups        = "drop"
  )

# Use Aug-Sep mean (all years) for scatter plots; shade incomplete-coverage points
p2 <- ggplot(year_summary, aes(x = temp_augsep, y = overshoot_pct, label = ReturnYear)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#4A90D9", fill = "#4A90D9",
              alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(size = n), colour = "#E07B54", alpha = 0.8) +
  geom_text(nudge_y = 1.5, size = 3, colour = "grey30") +
  scale_size_continuous(name = "n fish", range = c(2, 8)) +
  labs(
    title = "Annual overshoot rate vs. mean Aug–Sep temperature (all years)",
    x     = "Mean Aug–Sep forebay temperature (°C)",
    y     = "Wells Dam overshoot rate (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("fig_overshoot_vs_temp.png", p2, width = 6.5, height = 5, dpi = 180)

# ---- Plot 3: Annual C2 return rate vs. Aug–Sep mean temp --------------------
p3 <- ggplot(year_summary, aes(x = temp_augsep, y = return_pct, label = ReturnYear)) +
  geom_smooth(method = "lm", se = TRUE, colour = "#4A90D9", fill = "#4A90D9",
              alpha = 0.15, linewidth = 0.8) +
  geom_point(aes(size = n), colour = "#5B8A5E", alpha = 0.8) +
  geom_text(nudge_y = 1.5, size = 3, colour = "grey30") +
  scale_size_continuous(name = "n fish", range = c(2, 8)) +
  labs(
    title = "Annual return rate (C2) vs. mean Aug–Sep temperature (all years)",
    x     = "Mean Aug–Sep forebay temperature (°C)",
    y     = "C2 return rate (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("fig_return_vs_temp.png", p3, width = 6.5, height = 5, dpi = 180)

# ---- Plot 4: Coefficient plot for temperature effects -----------------------
coef_plot_df <- all_results %>%
  filter(parameter %in% temp_params,
         !grepl("interact", model)) %>%
  mutate(
    metric = case_when(
      grepl("max_ann",     model) ~ "Annual max",
      grepl("mean_augsep", model) ~ "Mean Aug–Sep (all yrs)",
      grepl("mean_augnov", model) ~ "Mean Aug–Nov (2012+)",
      grepl("dd_augnov",   model) ~ "Degree-days Aug–Nov (2012+)",
      grepl("days18",      model) ~ "Days >18°C Aug–Nov (2012+)",
      grepl("30d",         parameter) ~ "Arrival temp ±15d",
      grepl("arrival",     model) ~ "Arrival temp",
      TRUE                        ~ parameter
    ),
    outcome_label = recode(outcome,
      overshoot  = "Overshoot probability",
      return_c2  = "Return probability (C2)"
    )
  )

p4 <- ggplot(coef_plot_df,
             aes(x = mean_beta, y = reorder(metric, mean_beta),
                 colour = outcome_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = ci_lo_beta, xmax = ci_hi_beta),
                 linewidth = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    values = c("Overshoot probability" = "#E07B54",
               "Return probability (C2)" = "#5B8A5E"),
    name = NULL
  ) +
  labs(
    title = "Temperature effect estimates (standardised β, 95% CI)",
    x     = "β (log-odds per SD)",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave("fig_temp_coefs.png", p4, width = 7, height = 5, dpi = 180)

cat("\nAnalysis complete. Outputs saved to Entiat Wild Temperature/\n")
cat("  rrh_annual_temp_metrics.csv\n")
cat("  rrh_fish_temp_metrics.csv\n")
cat("  entiat_temperature_results.csv\n")
cat("  entiat_temperature_models.RData\n")
cat("  fig_temp_annual_metrics.png\n")
cat("  fig_overshoot_vs_temp.png\n")
cat("  fig_return_vs_temp.png\n")
cat("  fig_temp_coefs.png\n")
