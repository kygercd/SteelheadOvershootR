###############################################################################
# Rebuild entiat_analysis_with_detection_eff.csv
#
# Starts fresh from entiat_analysis_all_fish.csv (authoritative source for
# group assignments and fish attributes), carries over existing p_enl/p_lmr/
# p_okl values where available, and predicts efficiencies for any new fish.
# Run after entiat_return_bayesian_analysis.R.
###############################################################################

library(tidyverse)
library(lubridate)
library(brms)

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

load("gateway_bayesian_models.RData")

allf <- read_csv("entiat_analysis_all_fish.csv", show_col_types = FALSE)
eff  <- read_csv("entiat_analysis_with_detection_eff.csv", show_col_types = FALSE)

# Detection efficiency columns to carry over from existing file
eff_cols <- c("enl_q", "met_q", "oka_q", "return_date", "return_doy",
              "return_year", "year_c_enl", "year_c_lmr", "year_c_okl",
              "p_enl", "p_lmr", "p_okl")

# Start from allf — correct group assignments and attributes for all 404 fish
# Join on existing efficiency values where available
existing_eff <- eff %>% select(TagCode, any_of(eff_cols))
result <- allf %>% left_join(existing_eff, by = "TagCode")

new_tags <- result$TagCode[is.na(result$p_enl) & is.na(result$enl_q)]
cat(sprintf("Fish needing detection efficiency computed: %d / %d\n",
            length(new_tags), nrow(result)))

if (length(new_tags) == 0) {
  cat("All fish have existing detection efficiency values.\n")
  write_csv(result, "entiat_analysis_with_detection_eff.csv")
  cat(sprintf("Saved: entiat_analysis_with_detection_eff.csv (%d fish)\n", nrow(result)))
  quit(save = "no")
}

# =============================================================================
# SCALING PARAMETERS
# =============================================================================
enl_logq_mean <- mean(entiat$log_q);  enl_logq_sd <- sd(entiat$log_q)
met_logq_mean <- mean(methow$log_q);  met_logq_sd <- sd(methow$log_q)
oka_logq_mean <- mean(okanogan$log_q); oka_logq_sd <- sd(okanogan$log_q)
enl_base_year <- min(entiat$year)
scale_logq    <- function(lq, mu, s) (lq - mu) / s

# =============================================================================
# FLOW LOOKUP FOR NEW FISH
# =============================================================================
flow_enl <- flow_enl %>% mutate(date = as.Date(date))
flow_met <- flow_met %>% mutate(date = as.Date(date))
flow_oka <- flow_oka %>% mutate(date = as.Date(date))

new_fish <- result %>%
  filter(TagCode %in% new_tags) %>%
  mutate(anchor_date = as.Date(rrf_anchor)) %>%
  left_join(flow_enl %>% rename(enl_q_new = discharge), by = c("anchor_date" = "date")) %>%
  left_join(flow_met %>% rename(met_q_new = discharge), by = c("anchor_date" = "date")) %>%
  left_join(flow_oka %>% rename(oka_q_new = discharge), by = c("anchor_date" = "date")) %>%
  mutate(
    enl_q       = enl_q_new,
    met_q       = met_q_new,
    oka_q       = oka_q_new,
    return_date = anchor_date,
    return_doy  = yday(anchor_date),
    return_year = ReturnYear,
    year_c_enl  = ReturnYear - enl_base_year,
    year_c_lmr  = ReturnYear - min(methow$year),
    year_c_okl  = ReturnYear - min(okanogan$year)
  ) %>%
  select(-anchor_date, -enl_q_new, -met_q_new, -oka_q_new)

cat("\nFlow lookup for new fish:\n")
print(new_fish %>% select(TagCode, return_year, return_doy, enl_q, met_q, oka_q))

# =============================================================================
# PREDICT DETECTION EFFICIENCIES
# =============================================================================
predict_p_enl <- function(fish_row) {
  q <- fish_row$enl_q;  yr <- fish_row$year_c_enl
  if (is.na(q) || is.na(yr)) return(NA_real_)
  lq_s <- scale_logq(log(q), enl_logq_mean, enl_logq_sd)
  nd <- data.frame(log_q_s = lq_s, log_q2_s = lq_s^2, year = yr)
  mean(posterior_epred(fit_enl, newdata = nd,
                       allow_new_levels = TRUE, sample_new_levels = "gaussian"))
}
predict_p_lmr <- function(fish_row) {
  q <- fish_row$met_q;  yr <- fish_row$year_c_lmr
  if (is.na(q) || is.na(yr)) return(NA_real_)
  lq_s <- scale_logq(log(q), met_logq_mean, met_logq_sd)
  nd <- data.frame(log_q_s = lq_s, year = yr)
  mean(posterior_epred(fit_lmr, newdata = nd,
                       allow_new_levels = TRUE, sample_new_levels = "gaussian"))
}
predict_p_okl <- function(fish_row) {
  q <- fish_row$oka_q;  yr <- fish_row$year_c_okl
  if (is.na(q) || is.na(yr)) return(NA_real_)
  lq_s <- scale_logq(log(q), oka_logq_mean, oka_logq_sd)
  nd <- data.frame(log_q_s = lq_s, year = yr)
  mean(posterior_epred(fit_okl, newdata = nd,
                       allow_new_levels = TRUE, sample_new_levels = "gaussian"))
}

cat("\nPredicting detection efficiencies...\n")
new_fish <- new_fish %>%
  rowwise() %>%
  mutate(
    p_enl = predict_p_enl(cur_data()),
    p_lmr = predict_p_lmr(cur_data()),
    p_okl = predict_p_okl(cur_data())
  ) %>%
  ungroup()

cat("\nPredicted values:\n")
print(new_fish %>% select(TagCode, return_year, enl_q, p_enl, p_lmr, p_okl))

# =============================================================================
# MERGE BACK AND SAVE
# =============================================================================
result <- result %>%
  rows_update(new_fish %>% select(TagCode, any_of(eff_cols)), by = "TagCode")

cat(sprintf("\nFinal detection_eff: %d fish\n", nrow(result)))
cat(sprintf("  Group A: %d  |  Group B: %d\n",
            sum(result$group == "A"), sum(result$group == "B")))

write_csv(result, "entiat_analysis_with_detection_eff.csv")
cat("Saved: entiat_analysis_with_detection_eff.csv\n")
