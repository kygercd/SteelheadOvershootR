###############################################################################
# Sensitivity Analysis: Exclude Sparse Years (2022-2024, n=1 each)
# Companion to: run_full_hierarchical_analysis_with_harvest.R
#
# Refits all 6 hierarchical models on data excluding SpillYears 2022-2024
# (each n=1, all ReturnedDownstream=1) to assess whether these sparse,
# perfectly-separated years bias the harvest or spill effect estimates.
###############################################################################

library(brms)
library(tidyverse)
library(posterior)

set.seed(42)
setwd("/home/chas/SteelheadOvershootR")

cat("======================================================================\n")
cat("SENSITIVITY ANALYSIS: EXCLUDING 2022-2024 (SPARSE n=1 YEARS)\n")
cat("======================================================================\n\n")

analysis_df <- read_csv("analysis_fish_with_harvest_R.csv", show_col_types = FALSE)

cat(sprintf("Full dataset: %d fish, %d years\n", nrow(analysis_df),
            length(unique(analysis_df$SpillYear))))

excl_df <- analysis_df %>% filter(!SpillYear %in% c(2022, 2023, 2024))

cat(sprintf("Excluding 2022-2024: %d fish, %d years\n", nrow(excl_df),
            length(unique(excl_df$SpillYear))))
cat(sprintf("Removed: %d fish (all ReturnedDownstream=1)\n",
            nrow(analysis_df) - nrow(excl_df)))

excl_df$SpillYearFactor <- factor(excl_df$SpillYear)

model_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "sd")
)

spill_metrics <- c("FallSpillHours", "FallSpillDays", "FallTotalSpill",
                    "SpringSpillHours", "SpringSpillDays", "SpringTotalSpill")

run_model <- function(df, spill_metric, base_fit = NULL) {
  spill_values <- df[[spill_metric]]
  spill_mean <- mean(spill_values, na.rm = TRUE)
  spill_sd <- sd(spill_values, na.rm = TRUE)
  if (spill_sd == 0) spill_sd <- 1
  df$spill_std <- (spill_values - spill_mean) / spill_sd

  if (is.null(base_fit)) {
    fit <- brm(
      ReturnedDownstream ~ spill_std + HarvestExposure + (1 | SpillYearFactor),
      data = df,
      family = bernoulli(link = "logit"),
      prior = model_priors,
      iter = 4500, warmup = 1500, chains = 4, cores = 4, seed = 42,
      control = list(adapt_delta = 0.95),
      silent = 2, refresh = 0
    )
  } else {
    fit <- update(base_fit, newdata = df, recompile = FALSE,
                   iter = 4500, warmup = 1500, chains = 4, cores = 4, seed = 42,
                   control = list(adapt_delta = 0.95),
                   silent = 2, refresh = 0)
  }
  fit
}

results <- list()
base_fit <- NULL
for (m in spill_metrics) {
  cat(sprintf("\nFitting (excl. 2022-2024): %s...\n", m))
  fit <- run_model(excl_df, m, base_fit = base_fit)
  if (is.null(base_fit)) base_fit <- fit

  post <- as_draws_df(fit)
  beta_spill <- post$b_spill_std
  beta_harvest <- post$b_HarvestExposure
  sigma_alpha <- post$sd_SpillYearFactor__Intercept

  spill_hdi <- quantile(beta_spill, c(0.025, 0.975))
  harvest_hdi <- quantile(beta_harvest, c(0.025, 0.975))
  if (requireNamespace("bayestestR", quietly = TRUE)) {
    s <- bayestestR::hdi(beta_spill, ci = 0.95)
    spill_hdi <- c(s$CI_low, s$CI_high)
    h <- bayestestR::hdi(beta_harvest, ci = 0.95)
    harvest_hdi <- c(h$CI_low, h$CI_high)
  }

  results[[m]] <- tibble(
    Metric = m,
    Beta_Spill_Mean = mean(beta_spill),
    Spill_HDI_Lower = spill_hdi[1],
    Spill_HDI_Upper = spill_hdi[2],
    P_Spill_Positive = mean(beta_spill > 0),
    Beta_Harvest_Mean = mean(beta_harvest),
    Harvest_HDI_Lower = harvest_hdi[1],
    Harvest_HDI_Upper = harvest_hdi[2],
    Harvest_OR = exp(mean(beta_harvest)),
    Sigma_Alpha = mean(sigma_alpha),
    Rhat_max = max(rstan::summary(fit$fit)$summary[, "Rhat"], na.rm = TRUE)
  )
  cat(sprintf("  beta_spill = %.3f [%.3f, %.3f]; beta_harvest = %.3f [%.3f, %.3f]\n",
              mean(beta_spill), spill_hdi[1], spill_hdi[2],
              mean(beta_harvest), harvest_hdi[1], harvest_hdi[2]))
}

results_df <- bind_rows(results)
write_csv(results_df, "sensitivity_exclude_sparse_years_results.csv")
cat("\nSaved: sensitivity_exclude_sparse_years_results.csv\n")

# =============================================================================
# COMPARISON TO FULL-DATA RESULTS
# =============================================================================
cat("\n======================================================================\n")
cat("COMPARISON: FULL DATA (n=1321, 20 yrs) vs EXCLUDING 2022-2024 (n=1318, 17 yrs)\n")
cat("======================================================================\n")

full_df <- read_csv("full_hierarchical_results_with_harvest_R.csv", show_col_types = FALSE) %>%
  mutate(Metric = paste0(ifelse(Window == "Fall (Aug-Nov)", "Fall", "Spring"), Variable))

results_df2 <- results_df %>%
  mutate(Metric_clean = Metric)

comparison <- full_df %>%
  select(Metric, Beta_Spill_Mean_full = Beta_Spill_Mean,
         Spill_Lower_full = Spill_HDI_Lower, Spill_Upper_full = Spill_HDI_Upper,
         Beta_Harvest_Mean_full = Beta_Harvest_Mean,
         Harvest_Lower_full = Harvest_HDI_Lower, Harvest_Upper_full = Harvest_HDI_Upper) %>%
  inner_join(results_df %>%
               select(Metric, Beta_Spill_Mean_excl = Beta_Spill_Mean,
                      Spill_Lower_excl = Spill_HDI_Lower, Spill_Upper_excl = Spill_HDI_Upper,
                      Beta_Harvest_Mean_excl = Beta_Harvest_Mean,
                      Harvest_Lower_excl = Harvest_HDI_Lower, Harvest_Upper_excl = Harvest_HDI_Upper),
             by = "Metric")

print(comparison, width = Inf)
write_csv(comparison, "sensitivity_exclude_sparse_years_comparison.csv")
cat("\nSaved: sensitivity_exclude_sparse_years_comparison.csv\n")

cat("\n======================================================================\n")
cat("DONE\n")
cat("======================================================================\n")
