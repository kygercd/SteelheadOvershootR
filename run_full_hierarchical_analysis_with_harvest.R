###############################################################################
# Full Hierarchical Analysis for Steelhead Overshoot (2000-2024) - R Version
# Analyzes Fall (Aug-Nov) and Spring (Jan-Mar) spill effects.
# Includes individual fish covariate for harvest exposure.
#
# Harvest Exposure: Yes if fish has Rear Type Code = "H" (hatchery) AND
#                   was detected at Wells Dam in a year with Harvest_Status = "Open"
#
# Uses brms (Stan backend) to replicate the PyMC hierarchical model:
#   y_ij ~ Bernoulli(p_ij)
#   logit(p_ij) = alpha_j + beta_spill * Spill_j + beta_harvest * Harvest_i
#   alpha_j ~ Normal(mu_alpha, sigma_alpha)
#   mu_alpha ~ Normal(0, 2)
#   sigma_alpha ~ HalfNormal(1)
#   beta_spill ~ Normal(0, 1)
#   beta_harvest ~ Normal(0, 1)
###############################################################################

# Load required packages (install if needed)
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

cat("======================================================================\n")
cat("HIERARCHICAL LOGISTIC REGRESSION ANALYSIS (R/brms)\n")
cat("Fall (Aug-Nov) and Spring (Jan-Mar) Spill Windows\n")
cat("WITH HARVEST EXPOSURE COVARIATE\n")
cat("======================================================================\n\n")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================
cat("--- Loading Data ---\n")

# Set working directory
setwd("/home/chas/SteelheadOvershootR")

# Load overshoot fish data
overshoot_df <- read_csv("all_overshoot_fish_2001_2025.csv", show_col_types = FALSE)
overshoot_df <- overshoot_df %>%
  mutate(FirstWellsDate = as_datetime(FirstWellsDate),
         WellsYear = year(FirstWellsDate))

# Load rear type data
rear_type_df <- read_csv("OvershootRearType.csv", show_col_types = FALSE) %>%
  rename(TagCode = `Tag Code`) %>%
  select(TagCode, `Rear Type Code`) %>%
  distinct()

# Load harvest status data
harvest_df <- read_csv("Harvest_Open_Closed.csv", show_col_types = FALSE)

# Load spill data
spill_data <- read_csv("merged_spill_2000_2024.csv", show_col_types = FALSE) %>%
  mutate(DateTime = as_datetime(DateTime))

# Load site classifications
site_class_df <- read_csv("site_classifications.csv", show_col_types = FALSE)

cat(sprintf("Loaded %d overshoot fish\n", nrow(overshoot_df)))
cat(sprintf("Loaded %d fish with rear type info\n", nrow(rear_type_df)))
cat(sprintf("Loaded harvest status for %d years\n", nrow(harvest_df)))

# =============================================================================
# RECLASSIFY FISH FATE USING CORRECTED SITE CLASSIFICATIONS
# =============================================================================
cat("\n--- Reclassifying Fish Fate with Corrected Site Classifications ---\n")

# Create lookup for site to zone
site_to_zone <- setNames(site_class_df$Zone, site_class_df$Site)

# Function to classify fate based on last detection site
classify_fate <- function(last_site) {
  last_site_str <- as.character(last_site)
  # Check if site matches any classified site

  for (site in names(site_to_zone)) {
    if (grepl(site, last_site_str, fixed = TRUE)) {
      return(site_to_zone[site])
    }
  }
  # Default classification based on site name patterns
  if (grepl("Wells|WEA|WEL|WEJ", last_site_str)) {
    return("Wells Dam")
  }
  return("Unknown")
}

# Reclassify each fish
overshoot_df <- overshoot_df %>%
  mutate(ClassifiedZone = sapply(LastDetectionSite, classify_fate))

# Map to fate categories
fate_map <- c("Above" = "Above", "Below" = "Downstream",
              "Wells Dam" = "At Wells", "Unknown" = "Unknown")
overshoot_df <- overshoot_df %>%
  mutate(CorrectedFate = fate_map[ClassifiedZone])

# Show changes from original classification
cat("\nOriginal vs Corrected Fate comparison:\n")
print(table(Original = overshoot_df$Fate, Corrected = overshoot_df$CorrectedFate))

# Check NAU reclassifications
nau_fish <- overshoot_df %>% filter(grepl("NAU", LastDetectionSite))
cat(sprintf("\nFish with NAU as last detection site: %d\n", nrow(nau_fish)))
if (nrow(nau_fish) > 0) {
  cat("  Original fates: "); print(table(nau_fish$Fate))
  cat("  Corrected fates: "); print(table(nau_fish$CorrectedFate))
}

# Use corrected fate going forward
overshoot_df$Fate <- overshoot_df$CorrectedFate

# =============================================================================
# CREATE HARVEST EXPOSURE COVARIATE
# =============================================================================
cat("\n--- Creating Harvest Exposure Covariate ---\n")

# Merge rear type information
overshoot_df <- overshoot_df %>%
  left_join(rear_type_df, by = "TagCode")

n_with_rear <- sum(!is.na(overshoot_df$`Rear Type Code`))
cat(sprintf("Fish with rear type info: %d / %d\n", n_with_rear, nrow(overshoot_df)))

cat("\nRear Type distribution:\n")
print(table(overshoot_df$`Rear Type Code`, useNA = "ifany"))

# Merge harvest status based on Wells detection year
overshoot_df <- overshoot_df %>%
  left_join(harvest_df, by = c("WellsYear" = "Year"))

# Create harvest exposure covariate
# Yes (1) if: Rear Type Code = "H" AND Harvest_Status = "Open"
overshoot_df <- overshoot_df %>%
  mutate(HarvestExposure = as.integer(
    (`Rear Type Code` == "H") & (Harvest_Status == "Open") &
    !is.na(`Rear Type Code`) & !is.na(Harvest_Status)
  ))

cat("\nHarvest Exposure summary:\n")
overshoot_df %>%
  group_by(HarvestExposure) %>%
  summarise(Total_Fish = n(),
            Hatchery_Fish = sum(`Rear Type Code` == "H", na.rm = TRUE),
            .groups = "drop") %>%
  print()

# =============================================================================
# CALCULATE SPILL METRICS
# =============================================================================
cat("\n--- Calculating Spill Metrics ---\n")

spill_data <- spill_data %>%
  mutate(Year = year(DateTime),
         Month = month(DateTime),
         Day = day(DateTime))

spill_col <- "SPILL (KCFS)"

# Function to calculate spill metrics for a given month range
calc_spill_metrics <- function(spill_subset) {
  spill_subset %>%
    group_by(Year) %>%
    summarise(
      SpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
      SpillDays = n_distinct(paste(Month, Day)[.data[[spill_col]] > 0]),
      TotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

# Fall spill (Aug-Nov)
fall_spill <- spill_data %>% filter(Month >= 8, Month <= 11)
fall_metrics <- calc_spill_metrics(fall_spill) %>%
  mutate(SpillYear = Year + 1) %>%  # Fall spill assigned to next year

  rename(FallSpillHours = SpillHours,
         FallSpillDays = SpillDays,
         FallTotalSpill = TotalSpill) %>%
  select(SpillYear, FallSpillHours, FallSpillDays, FallTotalSpill)

# Spring spill (Jan-Mar)
spring_spill <- spill_data %>% filter(Month >= 1, Month <= 3)
spring_metrics <- calc_spill_metrics(spring_spill) %>%
  mutate(SpillYear = Year) %>%
  rename(SpringSpillHours = SpillHours,
         SpringSpillDays = SpillDays,
         SpringTotalSpill = TotalSpill) %>%
  select(SpillYear, SpringSpillHours, SpringSpillDays, SpringTotalSpill)

# Combine spill metrics
spill_metrics <- full_join(fall_metrics, spring_metrics, by = "SpillYear")
cat("Spill metrics calculated for years:", sort(spill_metrics$SpillYear), "\n")

# =============================================================================
# PREPARE ANALYSIS DATA
# =============================================================================
cat("\n--- Preparing Analysis Data ---\n")

# Filter to Wells vs Downstream
analysis_df <- overshoot_df %>%
  filter(Fate %in% c("At Wells", "Downstream"))
cat(sprintf("Fish for analysis (At Wells + Downstream): %d\n", nrow(analysis_df)))

# Merge with spill metrics
analysis_df <- analysis_df %>%
  left_join(spill_metrics, by = "SpillYear") %>%
  mutate(ReturnedDownstream = as.integer(Fate == "Downstream"))

# Drop rows without spill data
analysis_df <- analysis_df %>%
  filter(!is.na(FallTotalSpill), !is.na(SpringTotalSpill))
cat(sprintf("Fish with spill data: %d\n", nrow(analysis_df)))

# Year summary
year_summary <- analysis_df %>%
  group_by(SpillYear) %>%
  summarise(
    Downstream = sum(ReturnedDownstream),
    Total = n(),
    PctDownstream = mean(ReturnedDownstream),
    HarvestExposed = sum(HarvestExposure),
    PctHarvestExposed = mean(HarvestExposure),
    FallTotalSpill = first(FallTotalSpill),
    SpringTotalSpill = first(SpringTotalSpill),
    FallSpillHours = first(FallSpillHours),
    SpringSpillHours = first(SpringSpillHours),
    .groups = "drop"
  ) %>%
  mutate(AtWells = Total - Downstream)

cat("\nYear summary:\n")
print(year_summary %>%
        select(SpillYear, Downstream, AtWells, Total, PctDownstream,
               HarvestExposed, PctHarvestExposed))

# Save year summary
write_csv(year_summary, "full_year_summary_with_harvest_R.csv")

# Save analysis data
write_csv(analysis_df, "analysis_fish_with_harvest_R.csv")

# =============================================================================
# HIERARCHICAL MODELS WITH brms
# =============================================================================
cat("\n======================================================================\n")
cat("FITTING HIERARCHICAL MODELS WITH HARVEST EXPOSURE COVARIATE\n")
cat("======================================================================\n")

years <- sort(unique(analysis_df$SpillYear))
n_years <- length(years)

cat(sprintf("Number of years: %d (%d-%d)\n", n_years, min(years), max(years)))
cat(sprintf("Total fish: %d\n", nrow(analysis_df)))
cat(sprintf("Fish with harvest exposure: %d\n", sum(analysis_df$HarvestExposure)))

# Convert SpillYear to factor for random effects
analysis_df$SpillYearFactor <- factor(analysis_df$SpillYear)

# Define priors matching PyMC specification:
#   mu_alpha ~ Normal(0, 2)  -> Intercept prior
#   sigma_alpha ~ HalfNormal(1)  -> sd prior (half-normal(0,1))
#   beta_spill ~ Normal(0, 1)  -> coefficient prior
#   beta_harvest ~ Normal(0, 1)  -> coefficient prior
model_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "sd")
)

# Function to fit hierarchical model
run_hierarchical_model <- function(df, spill_metric, metric_name) {
  cat(sprintf("\nFitting %s...\n", metric_name))

  # Standardize spill values
  spill_values <- df[[spill_metric]]
  spill_mean <- mean(spill_values, na.rm = TRUE)
  spill_sd <- sd(spill_values, na.rm = TRUE)
  if (spill_sd == 0) spill_sd <- 1
  df$spill_std <- (spill_values - spill_mean) / spill_sd

  # Fit brms model
  # Model: logit(P(downstream)) = alpha_j + beta_spill * spill + beta_harvest * harvest
  # alpha_j ~ Normal(mu_alpha, sigma_alpha)
  fit <- brm(
    ReturnedDownstream ~ spill_std + HarvestExposure + (1 | SpillYearFactor),
    data = df,
    family = bernoulli(link = "logit"),
    prior = model_priors,
    iter = 4500,       # 3000 post-warmup + 1500 warmup = 4500 total
    warmup = 1500,
    chains = 4,
    cores = 4,
    seed = 42,
    control = list(adapt_delta = 0.95),
    silent = 2,
    refresh = 0
  )

  # Extract posterior samples
  post <- as_draws_df(fit)

  beta_spill_samples <- post$b_spill_std
  beta_harvest_samples <- post$b_HarvestExposure
  sigma_samples <- post$sd_SpillYearFactor__Intercept

  # Calculate HDI (95%)
  spill_hdi <- quantile(beta_spill_samples, probs = c(0.025, 0.975))
  harvest_hdi <- quantile(beta_harvest_samples, probs = c(0.025, 0.975))

  # For proper HDI, use bayestestR if available, otherwise use quantile-based
  if (requireNamespace("bayestestR", quietly = TRUE)) {
    spill_hdi_obj <- bayestestR::hdi(beta_spill_samples, ci = 0.95)
    spill_hdi <- c(spill_hdi_obj$CI_low, spill_hdi_obj$CI_high)
    harvest_hdi_obj <- bayestestR::hdi(beta_harvest_samples, ci = 0.95)
    harvest_hdi <- c(harvest_hdi_obj$CI_low, harvest_hdi_obj$CI_high)
  }

  # Extract year random effects
  year_effects <- ranef(fit)$SpillYearFactor[, , "Intercept"]

  # Print results
  cat(sprintf("  Spill: beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_spill_samples), spill_hdi[1], spill_hdi[2],
              mean(beta_spill_samples > 0)))
  cat(sprintf("  Harvest: beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_harvest_samples), harvest_hdi[1], harvest_hdi[2],
              mean(beta_harvest_samples > 0)))

  # Return results
  list(
    metric = metric_name,
    fit = fit,
    beta_spill_mean = mean(beta_spill_samples),
    beta_spill_sd = sd(beta_spill_samples),
    beta_spill_hdi = spill_hdi,
    p_spill_positive = mean(beta_spill_samples > 0),
    beta_harvest_mean = mean(beta_harvest_samples),
    beta_harvest_sd = sd(beta_harvest_samples),
    beta_harvest_hdi = harvest_hdi,
    p_harvest_positive = mean(beta_harvest_samples > 0),
    sigma_alpha_mean = mean(sigma_samples),
    alpha_means = year_effects[, "Estimate"],
    alpha_sds = year_effects[, "Est.Error"],
    beta_spill_samples = beta_spill_samples,
    beta_harvest_samples = beta_harvest_samples,
    spill_mean = spill_mean,
    spill_sd = spill_sd
  )
}

# --- Fall Spill (Aug-Nov) Models ---
cat("\n--- Fall Spill (Aug-Nov) Models with Harvest Exposure ---\n")
fall_results <- list()
for (metric in c("FallSpillHours", "FallSpillDays", "FallTotalSpill")) {
  fall_results[[metric]] <- run_hierarchical_model(analysis_df, metric, metric)
}

# --- Spring Spill (Jan-Mar) Models ---
cat("\n--- Spring Spill (Jan-Mar) Models with Harvest Exposure ---\n")
spring_results <- list()
for (metric in c("SpringSpillHours", "SpringSpillDays", "SpringTotalSpill")) {
  spring_results[[metric]] <- run_hierarchical_model(analysis_df, metric, metric)
}

# =============================================================================
# MODEL DIAGNOSTICS
# =============================================================================
cat("\n--- Model Diagnostics ---\n")

for (model_name in names(c(fall_results, spring_results))) {
  results <- c(fall_results, spring_results)
  fit <- results[[model_name]]$fit
  cat(sprintf("\n%s:\n", model_name))

  # Rhat and ESS
  diag_summary <- summary(fit)$fixed
  cat(sprintf("  Fixed effects Rhat range: [%.3f, %.3f]\n",
              min(diag_summary$Rhat), max(diag_summary$Rhat)))
  cat(sprintf("  Bulk ESS range: [%.0f, %.0f]\n",
              min(diag_summary$Bulk_ESS), max(diag_summary$Bulk_ESS)))
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
cat("\n--- Saving Results ---\n")

all_results <- bind_rows(
  # Fall results
  map_dfr(fall_results, function(r) {
    tibble(
      Window = "Fall (Aug-Nov)",
      Variable = gsub("Fall", "", r$metric),
      Beta_Spill_Mean = r$beta_spill_mean,
      Beta_Spill_SD = r$beta_spill_sd,
      Spill_HDI_Lower = r$beta_spill_hdi[1],
      Spill_HDI_Upper = r$beta_spill_hdi[2],
      P_Spill_Positive = r$p_spill_positive,
      Beta_Harvest_Mean = r$beta_harvest_mean,
      Beta_Harvest_SD = r$beta_harvest_sd,
      Harvest_HDI_Lower = r$beta_harvest_hdi[1],
      Harvest_HDI_Upper = r$beta_harvest_hdi[2],
      P_Harvest_Positive = r$p_harvest_positive,
      Sigma_Alpha = r$sigma_alpha_mean
    )
  }),
  # Spring results
  map_dfr(spring_results, function(r) {
    tibble(
      Window = "Spring (Jan-Mar)",
      Variable = gsub("Spring", "", r$metric),
      Beta_Spill_Mean = r$beta_spill_mean,
      Beta_Spill_SD = r$beta_spill_sd,
      Spill_HDI_Lower = r$beta_spill_hdi[1],
      Spill_HDI_Upper = r$beta_spill_hdi[2],
      P_Spill_Positive = r$p_spill_positive,
      Beta_Harvest_Mean = r$beta_harvest_mean,
      Beta_Harvest_SD = r$beta_harvest_sd,
      Harvest_HDI_Lower = r$beta_harvest_hdi[1],
      Harvest_HDI_Upper = r$beta_harvest_hdi[2],
      P_Harvest_Positive = r$p_harvest_positive,
      Sigma_Alpha = r$sigma_alpha_mean
    )
  })
)

write_csv(all_results, "full_hierarchical_results_with_harvest_R.csv")
cat("Saved: full_hierarchical_results_with_harvest_R.csv\n")

# =============================================================================
# VISUALIZATIONS
# =============================================================================
cat("\n--- Generating Visualizations ---\n")

# Figure 1: Data overview with harvest exposure
p1a <- ggplot(year_summary %>%
                pivot_longer(cols = c(Downstream, AtWells),
                             names_to = "Fate", values_to = "Count"),
              aes(x = SpillYear, y = Count, fill = Fate)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Downstream" = "#2196F3", "AtWells" = "#FF5722"),
                    labels = c("At Wells" = "Last at Wells",
                               "Downstream" = "Returned Downstream")) +
  labs(x = "Spill Year", y = "Number of Fish",
       title = sprintf("Fish Fate by Year (n=%d)", nrow(analysis_df))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1b <- ggplot(year_summary, aes(x = SpillYear, y = HarvestExposed)) +
  geom_col(fill = "#9C27B0", alpha = 0.7) +
  annotate("rect", xmin = 2015.5, xmax = 2023.5, ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.05) +
  labs(x = "Spill Year", y = "Number of Fish",
       title = "Harvest-Exposed Fish by Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1c <- analysis_df %>%
  group_by(HarvestExposure) %>%
  summarise(PctDownstream = mean(ReturnedDownstream) * 100,
            n = n(), .groups = "drop") %>%
  mutate(Label = ifelse(HarvestExposure == 0, "No Harvest\nExposure", "Harvest\nExposure")) %>%
  ggplot(aes(x = Label, y = PctDownstream, fill = factor(HarvestExposure))) +
  geom_col() +
  geom_text(aes(label = sprintf("n=%d", n)), vjust = -0.5) +
  scale_fill_manual(values = c("0" = "#4CAF50", "1" = "#E91E63"), guide = "none") +
  labs(x = "", y = "% Returned Downstream",
       title = "Downstream Return by Harvest Exposure") +
  theme_minimal()

p1d <- ggplot(year_summary, aes(x = FallTotalSpill/1000, y = PctDownstream * 100,
                                 size = Total)) +
  geom_point(alpha = 0.7, color = "#FF9800") +
  labs(x = "Fall Total Spill (thousand KCFS-hr)", y = "% Returned Downstream",
       title = "Fall Spill (Aug-Nov) vs Downstream Return") +
  theme_minimal() +
  guides(size = "none")

p1e <- ggplot(year_summary, aes(x = SpringTotalSpill/1000, y = PctDownstream * 100,
                                  size = Total)) +
  geom_point(alpha = 0.7, color = "#4CAF50") +
  labs(x = "Spring Total Spill (thousand KCFS-hr)", y = "% Returned Downstream",
       title = "Spring Spill (Jan-Mar) vs Downstream Return") +
  theme_minimal() +
  guides(size = "none")

p1f <- ggplot(year_summary, aes(x = SpillYear, y = Total,
                                  fill = ifelse(SpillYear >= 2016 & SpillYear <= 2023,
                                                "Closed", "Open"))) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = c("Closed" = "#E91E63", "Open" = "#4CAF50"),
                    name = "Harvest") +
  labs(x = "Spill Year", y = "Number of Fish",
       title = "Sample Size by Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- (p1a | p1b | p1c) / (p1d | p1e | p1f)
ggsave("full_data_overview_with_harvest_R.png", fig1, width = 16, height = 10, dpi = 150)
cat("Saved: full_data_overview_with_harvest_R.png\n")

# Figure 2: Spill effects - posterior distributions
fall_spill_posteriors <- bind_rows(
  tibble(samples = fall_results$FallSpillHours$beta_spill_samples, Metric = "SpillHours"),
  tibble(samples = fall_results$FallSpillDays$beta_spill_samples, Metric = "SpillDays"),
  tibble(samples = fall_results$FallTotalSpill$beta_spill_samples, Metric = "TotalSpill")
)

spring_spill_posteriors <- bind_rows(
  tibble(samples = spring_results$SpringSpillHours$beta_spill_samples, Metric = "SpillHours"),
  tibble(samples = spring_results$SpringSpillDays$beta_spill_samples, Metric = "SpillDays"),
  tibble(samples = spring_results$SpringTotalSpill$beta_spill_samples, Metric = "TotalSpill")
)

p2a <- ggplot(fall_spill_posteriors, aes(x = samples, fill = Metric)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("SpillHours" = "#FF9800", "SpillDays" = "#FFC107",
                                "TotalSpill" = "#FFE082")) +
  labs(x = "beta_spill (standardized)", y = "Density",
       title = sprintf("Fall Spill (Aug-Nov) Effect\nP(beta>0) = %.2f",
                       fall_results$FallTotalSpill$p_spill_positive)) +
  theme_minimal()

p2b <- ggplot(spring_spill_posteriors, aes(x = samples, fill = Metric)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("SpillHours" = "#4CAF50", "SpillDays" = "#8BC34A",
                                "TotalSpill" = "#CDDC39")) +
  labs(x = "beta_spill (standardized)", y = "Density",
       title = sprintf("Spring Spill (Jan-Mar) Effect\nP(beta>0) = %.2f",
                       spring_results$SpringTotalSpill$p_spill_positive)) +
  theme_minimal()

fig2 <- p2a | p2b
ggsave("full_spill_effects_with_harvest_R.png", fig2, width = 12, height = 5, dpi = 150)
cat("Saved: full_spill_effects_with_harvest_R.png\n")

# Figure 3: Harvest exposure effect posteriors
fall_harvest_posteriors <- bind_rows(
  tibble(samples = fall_results$FallSpillHours$beta_harvest_samples, Metric = "SpillHours"),
  tibble(samples = fall_results$FallSpillDays$beta_harvest_samples, Metric = "SpillDays"),
  tibble(samples = fall_results$FallTotalSpill$beta_harvest_samples, Metric = "TotalSpill")
)

spring_harvest_posteriors <- bind_rows(
  tibble(samples = spring_results$SpringSpillHours$beta_harvest_samples, Metric = "SpillHours"),
  tibble(samples = spring_results$SpringSpillDays$beta_harvest_samples, Metric = "SpillDays"),
  tibble(samples = spring_results$SpringTotalSpill$beta_harvest_samples, Metric = "TotalSpill")
)

p3a <- ggplot(fall_harvest_posteriors, aes(x = samples, fill = Metric)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("SpillHours" = "#FF9800", "SpillDays" = "#FFC107",
                                "TotalSpill" = "#FFE082")) +
  labs(x = "beta_harvest", y = "Density",
       title = sprintf("Harvest Exposure Effect (Fall Model)\nP(beta>0) = %.2f",
                       fall_results$FallTotalSpill$p_harvest_positive)) +
  theme_minimal()

p3b <- ggplot(spring_harvest_posteriors, aes(x = samples, fill = Metric)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("SpillHours" = "#4CAF50", "SpillDays" = "#8BC34A",
                                "TotalSpill" = "#CDDC39")) +
  labs(x = "beta_harvest", y = "Density",
       title = sprintf("Harvest Exposure Effect (Spring Model)\nP(beta>0) = %.2f",
                       spring_results$SpringTotalSpill$p_harvest_positive)) +
  theme_minimal()

fig3 <- p3a | p3b
ggsave("full_harvest_effects_R.png", fig3, width = 12, height = 5, dpi = 150)
cat("Saved: full_harvest_effects_R.png\n")

# Figure 4: Year random effects
fall_r <- fall_results$FallTotalSpill
spring_r <- spring_results$SpringTotalSpill

year_effects_df <- bind_rows(
  tibble(Year = years, Mean = fall_r$alpha_means, SD = fall_r$alpha_sds, Model = "Fall"),
  tibble(Year = years, Mean = spring_r$alpha_means, SD = spring_r$alpha_sds, Model = "Spring")
)

p4a <- ggplot(year_effects_df %>% filter(Model == "Fall"),
              aes(x = Year, y = Mean)) +
  geom_pointrange(aes(ymin = Mean - 1.96*SD, ymax = Mean + 1.96*SD),
                  color = "#FF9800", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Year", y = "Year Effect (logit scale)",
       title = "Fall Model: Year Random Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4b <- ggplot(year_effects_df %>% filter(Model == "Spring"),
              aes(x = Year, y = Mean)) +
  geom_pointrange(aes(ymin = Mean - 1.96*SD, ymax = Mean + 1.96*SD),
                  color = "#4CAF50", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Year", y = "Year Effect (logit scale)",
       title = "Spring Model: Year Random Effects") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig4 <- p4a | p4b
ggsave("full_year_effects_with_harvest_R.png", fig4, width = 14, height = 5, dpi = 150)
cat("Saved: full_year_effects_with_harvest_R.png\n")

# Figure 5: Trace plots for key parameters (diagnostics)
# Use the FallTotalSpill model as representative
cat("Generating trace plots...\n")
trace_plot <- mcmc_plot(fall_results$FallTotalSpill$fit,
                        pars = c("b_spill_std", "b_HarvestExposure",
                                 "sd_SpillYearFactor__Intercept"),
                        type = "trace")
ggsave("full_trace_plots_R.png", trace_plot, width = 10, height = 8, dpi = 150)
cat("Saved: full_trace_plots_R.png\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("SUMMARY\n")
cat("======================================================================\n")

fall_or <- exp(fall_results$FallTotalSpill$beta_harvest_mean)
spring_or <- exp(spring_results$SpringTotalSpill$beta_harvest_mean)

cat(sprintf("
Dataset: %d fish across %d years (%d-%d)
Harvest-exposed fish: %d (%.1f%%)

FALL SPILL MODEL (August 1 - November 30) - TotalSpill:
  Spill Effect:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f

  Harvest Exposure Effect:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

  Year Variability: sigma_alpha = %.3f

SPRING SPILL MODEL (January 1 - March 31) - TotalSpill:
  Spill Effect:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f

  Harvest Exposure Effect:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

  Year Variability: sigma_alpha = %.3f

INTERPRETATION:
- Harvest Exposure = 1 if fish is hatchery (Rear Type = H) AND detected at Wells
  in a year when harvest was open (2000-2015, 2024-2025)
- Harvest Exposure = 0 otherwise (wild fish OR hatchery in closed years 2016-2023)
- A positive beta_harvest indicates harvest-exposed fish are MORE likely to return downstream
- A negative beta_harvest indicates harvest-exposed fish are LESS likely to return downstream
",
nrow(analysis_df), n_years, min(years), max(years),
sum(analysis_df$HarvestExposure), mean(analysis_df$HarvestExposure)*100,
# Fall
fall_results$FallTotalSpill$beta_spill_mean,
fall_results$FallTotalSpill$beta_spill_hdi[1],
fall_results$FallTotalSpill$beta_spill_hdi[2],
fall_results$FallTotalSpill$p_spill_positive,
fall_results$FallTotalSpill$beta_harvest_mean,
fall_results$FallTotalSpill$beta_harvest_hdi[1],
fall_results$FallTotalSpill$beta_harvest_hdi[2],
fall_results$FallTotalSpill$p_harvest_positive,
fall_or,
fall_results$FallTotalSpill$sigma_alpha_mean,
# Spring
spring_results$SpringTotalSpill$beta_spill_mean,
spring_results$SpringTotalSpill$beta_spill_hdi[1],
spring_results$SpringTotalSpill$beta_spill_hdi[2],
spring_results$SpringTotalSpill$p_spill_positive,
spring_results$SpringTotalSpill$beta_harvest_mean,
spring_results$SpringTotalSpill$beta_harvest_hdi[1],
spring_results$SpringTotalSpill$beta_harvest_hdi[2],
spring_results$SpringTotalSpill$p_harvest_positive,
spring_or,
spring_results$SpringTotalSpill$sigma_alpha_mean
))

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")

# =============================================================================
# OPTIONAL: Save brms model objects for later inspection
# =============================================================================
save(fall_results, spring_results, analysis_df, year_summary,
     file = "hierarchical_models_with_harvest_R.RData")
cat("Saved model objects: hierarchical_models_with_harvest_R.RData\n")
