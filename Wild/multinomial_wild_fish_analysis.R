###############################################################################
# Multinomial Analysis for WILD Steelhead Fish - R Version
# Analyzes effects of spill and harvest status on fish fate:
#   - Last detected downstream (Below)
#   - Last detected in tributary upstream (Above Wells)
#   - Last detected at Wells Dam
#
# Uses only wild fish (Rear Type Code = 'W')
#
# Uses brms (Stan backend) for multinomial logistic regression:
#   y_ij ~ Categorical(p_ij)
#   log(p_ij[k] / p_ij[K]) = alpha_jk + beta_spill_k * Spill_j + beta_harvest_k * HarvestOpen_j
#   alpha_jk ~ Normal(mu_alpha_k, sigma_alpha_k)  (year random effects)
#   Reference category: Wells Dam
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
cat("MULTINOMIAL LOGISTIC REGRESSION ANALYSIS (R/brms)\n")
cat("WILD FISH ONLY\n")
cat("Effects of Spill and Harvest Status on Fish Fate\n")
cat("======================================================================\n\n")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================
cat("--- Loading Data ---\n")

# Set working directory to parent folder with data files
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

# =============================================================================
# CLASSIFY FISH FATE USING SITE CLASSIFICATIONS
# =============================================================================
cat("\n--- Classifying Fish Fate ---\n")

# Create lookup for site to zone
site_to_zone <- setNames(site_class_df$Zone, site_class_df$Site)

# Function to classify fate based on last detection site
classify_fate <- function(last_site) {
  last_site_str <- as.character(last_site)
  for (site in names(site_to_zone)) {
    if (grepl(site, last_site_str, fixed = TRUE)) {
      return(site_to_zone[site])
    }
  }
  if (grepl("Wells|WEA|WEL|WEJ", last_site_str)) {
    return("Wells Dam")
  }
  return("Unknown")
}

# Classify each fish
overshoot_df <- overshoot_df %>%
  mutate(ClassifiedZone = sapply(LastDetectionSite, classify_fate))

# Merge rear type
overshoot_df <- overshoot_df %>%
  left_join(rear_type_df, by = "TagCode")

# Merge harvest status
overshoot_df <- overshoot_df %>%
  left_join(harvest_df, by = c("WellsYear" = "Year"))

cat("\nRear Type distribution (all fish):\n")
print(table(overshoot_df$`Rear Type Code`, useNA = "ifany"))

# =============================================================================
# FILTER TO WILD FISH ONLY
# =============================================================================
cat("\n--- Filtering to Wild Fish Only ---\n")

wild_df <- overshoot_df %>%
  filter(`Rear Type Code` == "W")

cat(sprintf("Wild fish: %d\n", nrow(wild_df)))

cat("\nFate distribution for wild fish:\n")
print(table(wild_df$ClassifiedZone))

# =============================================================================
# CALCULATE SPILL METRICS
# =============================================================================
cat("\n--- Calculating Spill Metrics ---\n")

spill_data <- spill_data %>%
  mutate(Year = year(DateTime),
         Month = month(DateTime),
         Day = day(DateTime))

spill_col <- "SPILL (KCFS)"

# Function to calculate spill metrics
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
  mutate(SpillYear = Year + 1) %>%
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

# =============================================================================
# PREPARE ANALYSIS DATA
# =============================================================================
cat("\n--- Preparing Analysis Data ---\n")

# Create SpillYear
wild_df <- wild_df %>%
  mutate(SpillYear = WellsYear + 1)

# Filter to known fates only
analysis_df <- wild_df %>%
  filter(ClassifiedZone %in% c("Below", "Above", "Wells Dam"))

cat(sprintf("Wild fish with known fate: %d\n", nrow(analysis_df)))

# Merge with spill metrics
analysis_df <- analysis_df %>%
  left_join(spill_metrics, by = "SpillYear")

# Create outcome variable as factor with Wells Dam as reference
analysis_df <- analysis_df %>%
  mutate(
    Fate = factor(ClassifiedZone,
                  levels = c("Wells Dam", "Below", "Above")),  # Wells Dam as reference
    HarvestOpen = as.integer(Harvest_Status == "Open")
  )

# Drop rows with missing data
analysis_df <- analysis_df %>%
  filter(!is.na(FallTotalSpill), !is.na(SpringTotalSpill),
         !is.na(Fate), !is.na(HarvestOpen))

cat(sprintf("Wild fish with complete data: %d\n", nrow(analysis_df)))

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================
cat("\n======================================================================\n")
cat("SUMMARY STATISTICS - WILD FISH\n")
cat("======================================================================\n")

cat("\n--- Fate Distribution ---\n")
fate_table <- table(analysis_df$ClassifiedZone)
fate_pcts <- prop.table(fate_table) * 100
for (zone in c("Below", "Above", "Wells Dam")) {
  if (zone %in% names(fate_table)) {
    cat(sprintf("  %s: %d (%.1f%%)\n", zone, fate_table[zone], fate_pcts[zone]))
  }
}

cat("\n--- Fate by Harvest Status ---\n")
print(table(analysis_df$ClassifiedZone, analysis_df$Harvest_Status))

cat("\n--- Fate by Harvest Status (Percentages) ---\n")
crosstab_pct <- prop.table(table(analysis_df$ClassifiedZone, analysis_df$Harvest_Status),
                           margin = 2) * 100
print(round(crosstab_pct, 1))

# Year summary
year_summary <- analysis_df %>%
  group_by(SpillYear) %>%
  summarise(
    N_Fish = n(),
    Below = sum(ClassifiedZone == "Below"),
    Above = sum(ClassifiedZone == "Above"),
    WellsDam = sum(ClassifiedZone == "Wells Dam"),
    FallTotalSpill = first(FallTotalSpill),
    SpringTotalSpill = first(SpringTotalSpill),
    HarvestOpen = first(HarvestOpen),
    .groups = "drop"
  )

cat("\n--- Year Summary ---\n")
print(year_summary)

# Change to Wild directory for outputs
setwd("/home/chas/SteelheadOvershootR/Wild")

# Save year summary
write_csv(year_summary, "wild_fish_year_summary.csv")

# =============================================================================
# MULTINOMIAL LOGISTIC REGRESSION MODEL WITH brms
# =============================================================================
cat("\n======================================================================\n")
cat("FITTING MULTINOMIAL LOGISTIC REGRESSION MODELS\n")
cat("======================================================================\n")

years <- sort(unique(analysis_df$SpillYear))
n_years <- length(years)
n_categories <- 3

cat(sprintf("Number of years: %d (%d-%d)\n", n_years, min(years), max(years)))
cat(sprintf("Total wild fish: %d\n", nrow(analysis_df)))
cat("Outcome categories: Downstream (Below), Tributary (Above), Wells Dam (reference)\n")

# Convert SpillYear to factor for random effects
analysis_df$SpillYearFactor <- factor(analysis_df$SpillYear)

# Define priors for multinomial model
# For categorical family in brms, we get coefficients for each non-reference category
model_priors <- c(
  prior(normal(0, 2), class = "Intercept", dpar = "muBelow"),
  prior(normal(0, 2), class = "Intercept", dpar = "muAbove"),
  prior(normal(0, 1), class = "b", dpar = "muBelow"),
  prior(normal(0, 1), class = "b", dpar = "muAbove"),
  prior(normal(0, 1), class = "sd", dpar = "muBelow"),
  prior(normal(0, 1), class = "sd", dpar = "muAbove")
)

# Function to fit multinomial model
run_multinomial_model <- function(df, spill_metric, metric_name) {
  cat(sprintf("\nFitting %s...\n", metric_name))

  # Standardize spill values
  spill_values <- df[[spill_metric]]
  spill_mean <- mean(spill_values, na.rm = TRUE)
  spill_sd <- sd(spill_values, na.rm = TRUE)
  if (spill_sd == 0) spill_sd <- 1
  df$spill_std <- (spill_values - spill_mean) / spill_sd

  # Fit brms multinomial model
  # Model: log(P(Y=k)/P(Y=ref)) = alpha_jk + beta_spill_k * spill + beta_harvest_k * harvest
  # alpha_jk ~ Normal(mu_alpha_k, sigma_alpha_k)
  fit <- brm(
    Fate ~ spill_std + HarvestOpen + (1 | SpillYearFactor),
    data = df,
    family = categorical(link = "logit", refcat = "Wells Dam"),
    prior = model_priors,
    iter = 4500,
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

  # Get coefficient names for the two non-reference categories
  # For Below (Downstream)
  beta_spill_below <- post$b_muBelow_spill_std
  beta_harvest_below <- post$b_muBelow_HarvestOpen

  # For Above (Tributary)
  beta_spill_above <- post$b_muAbove_spill_std
  beta_harvest_above <- post$b_muAbove_HarvestOpen

  # Calculate HDIs
  if (requireNamespace("bayestestR", quietly = TRUE)) {
    spill_below_hdi <- bayestestR::hdi(beta_spill_below, ci = 0.95)
    spill_below_hdi <- c(spill_below_hdi$CI_low, spill_below_hdi$CI_high)

    spill_above_hdi <- bayestestR::hdi(beta_spill_above, ci = 0.95)
    spill_above_hdi <- c(spill_above_hdi$CI_low, spill_above_hdi$CI_high)

    harvest_below_hdi <- bayestestR::hdi(beta_harvest_below, ci = 0.95)
    harvest_below_hdi <- c(harvest_below_hdi$CI_low, harvest_below_hdi$CI_high)

    harvest_above_hdi <- bayestestR::hdi(beta_harvest_above, ci = 0.95)
    harvest_above_hdi <- c(harvest_above_hdi$CI_low, harvest_above_hdi$CI_high)
  } else {
    spill_below_hdi <- quantile(beta_spill_below, probs = c(0.025, 0.975))
    spill_above_hdi <- quantile(beta_spill_above, probs = c(0.025, 0.975))
    harvest_below_hdi <- quantile(beta_harvest_below, probs = c(0.025, 0.975))
    harvest_above_hdi <- quantile(beta_harvest_above, probs = c(0.025, 0.975))
  }

  # Print results
  cat("\nSpill Results (relative to Wells Dam):\n")
  cat(sprintf("  Downstream: beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_spill_below), spill_below_hdi[1], spill_below_hdi[2],
              mean(beta_spill_below > 0)))
  cat(sprintf("  Tributary:  beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_spill_above), spill_above_hdi[1], spill_above_hdi[2],
              mean(beta_spill_above > 0)))

  cat("\nHarvest Status Results (relative to Wells Dam):\n")
  cat(sprintf("  Downstream: beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_harvest_below), harvest_below_hdi[1], harvest_below_hdi[2],
              mean(beta_harvest_below > 0)))
  cat(sprintf("  Tributary:  beta = %.3f, 95%% HDI [%.3f, %.3f], P(beta>0) = %.3f\n",
              mean(beta_harvest_above), harvest_above_hdi[1], harvest_above_hdi[2],
              mean(beta_harvest_above > 0)))

  # Return results
  list(
    metric = metric_name,
    fit = fit,
    # Spill effects
    beta_spill_downstream_mean = mean(beta_spill_below),
    beta_spill_downstream_sd = sd(beta_spill_below),
    beta_spill_downstream_hdi = spill_below_hdi,
    p_spill_downstream_positive = mean(beta_spill_below > 0),

    beta_spill_tributary_mean = mean(beta_spill_above),
    beta_spill_tributary_sd = sd(beta_spill_above),
    beta_spill_tributary_hdi = spill_above_hdi,
    p_spill_tributary_positive = mean(beta_spill_above > 0),

    # Harvest effects
    beta_harvest_downstream_mean = mean(beta_harvest_below),
    beta_harvest_downstream_sd = sd(beta_harvest_below),
    beta_harvest_downstream_hdi = harvest_below_hdi,
    p_harvest_downstream_positive = mean(beta_harvest_below > 0),

    beta_harvest_tributary_mean = mean(beta_harvest_above),
    beta_harvest_tributary_sd = sd(beta_harvest_above),
    beta_harvest_tributary_hdi = harvest_above_hdi,
    p_harvest_tributary_positive = mean(beta_harvest_above > 0),

    # Samples for plotting
    beta_spill_downstream_samples = beta_spill_below,
    beta_spill_tributary_samples = beta_spill_above,
    beta_harvest_downstream_samples = beta_harvest_below,
    beta_harvest_tributary_samples = beta_harvest_above,

    spill_mean = spill_mean,
    spill_sd = spill_sd
  )
}

# --- Run Fall Spill Model ---
cat("\n--- Fall Spill (Aug-Nov) Model ---\n")
fall_result <- run_multinomial_model(analysis_df, "FallTotalSpill", "FallTotalSpill")

# --- Run Spring Spill Model ---
cat("\n--- Spring Spill (Jan-Mar) Model ---\n")
spring_result <- run_multinomial_model(analysis_df, "SpringTotalSpill", "SpringTotalSpill")

# =============================================================================
# MODEL DIAGNOSTICS
# =============================================================================
cat("\n--- Model Diagnostics ---\n")

for (model_name in c("Fall", "Spring")) {
  if (model_name == "Fall") {
    fit <- fall_result$fit
  } else {
    fit <- spring_result$fit
  }
  cat(sprintf("\n%s Model:\n", model_name))

  diag_summary <- summary(fit)$fixed
  cat(sprintf("  Fixed effects Rhat range: [%.3f, %.3f]\n",
              min(diag_summary$Rhat, na.rm = TRUE),
              max(diag_summary$Rhat, na.rm = TRUE)))
  cat(sprintf("  Bulk ESS range: [%.0f, %.0f]\n",
              min(diag_summary$Bulk_ESS, na.rm = TRUE),
              max(diag_summary$Bulk_ESS, na.rm = TRUE)))
}

# =============================================================================
# SAVE RESULTS
# =============================================================================
cat("\n--- Saving Results ---\n")

results_df <- bind_rows(
  tibble(
    Window = "Fall (Aug-Nov)",
    SpillMetric = fall_result$metric,
    Spill_Downstream_Beta = fall_result$beta_spill_downstream_mean,
    Spill_Downstream_SD = fall_result$beta_spill_downstream_sd,
    Spill_Downstream_HDI_Lower = fall_result$beta_spill_downstream_hdi[1],
    Spill_Downstream_HDI_Upper = fall_result$beta_spill_downstream_hdi[2],
    Spill_Downstream_P_Positive = fall_result$p_spill_downstream_positive,
    Spill_Tributary_Beta = fall_result$beta_spill_tributary_mean,
    Spill_Tributary_SD = fall_result$beta_spill_tributary_sd,
    Spill_Tributary_HDI_Lower = fall_result$beta_spill_tributary_hdi[1],
    Spill_Tributary_HDI_Upper = fall_result$beta_spill_tributary_hdi[2],
    Spill_Tributary_P_Positive = fall_result$p_spill_tributary_positive,
    Harvest_Downstream_Beta = fall_result$beta_harvest_downstream_mean,
    Harvest_Downstream_SD = fall_result$beta_harvest_downstream_sd,
    Harvest_Downstream_HDI_Lower = fall_result$beta_harvest_downstream_hdi[1],
    Harvest_Downstream_HDI_Upper = fall_result$beta_harvest_downstream_hdi[2],
    Harvest_Downstream_P_Positive = fall_result$p_harvest_downstream_positive,
    Harvest_Tributary_Beta = fall_result$beta_harvest_tributary_mean,
    Harvest_Tributary_SD = fall_result$beta_harvest_tributary_sd,
    Harvest_Tributary_HDI_Lower = fall_result$beta_harvest_tributary_hdi[1],
    Harvest_Tributary_HDI_Upper = fall_result$beta_harvest_tributary_hdi[2],
    Harvest_Tributary_P_Positive = fall_result$p_harvest_tributary_positive
  ),
  tibble(
    Window = "Spring (Jan-Mar)",
    SpillMetric = spring_result$metric,
    Spill_Downstream_Beta = spring_result$beta_spill_downstream_mean,
    Spill_Downstream_SD = spring_result$beta_spill_downstream_sd,
    Spill_Downstream_HDI_Lower = spring_result$beta_spill_downstream_hdi[1],
    Spill_Downstream_HDI_Upper = spring_result$beta_spill_downstream_hdi[2],
    Spill_Downstream_P_Positive = spring_result$p_spill_downstream_positive,
    Spill_Tributary_Beta = spring_result$beta_spill_tributary_mean,
    Spill_Tributary_SD = spring_result$beta_spill_tributary_sd,
    Spill_Tributary_HDI_Lower = spring_result$beta_spill_tributary_hdi[1],
    Spill_Tributary_HDI_Upper = spring_result$beta_spill_tributary_hdi[2],
    Spill_Tributary_P_Positive = spring_result$p_spill_tributary_positive,
    Harvest_Downstream_Beta = spring_result$beta_harvest_downstream_mean,
    Harvest_Downstream_SD = spring_result$beta_harvest_downstream_sd,
    Harvest_Downstream_HDI_Lower = spring_result$beta_harvest_downstream_hdi[1],
    Harvest_Downstream_HDI_Upper = spring_result$beta_harvest_downstream_hdi[2],
    Harvest_Downstream_P_Positive = spring_result$p_harvest_downstream_positive,
    Harvest_Tributary_Beta = spring_result$beta_harvest_tributary_mean,
    Harvest_Tributary_SD = spring_result$beta_harvest_tributary_sd,
    Harvest_Tributary_HDI_Lower = spring_result$beta_harvest_tributary_hdi[1],
    Harvest_Tributary_HDI_Upper = spring_result$beta_harvest_tributary_hdi[2],
    Harvest_Tributary_P_Positive = spring_result$p_harvest_tributary_positive
  )
)

write_csv(results_df, "wild_fish_multinomial_results.csv")
cat("Saved: wild_fish_multinomial_results.csv\n")

# Save analysis data
write_csv(analysis_df, "wild_fish_analysis_data.csv")
cat("Saved: wild_fish_analysis_data.csv\n")

# =============================================================================
# VISUALIZATIONS
# =============================================================================
cat("\n--- Generating Visualizations ---\n")

# Figure 1: Data overview
# Panel 1: Fate distribution
p1a <- analysis_df %>%
  count(ClassifiedZone) %>%
  mutate(ClassifiedZone = factor(ClassifiedZone,
                                  levels = c("Below", "Above", "Wells Dam"))) %>%
  ggplot(aes(x = ClassifiedZone, y = n, fill = ClassifiedZone)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = c("Below" = "#2196F3", "Above" = "#4CAF50",
                                "Wells Dam" = "#FF5722"),
                    guide = "none") +
  scale_x_discrete(labels = c("Below" = "Downstream\n(Below)",
                               "Above" = "Tributary\n(Above)",
                               "Wells Dam" = "Wells Dam")) +
  labs(x = "", y = "Number of Wild Fish",
       title = sprintf("Wild Fish Fate Distribution (n=%d)", nrow(analysis_df))) +
  theme_minimal()

# Panel 2: Fate by harvest status
p1b <- analysis_df %>%
  count(ClassifiedZone, HarvestOpen) %>%
  mutate(HarvestStatus = ifelse(HarvestOpen == 1, "Open", "Closed"),
         ClassifiedZone = factor(ClassifiedZone,
                                  levels = c("Below", "Above", "Wells Dam"))) %>%
  ggplot(aes(x = ClassifiedZone, y = n, fill = HarvestStatus)) +
  geom_col(position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("Open" = "#4CAF50", "Closed" = "#E91E63")) +
  scale_x_discrete(labels = c("Below" = "Downstream",
                               "Above" = "Tributary",
                               "Wells Dam" = "Wells Dam")) +
  labs(x = "", y = "Number of Wild Fish",
       title = "Wild Fish Fate by Harvest Status") +
  theme_minimal()

# Panel 3: Fate percentages by harvest status
p1c <- analysis_df %>%
  group_by(HarvestOpen, ClassifiedZone) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(HarvestOpen) %>%
  mutate(pct = n / sum(n) * 100,
         HarvestStatus = ifelse(HarvestOpen == 1, "Open", "Closed"),
         ClassifiedZone = factor(ClassifiedZone,
                                  levels = c("Below", "Above", "Wells Dam"))) %>%
  ggplot(aes(x = ClassifiedZone, y = pct, fill = HarvestStatus)) +
  geom_col(position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("Open" = "#4CAF50", "Closed" = "#E91E63")) +
  scale_x_discrete(labels = c("Below" = "Downstream",
                               "Above" = "Tributary",
                               "Wells Dam" = "Wells Dam")) +
  labs(x = "", y = "Percentage of Wild Fish",
       title = "Wild Fish Fate (%) by Harvest Status") +
  theme_minimal()

# Panel 4: Fall spill vs fate
yr_data <- year_summary %>%
  mutate(Below_pct = Below / N_Fish * 100,
         Above_pct = Above / N_Fish * 100)

p1d <- ggplot(yr_data) +
  geom_point(aes(x = FallTotalSpill/1000, y = Below_pct, size = N_Fish),
             alpha = 0.7, color = "#2196F3") +
  geom_point(aes(x = FallTotalSpill/1000, y = Above_pct, size = N_Fish),
             alpha = 0.7, color = "#4CAF50", shape = 15) +
  labs(x = "Fall Total Spill (thousand KCFS-hr)", y = "% of Wild Fish",
       title = "Fall Spill vs Fate (Wild Fish)") +
  theme_minimal() +
  guides(size = "none")

# Panel 5: Spring spill vs fate
p1e <- ggplot(yr_data) +
  geom_point(aes(x = SpringTotalSpill/1000, y = Below_pct, size = N_Fish),
             alpha = 0.7, color = "#2196F3") +
  geom_point(aes(x = SpringTotalSpill/1000, y = Above_pct, size = N_Fish),
             alpha = 0.7, color = "#4CAF50", shape = 15) +
  labs(x = "Spring Total Spill (thousand KCFS-hr)", y = "% of Wild Fish",
       title = "Spring Spill vs Fate (Wild Fish)") +
  theme_minimal() +
  guides(size = "none")

# Panel 6: Sample size by year
p1f <- ggplot(year_summary, aes(x = SpillYear, y = N_Fish,
                                 fill = factor(HarvestOpen))) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = c("0" = "#E91E63", "1" = "#4CAF50"),
                    labels = c("0" = "Closed", "1" = "Open"),
                    name = "Harvest") +
  labs(x = "Spill Year", y = "Number of Wild Fish",
       title = "Wild Fish Sample Size by Year") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig1 <- (p1a | p1b | p1c) / (p1d | p1e | p1f)
ggsave("wild_fish_data_overview.png", fig1, width = 16, height = 10, dpi = 150)
cat("Saved: wild_fish_data_overview.png\n")

# Figure 2: Model results - Spill and Harvest effects
# Fall spill effects
fall_posteriors <- bind_rows(
  tibble(samples = fall_result$beta_spill_downstream_samples,
         Category = "Downstream", Effect = "Spill"),
  tibble(samples = fall_result$beta_spill_tributary_samples,
         Category = "Tributary", Effect = "Spill"),
  tibble(samples = fall_result$beta_harvest_downstream_samples,
         Category = "Downstream", Effect = "Harvest"),
  tibble(samples = fall_result$beta_harvest_tributary_samples,
         Category = "Tributary", Effect = "Harvest")
)

spring_posteriors <- bind_rows(
  tibble(samples = spring_result$beta_spill_downstream_samples,
         Category = "Downstream", Effect = "Spill"),
  tibble(samples = spring_result$beta_spill_tributary_samples,
         Category = "Tributary", Effect = "Spill"),
  tibble(samples = spring_result$beta_harvest_downstream_samples,
         Category = "Downstream", Effect = "Harvest"),
  tibble(samples = spring_result$beta_harvest_tributary_samples,
         Category = "Tributary", Effect = "Harvest")
)

p2a <- fall_posteriors %>%
  filter(Effect == "Spill") %>%
  ggplot(aes(x = samples, fill = Category)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Downstream" = "#2196F3", "Tributary" = "#4CAF50")) +
  labs(x = "beta_spill (standardized)", y = "Density",
       title = sprintf("Fall Spill Effect on Fate (Wild Fish)\nDownstream: P(>0)=%.2f, Tributary: P(>0)=%.2f",
                       fall_result$p_spill_downstream_positive,
                       fall_result$p_spill_tributary_positive)) +
  theme_minimal()

p2b <- spring_posteriors %>%
  filter(Effect == "Spill") %>%
  ggplot(aes(x = samples, fill = Category)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Downstream" = "#2196F3", "Tributary" = "#4CAF50")) +
  labs(x = "beta_spill (standardized)", y = "Density",
       title = sprintf("Spring Spill Effect on Fate (Wild Fish)\nDownstream: P(>0)=%.2f, Tributary: P(>0)=%.2f",
                       spring_result$p_spill_downstream_positive,
                       spring_result$p_spill_tributary_positive)) +
  theme_minimal()

p2c <- fall_posteriors %>%
  filter(Effect == "Harvest") %>%
  ggplot(aes(x = samples, fill = Category)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Downstream" = "#2196F3", "Tributary" = "#4CAF50")) +
  labs(x = "beta_harvest (Harvest Open)", y = "Density",
       title = sprintf("Fall Model: Harvest Status Effect (Wild Fish)\nDownstream: P(>0)=%.2f, Tributary: P(>0)=%.2f",
                       fall_result$p_harvest_downstream_positive,
                       fall_result$p_harvest_tributary_positive)) +
  theme_minimal()

p2d <- spring_posteriors %>%
  filter(Effect == "Harvest") %>%
  ggplot(aes(x = samples, fill = Category)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Downstream" = "#2196F3", "Tributary" = "#4CAF50")) +
  labs(x = "beta_harvest (Harvest Open)", y = "Density",
       title = sprintf("Spring Model: Harvest Status Effect (Wild Fish)\nDownstream: P(>0)=%.2f, Tributary: P(>0)=%.2f",
                       spring_result$p_harvest_downstream_positive,
                       spring_result$p_harvest_tributary_positive)) +
  theme_minimal()

fig2 <- (p2a | p2b) / (p2c | p2d)
ggsave("wild_fish_model_results.png", fig2, width = 14, height = 10, dpi = 150)
cat("Saved: wild_fish_model_results.png\n")

# Figure 3: Trace plots for Fall model (diagnostics)
cat("Generating trace plots...\n")
trace_plot <- mcmc_plot(fall_result$fit,
                        pars = c("b_muBelow_spill_std", "b_muBelow_HarvestOpen",
                                 "b_muAbove_spill_std", "b_muAbove_HarvestOpen"),
                        type = "trace")
ggsave("wild_fish_trace_plots.png", trace_plot, width = 12, height = 10, dpi = 150)
cat("Saved: wild_fish_trace_plots.png\n")

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("SUMMARY - WILD FISH MULTINOMIAL ANALYSIS\n")
cat("======================================================================\n")

# Calculate odds ratios
or_spill_down_fall <- exp(fall_result$beta_spill_downstream_mean)
or_spill_trib_fall <- exp(fall_result$beta_spill_tributary_mean)
or_harvest_down_fall <- exp(fall_result$beta_harvest_downstream_mean)
or_harvest_trib_fall <- exp(fall_result$beta_harvest_tributary_mean)

or_spill_down_spring <- exp(spring_result$beta_spill_downstream_mean)
or_spill_trib_spring <- exp(spring_result$beta_spill_tributary_mean)
or_harvest_down_spring <- exp(spring_result$beta_harvest_downstream_mean)
or_harvest_trib_spring <- exp(spring_result$beta_harvest_tributary_mean)

cat(sprintf("
Dataset: %d WILD fish across %d years (%d-%d)

Fate Distribution:
  - Downstream (Below): %d (%.1f%%)
  - Tributary (Above): %d (%.1f%%)
  - Wells Dam: %d (%.1f%%)

Harvest Status:
  - Open years: %d years, %d fish
  - Closed years: %d years, %d fish

MODEL RESULTS (Reference Category: Wells Dam)
==============================================

FALL SPILL MODEL (August 1 - November 30):

  Spill Effect on Downstream vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio (per 1 SD increase in spill) = %.3f

  Spill Effect on Tributary vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio (per 1 SD increase in spill) = %.3f

  Harvest Open Effect on Downstream vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

  Harvest Open Effect on Tributary vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

SPRING SPILL MODEL (January 1 - March 31):

  Spill Effect on Downstream vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio (per 1 SD increase in spill) = %.3f

  Spill Effect on Tributary vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio (per 1 SD increase in spill) = %.3f

  Harvest Open Effect on Downstream vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

  Harvest Open Effect on Tributary vs Wells Dam:
    beta = %.3f
    95%% HDI [%.3f, %.3f]
    P(beta > 0) = %.3f
    Odds Ratio = %.3f

INTERPRETATION:
- Reference category is Wells Dam (fish last detected at Wells Dam)
- Positive beta: increases probability of that fate relative to Wells Dam
- Negative beta: decreases probability of that fate relative to Wells Dam
- Harvest Open = 1 when steelhead harvest was open in the return year (even though
  wild fish cannot legally be harvested, this tests for indirect effects)
",
nrow(analysis_df), n_years, min(years), max(years),
fate_table["Below"], fate_pcts["Below"],
fate_table["Above"], fate_pcts["Above"],
fate_table["Wells Dam"], fate_pcts["Wells Dam"],
sum(analysis_df$HarvestOpen == 1), sum(analysis_df$HarvestOpen == 1),
sum(analysis_df$HarvestOpen == 0), sum(analysis_df$HarvestOpen == 0),
# Fall spill
fall_result$beta_spill_downstream_mean,
fall_result$beta_spill_downstream_hdi[1], fall_result$beta_spill_downstream_hdi[2],
fall_result$p_spill_downstream_positive, or_spill_down_fall,
fall_result$beta_spill_tributary_mean,
fall_result$beta_spill_tributary_hdi[1], fall_result$beta_spill_tributary_hdi[2],
fall_result$p_spill_tributary_positive, or_spill_trib_fall,
# Fall harvest
fall_result$beta_harvest_downstream_mean,
fall_result$beta_harvest_downstream_hdi[1], fall_result$beta_harvest_downstream_hdi[2],
fall_result$p_harvest_downstream_positive, or_harvest_down_fall,
fall_result$beta_harvest_tributary_mean,
fall_result$beta_harvest_tributary_hdi[1], fall_result$beta_harvest_tributary_hdi[2],
fall_result$p_harvest_tributary_positive, or_harvest_trib_fall,
# Spring spill
spring_result$beta_spill_downstream_mean,
spring_result$beta_spill_downstream_hdi[1], spring_result$beta_spill_downstream_hdi[2],
spring_result$p_spill_downstream_positive, or_spill_down_spring,
spring_result$beta_spill_tributary_mean,
spring_result$beta_spill_tributary_hdi[1], spring_result$beta_spill_tributary_hdi[2],
spring_result$p_spill_tributary_positive, or_spill_trib_spring,
# Spring harvest
spring_result$beta_harvest_downstream_mean,
spring_result$beta_harvest_downstream_hdi[1], spring_result$beta_harvest_downstream_hdi[2],
spring_result$p_harvest_downstream_positive, or_harvest_down_spring,
spring_result$beta_harvest_tributary_mean,
spring_result$beta_harvest_tributary_hdi[1], spring_result$beta_harvest_tributary_hdi[2],
spring_result$p_harvest_tributary_positive, or_harvest_trib_spring
))

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")

# =============================================================================
# SAVE MODEL OBJECTS
# =============================================================================
save(fall_result, spring_result, analysis_df, year_summary,
     file = "wild_fish_multinomial_models.RData")
cat("Saved model objects: wild_fish_multinomial_models.RData\n")
