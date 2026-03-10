###############################################################################
# Interannual Spill Variability Analysis — Wells Dam 2005–2024
# Creates figures and summary statistics for the spill variability document
###############################################################################

library(ggplot2)
library(tidyverse)
library(patchwork)

set.seed(42)

# ---- Load data ---------------------------------------------------------------
df <- read.csv("full_year_summary_with_harvest_R.csv")

# Compute SpillDays from SpillHours (hours / 24)
df <- df %>%
  mutate(
    FallSpillDays  = FallSpillHours  / 24,
    SpringSpillDays = SpringSpillHours / 24,
    # Standardized versions (z-scores) as used in the model
    FallTotalSpill_z   = scale(FallTotalSpill)[,1],
    SpringTotalSpill_z = scale(SpringTotalSpill)[,1],
    FallSpillHours_z   = scale(FallSpillHours)[,1],
    SpringSpillHours_z = scale(SpringSpillHours)[,1],
    FallSpillDays_z    = scale(FallSpillDays)[,1],
    SpringSpillDays_z  = scale(SpringSpillDays)[,1],
    HarvestPeriod = ifelse(SpillYear <= 2015 | SpillYear >= 2024,
                           "Harvest Open", "Harvest Closed")
  )

# ---- Summary statistics ------------------------------------------------------
spill_stats <- df %>%
  summarise(
    Fall_TotalSpill_mean   = mean(FallTotalSpill),
    Fall_TotalSpill_sd     = sd(FallTotalSpill),
    Fall_TotalSpill_cv     = sd(FallTotalSpill)/mean(FallTotalSpill)*100,
    Fall_TotalSpill_min    = min(FallTotalSpill),
    Fall_TotalSpill_max    = max(FallTotalSpill),
    Fall_TotalSpill_min_yr = SpillYear[which.min(FallTotalSpill)],
    Fall_TotalSpill_max_yr = SpillYear[which.max(FallTotalSpill)],

    Spring_TotalSpill_mean   = mean(SpringTotalSpill),
    Spring_TotalSpill_sd     = sd(SpringTotalSpill),
    Spring_TotalSpill_cv     = sd(SpringTotalSpill)/mean(SpringTotalSpill)*100,
    Spring_TotalSpill_min    = min(SpringTotalSpill),
    Spring_TotalSpill_max    = max(SpringTotalSpill),
    Spring_TotalSpill_min_yr = SpillYear[which.min(SpringTotalSpill)],
    Spring_TotalSpill_max_yr = SpillYear[which.max(SpringTotalSpill)],

    Fall_Hours_mean   = mean(FallSpillHours),
    Fall_Hours_sd     = sd(FallSpillHours),
    Fall_Hours_cv     = sd(FallSpillHours)/mean(FallSpillHours)*100,
    Fall_Hours_min    = min(FallSpillHours),
    Fall_Hours_max    = max(FallSpillHours),
    Fall_Hours_min_yr = SpillYear[which.min(FallSpillHours)],
    Fall_Hours_max_yr = SpillYear[which.max(FallSpillHours)],

    Spring_Hours_mean   = mean(SpringSpillHours),
    Spring_Hours_sd     = sd(SpringSpillHours),
    Spring_Hours_cv     = sd(SpringSpillHours)/mean(SpringSpillHours)*100,
    Spring_Hours_min    = min(SpringSpillHours),
    Spring_Hours_max    = max(SpringSpillHours),
    Spring_Hours_min_yr = SpillYear[which.min(SpringSpillHours)],
    Spring_Hours_max_yr = SpillYear[which.max(SpringSpillHours)]
  )

cat("\n=== FALL TOTAL SPILL (KCFS-hours) ===\n")
cat(sprintf("  Mean: %.0f  SD: %.0f  CV: %.1f%%\n",
            spill_stats$Fall_TotalSpill_mean, spill_stats$Fall_TotalSpill_sd,
            spill_stats$Fall_TotalSpill_cv))
cat(sprintf("  Min: %.0f (%d)   Max: %.0f (%d)\n",
            spill_stats$Fall_TotalSpill_min, spill_stats$Fall_TotalSpill_min_yr,
            spill_stats$Fall_TotalSpill_max, spill_stats$Fall_TotalSpill_max_yr))

cat("\n=== SPRING TOTAL SPILL (KCFS-hours) ===\n")
cat(sprintf("  Mean: %.0f  SD: %.0f  CV: %.1f%%\n",
            spill_stats$Spring_TotalSpill_mean, spill_stats$Spring_TotalSpill_sd,
            spill_stats$Spring_TotalSpill_cv))
cat(sprintf("  Min: %.0f (%d)   Max: %.0f (%d)\n",
            spill_stats$Spring_TotalSpill_min, spill_stats$Spring_TotalSpill_min_yr,
            spill_stats$Spring_TotalSpill_max, spill_stats$Spring_TotalSpill_max_yr))

cat("\n=== FALL SPILL HOURS ===\n")
cat(sprintf("  Mean: %.0f  SD: %.0f  CV: %.1f%%\n",
            spill_stats$Fall_Hours_mean, spill_stats$Fall_Hours_sd,
            spill_stats$Fall_Hours_cv))
cat(sprintf("  Min: %.0f (%d)   Max: %.0f (%d)\n",
            spill_stats$Fall_Hours_min, spill_stats$Fall_Hours_min_yr,
            spill_stats$Fall_Hours_max, spill_stats$Fall_Hours_max_yr))

cat("\n=== SPRING SPILL HOURS ===\n")
cat(sprintf("  Mean: %.0f  SD: %.0f  CV: %.1f%%\n",
            spill_stats$Spring_Hours_mean, spill_stats$Spring_Hours_sd,
            spill_stats$Spring_Hours_cv))
cat(sprintf("  Min: %.0f (%d)   Max: %.0f (%d)\n",
            spill_stats$Spring_Hours_min, spill_stats$Spring_Hours_min_yr,
            spill_stats$Spring_Hours_max, spill_stats$Spring_Hours_max_yr))

# Correlation between fall and spring
r_total <- cor(df$FallTotalSpill, df$SpringTotalSpill)
r_hours <- cor(df$FallSpillHours, df$SpringSpillHours)
cat(sprintf("\n=== FALL-SPRING CORRELATIONS ===\n"))
cat(sprintf("  TotalSpill r = %.3f\n", r_total))
cat(sprintf("  SpillHours r = %.3f\n", r_hours))

# Trend over time
lm_fall   <- lm(FallTotalSpill   ~ SpillYear, data = df)
lm_spring <- lm(SpringTotalSpill ~ SpillYear, data = df)
cat(sprintf("\n=== TEMPORAL TRENDS ===\n"))
cat(sprintf("  Fall TotalSpill slope: %.0f KCFS-hr/yr  (p = %.3f)\n",
            coef(lm_fall)[2], summary(lm_fall)$coefficients[2,4]))
cat(sprintf("  Spring TotalSpill slope: %.0f KCFS-hr/yr  (p = %.3f)\n",
            coef(lm_spring)[2], summary(lm_spring)$coefficients[2,4]))

# ---- Shared theme/colors -----------------------------------------------------
harvest_colors <- c("Harvest Open" = "#E07B54", "Harvest Closed" = "#5B8DB8")

theme_spill <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    strip.background = element_rect(fill = "grey92")
  )

# ---- FIGURE 1: Bar charts — Fall and Spring Total Spill by year --------------
p1a <- ggplot(df, aes(x = factor(SpillYear), y = FallTotalSpill / 1000,
                       fill = HarvestPeriod)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(df$FallTotalSpill) / 1000,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = 1.5, y = mean(df$FallTotalSpill) / 1000 + 0.3,
           label = "Mean", size = 3.2, color = "black", hjust = 0) +
  scale_fill_manual(values = harvest_colors, name = "Harvest Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Fall Spill Volume (Aug–Nov)",
       x = NULL, y = "Total Spill (×1,000 KCFS-hrs)") +
  theme_spill

p1b <- ggplot(df, aes(x = factor(SpillYear), y = SpringTotalSpill / 1000,
                       fill = HarvestPeriod)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(df$SpringTotalSpill) / 1000,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  annotate("text", x = 1.5, y = mean(df$SpringTotalSpill) / 1000 + 0.5,
           label = "Mean", size = 3.2, color = "black", hjust = 0) +
  scale_fill_manual(values = harvest_colors, name = "Harvest Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "B. Spring Spill Volume (Jan–Mar)",
       x = NULL, y = "Total Spill (×1,000 KCFS-hrs)") +
  theme_spill

fig1 <- (p1a / p1b) +
  plot_annotation(
    title = "Interannual Variability in Spill Volume at Wells Dam (2005–2024)",
    subtitle = "Dashed line = 20-year mean. Color indicates harvest status for that spill year.",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("spill_variability_fig1_volume.png", fig1,
       width = 9, height = 7, dpi = 150, bg = "white")
cat("Saved fig1\n")

# ---- FIGURE 2: Spill Hours (duration) by year --------------------------------
p2a <- ggplot(df, aes(x = factor(SpillYear), y = FallSpillHours,
                       fill = HarvestPeriod)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(df$FallSpillHours),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = harvest_colors, name = "Harvest Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "A. Fall Spill Duration (Aug–Nov)",
       x = NULL, y = "Hours with Spill > 0") +
  theme_spill

p2b <- ggplot(df, aes(x = factor(SpillYear), y = SpringSpillHours,
                       fill = HarvestPeriod)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_hline(yintercept = mean(df$SpringSpillHours),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = harvest_colors, name = "Harvest Status") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "B. Spring Spill Duration (Jan–Mar)",
       x = NULL, y = "Hours with Spill > 0") +
  theme_spill

fig2 <- (p2a / p2b) +
  plot_annotation(
    title = "Interannual Variability in Spill Duration at Wells Dam (2005–2024)",
    subtitle = "Dashed line = 20-year mean. Color indicates harvest status for that spill year.",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("spill_variability_fig2_duration.png", fig2,
       width = 9, height = 7, dpi = 150, bg = "white")
cat("Saved fig2\n")

# ---- FIGURE 3: Fall vs Spring scatter, colored by harvest period -------------
p3 <- ggplot(df, aes(x = FallTotalSpill / 1000, y = SpringTotalSpill / 1000,
                      color = HarvestPeriod, label = SpillYear)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text(nudge_y = 0.5, size = 3, show.legend = FALSE) +
  geom_smooth(method = "lm", se = TRUE, aes(group = 1),
              color = "grey30", fill = "grey80", linewidth = 0.8, alpha = 0.3) +
  scale_color_manual(values = harvest_colors, name = "Harvest Status") +
  labs(
    title = "Fall vs. Spring Total Spill Volume at Wells Dam (2005–2024)",
    subtitle = sprintf("Pearson r = %.2f (n = 20 years). Shaded band = 95%% CI on linear fit.", r_total),
    x = "Fall Total Spill (×1,000 KCFS-hrs, Aug–Nov)",
    y = "Spring Total Spill (×1,000 KCFS-hrs, Jan–Mar)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12))

ggsave("spill_variability_fig3_scatter.png", p3,
       width = 7, height = 6, dpi = 150, bg = "white")
cat("Saved fig3\n")

# ---- FIGURE 4: Standardized spill values (model inputs) ----------------------
df_long_z <- df %>%
  select(SpillYear, HarvestPeriod,
         `Fall Total\nSpill` = FallTotalSpill_z,
         `Fall Spill\nHours`  = FallSpillHours_z,
         `Fall Spill\nDays`   = FallSpillDays_z,
         `Spring Total\nSpill` = SpringTotalSpill_z,
         `Spring Spill\nHours` = SpringSpillHours_z,
         `Spring Spill\nDays`  = SpringSpillDays_z) %>%
  pivot_longer(cols = -c(SpillYear, HarvestPeriod),
               names_to = "Metric", values_to = "Z_score") %>%
  mutate(Season = ifelse(grepl("Fall", Metric), "Fall (Aug–Nov)", "Spring (Jan–Mar)"),
         Metric = factor(Metric, levels = c("Fall Total\nSpill","Fall Spill\nHours",
                                             "Fall Spill\nDays","Spring Total\nSpill",
                                             "Spring Spill\nHours","Spring Spill\nDays")))

p4 <- ggplot(df_long_z, aes(x = factor(SpillYear), y = Z_score,
                              fill = HarvestPeriod)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_wrap(~ Metric, ncol = 3) +
  scale_fill_manual(values = harvest_colors, name = "Harvest Status") +
  scale_y_continuous(breaks = seq(-2, 3, 1)) +
  labs(
    title = "Standardized Spill Metrics Used as Model Inputs (2005–2024)",
    subtitle = "Values expressed as z-scores (mean = 0, SD = 1) as entered into Bayesian hierarchical models.\nPositive = above average; Negative = below average.",
    x = NULL, y = "Standardized Value (z-score)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9.5, color = "grey35"),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(size = 9.5)
  )

ggsave("spill_variability_fig4_standardized.png", p4,
       width = 10, height = 7, dpi = 150, bg = "white")
cat("Saved fig4\n")

# ---- FIGURE 5: Multi-metric comparison heatmap-style lollipop ----------------
# CV and range comparison across metrics
metric_summary <- tibble(
  Metric     = c("Fall\nTotal Spill","Fall\nSpill Hours","Fall\nSpill Days",
                 "Spring\nTotal Spill","Spring\nSpill Hours","Spring\nSpill Days"),
  Season     = c("Fall","Fall","Fall","Spring","Spring","Spring"),
  Mean       = c(mean(df$FallTotalSpill), mean(df$FallSpillHours), mean(df$FallSpillDays),
                 mean(df$SpringTotalSpill), mean(df$SpringSpillHours), mean(df$SpringSpillDays)),
  SD         = c(sd(df$FallTotalSpill), sd(df$FallSpillHours), sd(df$FallSpillDays),
                 sd(df$SpringTotalSpill), sd(df$SpringSpillHours), sd(df$SpringSpillDays)),
  CV_pct     = c(sd(df$FallTotalSpill)/mean(df$FallTotalSpill)*100,
                 sd(df$FallSpillHours)/mean(df$FallSpillHours)*100,
                 sd(df$FallSpillDays)/mean(df$FallSpillDays)*100,
                 sd(df$SpringTotalSpill)/mean(df$SpringTotalSpill)*100,
                 sd(df$SpringSpillHours)/mean(df$SpringSpillHours)*100,
                 sd(df$SpringSpillDays)/mean(df$SpringSpillDays)*100),
  Min        = c(min(df$FallTotalSpill), min(df$FallSpillHours), min(df$FallSpillDays),
                 min(df$SpringTotalSpill), min(df$SpringSpillHours), min(df$SpringSpillDays)),
  Max        = c(max(df$FallTotalSpill), max(df$FallSpillHours), max(df$FallSpillDays),
                 max(df$SpringTotalSpill), max(df$SpringSpillHours), max(df$SpringSpillDays))
) %>%
  mutate(Range_ratio = Max / Min,
         Metric = factor(Metric, levels = rev(Metric)))

p5 <- ggplot(metric_summary, aes(y = Metric, x = CV_pct, color = Season)) +
  geom_segment(aes(y = Metric, yend = Metric, x = 0, xend = CV_pct),
               linewidth = 1, color = "grey75") +
  geom_point(size = 5) +
  geom_text(aes(label = sprintf("%.0f%%", CV_pct)),
            nudge_x = 2.5, size = 3.5, color = "grey20") +
  scale_color_manual(values = c("Fall" = "#D4813A", "Spring" = "#4A86B8")) +
  scale_x_continuous(limits = c(0, 130), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Coefficient of Variation (CV) by Spill Metric",
    subtitle = "CV = SD / Mean × 100%. Higher CV = greater interannual variability.",
    x = "Coefficient of Variation (%)", y = NULL, color = "Season"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

ggsave("spill_variability_fig5_cv.png", p5,
       width = 7, height = 4.5, dpi = 150, bg = "white")
cat("Saved fig5\n")

# ---- FIGURE 6: Temporal trend lines ------------------------------------------
df_long_raw <- df %>%
  select(SpillYear, HarvestPeriod,
         `Fall Total Spill\n(KCFS-hrs)` = FallTotalSpill,
         `Spring Total Spill\n(KCFS-hrs)` = SpringTotalSpill,
         `Fall Spill Hours` = FallSpillHours,
         `Spring Spill Hours` = SpringSpillHours) %>%
  pivot_longer(cols = -c(SpillYear, HarvestPeriod),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Season = ifelse(grepl("Fall", Metric), "Fall (Aug–Nov)", "Spring (Jan–Mar)"))

p6 <- ggplot(df_long_raw %>% filter(grepl("Total Spill", Metric)),
             aes(x = SpillYear, y = Value / 1000, color = Season, shape = HarvestPeriod)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.7) +
  scale_color_manual(values = c("Fall (Aug–Nov)" = "#D4813A", "Spring (Jan–Mar)" = "#4A86B8")) +
  scale_shape_manual(values = c("Harvest Open" = 16, "Harvest Closed" = 17),
                     name = "Harvest Status") +
  scale_x_continuous(breaks = seq(2005, 2024, 2)) +
  labs(
    title = "Temporal Trends in Total Spill Volume at Wells Dam (2005–2024)",
    subtitle = "Dashed lines = linear trends. Triangles = harvest-closed years (2016–2023).",
    x = "Spill Year", y = "Total Spill (×1,000 KCFS-hrs)", color = "Season"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("spill_variability_fig6_trends.png", p6,
       width = 8, height = 5.5, dpi = 150, bg = "white")
cat("Saved fig6\n")

# ---- Print summary table for use in the report --------------------------------
cat("\n\n=== SUMMARY TABLE FOR REPORT ===\n")
cat(sprintf("%-30s %8s %8s %8s %6s %8s %8s\n",
            "Metric", "Mean", "SD", "Min", "Max", "CV(%)", "Range Ratio"))
cat(strrep("-", 80), "\n")
for (i in 1:nrow(metric_summary)) {
  cat(sprintf("%-30s %8.0f %8.0f %8.0f %6.0f %8.1f %8.1f\n",
              as.character(metric_summary$Metric[i]),
              metric_summary$Mean[i], metric_summary$SD[i],
              metric_summary$Min[i], metric_summary$Max[i],
              metric_summary$CV_pct[i], metric_summary$Range_ratio[i]))
}

cat("\nDone. All figures saved.\n")
