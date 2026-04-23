###############################################################################
# Generate all manuscript figures
###############################################################################
suppressPackageStartupMessages({
  library(tidyverse)
  library(brms)
  library(lubridate)
  library(patchwork)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

load("entiat_bayesian_models.RData")
load("entiat_season_interaction_models.RData")
load("entiat_cleaned_strict_models.RData")
d <- read_csv("entiat_analysis_with_detection_eff.csv", show_col_types = FALSE)

theme_ms <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"))

harvest_colors <- c("Open" = "#E69F00", "Closed" = "#56B4E9")

# =============================================================================
# FIGURE 1: Data Overview
# =============================================================================

# Panel A: Annual detection rate colored by harvest
annual <- d %>%
  group_by(ReturnYear, harvest_open) %>%
  summarise(n = n(), pct = mean(EntiatDetected) * 100, .groups = "drop") %>%
  mutate(Harvest = if_else(harvest_open == 1, "Open", "Closed"))

pA <- ggplot(annual, aes(x = ReturnYear, y = pct, fill = Harvest)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = harvest_colors) +
  scale_x_continuous(breaks = seq(2006, 2025, 3)) +
  labs(x = NULL, y = "Entiat detection rate (%)", fill = "Harvest",
       title = "A. Annual Entiat return rate") +
  theme_ms

# Panel B: Detection rate by group and season
season_grp <- d %>%
  group_by(group, ReturnSeason) %>%
  summarise(n = n(), pct = mean(EntiatDetected) * 100, .groups = "drop") %>%
  mutate(label = sprintf("%d%%\n(n=%d)", round(pct), n))

pB <- ggplot(season_grp, aes(x = ReturnSeason, y = pct, fill = group)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = label), position = position_dodge(0.7),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("A" = "#009E73", "B" = "#D55E00"),
                    labels = c("A (no Wells)", "B (Wells)")) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = NULL, y = "Entiat detection rate (%)", fill = "Group",
       title = "B. Detection rate by group and season") +
  theme_ms

# Panel C: Overall detection rate by group
overall <- d %>%
  group_by(group) %>%
  summarise(n = n(), pct = mean(EntiatDetected) * 100,
            se = sqrt(pct/100*(1-pct/100)/n) * 100, .groups = "drop") %>%
  mutate(Group = paste0("Group ", group, "\n(n=", n, ")"))

pC <- ggplot(overall, aes(x = Group, y = pct, fill = group)) +
  geom_col(width = 0.5) +
  geom_errorbar(aes(ymin = pct - 1.96*se, ymax = pct + 1.96*se), width = 0.15) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("A" = "#009E73", "B" = "#D55E00"), guide = "none") +
  scale_y_continuous(limits = c(0, 110)) +
  labs(x = NULL, y = "Entiat detection rate (%)",
       title = "C. Overall detection rate by group") +
  theme_ms

# Panel D: Annual sample size (A above axis, B below)
ann_grp <- d %>%
  group_by(ReturnYear, group) %>% summarise(n = n(), .groups = "drop") %>%
  mutate(n_signed = if_else(group == "A", n, -n))

pD <- ggplot(ann_grp, aes(x = ReturnYear, y = n_signed, fill = group)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("A" = "#009E73", "B" = "#D55E00"),
                    labels = c("A (no Wells)", "B (Wells)")) +
  scale_x_continuous(breaks = seq(2006, 2025, 3)) +
  scale_y_continuous(labels = abs) +
  labs(x = "Return year", y = "n fish", fill = "Group",
       title = "D. Annual sample size") +
  theme_ms

fig1 <- (pA + pB) / (pC + pD)
ggsave("fig1_data_overview.png", fig1, width = 10, height = 7, dpi = 150)
cat("Saved fig1_data_overview.png\n")

# =============================================================================
# FIGURE 2: Model 1 posteriors and OR
# =============================================================================

post1 <- as_draws_df(fit_model1)
b_w  <- post1$b_wells_interaction
b_h  <- post1$b_harvest_open

post_long <- bind_rows(
  tibble(param = "Wells interaction (Group B vs. A)", value = b_w),
  tibble(param = "Harvest open", value = b_h)
)

pF2a <- ggplot(post_long, aes(x = value, fill = param)) +
  geom_density(alpha = 0.7, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Wells interaction (Group B vs. A)" = "#0072B2",
                                "Harvest open" = "#E69F00")) +
  labs(x = "Log-odds coefficient (β)", y = "Density", fill = NULL,
       title = "A. Posterior distributions") +
  theme_ms + theme(legend.position = "bottom")

or_df <- tibble(
  Parameter = c("Wells interaction", "Harvest open"),
  OR  = c(exp(mean(b_w)), exp(mean(b_h))),
  lo  = c(exp(quantile(b_w, .025)), exp(quantile(b_h, .025))),
  hi  = c(exp(quantile(b_w, .975)), exp(quantile(b_h, .975)))
)

pF2b <- ggplot(or_df, aes(x = OR, y = Parameter)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, color = "#0072B2") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "Odds ratio (log scale)", y = NULL,
       title = "B. Odds ratio estimates with 95% CI") +
  theme_ms

fig2 <- pF2a + pF2b
ggsave("fig2_model1_posteriors.png", fig2, width = 10, height = 4, dpi = 150)
cat("Saved fig2_model1_posteriors.png\n")

# =============================================================================
# FIGURE 3: Season robustness (Model 1b) + detection-corrected comparison
# =============================================================================

load("entiat_corrected_models.RData")
post1b <- as_draws_df(fit_model1b)
b_w1b  <- post1b$b_wells_interaction
b_spr  <- post1b$b_spring_return
postc  <- as_draws_df(fit_corrected_m1)
b_wc   <- postc$b_wells_interaction

wells_comp <- bind_rows(
  tibble(Model = "M1: standard", value = b_w),
  tibble(Model = "M1b: + season covariate", value = b_w1b),
  tibble(Model = "M1c: detection-corrected", value = b_wc)
)

pF3a <- ggplot(wells_comp, aes(x = value, fill = Model)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("M1: standard" = "#E69F00",
                                "M1b: + season covariate" = "#CC79A7",
                                "M1c: detection-corrected" = "#0072B2")) +
  labs(x = "β (Wells interaction)", y = "Density", fill = NULL,
       title = "A. Wells effect across model variants") +
  theme_ms + theme(legend.position = "bottom")

spring_df <- tibble(value = b_spr)
pF3b <- ggplot(spring_df, aes(x = value)) +
  geom_density(fill = "#009E73", alpha = 0.7, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.5,
           label = sprintf("P(β < 0) = %.3f", mean(b_spr < 0))) +
  labs(x = "β (spring return season)", y = "Density",
       title = "B. Spring season coefficient (Model 1b)") +
  theme_ms

fig3 <- pF3a + pF3b
ggsave("fig3_robustness.png", fig3, width = 10, height = 4, dpi = 150)
cat("Saved fig3_robustness.png\n")

# =============================================================================
# FIGURE 4: Year random effects
# =============================================================================

re <- ranef(fit_model1)$ReturnYearFactor[,,"Intercept"] %>%
  as.data.frame() %>%
  rownames_to_column("Year") %>%
  mutate(Year = as.integer(Year)) %>%
  left_join(d %>% distinct(ReturnYear, harvest_open) %>%
              rename(Year = ReturnYear), by = "Year") %>%
  mutate(Harvest = if_else(harvest_open == 1, "Open", "Closed"))

fig4 <- ggplot(re, aes(x = Year, y = Estimate, color = Harvest)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = `Q2.5`, ymax = `Q97.5`), width = 0.3) +
  geom_point(size = 2.5) +
  scale_color_manual(values = harvest_colors) +
  scale_x_continuous(breaks = seq(2006, 2025, 3)) +
  labs(x = "Return year", y = "Year random effect (log-odds)",
       color = "Harvest", title = "Year-level random effects (Model 1)") +
  theme_ms

ggsave("fig4_year_effects.png", fig4, width = 9, height = 4, dpi = 150)
cat("Saved fig4_year_effects.png\n")

# =============================================================================
# FIGURE 5: Model 2 spill effects
# =============================================================================

post2  <- as_draws_df(fit_model2)
b2s    <- post2$b_SpringSpillHours_z
b2f    <- post2$b_FallSpillHours_z
b2h    <- post2$b_harvest_open

spill_post <- bind_rows(
  tibble(param = "Spring spill hours (std)", value = b2s),
  tibble(param = "Fall spill hours (std)",   value = b2f),
  tibble(param = "Harvest open",              value = b2h)
)

pF5a <- ggplot(spill_post, aes(x = value, fill = param)) +
  geom_density(alpha = 0.65, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Spring spill hours (std)" = "#0072B2",
                                "Fall spill hours (std)"   = "#D55E00",
                                "Harvest open"             = "#009E73")) +
  labs(x = "Log-odds coefficient (β)", y = "Density", fill = NULL,
       title = sprintf("A. Model 2 posteriors (n = %d Group B fish)", nrow(fit_model2$data))) +
  theme_ms + theme(legend.position = "bottom")

pF5b <- group_b_spill %>%
  group_by(ReturnYear, harvest_open) %>%
  summarise(pct = mean(EntiatDetected) * 100,
            FallSpillHours = first(FallSpillHours), .groups = "drop") %>%
  mutate(Harvest = if_else(harvest_open == 1, "Open", "Closed")) %>%
  ggplot(aes(x = FallSpillHours, y = pct, color = Harvest, label = ReturnYear)) +
  geom_text(size = 3) +
  scale_color_manual(values = harvest_colors) +
  labs(x = "Fall spill hours at Wells (Aug–Nov prior year)",
       y = "Group B Entiat detection rate (%)", color = "Harvest",
       title = "B. Fall spill vs. Group B detection rate by year") +
  theme_ms

fig5 <- pF5a + pF5b
ggsave("fig5_model2_spill.png", fig5, width = 11, height = 4.5, dpi = 150)
cat("Saved fig5_model2_spill.png\n")

# =============================================================================
# FIGURE 6: Strict vs. original outcome comparison
# =============================================================================

or_summary <- tribble(
  ~Model,                          ~OR,   ~lo,   ~hi,
  "M1: ENL detected\n(standard)",    exp(-1.348), exp(-1.881), exp(-0.834),
  "M1c: ENL detected\n(det-corr)",   exp(-1.836), exp(-2.737), exp(-1.067),
  "S1: Strict spawned\n(standard)",  exp(-0.593), exp(-1.016), exp(-0.182),
  "S2: Strict spawned\n(det-corr)",  exp(-0.723), exp(-1.251), exp(-0.212)
) %>%
  mutate(Model = factor(Model, levels = rev(Model)))

fig6 <- ggplot(or_summary, aes(x = OR, y = Model)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, color = "#0072B2") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1.0),
                labels = c("0.10", "0.20", "0.50", "1.00")) +
  labs(x = "Odds ratio — Wells Dam interaction (log scale)", y = NULL,
       title = "Wells Dam effect across outcome definitions") +
  theme_ms

ggsave("fig6_outcome_comparison.png", fig6, width = 7, height = 3.5, dpi = 150)
cat("Saved fig6_outcome_comparison.png\n")

cat("\nAll figures saved.\n")
