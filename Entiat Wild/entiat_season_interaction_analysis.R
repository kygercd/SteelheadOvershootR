###############################################################################
# Season × Wells Interaction Analysis: Wild Entiat Steelhead
#
# QUESTION
# --------
# Does the Wells Dam interaction effect on Entiat return probability differ
# between fall (Jul-Dec) and spring (Jan-May) returning fish?
#
# Raw rates:
#   Group A Fall:   74.0% (n=204)   Group A Spring: 68.6% (n=35)
#   Group B Fall:   57.7% (n=123)   Group B Spring: 64.7% (n=17)
#   Fall gap: -16.3 pp              Spring gap:  -3.9 pp
#
# MODELS
# ------
# Model 1e: EntiatDetected ~ wells_interaction * spring_return
#                           + harvest_open + (1|ReturnYear)
#   - Standard Bernoulli likelihood
#   - Interaction term tests whether Wells effect differs by season
#
# Model 1f: Y | vreal(p_enl) ~ wells_interaction * spring_return
#                              + harvest_open + (1|ReturnYear)
#   - Compound Bernoulli likelihood (detection-corrected)
#   - Tests if the seasonal interaction survives detection correction
#
# NOTE: Spring n is small (35 + 17 = 52 fish). Posterior for spring Wells
# effect will be wide. The interaction term is exploratory.
###############################################################################

library(brms)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(patchwork)
library(posterior)
library(tidybayes)

set.seed(42)
setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

cat("======================================================================\n")
cat("SEASON × WELLS INTERACTION ANALYSIS\n")
cat("======================================================================\n\n")

# =============================================================================
# LOAD DATA
# =============================================================================

fish <- read_csv("entiat_analysis_with_detection_eff.csv", show_col_types = FALSE) %>%
  mutate(
    p_enl_safe       = ifelse(is.na(p_enl), 1.0, p_enl),
    ReturnYearFactor = factor(ReturnYear),
    spring_return    = as.integer(ReturnSeason == "Spring"),
    wells_interaction = as.integer(group == "B")
  )

cat(sprintf("n = %d fish (Group A: %d, Group B: %d)\n",
            nrow(fish), sum(fish$group=="A"), sum(fish$group=="B")))

# =============================================================================
# DESCRIPTIVE STATISTICS
# =============================================================================

cat("\n--- Group × Season Detection Rates ---\n")

desc <- fish %>%
  mutate(Season = if_else(spring_return == 1, "Spring (Jan-May)", "Fall (Jul-Dec)")) %>%
  group_by(group, Season) %>%
  summarise(
    n          = n(),
    n_detected = sum(EntiatDetected),
    pct        = round(mean(EntiatDetected) * 100, 1),
    # Wilson score 95% CI
    lo95 = round(100 * (pct/100 + qnorm(0.025)^2/(2*n) -
                          qnorm(0.975) * sqrt(pct/100*(1-pct/100)/n + qnorm(0.025)^2/(4*n^2))) /
                   (1 + qnorm(0.025)^2/n), 1),
    hi95 = round(100 * (pct/100 + qnorm(0.025)^2/(2*n) +
                          qnorm(0.975) * sqrt(pct/100*(1-pct/100)/n + qnorm(0.025)^2/(4*n^2))) /
                   (1 + qnorm(0.025)^2/n), 1),
    .groups = "drop"
  )

print(desc)

# Gap by season
cat("\nWells effect by season (raw):\n")
fall_A   <- desc$pct[desc$group=="A" & grepl("Fall",   desc$Season)]
fall_B   <- desc$pct[desc$group=="B" & grepl("Fall",   desc$Season)]
spring_A <- desc$pct[desc$group=="A" & grepl("Spring", desc$Season)]
spring_B <- desc$pct[desc$group=="B" & grepl("Spring", desc$Season)]

cat(sprintf("  Fall gap:   Group A %.1f%% vs Group B %.1f%% = %.1f pp\n",
            fall_A, fall_B, fall_B - fall_A))
cat(sprintf("  Spring gap: Group A %.1f%% vs Group B %.1f%% = %.1f pp\n",
            spring_A, spring_B, spring_B - spring_A))

# Log-odds Wells effect by season
lo_fall_A   <- qlogis(fall_A/100);   lo_fall_B   <- qlogis(fall_B/100)
lo_spring_A <- qlogis(spring_A/100); lo_spring_B <- qlogis(spring_B/100)
cat(sprintf("\n  Fall Wells effect (log-odds):   %.3f\n", lo_fall_B   - lo_fall_A))
cat(sprintf("  Spring Wells effect (log-odds): %.3f\n", lo_spring_B - lo_spring_A))
cat(sprintf("  Interaction estimate (naive):   %.3f\n",
            (lo_spring_B - lo_spring_A) - (lo_fall_B - lo_fall_A)))

# =============================================================================
# PRIORS AND MCMC SETTINGS
# =============================================================================

model_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "sd")
)

MCMC_ITER   <- 4500
MCMC_WARMUP <- 1500
MCMC_CHAINS <- 4
MCMC_CORES  <- 4
ADAPT_DELTA <- 0.95

# =============================================================================
# MODEL 1e: STANDARD INTERACTION MODEL
# EntiatDetected ~ wells_interaction * spring_return + harvest_open + (1|Year)
# =============================================================================

cat("\n======================================================================\n")
cat("MODEL 1e: Standard — Wells × Season Interaction\n")
cat("======================================================================\n")

fit_1e <- brm(
  EntiatDetected ~ wells_interaction * spring_return + harvest_open
                 + (1 | ReturnYearFactor),
  data    = fish,
  family  = bernoulli(link = "logit"),
  prior   = model_priors,
  iter    = MCMC_ITER,  warmup = MCMC_WARMUP,
  chains  = MCMC_CHAINS, cores = MCMC_CORES,
  seed    = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2, refresh = 500
)

cat("\nModel 1e Summary:\n")
print(summary(fit_1e))

post_1e <- as_draws_df(fit_1e)

# Key parameters
b_wells_fall    <- post_1e$b_wells_interaction           # Wells effect in fall
b_spring_A      <- post_1e$b_spring_return               # Season effect, Group A
b_interaction   <- post_1e$`b_wells_interaction:spring_return`  # Interaction
b_wells_spring  <- b_wells_fall + b_interaction          # Wells effect in spring

ci_fall    <- quantile(b_wells_fall,   c(0.025, 0.975))
ci_spring  <- quantile(b_wells_spring, c(0.025, 0.975))
ci_interact <- quantile(b_interaction, c(0.025, 0.975))

cat(sprintf("\n--- Model 1e Key Results ---\n"))
cat(sprintf("Wells effect in FALL (reference):\n"))
cat(sprintf("  beta = %.3f, OR = %.3f, 95%% CI [%.3f, %.3f], P(<0) = %.4f\n",
            mean(b_wells_fall), exp(mean(b_wells_fall)),
            ci_fall[1], ci_fall[2], mean(b_wells_fall < 0)))
cat(sprintf("\nWells effect in SPRING (fall + interaction):\n"))
cat(sprintf("  beta = %.3f, OR = %.3f, 95%% CI [%.3f, %.3f], P(<0) = %.4f\n",
            mean(b_wells_spring), exp(mean(b_wells_spring)),
            ci_spring[1], ci_spring[2], mean(b_wells_spring < 0)))
cat(sprintf("\nInteraction term (spring vs fall difference in Wells effect):\n"))
cat(sprintf("  beta = %.3f, 95%% CI [%.3f, %.3f]\n",
            mean(b_interaction), ci_interact[1], ci_interact[2]))
cat(sprintf("  P(interaction > 0) = %.4f  [i.e. P(Wells effect smaller in spring)]\n",
            mean(b_interaction > 0)))

# Predicted probabilities (population-level, average year)
# Group A fall, Group A spring, Group B fall, Group B spring
alpha      <- post_1e$b_Intercept
b_harv     <- post_1e$b_harvest_open
harv_mean  <- mean(fish$harvest_open)

p_A_fall   <- plogis(alpha + b_harv * harv_mean)
p_A_spring <- plogis(alpha + b_spring_A + b_harv * harv_mean)
p_B_fall   <- plogis(alpha + b_wells_fall + b_harv * harv_mean)
p_B_spring <- plogis(alpha + b_wells_fall + b_spring_A + b_interaction + b_harv * harv_mean)

cat(sprintf("\nPredicted return probabilities (population-level):\n"))
cat(sprintf("  Group A Fall:   %.1f%% [%.1f, %.1f]\n",
            mean(p_A_fall)*100,
            quantile(p_A_fall, 0.025)*100, quantile(p_A_fall, 0.975)*100))
cat(sprintf("  Group A Spring: %.1f%% [%.1f, %.1f]\n",
            mean(p_A_spring)*100,
            quantile(p_A_spring, 0.025)*100, quantile(p_A_spring, 0.975)*100))
cat(sprintf("  Group B Fall:   %.1f%% [%.1f, %.1f]\n",
            mean(p_B_fall)*100,
            quantile(p_B_fall, 0.025)*100, quantile(p_B_fall, 0.975)*100))
cat(sprintf("  Group B Spring: %.1f%% [%.1f, %.1f]\n",
            mean(p_B_spring)*100,
            quantile(p_B_spring, 0.025)*100, quantile(p_B_spring, 0.975)*100))

diag_1e <- summary(fit_1e)$fixed
cat(sprintf("\nDiagnostics: Rhat [%.3f, %.3f], ESS [%.0f, %.0f]\n",
            min(diag_1e$Rhat), max(diag_1e$Rhat),
            min(diag_1e$Bulk_ESS), max(diag_1e$Bulk_ESS)))

# =============================================================================
# MODEL 1f: DETECTION-CORRECTED INTERACTION MODEL
# =============================================================================

cat("\n======================================================================\n")
cat("MODEL 1f: Detection-Corrected — Wells × Season Interaction\n")
cat("======================================================================\n")

# Reuse the custom family defined in entiat_detection_corrected_analysis.R
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

fit_1f <- brm(
  bf(EntiatDetected | vreal(p_enl_safe) ~
       wells_interaction * spring_return + harvest_open + (1 | ReturnYearFactor)),
  data     = fish,
  family   = det_corrected_family,
  stanvars = sv,
  prior    = model_priors,
  iter     = MCMC_ITER,  warmup = MCMC_WARMUP,
  chains   = MCMC_CHAINS, cores = MCMC_CORES,
  seed     = 42, control = list(adapt_delta = ADAPT_DELTA),
  silent   = 2, refresh = 500
)

cat("\nModel 1f Summary:\n")
print(summary(fit_1f))

post_1f <- as_draws_df(fit_1f)

b_wells_fall_f   <- post_1f$b_wells_interaction
b_spring_A_f     <- post_1f$b_spring_return
b_interaction_f  <- post_1f$`b_wells_interaction:spring_return`
b_wells_spring_f <- b_wells_fall_f + b_interaction_f

ci_fall_f     <- quantile(b_wells_fall_f,   c(0.025, 0.975))
ci_spring_f   <- quantile(b_wells_spring_f, c(0.025, 0.975))
ci_interact_f <- quantile(b_interaction_f,  c(0.025, 0.975))

cat(sprintf("\n--- Model 1f Key Results (Detection-Corrected) ---\n"))
cat(sprintf("Wells effect in FALL:\n"))
cat(sprintf("  beta = %.3f, OR = %.3f, 95%% CI [%.3f, %.3f], P(<0) = %.4f\n",
            mean(b_wells_fall_f), exp(mean(b_wells_fall_f)),
            ci_fall_f[1], ci_fall_f[2], mean(b_wells_fall_f < 0)))
cat(sprintf("Wells effect in SPRING:\n"))
cat(sprintf("  beta = %.3f, OR = %.3f, 95%% CI [%.3f, %.3f], P(<0) = %.4f\n",
            mean(b_wells_spring_f), exp(mean(b_wells_spring_f)),
            ci_spring_f[1], ci_spring_f[2], mean(b_wells_spring_f < 0)))
cat(sprintf("Interaction (spring vs fall):\n"))
cat(sprintf("  beta = %.3f, 95%% CI [%.3f, %.3f], P(>0) = %.4f\n",
            mean(b_interaction_f), ci_interact_f[1], ci_interact_f[2],
            mean(b_interaction_f > 0)))

# =============================================================================
# COMPARISON TABLE
# =============================================================================

cat("\n======================================================================\n")
cat("COMPARISON: Model 1e vs 1f — Wells Effect by Season\n")
cat("======================================================================\n")

comparison <- tibble(
  Season  = c("Fall", "Fall", "Spring", "Spring", "Interaction", "Interaction"),
  Model   = c("Standard (1e)", "Det-corrected (1f)",
              "Standard (1e)", "Det-corrected (1f)",
              "Standard (1e)", "Det-corrected (1f)"),
  Beta    = c(mean(b_wells_fall),    mean(b_wells_fall_f),
              mean(b_wells_spring),  mean(b_wells_spring_f),
              mean(b_interaction),   mean(b_interaction_f)),
  OR      = exp(c(mean(b_wells_fall),   mean(b_wells_fall_f),
                  mean(b_wells_spring), mean(b_wells_spring_f),
                  NA, NA)),
  CI_low  = c(ci_fall[1],    ci_fall_f[1],
              ci_spring[1],  ci_spring_f[1],
              ci_interact[1], ci_interact_f[1]),
  CI_high = c(ci_fall[2],    ci_fall_f[2],
              ci_spring[2],  ci_spring_f[2],
              ci_interact[2], ci_interact_f[2]),
  P_neg   = c(mean(b_wells_fall < 0),    mean(b_wells_fall_f < 0),
              mean(b_wells_spring < 0),  mean(b_wells_spring_f < 0),
              NA, NA)
)

print(comparison %>% mutate(across(c(Beta, OR, CI_low, CI_high, P_neg), ~ round(.x, 3))))

# =============================================================================
# SAVE RESULTS
# =============================================================================

save(fit_1e, fit_1f, fish, comparison,
     b_wells_fall, b_wells_spring, b_interaction,
     b_wells_fall_f, b_wells_spring_f, b_interaction_f,
     file = "entiat_season_interaction_models.RData")
cat("\nSaved: entiat_season_interaction_models.RData\n")

write_csv(comparison, "entiat_season_interaction_results.csv")
cat("Saved: entiat_season_interaction_results.csv\n")

# =============================================================================
# VISUALIZATIONS
# =============================================================================

cat("\n--- Generating Figures ---\n")

harvest_colors <- c("Open" = "#E07B54", "Closed" = "#5B8DB8")
group_colors   <- c("A" = "#4CAF50", "B" = "#FF5722")
season_colors  <- c("Fall (Jul-Dec)" = "#E07B54", "Spring (Jan-May)" = "#5B8DB8")

# ---- Figure 1: Raw rates by group and season with CIs ----

desc_plot <- desc %>%
  mutate(
    Group  = if_else(group == "A", "Group A\n(no Wells)", "Group B\n(Wells interaction)"),
    label  = sprintf("%.1f%%\n(n=%d)", pct, n)
  )

p1 <- ggplot(desc_plot, aes(x = Group, y = pct, fill = Season)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                position = position_dodge(0.7), width = 0.25, linewidth = 0.7) +
  geom_text(aes(label = label, y = hi95 + 2),
            position = position_dodge(0.7), size = 3.2, vjust = 0) +
  scale_fill_manual(values = season_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  labs(x = "", y = "% Detected in Entiat River",
       title = "Entiat Return Rate by Group and Season",
       subtitle = "Error bars = 95% Wilson confidence interval",
       fill = "Return season") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ---- Figure 2: Posterior Wells effect by season ----

wells_post_df <- bind_rows(
  tibble(beta = b_wells_fall,    season = "Fall\n(standard)",   model = "Standard"),
  tibble(beta = b_wells_spring,  season = "Spring\n(standard)", model = "Standard"),
  tibble(beta = b_wells_fall_f,  season = "Fall\n(corrected)",  model = "Det-corrected"),
  tibble(beta = b_wells_spring_f,season = "Spring\n(corrected)",model = "Det-corrected")
) %>%
  mutate(
    season_clean = gsub("\n.*","", season),
    model_f = factor(model, levels = c("Standard","Det-corrected"))
  )

p2 <- ggplot(wells_post_df, aes(x = beta, fill = season_clean, linetype = model_f)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Fall" = "#E07B54", "Spring" = "#5B8DB8"),
                    name = "Return season") +
  scale_linetype_manual(values = c("Standard" = "solid", "Det-corrected" = "dashed"),
                        name = "Model") +
  labs(x = "Wells Interaction Coefficient (log-odds)",
       y = "Posterior Density",
       title = "Posterior Wells Effect by Season",
       subtitle = "Solid = standard likelihood; Dashed = detection-corrected") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ---- Figure 3: Interaction term posterior ----

interact_df <- bind_rows(
  tibble(beta = b_interaction,   model = "Standard (1e)"),
  tibble(beta = b_interaction_f, model = "Det-corrected (1f)")
)

# Add P(>0) annotation
p_pos_std  <- mean(b_interaction   > 0)
p_pos_corr <- mean(b_interaction_f > 0)

p3 <- ggplot(interact_df, aes(x = beta, fill = model)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("Standard (1e)" = "#9C27B0",
                               "Det-corrected (1f)" = "#009688"),
                    guide = "none") +
  facet_wrap(~model, ncol = 1,
             labeller = labeller(model = c(
               "Standard (1e)"      = sprintf("Standard (1e)  P(interaction>0) = %.3f", p_pos_std),
               "Det-corrected (1f)" = sprintf("Det-corrected (1f)  P(interaction>0) = %.3f", p_pos_corr)
             ))) +
  labs(x = "Interaction: Spring − Fall Wells Effect (log-odds)",
       y = "Posterior Density",
       title = "Season × Wells Interaction Posterior",
       subtitle = "Positive = Wells effect attenuated in spring relative to fall") +
  theme_minimal(base_size = 12)

# ---- Figure 4: Predicted probabilities with CIs ----

pred_df <- tibble(
  Group  = rep(c("Group A\n(no Wells)", "Group B\n(Wells)"), each = 2),
  Season = rep(c("Fall", "Spring"), 2),
  Mean   = c(mean(p_A_fall), mean(p_A_spring), mean(p_B_fall), mean(p_B_spring)) * 100,
  Lo     = c(quantile(p_A_fall,   0.025), quantile(p_A_spring, 0.025),
             quantile(p_B_fall,   0.025), quantile(p_B_spring, 0.025)) * 100,
  Hi     = c(quantile(p_A_fall,   0.975), quantile(p_A_spring, 0.975),
             quantile(p_B_fall,   0.975), quantile(p_B_spring, 0.975)) * 100,
  Type   = "Standard"
)

p4 <- ggplot(pred_df, aes(x = Season, y = Mean, color = Group, group = Group)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi, fill = Group), alpha = 0.15, color = NA) +
  geom_point(data = desc_plot %>%
               mutate(Season  = if_else(grepl("Fall",   Season), "Fall", "Spring"),
                      Group   = if_else(group == "A", "Group A\n(no Wells)", "Group B\n(Wells)")),
             aes(x = Season, y = pct), shape = 1, size = 4, stroke = 1.2) +
  scale_color_manual(values = c("Group A\n(no Wells)" = "#4CAF50",
                                "Group B\n(Wells)"     = "#FF5722")) +
  scale_fill_manual(values  = c("Group A\n(no Wells)" = "#4CAF50",
                                "Group B\n(Wells)"     = "#FF5722"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Return Season", y = "Predicted % Detected in Entiat River",
       title = "Model 1e: Predicted Return Probability by Group × Season",
       subtitle = "Shaded = 95% posterior CI (Model 1e); Open circles = observed rates",
       color = "") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ---- Figure 5: Annual rates by group colored by season ----

annual_df <- fish %>%
  mutate(Season = if_else(spring_return == 1, "Spring", "Fall")) %>%
  group_by(ReturnYear, group, Season) %>%
  summarise(n = n(), pct = mean(EntiatDetected) * 100, .groups = "drop") %>%
  filter(n >= 3)  # only show years with >= 3 fish per cell

p5 <- ggplot(annual_df, aes(x = ReturnYear, y = pct,
                             color = group, shape = Season, size = n)) +
  geom_point(alpha = 0.85) +
  geom_smooth(aes(group = interaction(group, Season)),
              method = "loess", span = 0.9, se = FALSE, linewidth = 0.8, alpha = 0.6) +
  scale_color_manual(values = group_colors,
                     labels = c("A" = "Group A (no Wells)", "B" = "Group B (Wells)"),
                     name = "") +
  scale_shape_manual(values = c("Fall" = 16, "Spring" = 17), name = "Season") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  scale_y_continuous(limits = c(0, 105)) +
  labs(x = "Return Year", y = "% Detected in Entiat River",
       title = "Annual Entiat Return Rates by Group and Season",
       subtitle = "Point size proportional to n; years with <3 fish excluded") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Combine and save ----

fig_main <- (p1 | p4) / (p2 | p3)
ggsave("entiat_season_interaction.png", fig_main, width = 14, height = 11, dpi = 150)
cat("Saved: entiat_season_interaction.png\n")

ggsave("entiat_season_annual.png", p5, width = 11, height = 5, dpi = 150)
cat("Saved: entiat_season_annual.png\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n======================================================================\n")
cat("FINAL SUMMARY\n")
cat("======================================================================\n\n")

cat("RAW DATA:\n")
cat(sprintf("  Fall:   Group A %.1f%% vs Group B %.1f%% (gap = %.1f pp, n=%d+%d)\n",
            fall_A, fall_B, fall_B-fall_A,
            desc$n[desc$group=="A" & grepl("Fall",desc$Season)],
            desc$n[desc$group=="B" & grepl("Fall",desc$Season)]))
cat(sprintf("  Spring: Group A %.1f%% vs Group B %.1f%% (gap = %.1f pp, n=%d+%d)\n",
            spring_A, spring_B, spring_B-spring_A,
            desc$n[desc$group=="A" & grepl("Spring",desc$Season)],
            desc$n[desc$group=="B" & grepl("Spring",desc$Season)]))

cat("\nMODEL 1e (Standard):\n")
cat(sprintf("  Wells effect FALL:   beta=%.3f, OR=%.2f, P(<0)=%.3f\n",
            mean(b_wells_fall), exp(mean(b_wells_fall)), mean(b_wells_fall<0)))
cat(sprintf("  Wells effect SPRING: beta=%.3f, OR=%.2f, P(<0)=%.3f\n",
            mean(b_wells_spring), exp(mean(b_wells_spring)), mean(b_wells_spring<0)))
cat(sprintf("  Interaction (spring-fall diff): beta=%.3f, P(>0)=%.3f\n",
            mean(b_interaction), mean(b_interaction>0)))

cat("\nMODEL 1f (Detection-Corrected):\n")
cat(sprintf("  Wells effect FALL:   beta=%.3f, OR=%.2f, P(<0)=%.3f\n",
            mean(b_wells_fall_f), exp(mean(b_wells_fall_f)), mean(b_wells_fall_f<0)))
cat(sprintf("  Wells effect SPRING: beta=%.3f, OR=%.2f, P(<0)=%.3f\n",
            mean(b_wells_spring_f), exp(mean(b_wells_spring_f)), mean(b_wells_spring_f<0)))
cat(sprintf("  Interaction (spring-fall diff): beta=%.3f, P(>0)=%.3f\n",
            mean(b_interaction_f), mean(b_interaction_f>0)))

cat(sprintf("\nCAUTION: Spring n is small (Group A: %d, Group B: %d).\n",
            desc$n[desc$group=="A" & grepl("Spring",desc$Season)],
            desc$n[desc$group=="B" & grepl("Spring",desc$Season)]))
cat("  The interaction posterior will be wide. Treat as exploratory.\n")

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
