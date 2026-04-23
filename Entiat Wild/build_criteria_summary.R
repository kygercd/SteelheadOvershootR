# Build per-fish summary table across all outcome criteria discussed
#
# Criteria:
#   1. Original broad    — any Entiat detection post-anchor (EntiatDetected)
#   2. Original strict   — fall: ENL + upstream; spring: ENL alone
#                          (all fish fall-anchored, so effectively all need ENL + upstream)
#   3. Reclassified      — reclassified by actual Entiat entry season;
#                          spring-entry: ENL alone; fall-entry: ENL + upstream;
#                          Group B fall-entry ENL-only recoded to 1
#   4. ENL sensitivity   — reclassified base + credits Group A fall-entry ENL-only
#                          fish with confirmed upstream antenna movement

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(lubridate)
})

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

load("entiat_cleaned_strict_models.RData")    # fish_clean, fit_s1c
load("entiat_reclassified_strict_models.RData") # fish_r
load("entiat_enl_sensitivity_models.RData")   # fish_s

det <- read_csv("EntiatWildSteel.csv", show_col_types = FALSE) |>
  rename(tag = `Tag Code`, site = `Site Name`) |>
  mutate(first_dt = mdy_hm(`First Time Value`))

enl_site     <- "ENL - Lower Entiat River"
upstream_ent <- c("ENM - Middle Entiat River","ENA - Upper Entiat River",
                  "ENS - Entiat River","ENT - Entiat River Trap",
                  "ENA - Upper Entiat River at rkm 17.1",
                  "ENF - Upper Entiat River at rkm 40.6",
                  "ENS - Upper Entiat River at rkm 35.7",
                  "EHL - Entiat NFH Adult Ladder",
                  "MAD - Mad River, Entiat River Basin",
                  "TLT - Tillicum Creek Temporary Array",
                  "RCT - Roaring Creek Temporary Array")
entiat_sites <- c(enl_site, upstream_ent)

# ── Per-fish Entiat detection detail ──────────────────────────────────────────
entiat_detail <- lapply(fish_clean$TagCode, function(tc) {
  fm      <- fish_clean |> filter(TagCode == tc)
  anchor  <- as.POSIXct(fm$rrf_anchor, tz = "UTC")
  win_end <- anchor %m+% months(18)

  d <- det |> filter(tag == tc, first_dt >= anchor, first_dt <= win_end,
                     site %in% entiat_sites)

  has_enl <- any(d$site == enl_site)
  has_up  <- any(d$site %in% upstream_ent)

  first_ent <- det |> filter(tag == tc, first_dt > anchor,
                              site %in% entiat_sites) |>
    slice_min(first_dt, n = 1)

  first_enl_month <- if (has_enl) {
    month(min(d$first_dt[d$site == enl_site]))
  } else NA_integer_

  first_ent_month <- if (nrow(first_ent) > 0) month(first_ent$first_dt) else NA_integer_

  actual_entry_season <- case_when(
    nrow(first_ent) == 0             ~ "never",
    !is.na(first_enl_month) & first_enl_month %in% 1:6 ~ "spring",
    !is.na(first_enl_month)          ~ "fall",
    !is.na(first_ent_month) & first_ent_month %in% 1:6 ~ "spring",
    TRUE                             ~ "fall"
  )

  data.frame(TagCode = tc, has_enl = has_enl, has_upstream = has_up,
             actual_entry_season = actual_entry_season)
}) |> bind_rows()

# ── Load directionality ────────────────────────────────────────────────────────
dir_df <- read_csv("~/ptagis_data/enl_directionality_ga_enl_only_29fish.csv",
                   show_col_types = FALSE) |>
  select(TagCode, direction_class, n_obs_enl = n_obs)

# Add the one missing spring fish manually
dir_df <- bind_rows(
  dir_df,
  data.frame(TagCode = "3D9.1C2DE55E7C",
             direction_class = "Upper/Upstream only (single obs)",
             n_obs_enl = 1L)
)

# ── Assemble table ─────────────────────────────────────────────────────────────
tbl <- fish_clean |>
  select(TagCode, group, ReturnYear, rrf_anchor, EntiatDetected,
         EntiatSpawned_orig = EntiatSpawned) |>
  left_join(entiat_detail, by = "TagCode") |>
  left_join(fish_r  |> select(TagCode, spring_return_r, EntiatSpawned_r),
            by = "TagCode") |>
  left_join(fish_s  |> select(TagCode, EntiatSpawned_s),
            by = "TagCode") |>
  left_join(dir_df, by = "TagCode") |>
  mutate(
    rrf_anchor_season = "fall",   # all fish anchored Jul–Dec by design
    enl_only = has_enl & !has_upstream,
    entry_season_reclassified = if_else(spring_return_r == 1, "spring", "fall"),
    direction_class = if_else(
      enl_only, direction_class, NA_character_
    )
  ) |>
  select(
    TagCode, group, ReturnYear,
    rrf_anchor_season,
    actual_entry_season,
    has_enl, has_upstream, enl_only,
    enl_direction = direction_class,
    # Outcome under each criterion
    criterion_1_broad     = EntiatDetected,
    criterion_2_orig_strict = EntiatSpawned_orig,
    criterion_3_reclassified = EntiatSpawned_r,
    criterion_4_enl_sens  = EntiatSpawned_s
  ) |>
  arrange(group, actual_entry_season, TagCode)

# ── Group-level summary ────────────────────────────────────────────────────────
cat("=== PER-CRITERION RETURN RATES ===\n\n")

summarise_criterion <- function(col) {
  tbl |> group_by(group) |>
    summarise(n = n(), spawned = sum(.data[[col]], na.rm=TRUE),
              pct = round(mean(.data[[col]], na.rm=TRUE)*100, 1),
              .groups = "drop") |>
    mutate(criterion = col)
}

bind_rows(
  summarise_criterion("criterion_1_broad"),
  summarise_criterion("criterion_2_orig_strict"),
  summarise_criterion("criterion_3_reclassified"),
  summarise_criterion("criterion_4_enl_sens")
) |> select(criterion, group, n, spawned, pct) |> print(n=Inf)

cat("\n=== ENL-ONLY FISH BREAKDOWN ===\n")
tbl |> filter(enl_only) |>
  group_by(group, actual_entry_season) |>
  summarise(n=n(),
            c2=sum(criterion_2_orig_strict),
            c3=sum(criterion_3_reclassified),
            c4=sum(criterion_4_enl_sens),
            .groups="drop") |>
  print()

# ── Save ──────────────────────────────────────────────────────────────────────
out <- "entiat_outcome_criteria_summary.csv"
write_csv(tbl, out)
cat(sprintf("\nSaved: %s  (%d rows, %d cols)\n", out, nrow(tbl), ncol(tbl)))
cat("Columns:", paste(names(tbl), collapse=", "), "\n")
