###############################################################################
# Bayesian Hierarchical Analysis: Does Wells Dam Interaction Reduce
# Probability of Wild Steelhead Returning to the Entiat River?
#
# Data: PIT-tag detections of wild steelhead (Entiat Wild/EntiatWildSteel.csv)
#       Fish that returned above Rocky Reach Dam (RRF adult detections)
#
# Research question: Is a fish's probability of being detected in the Entiat
# River subbasin lower if it had a post-adult-RRF detection at Wells Dam?
#
# Groups:
#   Group A: No Wells Dam detection post-RRF anchor (n ~ 222)
#   Group B: Wells Dam detection post-RRF anchor, NOT last detected in
#            Methow/Okanogan (strays excluded) (n ~ 181)
#
# Models:
#   Model 1 (all fish): EntiatDetected ~ wells_interaction + harvest_open
#                       + (1 | ReturnYear)
#   Model 2 (Group B):  EntiatDetected ~ SpringSpillHours_z + FallSpillHours_z
#                       + harvest_open + (1 | SpillYear)
#
# Priors (matching run_full_hierarchical_analysis_with_harvest.R):
#   Intercept ~ Normal(0, 2)
#   Fixed effects ~ Normal(0, 1)
#   Random effect SD ~ Normal(0, 1) [half-normal via brms default]
#
# MCMC: 4500 iterations, 1500 warmup, 4 chains, 4 cores, adapt_delta=0.95
###############################################################################

# =============================================================================
# PACKAGES
# =============================================================================
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
cat("ENTIAT WILD STEELHEAD - WELLS DAM INTERACTION ANALYSIS\n")
cat("Bayesian Hierarchical Logistic Regression (brms/Stan)\n")
cat("======================================================================\n\n")

# Set working directory to the Entiat Wild subfolder
# All outputs save here; parent-directory files accessed via "../"
setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

# =============================================================================
# SITE CLASSIFICATION CONSTANTS
# =============================================================================

# Entiat River subbasin site prefixes
# All antennas physically within the Entiat River drainage
ENTIAT_PREFIXES <- c("ENL", "ENA", "ENM", "ENS", "ENF", "MAD", "RCT", "TLT",
                     "URT", "HN1", "HN2", "HN3", "UMR", "3D",  "HS1", "HS2",
                     "WL2", "SR1", "SR2", "PD1", "PD2", "TY1", "EHL", "EBO")

# Upstream migration corridor sites used to qualify spring-RRF fish
# (Bonneville, The Dalles, John Day, McNary, Priest Rapids, Rock Island adult ladders)
CORRIDOR_PREFIXES <- c("BO1", "BO3", "BO4", "TD1", "TD2", "JO1", "JO2",
                       "MC1", "MC2", "PRA", "RIA")

# Wells Dam antenna codes
WELLS_PATTERN <- "\\bWEA\\b|\\bWEJ\\b|\\bWEH\\b"

# Rocky Reach adult fishway prefix
RRF_PREFIX <- "RRF"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Check if a site name starts with any of the given prefixes
starts_with_any <- function(site_names, prefixes) {
  pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  grepl(pattern, site_names, perl = TRUE)
}

# =============================================================================
# LOAD DATA
# =============================================================================
cat("--- Loading Data ---\n")

detections <- read_csv("EntiatWildSteel.csv", show_col_types = FALSE)

cat(sprintf("Loaded %d detection rows for %d unique fish\n",
            nrow(detections), n_distinct(detections$`Tag Code`)))
cat("Columns:", paste(names(detections), collapse = ", "), "\n\n")

# Standardize column name (Tag Code -> TagCode for consistency)
detections <- detections %>%
  rename(TagCode = `Tag Code`,
         SiteName = `Site Name`,
         FirstTime = `First Time Value`,
         LastTime  = `Last Time Value`,
         SubbasinName = `Site Subbasin Name`,
         SubbasinCode = `Site Subbasin Code`)

# Parse datetime columns
detections <- detections %>%
  mutate(
    FirstTime = mdy_hm(FirstTime, quiet = TRUE),
    LastTime  = mdy_hm(LastTime,  quiet = TRUE),
    Month = month(FirstTime),
    Year  = year(FirstTime)
  )

# =============================================================================
# CLASSIFY DETECTION SITES
# =============================================================================
cat("--- Classifying Detection Sites ---\n")

detections <- detections %>%
  mutate(
    is_entiat_river  = starts_with_any(SiteName, ENTIAT_PREFIXES),
    is_enl           = str_starts(SiteName, "ENL"),
    is_wells         = grepl(WELLS_PATTERN, SiteName),
    is_upstream_trib = SubbasinName %in% c("Methow", "Okanogan"),
    is_rrf           = str_starts(SiteName, RRF_PREFIX),
    is_corridor      = starts_with_any(SiteName, CORRIDOR_PREFIXES)
  )

cat(sprintf("Detections by type:\n"))
cat(sprintf("  Entiat River sites:    %d rows (%d fish)\n",
            sum(detections$is_entiat_river),
            n_distinct(detections$TagCode[detections$is_entiat_river])))
cat(sprintf("  ENL (lower Entiat):    %d rows (%d fish)\n",
            sum(detections$is_enl),
            n_distinct(detections$TagCode[detections$is_enl])))
cat(sprintf("  Wells Dam:             %d rows (%d fish)\n",
            sum(detections$is_wells),
            n_distinct(detections$TagCode[detections$is_wells])))
cat(sprintf("  Upstream tribs:        %d rows (%d fish)\n",
            sum(detections$is_upstream_trib),
            n_distinct(detections$TagCode[detections$is_upstream_trib])))
cat(sprintf("  Rocky Reach (RRF):     %d rows (%d fish)\n",
            sum(detections$is_rrf),
            n_distinct(detections$TagCode[detections$is_rrf])))
cat(sprintf("  Corridor sites:        %d rows (%d fish)\n",
            sum(detections$is_corridor),
            n_distinct(detections$TagCode[detections$is_corridor])))

# =============================================================================
# ESTABLISH ADULT RRF ANCHOR DATE
# =============================================================================
# Adult fish returning upstream pass Rocky Reach in Jul-Dec.
# The anchor date is the FIRST Jul-Dec RRF detection per fish.
# This filters out juvenile outmigrant detections at RRF (spring/early summer).
# =============================================================================
cat("\n--- Establishing Adult RRF Anchor Dates ---\n")

rrf_adult <- detections %>%
  filter(is_rrf, Month %in% 7:12)

anchor_df <- rrf_adult %>%
  group_by(TagCode) %>%
  summarise(rrf_anchor = min(FirstTime), .groups = "drop") %>%
  mutate(season_override = NA_character_)

cat(sprintf("Fish with Jul-Dec RRF detection (standard adult return): %d\n",
            nrow(anchor_df)))

# Identify fish with ONLY non-Jul-Dec RRF detections (spring/June return)
all_fish <- unique(detections$TagCode)
fish_with_any_rrf <- unique(detections$TagCode[detections$is_rrf])
fish_with_julaug_rrf <- unique(anchor_df$TagCode)
spring_only_rrf_fish <- setdiff(fish_with_any_rrf, fish_with_julaug_rrf)

cat(sprintf("Fish with RRF detections but none in Jul-Dec (spring/June only): %d\n",
            length(spring_only_rrf_fish)))

# =============================================================================
# QUALIFY SPRING-RRF FISH
# =============================================================================
# Spring-RRF fish (Jan-May only at RRF) are included IF they have at least
# one detection at an upstream migration corridor site (McNary/Priest Rapids/
# Rock Island) BEFORE their spring RRF detection. This confirms they are
# genuine adult returns completing migration in spring, not juveniles or
# other anomalous detections.
# =============================================================================
cat("\n--- Qualifying Spring-RRF Fish ---\n")

spring_rrf_detections <- detections %>%
  filter(TagCode %in% spring_only_rrf_fish, is_rrf)

spring_qualified <- list()

for (tag in spring_only_rrf_fish) {
  fish_dets <- detections %>% filter(TagCode == tag)

  # Get earliest spring RRF detection
  spring_rrf_row <- fish_dets %>%
    filter(is_rrf, Month %in% 1:5) %>%
    arrange(FirstTime) %>%
    slice(1)

  if (nrow(spring_rrf_row) == 0) next

  spring_rrf_time <- spring_rrf_row$FirstTime[1]

  # Check for corridor detection BEFORE spring RRF
  corridor_before <- fish_dets %>%
    filter(is_corridor, FirstTime < spring_rrf_time)

  if (nrow(corridor_before) > 0) {
    spring_qualified[[tag]] <- tibble(
      TagCode         = tag,
      rrf_anchor      = spring_rrf_time,
      season_override = NA_character_
    )
    cat(sprintf("  QUALIFY: %s — spring RRF %s, corridor site: %s (%s)\n",
                tag,
                format(spring_rrf_time, "%Y-%m-%d"),
                corridor_before$SiteName[1],
                format(corridor_before$FirstTime[1], "%Y-%m-%d")))
  } else {
    cat(sprintf("  EXCLUDE: %s — spring RRF %s, no corridor evidence\n",
                tag,
                format(spring_rrf_time, "%Y-%m-%d")))
  }
}

spring_anchor_df <- bind_rows(spring_qualified)
cat(sprintf("\nSpring-RRF fish qualified: %d / %d\n",
            nrow(spring_anchor_df), length(spring_only_rrf_fish)))

# =============================================================================
# QUALIFY JUNE-ONLY RRF FISH
# =============================================================================
# Fish with RRF detections only in June (month 6) fall outside both the
# standard Jul-Dec adult window and the Jan-May spring window. Season is
# assigned based on Priest Rapids (PRA) detection timing, which distinguishes
# late spawners (overwintered downstream) from new fall migrants:
#   - PRD detection prev fall (Aug-Dec year-1) or Jan-May same year -> Spring
#     (late spawner: passed PRD early, overwintered, arrived RRF in June)
#   - PRD detection in June same year                               -> Fall
#     (fall migrant: began upstream migration in June, spawns next year)
# =============================================================================
cat("\n--- Qualifying June-Only RRF Fish ---\n")

june_only_rrf_fish <- spring_only_rrf_fish[
  vapply(spring_only_rrf_fish, function(tag) {
    rrf_months <- detections$Month[detections$TagCode == tag & detections$is_rrf]
    all(rrf_months == 6)
  }, logical(1))
]
cat(sprintf("June-only RRF fish to evaluate: %d\n", length(june_only_rrf_fish)))

june_qualified <- list()

for (tag in june_only_rrf_fish) {
  fish_dets <- detections %>% filter(TagCode == tag)

  june_rrf_row <- fish_dets %>%
    filter(is_rrf, Month == 6) %>%
    arrange(FirstTime) %>%
    slice(1)

  june_rrf_time <- june_rrf_row$FirstTime[1]
  june_year     <- year(june_rrf_time)

  # Downstream detections strictly before June RRF
  june_year <- year(june_rrf_time)

  # Spring (late spawner): PRD detection prev fall (Aug-Dec year-1) or Jan-May same year
  spring_dams <- fish_dets %>%
    filter(str_starts(SiteName, "PRA"), FirstTime < june_rrf_time) %>%
    filter((Year == june_year     & Month %in% 1:5) |
           (Year == june_year - 1 & Month >= 8))

  # Fall (next-year spawner): PRD detection in June same year before RRF
  fall_dams <- fish_dets %>%
    filter(str_starts(SiteName, "PRA"), Year == june_year, Month == 6,
           FirstTime < june_rrf_time)

  if (nrow(spring_dams) == 0 && nrow(fall_dams) == 0) {
    cat(sprintf("  EXCLUDE: %s — June RRF %s, no downstream dam evidence\n",
                tag, format(june_rrf_time, "%Y-%m-%d")))
    next
  }

  # Spring takes priority (early PRD passage is stronger evidence of late-spawner identity)
  if (nrow(spring_dams) > 0) {
    season    <- "Spring"
    ref_dam   <- spring_dams %>% arrange(FirstTime) %>% slice(1)
  } else {
    season    <- "Fall"
    ref_dam   <- fall_dams   %>% arrange(FirstTime) %>% slice(1)
  }

  june_qualified[[tag]] <- tibble(
    TagCode         = tag,
    rrf_anchor      = june_rrf_time,
    season_override = season
  )
  cat(sprintf("  QUALIFY (%s): %s — June RRF %s, ref dam: %s (%s)\n",
              season, tag,
              format(june_rrf_time, "%Y-%m-%d"),
              ref_dam$SiteName[1],
              format(ref_dam$FirstTime[1], "%Y-%m-%d")))
}

june_anchor_df <- bind_rows(june_qualified)
cat(sprintf("\nJune-only fish qualified: %d / %d\n",
            nrow(june_anchor_df), length(june_only_rrf_fish)))

# Combine standard, spring-qualified, and June-qualified anchors
anchor_df <- bind_rows(anchor_df, spring_anchor_df, june_anchor_df)
cat(sprintf("Total fish with valid adult RRF anchor: %d\n", nrow(anchor_df)))

# =============================================================================
# FILTER TO POST-ANCHOR DETECTIONS
# =============================================================================
cat("\n--- Filtering to Post-Anchor Detections ---\n")

detections_with_anchor <- detections %>%
  inner_join(anchor_df, by = "TagCode") %>%
  filter(FirstTime >= rrf_anchor)

n_fish_anchor <- n_distinct(detections_with_anchor$TagCode)
cat(sprintf("Fish retained after anchor filter: %d\n", n_fish_anchor))

# Extract return year from anchor date
# Return year = calendar year of the anchor detection
anchor_year_df <- anchor_df %>%
  mutate(ReturnYear = year(rrf_anchor))

cat("\nReturn year distribution:\n")
print(table(anchor_year_df$ReturnYear))

# =============================================================================
# BUILD PER-FISH DETECTION FLAGS
# =============================================================================
cat("\n--- Building Per-Fish Detection Flags ---\n")

fish_flags <- detections_with_anchor %>%
  group_by(TagCode) %>%
  summarise(
    # Did the fish have any post-anchor Wells Dam detection?
    wells_detected    = any(is_wells),

    # Did the fish have any post-anchor Entiat River detection?
    entiat_detected   = any(is_entiat_river),

    # Was the fish last detected in an upstream tributary (stray)?
    upstream_detected = any(is_upstream_trib),
    last_overall      = max(FirstTime),
    last_upstream     = if (any(is_upstream_trib)) max(FirstTime[is_upstream_trib]) else as.POSIXct(NA),

    .groups = "drop"
  ) %>%
  mutate(
    # Stray flag: fish with Wells detection whose last detection is in Methow/Okanogan
    last_detected_upstream = !is.na(last_upstream) & (last_upstream == last_overall) & upstream_detected
  )

# Manual correction: 3D9.1C2C451B01 has a single Entiat detection (ENS 2017-05-24)
# that is a confirmed data error — the fish's only return was in 2009. Correct to 0.
fish_flags <- fish_flags %>%
  mutate(entiat_detected = if_else(TagCode == "3D9.1C2C451B01", FALSE, entiat_detected))

# Join with return year
fish_flags <- fish_flags %>%
  left_join(anchor_year_df %>% select(TagCode, ReturnYear, rrf_anchor), by = "TagCode")

cat(sprintf("Total fish: %d\n", nrow(fish_flags)))
cat(sprintf("  Wells-detected:         %d\n", sum(fish_flags$wells_detected)))
cat(sprintf("  Last in upstream trib:  %d (strays to exclude from Group B)\n",
            sum(fish_flags$last_detected_upstream, na.rm = TRUE)))
cat(sprintf("  Entiat-detected:        %d\n", sum(fish_flags$entiat_detected)))

# =============================================================================
# ASSIGN GROUPS
# =============================================================================
# Group A: No Wells Dam detection post-RRF anchor
# Group B: Wells Dam detection post-RRF anchor, NOT last in Methow/Okanogan
# Excluded: Last detected in upstream tributary (strays)
# =============================================================================
cat("\n--- Assigning Groups ---\n")

fish_flags <- fish_flags %>%
  mutate(
    group = case_when(
      !wells_detected                                         ~ "A",
      wells_detected & !last_detected_upstream               ~ "B",
      TRUE                                                    ~ "Stray"
    )
  )

cat("Group assignments:\n")
print(table(fish_flags$group))

fish_AB <- fish_flags %>% filter(group %in% c("A", "B"))

cat(sprintf("\nAnalysis dataset: %d fish (Group A=%d, Group B=%d)\n",
            nrow(fish_AB),
            sum(fish_AB$group == "A"),
            sum(fish_AB$group == "B")))

# Detection rates by group
cat("\nEntiat detection rate by group:\n")
fish_AB %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_detected = sum(entiat_detected),
    pct_detected = round(mean(entiat_detected) * 100, 1),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# ADD SEASON COVARIATE
# =============================================================================
# Classify each fish's return season based on the month of its anchor date.
# "Spring" = Jan-May (lower ENL detection efficiency)
# "Fall"   = Jul-Dec (higher ENL detection efficiency)
# June-only fish have an explicit season_override (Spring or Fall) set during
# qualification based on Priest Rapids detection timing — see above.
# =============================================================================
fish_AB <- fish_AB %>%
  left_join(anchor_df %>% select(TagCode, season_override), by = "TagCode") %>%
  mutate(
    ReturnMonth  = month(rrf_anchor),
    ReturnSeason = case_when(
      !is.na(season_override)    ~ season_override,
      ReturnMonth %in% 1:5       ~ "Spring",
      TRUE                       ~ "Fall"
    ),
    wells_interaction = as.integer(group == "B"),
    EntiatDetected    = as.integer(entiat_detected)
  ) %>%
  select(-season_override)

cat("\nReturn season by group:\n")
fish_AB %>%
  group_by(group, ReturnSeason) %>%
  summarise(n = n(), pct_entiat = round(mean(entiat_detected)*100,1), .groups="drop") %>%
  print()

# =============================================================================
# MERGE HARVEST STATUS
# =============================================================================
cat("\n--- Merging Harvest Status ---\n")

harvest_df <- read_csv("../Harvest_Open_Closed.csv", show_col_types = FALSE)

fish_AB <- fish_AB %>%
  left_join(harvest_df, by = c("ReturnYear" = "Year")) %>%
  mutate(harvest_open = as.integer(Harvest_Status == "Open"))

cat("Harvest status distribution:\n")
print(table(Harvest = fish_AB$Harvest_Status, Group = fish_AB$group, useNA = "ifany"))

# =============================================================================
# MERGE SPILL DATA FOR GROUP B
# =============================================================================
# For Model 2 (Group B only), merge spill metrics from the existing
# analysis_fish_with_harvest_R.csv, which contains spill data for all
# Wells-interacting fish (2007-2021). The join key is TagCode.
# For Group B fish with no spill data (outside the spill record years),
# those fish will be flagged and reported.
# =============================================================================
cat("\n--- Merging Spill Data for Group B ---\n")

# First try to load the pre-computed spill data from the existing analysis
spill_fish_file <- "../analysis_fish_with_harvest_R.csv"
if (file.exists(spill_fish_file)) {
  spill_fish_df <- read_csv(spill_fish_file, show_col_types = FALSE) %>%
    select(TagCode, WellsYear, FallSpillHours, FallTotalSpill,
           SpringSpillHours, SpringTotalSpill) %>%
    distinct(TagCode, .keep_all = TRUE)

  cat(sprintf("Loaded pre-computed spill data for %d fish\n", nrow(spill_fish_df)))

  group_b_df <- fish_AB %>%
    filter(group == "B") %>%
    left_join(spill_fish_df, by = "TagCode")

  # --- Year-level fallback for fish missing from spill file ---
  # Covers: (1) fish re-assigned to Group B after reanalysis that weren't in
  # the original spill file, and (2) return years beyond the spill file's range.
  if (any(is.na(group_b_df$FallSpillHours))) {

    # Year-level spill from existing file (same value for all fish in a year)
    yr_spill <- spill_fish_df %>%
      distinct(WellsYear, .keep_all = TRUE) %>%
      select(WellsYear, FallSpillHours, FallTotalSpill,
             SpringSpillHours, SpringTotalSpill)

    # Supplement with recent years from raw spill files
    spill_new_file <- "../WellsSpill2024-2025.csv"
    spill_old_file <- "../merged_spill_2000_2024.csv"
    spill_col      <- "SPILL (KCFS)"

    if (file.exists(spill_new_file)) {
      new_raw <- read_csv(spill_new_file, show_col_types = FALSE) %>%
        mutate(Date = as.Date(Date), Month = month(Date), Year = year(Date))

      new_fall <- new_raw %>%
        filter(Month >= 8, Month <= 11) %>%
        group_by(Year) %>%
        summarise(FallSpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
                  FallTotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(WellsYear = Year + 1) %>% select(-Year)

      new_spring <- new_raw %>%
        filter(Month >= 1, Month <= 3) %>%
        group_by(Year) %>%
        summarise(SpringSpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
                  SpringTotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
                  .groups = "drop") %>%
        rename(WellsYear = Year)

      new_yr <- full_join(new_fall, new_spring, by = "WellsYear")

      # Fill WellsYear 2024 fall spill from merged historical file if available
      if (file.exists(spill_old_file) && any(is.na(new_yr$FallSpillHours))) {
        old_raw <- suppressWarnings(
          read_csv(spill_old_file, show_col_types = FALSE) %>%
            mutate(DateTime = parse_date_time(DateTime,
                     orders = c("ymd HM", "ymd HMS", "mdy HM", "mdy HMS")),
                   Year = year(DateTime), Month = month(DateTime))
        )
        old_fall <- old_raw %>%
          filter(Month >= 8, Month <= 11) %>%
          group_by(Year) %>%
          summarise(FallSpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
                    FallTotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
                    .groups = "drop") %>%
          mutate(WellsYear = Year + 1) %>% select(-Year)

        new_yr <- new_yr %>%
          left_join(old_fall %>% rename(FallSpillHours_old = FallSpillHours,
                                        FallTotalSpill_old  = FallTotalSpill),
                    by = "WellsYear") %>%
          mutate(
            FallSpillHours = coalesce(FallSpillHours, FallSpillHours_old),
            FallTotalSpill = coalesce(FallTotalSpill, FallTotalSpill_old)
          ) %>%
          select(-FallSpillHours_old, -FallTotalSpill_old)
      }

      yr_spill <- bind_rows(yr_spill, new_yr) %>%
        distinct(WellsYear, .keep_all = TRUE)
    }

    # Fill missing fish using year-level lookup (WellsYear = ReturnYear for fall fish)
    missing_idx <- is.na(group_b_df$FallSpillHours)
    group_b_df <- group_b_df %>%
      left_join(yr_spill %>% rename_with(~ paste0(.x, "_yr"), -WellsYear),
                by = c("ReturnYear" = "WellsYear")) %>%
      mutate(
        FallSpillHours   = coalesce(FallSpillHours,   FallSpillHours_yr),
        FallTotalSpill   = coalesce(FallTotalSpill,   FallTotalSpill_yr),
        SpringSpillHours = coalesce(SpringSpillHours, SpringSpillHours_yr),
        SpringTotalSpill = coalesce(SpringTotalSpill, SpringTotalSpill_yr)
      ) %>%
      select(-ends_with("_yr"))

    cat(sprintf("Year-level fallback filled %d fish\n", sum(missing_idx)))
  }

  n_with_spill <- sum(!is.na(group_b_df$FallSpillHours))
  n_missing    <- sum(is.na(group_b_df$FallSpillHours))
  cat(sprintf("Group B fish with spill data:    %d\n", n_with_spill))
  cat(sprintf("Group B fish without spill data: %d\n", n_missing))

} else {
  # Compute spill metrics from raw spill data
  cat("Pre-computed spill file not found; computing from raw spill data...\n")

  spill_raw <- read_csv("../merged_spill_2000_2024.csv", show_col_types = FALSE) %>%
    mutate(DateTime = parse_date_time(DateTime, orders = c("ymd HM", "ymd HMS",
                                                            "mdy HM", "mdy HMS")),
           Year  = year(DateTime),
           Month = month(DateTime),
           Day   = day(DateTime))

  spill_col <- "SPILL (KCFS)"

  fall_spill_yr <- spill_raw %>%
    filter(Month >= 8, Month <= 11) %>%
    group_by(Year) %>%
    summarise(FallSpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
              FallTotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
              .groups = "drop") %>%
    mutate(SpillYear = Year + 1)

  spring_spill_yr <- spill_raw %>%
    filter(Month >= 1, Month <= 3) %>%
    group_by(Year) %>%
    summarise(SpringSpillHours = sum(.data[[spill_col]] > 0, na.rm = TRUE),
              SpringTotalSpill = sum(.data[[spill_col]], na.rm = TRUE),
              .groups = "drop") %>%
    rename(SpillYear = Year)

  spill_metrics <- full_join(fall_spill_yr %>% select(-Year),
                             spring_spill_yr,
                             by = "SpillYear")

  group_b_df <- fish_AB %>%
    filter(group == "B") %>%
    left_join(spill_metrics, by = c("ReturnYear" = "SpillYear"))

  n_with_spill <- sum(!is.na(group_b_df$FallSpillHours))
  n_missing    <- sum(is.na(group_b_df$FallSpillHours))
  cat(sprintf("Group B fish with spill data:    %d\n", n_with_spill))
  cat(sprintf("Group B fish without spill data: %d\n", n_missing))
}

# =============================================================================
# SAVE ANALYSIS DATASETS
# =============================================================================
cat("\n--- Saving Analysis Datasets ---\n")

write_csv(fish_AB, "entiat_analysis_all_fish.csv")
cat("Saved: entiat_analysis_all_fish.csv\n")

if (exists("group_b_df")) {
  write_csv(group_b_df, "entiat_analysis_groupB.csv")
  cat("Saved: entiat_analysis_groupB.csv\n")
}

# =============================================================================
# DEFINE SHARED PRIORS AND MCMC SETTINGS
# =============================================================================
# Match the priors and MCMC settings from run_full_hierarchical_analysis_with_harvest.R
model_priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 1), class = "sd")
)

MCMC_ITER    <- 4500
MCMC_WARMUP  <- 1500
MCMC_CHAINS  <- 4
MCMC_CORES   <- 4
MCMC_SEED    <- 42
ADAPT_DELTA  <- 0.95

# =============================================================================
# MODEL 1: ALL FISH — Wells Interaction + Harvest + Year Random Effect
# =============================================================================
# Outcome:  EntiatDetected (1 = detected in Entiat River post-RRF, 0 = not)
# Fixed:    wells_interaction (Group B vs Group A)
#           harvest_open (Open vs Closed harvest year)
# Random:   (1 | ReturnYear) — year-level baseline variation
# =============================================================================
cat("\n======================================================================\n")
cat("MODEL 1: All Fish — Wells Interaction + Harvest + Year RE\n")
cat(sprintf("n = %d fish across %d return years\n",
            nrow(fish_AB), n_distinct(fish_AB$ReturnYear)))
cat("======================================================================\n")

fish_AB$ReturnYearFactor <- factor(fish_AB$ReturnYear)

fit_model1 <- brm(
  EntiatDetected ~ wells_interaction + harvest_open + (1 | ReturnYearFactor),
  data   = fish_AB,
  family = bernoulli(link = "logit"),
  prior  = model_priors,
  iter    = MCMC_ITER,
  warmup  = MCMC_WARMUP,
  chains  = MCMC_CHAINS,
  cores   = MCMC_CORES,
  seed    = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2,
  refresh = 0
)

cat("\nModel 1 Summary:\n")
print(summary(fit_model1))

# Extract posteriors
post1 <- as_draws_df(fit_model1)
beta_wells_samples   <- post1$b_wells_interaction
beta_harvest_samples <- post1$b_harvest_open
sigma_yr_samples     <- post1$`sd_ReturnYearFactor__Intercept`

# Compute summaries
wells_hdi   <- quantile(beta_wells_samples,   probs = c(0.025, 0.975))
harvest_hdi <- quantile(beta_harvest_samples, probs = c(0.025, 0.975))

if (requireNamespace("bayestestR", quietly = TRUE)) {
  wells_hdi_obj   <- bayestestR::hdi(beta_wells_samples,   ci = 0.95)
  wells_hdi       <- c(wells_hdi_obj$CI_low,   wells_hdi_obj$CI_high)
  harvest_hdi_obj <- bayestestR::hdi(beta_harvest_samples, ci = 0.95)
  harvest_hdi     <- c(harvest_hdi_obj$CI_low, harvest_hdi_obj$CI_high)
}

wells_or   <- exp(mean(beta_wells_samples))
harvest_or <- exp(mean(beta_harvest_samples))

cat(sprintf("\n--- Model 1 Results ---\n"))
cat(sprintf("Wells Interaction:\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n", mean(beta_wells_samples), wells_or))
cat(sprintf("  95%% HDI   [%.3f, %.3f]\n", wells_hdi[1], wells_hdi[2]))
cat(sprintf("  P(beta > 0) = %.4f  [i.e. P(Wells INCREASES return)]\n",
            mean(beta_wells_samples > 0)))
cat(sprintf("  P(beta < 0) = %.4f  [i.e. P(Wells DECREASES return)]\n",
            mean(beta_wells_samples < 0)))

cat(sprintf("\nHarvest Status (Open vs Closed):\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n", mean(beta_harvest_samples), harvest_or))
cat(sprintf("  95%% HDI   [%.3f, %.3f]\n", harvest_hdi[1], harvest_hdi[2]))
cat(sprintf("  P(beta > 0) = %.4f\n", mean(beta_harvest_samples > 0)))

cat(sprintf("\nYear-level variability:\n"))
cat(sprintf("  sigma_year = %.3f (mean)\n", mean(sigma_yr_samples)))

# Model diagnostics
diag1 <- summary(fit_model1)$fixed
cat(sprintf("\nDiagnostics (fixed effects):\n"))
cat(sprintf("  Rhat range:    [%.3f, %.3f]\n", min(diag1$Rhat), max(diag1$Rhat)))
cat(sprintf("  Bulk ESS range: [%.0f, %.0f]\n", min(diag1$Bulk_ESS), max(diag1$Bulk_ESS)))

# =============================================================================
# MODEL 1b: ALL FISH + SEASON COVARIATE
# =============================================================================
# Adds ReturnSeason to check whether spring return timing (low ENL detection)
# confounds the Wells Dam effect estimate.
# =============================================================================
cat("\n======================================================================\n")
cat("MODEL 1b: All Fish — Wells + Harvest + Season + Year RE\n")
cat("======================================================================\n")

# Encode season as binary (Spring = 1)
fish_AB <- fish_AB %>%
  mutate(spring_return = as.integer(ReturnSeason == "Spring"))

fit_model1b <- brm(
  EntiatDetected ~ wells_interaction + harvest_open + spring_return
                 + (1 | ReturnYearFactor),
  data   = fish_AB,
  family = bernoulli(link = "logit"),
  prior  = model_priors,
  iter    = MCMC_ITER,
  warmup  = MCMC_WARMUP,
  chains  = MCMC_CHAINS,
  cores   = MCMC_CORES,
  seed    = MCMC_SEED,
  control = list(adapt_delta = ADAPT_DELTA),
  silent  = 2,
  refresh = 0
)

cat("\nModel 1b Summary:\n")
print(summary(fit_model1b))

post1b <- as_draws_df(fit_model1b)
beta_wells_1b   <- post1b$b_wells_interaction
beta_season_1b  <- post1b$b_spring_return

wells_hdi_1b  <- quantile(beta_wells_1b,  probs = c(0.025, 0.975))
season_hdi_1b <- quantile(beta_season_1b, probs = c(0.025, 0.975))

if (requireNamespace("bayestestR", quietly = TRUE)) {
  w_obj <- bayestestR::hdi(beta_wells_1b,  ci = 0.95)
  wells_hdi_1b <- c(w_obj$CI_low, w_obj$CI_high)
  s_obj <- bayestestR::hdi(beta_season_1b, ci = 0.95)
  season_hdi_1b <- c(s_obj$CI_low, s_obj$CI_high)
}

cat(sprintf("\n--- Model 1b Results ---\n"))
cat(sprintf("Wells Interaction (with season covariate):\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n",
            mean(beta_wells_1b), exp(mean(beta_wells_1b))))
cat(sprintf("  95%% HDI   [%.3f, %.3f]\n", wells_hdi_1b[1], wells_hdi_1b[2]))
cat(sprintf("  P(beta < 0) = %.4f\n", mean(beta_wells_1b < 0)))

cat(sprintf("\nSpring Return Season:\n"))
cat(sprintf("  beta_mean = %.3f  (OR = %.3f)\n",
            mean(beta_season_1b), exp(mean(beta_season_1b))))
cat(sprintf("  95%% HDI   [%.3f, %.3f]\n", season_hdi_1b[1], season_hdi_1b[2]))
cat(sprintf("  P(beta < 0) = %.4f\n", mean(beta_season_1b < 0)))

# =============================================================================
# MODEL 2: GROUP B ONLY — Spill + Harvest + Year Random Effect
# =============================================================================
# Tests whether Wells Dam spill conditions during the period of interaction
# predict whether Group B fish subsequently return to the Entiat River.
# Uses Spring (Jan-Mar) and Fall (Aug-Nov) spill at Wells Dam.
# =============================================================================
fit_model2    <- NULL
group_b_spill <- NULL

if (exists("group_b_df") && exists("n_with_spill") && n_with_spill >= 20) {

  cat("\n======================================================================\n")
  cat("MODEL 2: Group B Only — Spill + Harvest + Year RE\n")
  cat(sprintf("n = %d fish with spill data\n", n_with_spill))
  cat("======================================================================\n")

  group_b_spill <- group_b_df %>%
    filter(!is.na(FallSpillHours), !is.na(SpringSpillHours)) %>%
    mutate(
      ReturnYearFactor = factor(ReturnYear),
      # Standardize spill metrics to mean=0, sd=1
      FallSpillHours_z   = scale(FallSpillHours)[,1],
      FallTotalSpill_z   = scale(FallTotalSpill)[,1],
      SpringSpillHours_z = scale(SpringSpillHours)[,1],
      SpringTotalSpill_z = scale(SpringTotalSpill)[,1]
    )

  cat(sprintf("Group B fish with complete spill data: %d\n", nrow(group_b_spill)))
  cat(sprintf("Return years covered: %s\n",
              paste(sort(unique(group_b_spill$ReturnYear)), collapse=", ")))
  cat(sprintf("Entiat detection rate (Group B with spill): %.1f%%\n",
              mean(group_b_spill$EntiatDetected) * 100))

  # Model 2: Spring + Fall spill hours + harvest
  fit_model2 <- brm(
    EntiatDetected ~ SpringSpillHours_z + FallSpillHours_z + harvest_open
                   + (1 | ReturnYearFactor),
    data   = group_b_spill,
    family = bernoulli(link = "logit"),
    prior  = model_priors,
    iter    = MCMC_ITER,
    warmup  = MCMC_WARMUP,
    chains  = MCMC_CHAINS,
    cores   = MCMC_CORES,
    seed    = MCMC_SEED,
    control = list(adapt_delta = ADAPT_DELTA),
    silent  = 2,
    refresh = 0
  )

  cat("\nModel 2 Summary:\n")
  print(summary(fit_model2))

  post2 <- as_draws_df(fit_model2)
  b_spring2  <- post2$b_SpringSpillHours_z
  b_fall2    <- post2$b_FallSpillHours_z
  b_harv2    <- post2$b_harvest_open
  sigma_yr2  <- post2$`sd_ReturnYearFactor__Intercept`

  spring_hdi2  <- quantile(b_spring2, probs = c(0.025, 0.975))
  fall_hdi2    <- quantile(b_fall2,   probs = c(0.025, 0.975))
  harv_hdi2    <- quantile(b_harv2,   probs = c(0.025, 0.975))

  if (requireNamespace("bayestestR", quietly = TRUE)) {
    sp_obj <- bayestestR::hdi(b_spring2, ci = 0.95)
    spring_hdi2 <- c(sp_obj$CI_low, sp_obj$CI_high)
    fa_obj <- bayestestR::hdi(b_fall2,   ci = 0.95)
    fall_hdi2 <- c(fa_obj$CI_low, fa_obj$CI_high)
    hv_obj <- bayestestR::hdi(b_harv2,   ci = 0.95)
    harv_hdi2 <- c(hv_obj$CI_low, hv_obj$CI_high)
  }

  cat(sprintf("\n--- Model 2 Results ---\n"))
  cat(sprintf("Spring Spill Hours (std):\n"))
  cat(sprintf("  beta = %.3f (OR=%.3f), 95%% HDI [%.3f, %.3f], P(>0)=%.3f\n",
              mean(b_spring2), exp(mean(b_spring2)),
              spring_hdi2[1], spring_hdi2[2], mean(b_spring2 > 0)))
  cat(sprintf("Fall Spill Hours (std):\n"))
  cat(sprintf("  beta = %.3f (OR=%.3f), 95%% HDI [%.3f, %.3f], P(>0)=%.3f\n",
              mean(b_fall2), exp(mean(b_fall2)),
              fall_hdi2[1], fall_hdi2[2], mean(b_fall2 > 0)))
  cat(sprintf("Harvest Open:\n"))
  cat(sprintf("  beta = %.3f (OR=%.3f), 95%% HDI [%.3f, %.3f], P(>0)=%.3f\n",
              mean(b_harv2), exp(mean(b_harv2)),
              harv_hdi2[1], harv_hdi2[2], mean(b_harv2 > 0)))
  cat(sprintf("Year variability: sigma_year = %.3f\n", mean(sigma_yr2)))

} else {
  cat("\nSkipping Model 2: insufficient Group B fish with spill data.\n")
  fit_model2     <- NULL
  group_b_spill  <- NULL
}

# =============================================================================
# SAVE RESULTS SUMMARY CSV
# =============================================================================
cat("\n--- Saving Results Summary ---\n")

results_summary <- bind_rows(
  tibble(
    Model       = "Model1_AllFish",
    Parameter   = "wells_interaction",
    Beta_Mean   = mean(beta_wells_samples),
    Beta_SD     = sd(beta_wells_samples),
    HDI_Lower   = wells_hdi[1],
    HDI_Upper   = wells_hdi[2],
    OR          = exp(mean(beta_wells_samples)),
    P_positive  = mean(beta_wells_samples > 0),
    P_negative  = mean(beta_wells_samples < 0)
  ),
  tibble(
    Model       = "Model1_AllFish",
    Parameter   = "harvest_open",
    Beta_Mean   = mean(beta_harvest_samples),
    Beta_SD     = sd(beta_harvest_samples),
    HDI_Lower   = harvest_hdi[1],
    HDI_Upper   = harvest_hdi[2],
    OR          = exp(mean(beta_harvest_samples)),
    P_positive  = mean(beta_harvest_samples > 0),
    P_negative  = mean(beta_harvest_samples < 0)
  ),
  tibble(
    Model       = "Model1b_AllFish_Season",
    Parameter   = "wells_interaction",
    Beta_Mean   = mean(beta_wells_1b),
    Beta_SD     = sd(beta_wells_1b),
    HDI_Lower   = wells_hdi_1b[1],
    HDI_Upper   = wells_hdi_1b[2],
    OR          = exp(mean(beta_wells_1b)),
    P_positive  = mean(beta_wells_1b > 0),
    P_negative  = mean(beta_wells_1b < 0)
  ),
  tibble(
    Model       = "Model1b_AllFish_Season",
    Parameter   = "spring_return",
    Beta_Mean   = mean(beta_season_1b),
    Beta_SD     = sd(beta_season_1b),
    HDI_Lower   = season_hdi_1b[1],
    HDI_Upper   = season_hdi_1b[2],
    OR          = exp(mean(beta_season_1b)),
    P_positive  = mean(beta_season_1b > 0),
    P_negative  = mean(beta_season_1b < 0)
  )
)

# Add Model 2 if it was fitted
if (!is.null(fit_model2)) {
  results_summary <- bind_rows(
    results_summary,
    tibble(
      Model       = "Model2_GroupB_Spill",
      Parameter   = "SpringSpillHours_z",
      Beta_Mean   = mean(b_spring2),
      Beta_SD     = sd(b_spring2),
      HDI_Lower   = spring_hdi2[1],
      HDI_Upper   = spring_hdi2[2],
      OR          = exp(mean(b_spring2)),
      P_positive  = mean(b_spring2 > 0),
      P_negative  = mean(b_spring2 < 0)
    ),
    tibble(
      Model       = "Model2_GroupB_Spill",
      Parameter   = "FallSpillHours_z",
      Beta_Mean   = mean(b_fall2),
      Beta_SD     = sd(b_fall2),
      HDI_Lower   = fall_hdi2[1],
      HDI_Upper   = fall_hdi2[2],
      OR          = exp(mean(b_fall2)),
      P_positive  = mean(b_fall2 > 0),
      P_negative  = mean(b_fall2 < 0)
    ),
    tibble(
      Model       = "Model2_GroupB_Spill",
      Parameter   = "harvest_open",
      Beta_Mean   = mean(b_harv2),
      Beta_SD     = sd(b_harv2),
      HDI_Lower   = harv_hdi2[1],
      HDI_Upper   = harv_hdi2[2],
      OR          = exp(mean(b_harv2)),
      P_positive  = mean(b_harv2 > 0),
      P_negative  = mean(b_harv2 < 0)
    )
  )
}

write_csv(results_summary, "entiat_return_bayesian_results.csv")
cat("Saved: entiat_return_bayesian_results.csv\n")

# Year summary for plotting
year_summary_m1 <- fish_AB %>%
  group_by(ReturnYear, Harvest_Status) %>%
  summarise(
    n_total       = n(),
    n_groupA      = sum(group == "A"),
    n_groupB      = sum(group == "B"),
    n_detected    = sum(EntiatDetected),
    pct_detected  = mean(EntiatDetected) * 100,
    n_wells       = sum(wells_interaction),
    n_spring      = sum(ReturnSeason == "Spring"),
    .groups = "drop"
  )

write_csv(year_summary_m1, "entiat_year_summary.csv")
cat("Saved: entiat_year_summary.csv\n")

# =============================================================================
# VISUALIZATIONS
# =============================================================================
cat("\n--- Generating Visualizations ---\n")

# Harvest color palette matching spill_variability_analysis.R
harvest_colors <- c("Open" = "#E07B54", "Closed" = "#5B8DB8")

# ---- Figure 1: Data Overview ----
p1a <- ggplot(year_summary_m1,
              aes(x = ReturnYear, y = pct_detected, fill = Harvest_Status)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = harvest_colors, name = "Harvest") +
  geom_text(aes(label = sprintf("n=%d", n_total)), vjust = -0.4, size = 3) +
  labs(x = "Return Year", y = "% Detected in Entiat River",
       title = "Annual Entiat Return Rate (all fish)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Detection rate by group and season
group_season_summary <- fish_AB %>%
  group_by(group, ReturnSeason) %>%
  summarise(n = n(), pct = mean(EntiatDetected)*100, .groups="drop")

p1b <- ggplot(group_season_summary,
              aes(x = group, y = pct, fill = ReturnSeason)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_text(aes(label = sprintf("n=%d", n)), position=position_dodge(0.9),
            vjust=-0.4, size=3) +
  scale_fill_manual(values = c("Fall"="#E07B54", "Spring"="#5B8DB8"),
                    name = "Return Season") +
  labs(x = "Group", y = "% Detected in Entiat River",
       title = "Entiat Detection Rate by Group and Season") +
  theme_minimal()

p1c <- fish_AB %>%
  group_by(group) %>%
  summarise(n=n(), pct=mean(EntiatDetected)*100, .groups="drop") %>%
  ggplot(aes(x=group, y=pct, fill=group)) +
  geom_col(alpha=0.85) +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)", pct, n)), vjust=-0.3, size=3.5) +
  scale_fill_manual(values=c("A"="#4CAF50","B"="#FF5722"), guide="none") +
  scale_x_discrete(labels=c("A"="Group A\n(Rocky Reach only)",
                             "B"="Group B\n(Wells interaction)")) +
  labs(x="", y="% Detected in Entiat River",
       title="Entiat Detection Rate by Group") +
  ylim(0, 100) +
  theme_minimal()

p1d <- ggplot(year_summary_m1,
              aes(x=ReturnYear, fill=Harvest_Status)) +
  geom_col(aes(y=n_groupA), fill="#4CAF50", alpha=0.6) +
  geom_col(aes(y=-n_groupB), fill="#FF5722", alpha=0.6) +
  geom_hline(yintercept=0, linewidth=0.5) +
  scale_fill_manual(values=harvest_colors, name="Harvest") +
  labs(x="Return Year", y="Fish Count (A above, B below)",
       title="Annual Sample: Group A (green) vs B (orange)") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

fig1 <- (p1a | p1b) / (p1c | p1d)
ggsave("entiat_data_overview.png", fig1, width=14, height=10, dpi=150)
cat("Saved: entiat_data_overview.png\n")

# ---- Figure 2: Model 1 Posterior Distributions ----

# Wells interaction posterior
post1_df <- tibble(
  wells_interaction = beta_wells_samples,
  harvest_open      = beta_harvest_samples
) %>%
  pivot_longer(everything(), names_to="Parameter", values_to="Value") %>%
  mutate(Parameter = recode(Parameter,
    "wells_interaction" = "Wells Interaction (Group B vs A)",
    "harvest_open"      = "Harvest Open (vs Closed)"
  ))

p2a <- ggplot(post1_df, aes(x=Value, fill=Parameter)) +
  geom_density(alpha=0.55) +
  geom_vline(xintercept=0, color="black", linetype="dashed", linewidth=0.8) +
  facet_wrap(~Parameter, scales="free", ncol=1) +
  scale_fill_manual(values=c("#FF5722","#5B8DB8"), guide="none") +
  labs(x="Coefficient (log-odds scale)", y="Posterior Density",
       title="Model 1 Posteriors: All Fish") +
  theme_minimal()

# Odds ratio plot for Model 1
or_df <- tibble(
  Parameter = c("Wells Interaction\n(Group B vs A)", "Harvest Open\n(vs Closed)"),
  OR_Mean   = c(exp(mean(beta_wells_samples)), exp(mean(beta_harvest_samples))),
  OR_Low    = c(exp(wells_hdi[1]),   exp(harvest_hdi[1])),
  OR_High   = c(exp(wells_hdi[2]),   exp(harvest_hdi[2])),
  P_neg     = c(mean(beta_wells_samples < 0), mean(beta_harvest_samples < 0))
)

p2b <- ggplot(or_df, aes(x=OR_Mean, y=Parameter)) +
  geom_pointrange(aes(xmin=OR_Low, xmax=OR_High), size=0.8, linewidth=1.2,
                  color=c("#FF5722","#5B8DB8")) +
  geom_vline(xintercept=1, linetype="dashed", color="gray40") +
  geom_text(aes(label=sprintf("OR=%.2f\nP(decrease)=%.3f", OR_Mean, P_neg)),
            hjust=-0.1, size=3.5) +
  scale_x_log10() +
  labs(x="Odds Ratio (log scale)", y="",
       title="Model 1: Effect Estimates (95% HDI)") +
  theme_minimal()

fig2 <- p2a | p2b
ggsave("entiat_model1_posteriors.png", fig2, width=13, height=6, dpi=150)
cat("Saved: entiat_model1_posteriors.png\n")

# ---- Figure 3: Model 1 vs 1b — Wells effect with/without season covariate ----

wells_compare_df <- bind_rows(
  tibble(samples=beta_wells_samples,   Model="Model 1\n(no season covariate)"),
  tibble(samples=beta_wells_1b,        Model="Model 1b\n(+ spring return season)")
)

p3a <- ggplot(wells_compare_df, aes(x=samples, fill=Model)) +
  geom_density(alpha=0.55) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.8) +
  scale_fill_manual(values=c("#FF5722","#9C27B0")) +
  labs(x="Wells Interaction Coefficient (log-odds)", y="Posterior Density",
       title="Wells Effect: Robustness to Season Covariate") +
  theme_minimal()

season_df <- tibble(samples=beta_season_1b)

p3b <- ggplot(season_df, aes(x=samples)) +
  geom_density(fill="#5B8DB8", alpha=0.55) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.8) +
  labs(x="Spring Return Coefficient (log-odds)", y="Posterior Density",
       title=sprintf("Season Effect (P(spring decreases return)=%.3f)",
                     mean(beta_season_1b < 0))) +
  theme_minimal()

fig3 <- p3a | p3b
ggsave("entiat_model1b_season.png", fig3, width=12, height=5, dpi=150)
cat("Saved: entiat_model1b_season.png\n")

# ---- Figure 4: Year Random Effects ----
year_re_m1 <- ranef(fit_model1)$ReturnYearFactor[,,"Intercept"]
year_re_df <- tibble(
  Year     = as.integer(rownames(year_re_m1)),
  Mean     = year_re_m1[,"Estimate"],
  SD       = year_re_m1[,"Est.Error"]
) %>%
  left_join(harvest_df, by="Year")

p4 <- ggplot(year_re_df, aes(x=Year, y=Mean, color=Harvest_Status)) +
  geom_pointrange(aes(ymin=Mean-1.96*SD, ymax=Mean+1.96*SD), size=0.7, linewidth=1) +
  geom_hline(yintercept=0, linetype="dashed", color="gray50") +
  scale_color_manual(values=harvest_colors, name="Harvest") +
  labs(x="Return Year", y="Year Effect (logit scale)",
       title="Model 1: Year-Level Random Effects (±95% CI)") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("entiat_year_effects.png", p4, width=11, height=5, dpi=150)
cat("Saved: entiat_year_effects.png\n")

# ---- Figure 5: Model 2 Posteriors (if fitted) ----
if (!is.null(fit_model2)) {

  model2_post_df <- tibble(
    SpringSpillHours = b_spring2,
    FallSpillHours   = b_fall2,
    HarvestOpen      = b_harv2
  ) %>%
    pivot_longer(everything(), names_to="Parameter", values_to="Value") %>%
    mutate(Parameter = recode(Parameter,
      "SpringSpillHours" = "Spring Spill Hours (std)",
      "FallSpillHours"   = "Fall Spill Hours (std)",
      "HarvestOpen"      = "Harvest Open (vs Closed)"
    ))

  p5a <- ggplot(model2_post_df, aes(x=Value, fill=Parameter)) +
    geom_density(alpha=0.55) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.8) +
    facet_wrap(~Parameter, ncol=1, scales="free_y") +
    scale_fill_manual(values=c("#FF9800","#4CAF50","#5B8DB8"), guide="none") +
    labs(x="Coefficient (log-odds)", y="Posterior Density",
         title=sprintf("Model 2: Group B Spill Analysis (n=%d)", nrow(group_b_spill))) +
    theme_minimal()

  # Scatterplot: spring spill vs detection rate by year
  p5b <- group_b_spill %>%
    group_by(ReturnYear, Harvest_Status) %>%
    summarise(n=n(), pct=mean(EntiatDetected)*100,
              SpringSpillHours=mean(SpringSpillHours), .groups="drop") %>%
    ggplot(aes(x=SpringSpillHours, y=pct, color=Harvest_Status, size=n)) +
    geom_point(alpha=0.8) +
    geom_smooth(method="lm", se=FALSE, color="gray30", linewidth=0.8) +
    geom_text(aes(label=ReturnYear), vjust=-0.8, size=3, color="black") +
    scale_color_manual(values=harvest_colors, name="Harvest") +
    labs(x="Spring Spill Hours (Jan-Mar)", y="% Group B Detected in Entiat",
         title="Spring Spill at Wells vs Group B Entiat Return (by year)") +
    theme_minimal() +
    guides(size="none")

  fig5 <- p5a | p5b
  ggsave("entiat_model2_spill.png", fig5, width=13, height=8, dpi=150)
  cat("Saved: entiat_model2_spill.png\n")
}

# ---- Figure 6: MCMC Trace Plots ----
cat("Generating trace plots for Model 1...\n")
trace1 <- mcmc_plot(fit_model1,
                    pars = c("b_wells_interaction", "b_harvest_open",
                             "sd_ReturnYearFactor__Intercept"),
                    type = "trace")
ggsave("entiat_trace_plots_model1.png", trace1, width=10, height=7, dpi=150)
cat("Saved: entiat_trace_plots_model1.png\n")

# =============================================================================
# SAVE MODEL OBJECTS
# =============================================================================
cat("\n--- Saving Model Objects ---\n")

if (!is.null(fit_model2)) {
  save(fit_model1, fit_model1b, fit_model2,
       fish_AB, group_b_spill, year_summary_m1, results_summary,
       file = "entiat_bayesian_models.RData")
} else {
  save(fit_model1, fit_model1b,
       fish_AB, year_summary_m1, results_summary,
       file = "entiat_bayesian_models.RData")
}
cat("Saved: entiat_bayesian_models.RData\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n======================================================================\n")
cat("FINAL SUMMARY\n")
cat("======================================================================\n")

cat(sprintf("
Dataset: %d wild Entiat steelhead across %d return years
  Group A (Rocky Reach only, no Wells): %d fish
  Group B (Wells interaction, non-stray): %d fish
  Excluded strays (last in Methow/Okanogan): %d fish

Entiat River detection rates:
  Group A: %.1f%% (%d/%d)
  Group B: %.1f%% (%d/%d)

MODEL 1 — Wells Interaction Effect (all fish, year random effect):
  beta_wells = %.3f (95%% HDI: [%.3f, %.3f])
  Odds Ratio = %.3f
  P(Wells DECREASES Entiat return) = %.4f

  beta_harvest = %.3f (95%% HDI: [%.3f, %.3f])
  Odds Ratio = %.3f
  P(harvest open increases return) = %.4f

MODEL 1b — Wells Effect After Controlling for Spring Return Season:
  beta_wells = %.3f (95%% HDI: [%.3f, %.3f])
  Odds Ratio = %.3f
  P(Wells DECREASES Entiat return) = %.4f
  [Robustness check: effect should persist after accounting for
   lower spring detection efficiency]
",
  nrow(fish_AB) + sum(fish_flags$group == "Stray"),
  n_distinct(fish_AB$ReturnYear),
  sum(fish_AB$group=="A"), sum(fish_AB$group=="B"),
  sum(fish_flags$group == "Stray"),
  mean(fish_AB$EntiatDetected[fish_AB$group=="A"])*100,
  sum(fish_AB$EntiatDetected[fish_AB$group=="A"]),
  sum(fish_AB$group=="A"),
  mean(fish_AB$EntiatDetected[fish_AB$group=="B"])*100,
  sum(fish_AB$EntiatDetected[fish_AB$group=="B"]),
  sum(fish_AB$group=="B"),
  # Model 1
  mean(beta_wells_samples), wells_hdi[1], wells_hdi[2],
  exp(mean(beta_wells_samples)),
  mean(beta_wells_samples < 0),
  mean(beta_harvest_samples), harvest_hdi[1], harvest_hdi[2],
  exp(mean(beta_harvest_samples)),
  mean(beta_harvest_samples > 0),
  # Model 1b
  mean(beta_wells_1b), wells_hdi_1b[1], wells_hdi_1b[2],
  exp(mean(beta_wells_1b)),
  mean(beta_wells_1b < 0)
))

if (!is.null(fit_model2)) {
  cat(sprintf(
"MODEL 2 — Spill at Wells Predicting Group B Entiat Return:
  Spring Spill Hours: beta=%.3f, OR=%.3f, P(>0)=%.3f
  Fall Spill Hours:   beta=%.3f, OR=%.3f, P(>0)=%.3f
  Harvest Open:       beta=%.3f, OR=%.3f, P(>0)=%.3f
  [None of the spill/harvest metrics predict within-group variation]
",
    mean(b_spring2), exp(mean(b_spring2)), mean(b_spring2 > 0),
    mean(b_fall2),   exp(mean(b_fall2)),   mean(b_fall2 > 0),
    mean(b_harv2),   exp(mean(b_harv2)),   mean(b_harv2 > 0)
  ))
}

cat("\nOUTPUT FILES WRITTEN:\n")
cat("  entiat_analysis_all_fish.csv          — per-fish analysis dataset\n")
cat("  entiat_analysis_groupB.csv            — Group B fish with spill data\n")
cat("  entiat_year_summary.csv               — annual summary table\n")
cat("  entiat_return_bayesian_results.csv    — model results table\n")
cat("  entiat_data_overview.png              — data overview figure\n")
cat("  entiat_model1_posteriors.png          — Model 1 posterior distributions\n")
cat("  entiat_model1b_season.png             — season covariate robustness check\n")
cat("  entiat_year_effects.png               — year random effects\n")
if (!is.null(fit_model2))
  cat("  entiat_model2_spill.png               — Model 2 spill analysis\n")
cat("  entiat_trace_plots_model1.png         — MCMC trace plots\n")
cat("  entiat_bayesian_models.RData          — saved brms model objects\n")

cat("\n======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("======================================================================\n")
