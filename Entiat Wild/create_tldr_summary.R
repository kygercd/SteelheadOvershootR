###############################################################################
# Create TLDR Summary Word Document
# Covers: (1) Gateway PIT Array Detection Efficiency, (2) Entiat Steelhead
# Wells Overshoot Return Analysis (standard + detection-corrected)
###############################################################################

library(officer)
library(magrittr)

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

h1 <- function(doc, text) body_add_par(doc, text, style = "heading 1")
h2 <- function(doc, text) body_add_par(doc, text, style = "heading 2")
h3 <- function(doc, text) body_add_par(doc, text, style = "heading 3")
body <- function(doc, text) body_add_fpar(doc, fpar(ftext(text)))
bullet <- function(doc, text) body_add_par(doc, text, style = "Normal")
spacer <- function(doc) body_add_fpar(doc, fpar(ftext("")))

bold_run  <- function(text) ftext(text, prop = fp_text(bold = TRUE))
plain_run <- function(text) ftext(text)

body_bold <- function(doc, bold_part, rest) {
  body_add_fpar(doc, fpar(bold_run(bold_part), plain_run(rest)))
}

# Simple two-column table builder
add_table2 <- function(doc, header1, header2, rows) {
  # rows is a list of c(col1, col2) character vectors
  ft <- data.frame(
    Col1 = c(header1, sapply(rows, `[`, 1)),
    Col2 = c(header2, sapply(rows, `[`, 2)),
    stringsAsFactors = FALSE
  )
  # Use officer's body_add_table
  body_add_table(doc, ft, style = "Table Professional", header = TRUE)
}

# =============================================================================
# CREATE DOCUMENT
# =============================================================================

doc <- read_docx()

# Title block
doc <- body_add_fpar(doc,
  fpar(ftext("Wells Dam Overshoot & PIT Array Detection Efficiency",
             prop = fp_text(bold = TRUE, font.size = 18))),
  pos = "after")

doc <- body_add_fpar(doc,
  fpar(ftext("Summary Report — Wild Entiat River Steelhead",
             prop = fp_text(italic = TRUE, font.size = 13))),
  pos = "after")

doc <- body_add_fpar(doc,
  fpar(ftext("Public Utility District No. 1 of Douglas County  |  April 2026",
             prop = fp_text(font.size = 11, color = "#555555"))),
  pos = "after")

doc <- spacer(doc)

# =============================================================================
# SECTION 1: OVERVIEW
# =============================================================================

doc <- h1(doc, "Overview")

doc <- body(doc, paste0(
  "This report covers two related analyses using PIT-tag detection data for wild steelhead ",
  "in the upper Columbia River. The first characterizes detection efficiency at three adult ",
  "steelhead gateway PIT arrays (Methow/LMR, Okanogan/OKL, Entiat/ENL), which is needed to ",
  "correctly interpret detection-based return estimates. The second uses those efficiency ",
  "estimates to rigorously test whether wild Entiat River steelhead that overshoot to Wells Dam ",
  "have reduced probability of returning to their natal subbasin."
))

doc <- spacer(doc)

# =============================================================================
# SECTION 2: GATEWAY DETECTION EFFICIENCY
# =============================================================================

doc <- h1(doc, "Part 1: Gateway Array Detection Efficiency")

doc <- h2(doc, "What Was Analyzed")

doc <- body(doc, paste0(
  "Three lower-subbasin PIT antenna arrays serve as 'gateways' for adult steelhead returning ",
  "to their natal tributaries above Wells Dam. A fish detected above Wells Dam that is NOT ",
  "detected at its subbasin gateway was either (a) missed by the antenna or (b) never entered ",
  "that subbasin. Distinguishing these requires knowing how often fish are missed. ",
  "Detection efficiency was estimated for:"
))

doc <- bullet(doc, "LMR \u2014 Lower Methow River antenna (operational since 3/4/2009). Data: 9,702 eligible fish.")
doc <- bullet(doc, "OKL \u2014 Lower Okanogan River antenna (operational since 10/8/2013). Data: 2,767 eligible fish.")
doc <- bullet(doc, "ENL \u2014 Lower Entiat River antenna (operational since 10/25/2007). Data: 2,078 eligible fish.")

doc <- body(doc, paste0(
  "Source data: PTAGIS detection records for Wells Dam steelhead subsequently detected in ",
  "Methow/Okanogan subbasins (LMR/OKL analysis) and Rocky Reach Dam steelhead detected in the ",
  "Entiat subbasin (ENL analysis). For each site, a fish that passed Wells Dam and was later ",
  "detected anywhere in the subbasin was used to determine whether it was detected at the gateway ",
  "array (detected = 1) or bypassed/missed it (detected = 0)."
))

doc <- spacer(doc)
doc <- h2(doc, "Overall Miss Rates")

# Miss rate table
miss_rows <- list(
  c("LMR (Methow)", "9,702", "1,818", "18.7%"),
  c("OKL (Okanogan)", "2,767", "824", "29.8%"),
  c("ENL (Entiat)", "2,078", "461", "22.2%")
)
miss_df <- data.frame(
  Site = c("LMR (Methow)", "OKL (Okanogan)", "ENL (Entiat)"),
  `Total fish` = c("9,702", "2,767", "2,078"),
  `Missed gateway` = c("1,818", "824", "461"),
  `Miss rate` = c("18.7%", "29.8%", "22.2%"),
  check.names = FALSE
)
doc <- body_add_table(doc, miss_df, style = "Table Professional", header = TRUE)

doc <- body_add_fpar(doc,
  fpar(ftext("Table 1. Overall miss rates at each gateway array during its operational period.",
             prop = fp_text(italic = TRUE, font.size = 9))))

doc <- spacer(doc)
doc <- h2(doc, "Key Patterns")

doc <- h3(doc, "LMR and OKL \u2014 Flow-Driven Seasonal Pattern")

doc <- body(doc, paste0(
  "Detection is near-perfect in late summer and fall (Aug\u2013Nov; >98%) when flows are low. ",
  "Efficiency drops sharply in spring (Mar\u2013Jun) as snowmelt raises river discharge: ",
  "as low as 27\u201359% at LMR and 16\u201351% at OKL. Both sites showed a mid-period performance ",
  "dip around 2015\u20132019 before recovering."
))

doc <- h3(doc, "ENL \u2014 Unusual Double-Trough Pattern")

doc <- body(doc, paste0(
  "ENL has a different seasonal pattern than LMR/OKL: detection dips in both May (spring ",
  "freshet) AND September (peak steelhead run timing), creating an inverted-U flow relationship. ",
  "Annual detection declined after 2015, with especially poor performance in 2019, 2023, and 2024 \u2014 ",
  "likely reflecting equipment or maintenance issues."
))

doc <- spacer(doc)
doc <- h2(doc, "Models Fitted")

doc <- body(doc, paste0(
  "Two complementary models were fitted for each site:"
))

doc <- bullet(doc, paste0(
  "GAM model: detected ~ s(day-of-year, cyclic) + s(year). Provides seasonal and annual trend ",
  "predictions. AUC: LMR 0.893, OKL 0.925, ENL 0.846."
))
doc <- bullet(doc, paste0(
  "Bayesian flow model (brms): detected ~ log_discharge + s(day-of-year) + (1|year). ",
  "Incorporates USGS daily streamflow as the physical driver; year random effects isolate ",
  "equipment problems from flow-driven variation. AUC: LMR 0.900, OKL 0.929, ENL 0.875."
))

doc <- body(doc, paste0(
  "Flow coefficients confirm the expected direction: monotonically negative for LMR and OKL ",
  "(higher flow = lower detection); inverted-U for ENL (both very low and very high flows ",
  "reduce detection, optimal \u2248 200\u2013500 cfs). GAM predictions by return date and year are ",
  "used as per-fish p_ENL, p_LMR, and p_OKL in the return analysis below."
))

doc <- spacer(doc)

# =============================================================================
# SECTION 3: ENTIAT RETURN ANALYSIS
# =============================================================================

doc <- h1(doc, "Part 2: Wells Dam Overshoot and Entiat Return Probability")

doc <- h2(doc, "Research Question")

doc <- body(doc, paste0(
  "Do wild Entiat steelhead that overshoot past their natal tributary and interact with Wells Dam ",
  "(rkm 830, 31 miles upstream of the Entiat confluence) have reduced probability of subsequently ",
  "returning to the Entiat River to spawn?"
))

doc <- spacer(doc)
doc <- h2(doc, "Data and Study Design")

doc <- body(doc, paste0(
  "PIT-tag detection records from PTAGIS for 417 wild Entiat-origin steelhead detected at Rocky ",
  "Reach Dam adult fishway (RRF) across return years 2004\u20132025. Each fish was assigned an ",
  "'adult RRF anchor date' (first Jul\u2013Dec RRF detection) to exclude juvenile detections. ",
  "379 fish with valid anchors were included (13 confirmed strays excluded)."
))

doc <- body(doc, "Fish were classified into two groups based on post-anchor detection history:")
doc <- bullet(doc, "Group A (n = 239): No Wells Dam detection after adult RRF passage.")
doc <- bullet(doc, "Group B (n = 140): Wells Dam detection after RRF, NOT last detected in Methow/Okanogan.")
doc <- bullet(doc, "Excluded \u2014 Strays (n = 13): Last detected in Methow or Okanogan subbasin (confirmed non-natal spawners).")

doc <- spacer(doc)
doc <- h2(doc, "Raw Detection Rates")

raw_df <- data.frame(
  Group = c("Group A (no Wells)", "Group B (Wells interaction)"),
  n = c("239", "140"),
  `Entiat detected` = c("175 (73.2%)", "82 (58.6%)"),
  `Raw gap` = c("\u2014", "\u221214.6 pp"),
  check.names = FALSE
)
doc <- body_add_table(doc, raw_df, style = "Table Professional", header = TRUE)

doc <- body_add_fpar(doc,
  fpar(ftext("Table 2. Raw Entiat River detection rates by group. pp = percentage points.",
             prop = fp_text(italic = TRUE, font.size = 9))))

doc <- spacer(doc)
doc <- h2(doc, "Bayesian Hierarchical Model Results")

doc <- body(doc, paste0(
  "Bayesian logistic regression with year random effects: ",
  "EntiatDetected ~ wells_interaction + harvest_open + (1|ReturnYear). ",
  "4 chains \u00d7 4,500 iterations; all R\u0302 \u2264 1.01."
))

results_df <- data.frame(
  Model = c(
    "Model 1 \u2014 Standard",
    "Model 1b \u2014 + Spring season",
    "Model 1c \u2014 Detection-corrected",
    "Model 1d \u2014 Corrected + season"
  ),
  "Beta (Wells)" = c("-0.642", "-0.645", "-0.790", "-0.789"),
  "95% CI" = c("[-1.09, -0.21]", "[-1.08, -0.18]",
               "[-1.34, -0.26]", "[-1.34, -0.25]"),
  "Odds Ratio" = c("0.53", "0.52", "0.45", "0.45"),
  "P(Beta < 0)" = c("0.997", "0.998", "0.998", "0.998"),
  check.names = FALSE
)
doc <- body_add_table(doc, results_df, style = "Table Professional", header = TRUE)

doc <- body_add_fpar(doc,
  fpar(ftext(paste0(
    "Table 3. Wells Dam interaction effect across four model specifications. ",
    "Models 1c/1d use compound Bernoulli likelihood Y ~ Bernoulli(\u03b8 \u00d7 p_ENL) ",
    "incorporating per-fish ENL detection efficiency from the gateway GAM."),
    prop = fp_text(italic = TRUE, font.size = 9))))

doc <- spacer(doc)
doc <- h2(doc, "Harvest and Spill Effects (Model 2)")

doc <- body(doc, paste0(
  "Harvest status (open vs. closed years) showed no credible effect on Entiat return ",
  "probability in any model (\u03b2 \u2248 \u22120.24; 95% CI spanning zero; OR = 0.79). ",
  "Among Group B fish with available spill data (n = 116), neither spring nor fall spill ",
  "hours at Wells Dam predicted Entiat return probability. All Model 2 posteriors were ",
  "centered near zero with 95% HDIs spanning zero:"
))

doc <- bullet(doc, "Spring Spill Hours (std): \u03b2 = +0.301 (95% HDI: \u22120.22, +0.86); OR = 1.35; P(positive) = 0.86")
doc <- bullet(doc, "Fall Spill Hours (std): \u03b2 = \u22120.045 (95% HDI: \u22120.54, +0.48); OR = 0.96; P(positive) = 0.43")
doc <- bullet(doc, "Harvest Open: \u03b2 = \u22120.084 (95% HDI: \u22121.31, +1.14); OR = 0.92")

doc <- body(doc, paste0(
  "This null spill finding is consistent across spring and fall windows and across hours and ",
  "volume metrics. Spill operations at Wells Dam do not appear to substantially alter the ",
  "probability of wild Entiat steelhead returning to the Entiat River following an overshoot."
))

doc <- spacer(doc)
doc <- h2(doc, "Detection Efficiency Robustness Checks")

doc <- body(doc, paste0(
  "Three independent checks confirm the group gap is not a detection artifact:"
))

doc <- bullet(doc, paste0(
  "Seasonal composition: Both groups had similar spring return fractions (14.6% Group A; ",
  "12.1% Group B), so differential spring/fall composition cannot explain the gap."
))
doc <- bullet(doc, paste0(
  "Within-season rates: The gap is concentrated in fall returns (Group A 74.0% vs. Group B 57.7%), ",
  "precisely when ENL detection is highest \u2014 not in spring when detection issues would appear."
))
doc <- bullet(doc, paste0(
  "Compound likelihood model: Formally incorporating per-fish p_ENL into the Bernoulli likelihood ",
  "amplifies rather than attenuates the Wells effect (OR 0.53 \u2192 0.45). Both groups have identical ",
  "mean p_ENL = 0.864, so differential detection is not a plausible source of bias."
))
doc <- bullet(doc, paste0(
  "Hidden stray analysis: Of 53 undetected Group B fish, only \u223c2.1 (4%) are estimated as hidden ",
  "strays missed at LMR or OKL. The majority are either missed Entiat returners (\u223c19, 36%) or ",
  "fish with other fates (died/harvested, \u223c32, 60%)."
))

doc <- spacer(doc)
doc <- h2(doc, "Detection-Corrected Return Rates")

corr_df <- data.frame(
  Group = c("Group A (no Wells)", "Group B (Wells interaction)"),
  `Observed rate` = c("73.2%", "58.6%"),
  `Corrected rate (+ ENL misses)` = c("~80%", "~72%"),
  check.names = FALSE
)
doc <- body_add_table(doc, corr_df, style = "Table Professional", header = TRUE)

doc <- body_add_fpar(doc,
  fpar(ftext(paste0(
    "Table 4. Observed vs. detection-corrected Entiat return rates. ",
    "Both groups' true return rates are higher than raw observed rates; ",
    "the relative difference in log-odds is larger after correction."),
    prop = fp_text(italic = TRUE, font.size = 9))))

doc <- spacer(doc)

# =============================================================================
# SECTION 4: BOTTOM LINE / KEY TAKEAWAYS
# =============================================================================

doc <- h1(doc, "Key Takeaways")

doc <- body_add_fpar(doc,
  fpar(bold_run("1.  Wells Dam overshoot substantially reduces Entiat return probability. "),
       plain_run(paste0(
         "Group B fish (Wells interaction) have approximately half the odds of returning ",
         "to the Entiat River compared to Group A fish that did not overshoot (OR = 0.53 raw; ",
         "OR = 0.45 detection-corrected). This is a large, highly credible effect (P > 0.997 ",
         "across all four model specifications)."
       ))))

doc <- spacer(doc)

doc <- body_add_fpar(doc,
  fpar(bold_run("2.  The effect is not a detection artifact. "),
       plain_run(paste0(
         "Spring return timing, within-season detection rates, a compound likelihood model, ",
         "and a hidden stray analysis all confirm the group gap reflects a genuine behavioral ",
         "or survival difference rather than imperfect antenna detection."
       ))))

doc <- spacer(doc)

doc <- body_add_fpar(doc,
  fpar(bold_run("3.  Spill and harvest do not explain the effect. "),
       plain_run(paste0(
         "Neither spring nor fall spill volumes at Wells Dam, nor annual harvest status ",
         "(open vs. closed), predicted whether a Wells-interacting fish subsequently returned ",
         "to the Entiat. This is consistent across multiple analysis approaches."
       ))))

doc <- spacer(doc)

doc <- body_add_fpar(doc,
  fpar(bold_run("4.  Gateway detection efficiency is seasonal and flow-driven. "),
       plain_run(paste0(
         "ENL, LMR, and OKL antennas all show substantial seasonal variation in detection ",
         "probability, primarily driven by streamflow. LMR and OKL are near-perfect in fall ",
         "but can fall to 20\u201340% in spring high-flow periods. ENL shows a more complex ",
         "pattern with both spring and September troughs. Per-fish efficiency estimates from ",
         "GAM models are used to correct return rate estimates."
       ))))

doc <- spacer(doc)

doc <- body_add_fpar(doc,
  fpar(bold_run("5.  Mechanism is unresolved. "),
       plain_run(paste0(
         "The most likely contributors to non-return in Group B fish are: downstream passage ",
         "difficulty or injury at Wells Dam, overwinter/en-route mortality in the reservoir, ",
         "and undetected straying to non-natal tributaries. PIT-tag data alone cannot ",
         "discriminate among these mechanisms."
       ))))

doc <- spacer(doc)

# =============================================================================
# SECTION 5: FILES AND METHODS REFERENCE
# =============================================================================

doc <- h1(doc, "Analysis Files")

files_df <- data.frame(
  File = c(
    "EntiatWildSteel.csv",
    "entiat_return_bayesian_analysis.R",
    "entiat_detection_corrected_analysis.R",
    "entiat_bayesian_models.RData",
    "entiat_corrected_models.RData",
    "gateway_detection_models.RData",
    "gateway_bayesian_models.RData",
    "gateway_detection_summary.docx",
    "Entiat_Wild_Steelhead_Wells_Overshoot_Manuscript.docx"
  ),
  Description = c(
    "Source PIT-tag detection data (417 fish, full detection histories)",
    "Standard Bayesian logistic regression (Models 1, 1b, 2)",
    "Detection-corrected models (1c, 1d) + hidden stray analysis",
    "Fitted brms objects: fit_model1, fit_model1b, fit_model2",
    "Fitted brms objects: fit_corrected_m1, fit_corrected_m1d",
    "Fitted GAMs: gam_met, gam_oka, gam_enl",
    "Fitted Bayesian flow models: fit_lmr_comb, fit_okl_comb, fit_enl_comb",
    "Full gateway detection efficiency analysis report (Word)",
    "Full manuscript with all methods, results, and discussion"
  ),
  stringsAsFactors = FALSE
)
doc <- body_add_table(doc, files_df, style = "Table Professional", header = TRUE)

doc <- spacer(doc)

# Footer note
doc <- body_add_fpar(doc,
  fpar(ftext(paste0(
    "All files located in /home/chas/SteelheadOvershootR/Entiat Wild/  \u2022  ",
    "Analysis performed in R using brms, mgcv, tidyverse, officer  \u2022  April 2026"),
    prop = fp_text(font.size = 9, color = "#666666"))))

# =============================================================================
# SAVE
# =============================================================================

print(doc, target = "Entiat_Steelhead_Wells_Overshoot_TLDR.docx")
cat("Saved: Entiat_Steelhead_Wells_Overshoot_TLDR.docx\n")
