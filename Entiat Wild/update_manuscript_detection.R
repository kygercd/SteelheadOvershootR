###############################################################################
# Update Entiat_Wild_Steelhead_Wells_Overshoot_Manuscript.docx
# with detection-corrected model results (Models 1c, 1d) and
# hidden stray analysis
###############################################################################

library(officer)

setwd("/home/chas/SteelheadOvershootR/Entiat Wild")

doc <- read_docx("Entiat_Wild_Steelhead_Wells_Overshoot_Manuscript.docx")

# =============================================================================
# HELPER: add a body paragraph using fpar (no named style required)
# =============================================================================
add_body <- function(doc, text, pos = "after") {
  body_add_fpar(doc, fpar(ftext(text)), pos = pos)
}

add_heading3 <- function(doc, text, pos = "after") {
  body_add_par(doc, value = text, style = "Heading 3", pos = pos)
}

# =============================================================================
# CHANGE 1: ABSTRACT
# Add detection-corrected result sentence after the spring season sentence
# =============================================================================

doc <- body_replace_all_text(
  doc,
  old_value = paste0("ruling out lower spring detection efficiency at Entiat antennas as a ",
                     "confounding explanation. Among Group B fish"),
  new_value = paste0(
    "ruling out lower spring detection efficiency at Entiat antennas as a confounding explanation. ",
    "A compound-detection-likelihood model incorporating per-fish ENL detection efficiency ",
    "(Y \u223c Bernoulli(\u03b8 \u00d7 p_ENL); mean p_ENL = 0.86) yielded a stronger estimate ",
    "(OR = 0.45; 95% CI: 0.26\u20130.77; P(\u03b2 < 0) = 0.998), and a hidden stray analysis ",
    "estimated fewer than three undetected strays among 53 undetected Group B fish, ",
    "confirming the effect is not attributable to imperfect detection or hidden straying. ",
    "Among Group B fish"
  ),
  only_at_cursor = FALSE
)

# =============================================================================
# CHANGE 2: INTRODUCTION — Detection Efficiency Considerations
# Extend paragraph with 4th approach (compound likelihood)
# =============================================================================

doc <- body_replace_all_text(
  doc,
  old_value = "fitting a robustness model that includes spring return timing as an explicit covariate.",
  new_value = paste0(
    "fitting a robustness model that includes spring return timing as an explicit covariate; and (4) ",
    "fitting a compound Bernoulli likelihood model, Y_i \u223c Bernoulli(\u03b8_i \u00d7 p_ENL_i), in which ",
    "per-fish detection efficiency p_ENL_i was estimated from a seasonal-by-annual GAM of ENL ",
    "detection data, allowing undetected fish to contribute to the likelihood either as ",
    "non-returners or as missed returners in proportion to their known detection probability."
  ),
  only_at_cursor = FALSE
)

# =============================================================================
# CHANGE 3: METHODS — Add Models 1c and 1d after the Model 2 list item
# =============================================================================

doc <- cursor_reach(doc, keyword = "Model 2 \\(Group B only\\)")

doc <- add_body(
  doc,
  paste0(
    "Model 1c (detection-corrected): Y_i | vreal(p_ENL_i) \u223c wells_interaction + ",
    "harvest_open + (1 | ReturnYearFactor), fitted with compound Bernoulli likelihood ",
    "P(Y_i = 1) = \u03b8_i \u00d7 p_ENL_i implemented via a brms custom family. ",
    "Per-fish p_ENL_i estimated from ENL gateway GAM at each fish's return date and year; ",
    "11 pre-operational returns used p_ENL_i = 1.0."
  )
)

doc <- add_body(
  doc,
  "Model 1d (detection-corrected + season): adds spring_return to Model 1c.",
  pos = "after"
)

# =============================================================================
# CHANGE 4: RESULTS — New subsection before "Year Random Effects"
# =============================================================================

doc <- cursor_reach(doc, keyword = "Year Random Effects")

# Add a blank spacer first
doc <- add_body(doc, "", pos = "before")

# Add fate paragraph (reverse order since we're inserting "before")
doc <- add_body(
  doc,
  paste0(
    "Detection-corrected return rates (Group A: \u223c80%; Group B: \u223c72%) are both higher than raw ",
    "observed rates (73.2% and 58.6% respectively), reflecting ENL detection misses. We additionally ",
    "estimated expected fates for the 53 undetected Group B fish using per-fish p_LMR and p_OKL ",
    "from Methow and Okanogan gateway detection GAMs: \u223c19.2 fish (36%) were expected to be true ",
    "Entiat returners missed at ENL; \u223c2.1 fish (4%) hidden strays missed at LMR or OKL; and ",
    "\u223c31.7 fish (60%) likely died, were harvested, or had other unobservable fates. The negligible ",
    "hidden stray count confirms that misclassification of undetected strays is not driving the ",
    "Wells interaction result."
  ),
  pos = "before"
)

doc <- add_body(
  doc,
  paste0(
    "The detection-corrected model (Model 1c) yielded \u03b2_wells = \u22120.790 ",
    "(95% CI: \u22121.342, \u22120.260; OR = 0.45; P(\u03b2 < 0) = 0.998), and Model 1d (adding spring ",
    "season) yielded \u03b2_wells = \u22120.789 (95% CI: \u22121.344, \u22120.252; OR = 0.45; P(\u03b2 < 0) = 0.998). ",
    "Both are \u223c0.15 log-odds units more negative than Model 1, because the compound likelihood ",
    "infers that Group A's higher baseline return rate makes its undetected fish relatively more ",
    "likely to be missed returners, slightly widening the inferred gap between groups. ",
    "Both groups had virtually identical mean ENL detection efficiencies (Group A: 0.864; Group B: 0.864), ",
    "so differential per-group detection probability is not a plausible source of bias. ",
    "Convergence diagnostics were excellent (R\u0302 \u2264 1.001; bulk ESS > 12,000 for all fixed effects)."
  ),
  pos = "before"
)

doc <- add_body(
  doc,
  paste0(
    "To formally account for imperfect ENL detection, we fitted a compound Bernoulli likelihood ",
    "model: Y_i \u223c Bernoulli(\u03b8_i \u00d7 p_ENL_i), where \u03b8_i = inv_logit(\u03b7_i) is the true return ",
    "probability and p_ENL_i is per-fish ENL detection efficiency from the ENL gateway GAM ",
    "(return date and year). Detection efficiencies ranged from 0.25 to 1.0 (mean = 0.860; ",
    "11 pre-operational returns set to p_ENL = 1.0)."
  ),
  pos = "before"
)

doc <- add_heading3(doc, "Detection-Corrected Models (Models 1c and 1d)", pos = "before")

# =============================================================================
# CHANGE 5: TABLE 3 caption — update "three models" to "all models" and add note
# =============================================================================

doc <- body_replace_all_text(
  doc,
  old_value = "Table 3. Posterior summaries for all three models.",
  new_value  = "Table 3. Posterior summaries for all models.",
  only_at_cursor = FALSE
)

doc <- body_replace_all_text(
  doc,
  old_value = "Bold P(decrease) values indicate strong directional evidence (\u22650.95).",
  new_value  = paste0(
    "Bold P(decrease) values indicate strong directional evidence (\u22650.95). ",
    "Detection-corrected models (1c, 1d; see Results): \u03b2_wells = \u22120.790/\u22120.789, ",
    "OR = 0.45, P(decrease) = 0.998."
  ),
  only_at_cursor = FALSE
)

# =============================================================================
# CHANGE 6: DISCUSSION — Add 4th line of evidence; update "Taken together"
# =============================================================================

doc <- body_replace_all_text(
  doc,
  old_value = paste0(
    "Taken together, these results firmly establish that the group gap reflects a genuine ",
    "behavioral or survival difference rather than a detection artifact."
  ),
  new_value = paste0(
    "Fourth, a compound Bernoulli likelihood model (Model 1c) incorporating per-fish ENL detection ",
    "efficiency yielded a stronger Wells effect estimate (\u03b2 = \u22120.790, OR = 0.45; P(\u03b2 < 0) = 0.998) ",
    "than the standard model, demonstrating that accounting for detection imperfection amplifies ",
    "rather than attenuates the inferred effect. A hidden stray analysis estimated only \u223c2.1 ",
    "hidden strays (4% of 53 undetected Group B fish) missed at LMR or OKL, confirming that ",
    "contamination with undetected non-natal spawners is negligible. ",
    "Taken together, these four lines of evidence firmly establish that the group gap reflects a ",
    "genuine behavioral or survival difference rather than a detection artifact."
  ),
  only_at_cursor = FALSE
)

# Update the detection gaps bullet point
doc <- body_replace_all_text(
  doc,
  old_value = "Correcting for these rates narrows but does not eliminate the group gap.",
  new_value = paste0(
    "Formal correction via compound Bernoulli likelihood (Model 1c) yielded OR = 0.45, ",
    "indicating the true effect is larger than the naive estimate suggests."
  ),
  only_at_cursor = FALSE
)

# =============================================================================
# CHANGE 7: CONCLUSIONS — Add reference to detection-corrected result
# =============================================================================

doc <- body_replace_all_text(
  doc,
  old_value = paste0(
    "This effect is robust to seasonal detection efficiency confounds and consistent across ",
    "dam operating conditions."
  ),
  new_value = paste0(
    "A compound-detection-likelihood model corroborated this finding (OR = 0.45; P(\u03b2 < 0) = 0.998), ",
    "and a hidden stray analysis found fewer than three hidden strays among 53 undetected Group B fish. ",
    "This effect is therefore robust to both seasonal detection efficiency confounds and potential ",
    "hidden straying, and consistent across dam operating conditions."
  ),
  only_at_cursor = FALSE
)

# =============================================================================
# SAVE
# =============================================================================

print(doc, target = "Entiat_Wild_Steelhead_Wells_Overshoot_Manuscript.docx")

cat("\nManuscript updated and saved successfully.\n")
cat("Changes made:\n")
cat("  1. Abstract: added detection-corrected result sentence\n")
cat("  2. Introduction: extended detection efficiency paragraph with 4th approach\n")
cat("  3. Methods: added Models 1c and 1d descriptions after Model 2\n")
cat("  4. Results: added 'Detection-Corrected Models' subsection before Year RE section\n")
cat("  5. Table 3 caption: updated heading and added detection-corrected note\n")
cat("  6. Discussion: added 4th line of evidence; updated detection gap bullet\n")
cat("  7. Conclusions: added detection-corrected robustness statement\n")
