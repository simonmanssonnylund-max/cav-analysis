###############################################################################
# CAV PREVALENCE FIGURE WITH TIME SINCE LAST CAV ASCERTAINMENT
#
# Existing objects required in the R session:
#   test          - final MSM-ready long-format dataset
#   prev_sahl_df  - observed and model-predicted prevalence by time and state
#
# Expected columns:
#   test:         lopnr, years, state, unit
#                 optional: artificial, obs_yrs
#   prev_sahl_df: time, state_name, observed, modelled
#
# Assumptions:
#   state 1, 2, 3 = observed/conclusive CAV states
#   state 4       = exact death/retransplantation
#   state 99      = alive; CAV state unknown (including covariate-update rows)
#   baseline at time 0 is assumed to establish CAV-free status
#
# The grey median/IQR describes time since the most recent REAL, CONCLUSIVE
# post-transplant CAV assessment (or baseline if no such assessment exists).
# It is computed in the living, unretransplanted, under-follow-up population.
# The prevalence curves can also include the absorbing death/reTx state;
# consequently the recency and prevalence denominators are not identical.
###############################################################################


# =============================================================================
# 1. PACKAGES AND SETTINGS
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)

TIME_MAX <- 7
TIME_STEP <- 0.05
SAHL_LEVEL <- "Göteborg"             # Change if the unit level differs
PLOT_FONT <- "Aptos"

OUT_PNG <- "Figure_CAV_prevalence_with_ascertainment_recency.png"
OUT_PDF <- "Figure_CAV_prevalence_with_ascertainment_recency.pdf"

required_test <- c("lopnr", "years", "state", "unit")
required_prev <- c("time", "state_name", "observed", "modelled")

stopifnot(
  exists("test"),
  exists("prev_sahl_df"),
  all(required_test %in% names(test)),
  all(required_prev %in% names(prev_sahl_df))
)

# Four facet labels from the original figure.
state_labels <- c(
  "max(FYND) = 1",
  "max(FYND) = 2",
  "max(FYND) >= 3",
  "Death/Re-HTx"
)

if (!setequal(as.character(unique(prev_sahl_df$state_name)), state_labels)) {
  stop(
    "Unexpected labels in prev_sahl_df$state_name. Check state_labels ",
    "against unique(prev_sahl_df$state_name)."
  )
}

# Keep the facet order of the original prevalence data.
prev_sahl_df <- prev_sahl_df %>%
  mutate(state_name = factor(as.character(state_name), levels = state_labels))


# =============================================================================
# 2. PREPARE FINAL MSM DATASET AND IDENTIFY ARTIFICIAL OBSERVATIONS
# =============================================================================

# Critically, use the original observed state, NOT statemax/statemax_msm:
# carrying a state forward does not constitute a new CAV assessment.
dat <- test %>%
  ungroup() %>%
  mutate(
    lopnr = as.character(lopnr),
    years = as.numeric(years),
    state = as.integer(as.character(state)),
    is_artificial = if ("artificial" %in% names(test)) {
      coalesce(as.logical(artificial), FALSE)
    } else {
      FALSE
    }
  )

# If the original dataset contains follow-up end times, use them.
# Otherwise, assume each patient has a terminal row at their true end of
# follow-up (NOT simply their most recent angiography).
if ("obs_yrs" %in% names(dat)) {
  dat <- dat %>% mutate(obs_yrs = as.numeric(obs_yrs))
} else {
  warning(
    "No obs_yrs in test: using each patient's maximum years value. ",
    "This is valid only if test includes the true end-of-follow-up row."
  )
  dat <- dat %>%
    group_by(lopnr) %>%
    mutate(obs_yrs = max(years, na.rm = TRUE)) %>%
    ungroup()
}

# Restrict to the same site used in the prevalence figure.
dat_sahl <- dat %>%
  filter(as.character(unit) == SAHL_LEVEL)

if (nrow(dat_sahl) == 0L) {
  stop("No patients found for SAHL_LEVEL = '", SAHL_LEVEL, "'.")
}


# =============================================================================
# 3. RECONSTRUCT TIMES OF CONCLUSIVE CAV ASSESSMENTS
# =============================================================================

# Include ALL genuine, conclusive assessments, including those that do not
# change a patient's previously recorded maximum CAV grade.
# Exclude state=99 (no new knowledge), state=4 (death), artificial rows,
# and the assumed baseline state (handled separately as time zero).
assessment_times <- dat_sahl %>%
  filter(
    years > 0,
    years <= TIME_MAX,
    state %in% 1:3,
    !is_artificial
  ) %>%
  transmute(
    lopnr,
    exam_year = years
  ) %>%
  distinct() %>%
  arrange(lopnr, exam_year)


# =============================================================================
# 4. RECONSTRUCT FOLLOW-UP AND EXACT DEATH / RETRANSPLANTATION TIMES
# =============================================================================

recency_patients <- dat_sahl %>%
  group_by(lopnr) %>%
  summarise(
    obs_yrs = max(obs_yrs, na.rm = TRUE),
    death_time = if (any(state == 4L, na.rm = TRUE)) {
      min(years[state == 4L], na.rm = TRUE)
    } else {
      Inf
    },
    .groups = "drop"
  ) %>%
  mutate(end_time = pmin(obs_yrs, TIME_MAX))

if (
  any(!is.finite(recency_patients$end_time)) ||
  any(recency_patients$end_time < 0) ||
  anyNA(recency_patients$lopnr)
) {
  stop("Missing or invalid patient-level follow-up times.")
}

# At any time t, 'at risk' means:
#   * alive and not retransplanted at t, and
#   * still under observation at t.
# This is for the recency calculation, not the prevalence denominator.


# =============================================================================
# 5. TIME SINCE LAST CAV ASCERTAINMENT AT EACH FOLLOW-UP TIME
# =============================================================================

recency_patient_times <- tidyr::crossing(
  recency_patients,
  time = seq(0, TIME_MAX, by = TIME_STEP)
) %>%
  filter(
    time <= end_time,
    time < death_time
  ) %>%
  left_join(
    assessment_times,
    by = join_by(lopnr, closest(time >= exam_year))
  ) %>%
  mutate(
    # A patient without a follow-up assessment retains baseline time zero.
    last_assessment = coalesce(exam_year, 0),
    years_since_assessment = time - last_assessment,
    baseline_only = is.na(exam_year)
  )

stopifnot(
  all(recency_patient_times$years_since_assessment >= -1e-10),
  all(recency_patient_times$years_since_assessment <= TIME_MAX + 1e-10)
)


# =============================================================================
# 6. SUMMARISE MEDIAN, IQR, AND PATIENTS WITHOUT FOLLOW-UP ANGIOGRAPHY
# =============================================================================

recency_summary <- recency_patient_times %>%
  group_by(time) %>%
  summarise(
    n_at_risk = n(),
    median_recency = median(years_since_assessment),
    q25 = unname(quantile(years_since_assessment, probs = 0.25)),
    q75 = unname(quantile(years_since_assessment, probs = 0.75)),
    pct_baseline_only = 100 * mean(baseline_only),
    .groups = "drop"
  )

# Check the distribution at clinically meaningful time points.
recency_key_times <- recency_summary %>%
  filter(round(time, 2) %in% c(1, 3, 6, 7)) %>%
  mutate(across(
    c(median_recency, q25, q75, pct_baseline_only),
    ~ round(.x, 2)
  ))

print(recency_key_times)


# =============================================================================
# 7. PREPARE REPEATED, FACET-SPECIFIC RECENCY SERIES
# =============================================================================

# In the original figure the prevalence y-axis ranges differ by state:
#   max(FYND)=1, Death/Re-HTx: 0-100%
#   max(FYND)=2, max(FYND)>=3: 0-10%
# A separate conversion is therefore required in each panel.
recency_plot <- tidyr::crossing(
  recency_summary,
  state_name = factor(state_labels, levels = state_labels)
) %>%
  mutate(
    ymax = case_when(
      as.character(state_name) %in%
        c("max(FYND) = 2", "max(FYND) >= 3") ~ 10,
      TRUE ~ 100
    ),
    median_scaled = median_recency * ymax / TIME_MAX,
    q25_scaled = q25 * ymax / TIME_MAX,
    q75_scaled = q75 * ymax / TIME_MAX
  )

# This function sets the left prevalence axis and a right recency axis.
# Both axes are linearly aligned within each individual facet.
make_yscale <- function(ymax) {
  scale_y_continuous(
    limits = c(0, ymax),
    breaks = seq(0, ymax, by = if (ymax == 10) 1 else 10),
    expand = expansion(mult = c(0, 0.02)),
    oob = scales::squish,
    sec.axis = sec_axis(
      transform = function(y) y * TIME_MAX / ymax,
      breaks = 0:7,
      name = "Years since last CAV ascertainment"
    )
  )
}

# Supplement the figure caption with year-7 data completeness information.
year7 <- recency_summary %>% filter(abs(time - TIME_MAX) < 1e-8)

year7_note <- if (nrow(year7) == 1L) {
  sprintf(
    "At year 7, %.1f%% of %d at-risk patients had no conclusive post-transplant CAV assessment.",
    year7$pct_baseline_only,
    year7$n_at_risk
  )
} else {
  ""
}


# =============================================================================
# 8. CREATE THE PUBLICATION-READY TWO-BY-TWO FIGURE
# =============================================================================

p_prev_sahl <- ggplot(
  prev_sahl_df,
  aes(x = time)
) +

  # 8.1 Grey IQR band, drawn behind the prevalence trajectories.
  geom_ribbon(
    data = recency_plot,
    aes(x = time, ymin = q25_scaled, ymax = q75_scaled),
    inherit.aes = FALSE,
    fill = "grey60",
    alpha = 0.18
  ) +

  # 8.2 Median time since latest CAV ascertainment.
  geom_line(
    data = recency_plot,
    aes(x = time, y = median_scaled),
    inherit.aes = FALSE,
    colour = "grey35",
    linewidth = 0.6
  ) +

  # 8.3 Observed prevalence (original data and aesthetics).
  geom_line(
    aes(
      y = observed,
      colour = "Observed",
      linetype = "Observed"
    ),
    linewidth = 0.9
  ) +

  # 8.4 Model-predicted prevalence (original data and aesthetics).
  geom_line(
    aes(
      y = modelled,
      colour = "Model-predicted",
      linetype = "Model-predicted"
    ),
    linewidth = 0.9
  ) +

  # 8.5 Four panels, with the original state-specific prevalence scales.
  facet_wrap(
    ~ state_name,
    ncol = 2,
    scales = "free_y"
  ) +

  ggh4x::facetted_pos_scales(
    y = list(
      state_name == "max(FYND) = 1" ~ make_yscale(100),
      state_name == "max(FYND) = 2" ~ make_yscale(10),
      state_name == "max(FYND) >= 3" ~ make_yscale(10),
      state_name == "Death/Re-HTx" ~ make_yscale(100)
    )
  ) +

  # 8.6 Common time axis.
  scale_x_continuous(
    limits = c(0, TIME_MAX),
    breaks = 0:7,
    expand = expansion(mult = c(0, 0.01))
  ) +

  # 8.7 Prevalence curves: colours and line types as in the original figure.
  scale_colour_manual(
    values = c(
      "Observed" = "#9e8b7f",
      "Model-predicted" = "#651220"
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Observed" = "solid",
      "Model-predicted" = "dashed"
    )
  ) +

  # 8.8 Labels and explanation of the grey overlay.
  labs(
    x = "Years after heart transplantation",
    y = "Probability in state (%)",
    colour = NULL,
    linetype = NULL,
    caption = paste(
      "Grey line: median time since last CAV state ascertainment;",
      "grey band: interquartile range.",
      "Baseline CAV-free status is assumed at transplantation.",
      year7_note
    )
  ) +

  # 8.9 Theme, gridlines, facet-title boxes, and margins.
  theme_minimal(
    base_size = 11,
    base_family = PLOT_FONT
  ) +

  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(
      linewidth = 0.25,
      colour = "grey85"
    ),
    panel.grid.major.y = element_line(
      linewidth = 0.25,
      colour = "grey85"
    ),

    panel.border = element_rect(
      fill = NA,
      colour = "grey40",
      linewidth = 0.5
    ),

    strip.background = element_rect(
      fill = "white",
      colour = "grey40",
      linewidth = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 11,
      margin = margin(t = 5, b = 5)
    ),

    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(1.5, "cm"),

    plot.caption = element_text(
      size = 8.5,
      hjust = 0
    ),

    plot.margin = margin(t = 10, r = 12, b = 8, l = 10)
  )

print(p_prev_sahl)


# =============================================================================
# 9. EXPORT FIGURE AND OPTIONAL DIAGNOSTIC TABLE
# =============================================================================

ggsave(
  filename = OUT_PNG,
  plot = p_prev_sahl,
  width = 12,
  height = 8,
  dpi = 600,
  bg = "white"
)

# Prefer Cairo to embed/render Aptos reliably when Cairo is available.
if (isTRUE(capabilities("cairo"))) {
  ggsave(
    filename = OUT_PDF,
    plot = p_prev_sahl,
    width = 12,
    height = 8,
    device = grDevices::cairo_pdf
  )
} else {
  message("Cairo PDF unavailable; PNG was saved. You can export PDF manually.")
}

# Optional: save summaries for table/supplementary reporting.
# write.csv(recency_summary, "CAV_ascertainment_recency_by_time.csv",
#           row.names = FALSE)

message("Saved: ", OUT_PNG)
if (isTRUE(capabilities("cairo"))) message("Saved: ", OUT_PDF)

###############################################################################
# NOTE ON INTERPRETATION
#
# The grey line/band characterises recency of known CAV assessments, including
# the assumed initial state at transplantation. A large upper quartile can
# indicate prolonged periods without reassessment, but the association with
# an observed-versus-predicted discrepancy is descriptive, not causal.
#
# If test has neither obs_yrs nor a true terminal row for every patient,
# reconstruct administrative end-of-follow-up dates first. Otherwise the
# at-risk denominator, and hence the recency distribution, may be biased.
###############################################################################
