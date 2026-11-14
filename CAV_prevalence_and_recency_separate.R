###############################################################################
# CAV: OBSERVED VS MODEL-PREDICTED PREVALENCE + ASCERTAINMENT RECENCY
#
# Outputs
#   1. Figure_CAV_prevalence.{png,pdf}          - original 2-by-2 CAV plot
#   2. Figure_CAV_recency.{png,pdf}             - separate recency plot
#   3. Figure_CAV_prevalence_and_recency.{png,pdf} (if patchwork installed)
#   4. Table_CAV_recency_summary.csv
#
# Objects required in the R session
#   test         : final MSM-ready longitudinal data
#   prev_sahl_df : observed and model-predicted prevalence
#
# Expected columns:
#   test:          lopnr, years, state, unit
#                  optional: artificial, obs_yrs
#   prev_sahl_df:  time, state_name, observed, modelled
#
# IMPORTANT
#   * test$state must contain the ORIGINAL observed CAV state, not a forward-
#     filled state. Do NOT use statemax to identify examinations.
#   * state 1:3 = conclusive CAV observations; 4 = death/retransplantation;
#     99 = alive with unknown CAV state / artificial covariate-update row.
#   * The initial CAV-free state at year 0 is ASSUMED, not an angiography.
#   * Recency summaries are among patients alive, not retransplanted, and
#     under follow-up at the respective time; prevalence denominators differ.
#   * A patient with no post-transplant conclusive examination has recency = t.
###############################################################################


# =============================================================================
# 1. PACKAGES, CONSTANTS AND INPUT CHECKS
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)

TIME_MAX <- 7
TIME_STEP <- 0.05
SAHL_LEVEL <- "Göteborg"
PLOT_FONT <- "Aptos"
SHOW_P10_P90 <- TRUE           # FALSE = show median and IQR only

stopifnot(
  exists("test"),
  exists("prev_sahl_df"),
  all(c("lopnr", "years", "state", "unit") %in% names(test)),
  all(c("time", "state_name", "observed", "modelled") %in%
        names(prev_sahl_df))
)

STATE_LABELS <- c(
  "max(FYND) = 1",
  "max(FYND) = 2",
  "max(FYND) >= 3",
  "Death/Re-HTx"
)

if (!setequal(unique(as.character(prev_sahl_df$state_name)), STATE_LABELS)) {
  stop(
    "Unexpected facet labels. Inspect unique(prev_sahl_df$state_name), ",
    "then update STATE_LABELS in Section 1."
  )
}

prev_sahl_df <- prev_sahl_df %>%
  mutate(state_name = factor(as.character(state_name), levels = STATE_LABELS))


# =============================================================================
# 2. PREPARE THE FINAL LONGITUDINAL DATA
# =============================================================================

# Raw state MUST identify true observations. Do not replace with statemax.
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

# If obs_yrs is unavailable, the maximum observed time is a valid substitute
# ONLY when there is a terminal follow-up row for every patient.
if ("obs_yrs" %in% names(dat)) {
  dat <- dat %>% mutate(obs_yrs = as.numeric(obs_yrs))
} else {
  warning(
    "test has no obs_yrs. Using max(years) per patient. This is valid ",
    "ONLY when each patient has a true terminal follow-up row."
  )
  dat <- dat %>%
    group_by(lopnr) %>%
    mutate(obs_yrs = max(years, na.rm = TRUE)) %>%
    ungroup()
}

# Same center as the published Sahlgrenska prevalence plot.
dat_sahl <- dat %>%
  filter(as.character(unit) == SAHL_LEVEL)

if (nrow(dat_sahl) == 0L) {
  stop("No patients found for SAHL_LEVEL = '", SAHL_LEVEL, "'.")
}

if (!"artificial" %in% names(test)) {
  message(
    "No artificial indicator: assuming all inserted non-examination rows ",
    "remain coded state 99 in test$state."
  )
}


# =============================================================================
# 3. RECONSTRUCT THE TIMES OF CONCLUSIVE CAV ASSESSMENTS
# =============================================================================

# Keep ALL genuine assessments (also those with unchanged CAV grade).
# Exclude unknown grades, exact absorbing events, synthetic observations,
# and the assumed CAV-free baseline at year zero.
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
# 4. RECONSTRUCT FOLLOW-UP AND EXACT ABSORBING EVENT TIMES
# =============================================================================

recency_patients <- dat_sahl %>%
  group_by(lopnr) %>%
  summarise(
    obs_yrs = if (all(is.na(obs_yrs))) NA_real_ else
      max(obs_yrs, na.rm = TRUE),
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


# =============================================================================
# 5. RECENCY FOR EVERY PATIENT STILL AT RISK AT EACH TIME POINT
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
    # Year 0 is an assumed CAV-free starting state, not an angiography.
    last_assessment = coalesce(exam_year, 0),
    years_since_assessment = time - last_assessment,
    baseline_only = is.na(exam_year)
  )

stopifnot(
  all(recency_patient_times$years_since_assessment >= -1e-10),
  all(recency_patient_times$years_since_assessment <= TIME_MAX + 1e-10)
)


# =============================================================================
# 6. RECENCY SUMMARY: MEDIAN, IQR, 10TH-90TH AND DATA COMPLETENESS
# =============================================================================

recency_summary <- recency_patient_times %>%
  group_by(time) %>%
  summarise(
    n_at_risk = n(),
    median_recency = median(years_since_assessment),
    p10 = unname(quantile(years_since_assessment, 0.10)),
    q25 = unname(quantile(years_since_assessment, 0.25)),
    q75 = unname(quantile(years_since_assessment, 0.75)),
    p90 = unname(quantile(years_since_assessment, 0.90)),
    pct_baseline_only = 100 * mean(baseline_only),
    .groups = "drop"
  )

# Print the clinically relevant time points, including size of the risk set.
recency_key_times <- recency_summary %>%
  filter(round(time, 2) %in% c(1, 3, 6, 7)) %>%
  mutate(across(
    c(median_recency, p10, q25, q75, p90, pct_baseline_only),
    ~ round(.x, 2)
  ))

print(recency_key_times)
write.csv(recency_summary, "Table_CAV_recency_summary.csv", row.names = FALSE)


# =============================================================================
# 7. ORIGINAL PREVALENCE FIGURE: TWO BY TWO, ORIGINAL COLORS
# =============================================================================

# Keep the original free y-ranges, WITHOUT any secondary y-axes.
make_prevalence_scale <- function(ymax) {
  scale_y_continuous(
    limits = c(0, ymax),
    breaks = seq(0, ymax, by = if (ymax == 10) 1 else 10),
    expand = expansion(mult = c(0, 0.02)),
    oob = scales::squish
  )
}

p_prev_sahl <- ggplot(prev_sahl_df, aes(x = time)) +
  geom_line(
    aes(y = observed, colour = "Observed", linetype = "Observed"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = modelled,
        colour = "Model-predicted", linetype = "Model-predicted"),
    linewidth = 0.9
  ) +
  facet_wrap(~ state_name, ncol = 2, scales = "free_y") +
  ggh4x::facetted_pos_scales(
    y = list(
      state_name == "max(FYND) = 1" ~ make_prevalence_scale(100),
      state_name == "max(FYND) = 2" ~ make_prevalence_scale(10),
      state_name == "max(FYND) >= 3" ~ make_prevalence_scale(10),
      state_name == "Death/Re-HTx" ~ make_prevalence_scale(100)
    )
  ) +
  scale_x_continuous(
    limits = c(0, TIME_MAX),
    breaks = 0:7,
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_colour_manual(
    values = c(
      "Observed" = "#9eb8b7",            # original muted blue-green
      "Model-predicted" = "#651220"      # original burgundy
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Observed" = "solid",
      "Model-predicted" = "dashed"
    )
  ) +
  labs(
    x = "Years after heart transplantation",
    y = "Probability in state (%)",
    colour = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 11, base_family = PLOT_FONT) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey85"),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey85"),
    panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.5),
    strip.background = element_rect(
      fill = "white", colour = "grey40", linewidth = 0.5
    ),
    strip.text = element_text(
      face = "bold", size = 11, margin = margin(t = 5, b = 5)
    ),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.width = grid::unit(1.5, "cm"),
    plot.margin = margin(t = 10, r = 10, b = 8, l = 10)
  )


# =============================================================================
# 8. SEPARATE RECENCY FIGURE: MEDIAN + IQR + OPTIONAL 10TH-90TH
# =============================================================================

p_recency <- ggplot(recency_summary, aes(x = time))

# Outer band communicates the tails without dominating the figure.
if (SHOW_P10_P90) {
  p_recency <- p_recency +
    geom_ribbon(
      aes(ymin = p10, ymax = p90, fill = "10th–90th percentile"),
      alpha = 0.45,
      colour = NA
    )
}

# IQR is more prominent than the outer 10th–90th percentile band.
p_recency <- p_recency +
  geom_ribbon(
    aes(ymin = q25, ymax = q75, fill = "25th–75th percentile (IQR)"),
    alpha = 0.65,
    colour = NA
  ) +
  geom_line(
    aes(y = median_recency, colour = "Median"),
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, TIME_MAX),
    breaks = 0:7,
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, TIME_MAX),
    breaks = 0:7,
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "10th–90th percentile" = "#D0D0D0",
      "25th–75th percentile (IQR)" = "#85959B"
    ),
    breaks = if (SHOW_P10_P90) {
      c("25th–75th percentile (IQR)", "10th–90th percentile")
    } else {
      "25th–75th percentile (IQR)"
    },
    name = NULL
  ) +
  scale_colour_manual(values = c("Median" = "#30383C"), name = NULL) +
  labs(
    x = "Years after heart transplantation",
    y = "Years since last CAV state ascertainment",
    caption = paste(
      "Recency is evaluated among patients alive, not retransplanted and still under follow-up.",
      "Time zero is the assumed CAV-free starting state; only conclusive follow-up",
      "angiographies reset ascertainment recency."
    )
  ) +
  guides(
    colour = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  theme_minimal(base_size = 11, base_family = PLOT_FONT) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey85"),
    panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey85"),
    panel.border = element_rect(fill = NA, colour = "grey40", linewidth = 0.5),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    plot.caption = element_text(size = 9, hjust = 0),
    plot.margin = margin(t = 5, r = 10, b = 8, l = 10)
  )


# =============================================================================
# 9. DISPLAY AND SAVE THE TWO FIGURES INDIVIDUALLY
# =============================================================================

print(p_prev_sahl)
print(p_recency)

save_figure <- function(plot, filename, width, height) {
  ggsave(
    paste0(filename, ".png"),
    plot, width = width, height = height,
    dpi = 600, bg = "white"
  )
  if (isTRUE(capabilities("cairo"))) {
    ggsave(
      paste0(filename, ".pdf"),
      plot, width = width, height = height,
      device = grDevices::cairo_pdf
    )
  } else {
    message("Cairo unavailable; saved PNG only for: ", filename)
  }
}

save_figure(p_prev_sahl, "Figure_CAV_prevalence", width = 9, height = 7)
save_figure(p_recency, "Figure_CAV_recency", width = 9, height = 4.5)


# =============================================================================
# 10. OPTIONAL VERTICAL COMPOSITE: PREVALENCE ABOVE, RECENCY BELOW
# =============================================================================

# The recency plot is an independent, full-width panel below the 2-by-2 plot.
# There are NO secondary axes and no repeated recency lines across the facets.
# The two component plots above remain available and are saved separately.
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_combined <- patchwork::wrap_plots(
    p_prev_sahl,
    p_recency,
    ncol = 1,
    heights = c(2.1, 1)
  )

  print(p_combined)
  save_figure(
    p_combined,
    "Figure_CAV_prevalence_and_recency",
    width = 10,
    height = 12
  )
} else {
  message(
    "For a combined vertical layout, install patchwork with ",
    "install.packages('patchwork'). Individual figures have been saved."
  )
}


# =============================================================================
# 11. INTERPRETATION AND ADDITIONAL CHECKS
# =============================================================================

# The 90th percentile may equal the full time since transplantation if >=10%
# of the at-risk patients have NEVER had a conclusive post-transplant exam.
# Inspect pct_baseline_only and n_at_risk when interpreting rising percentiles.
# A widening gap between observed and predicted prevalence as recency increases
# is suggestive but does not establish that delayed ascertainment caused it.
#
# If test lacks obs_yrs AND terminal follow-up rows, end-of-follow-up times must
# be reconstructed before interpreting the risk-set summaries.
# If test$state contains forward-filled values on artificial rows, reconstruct
# original assessment indicators before using this script.
###############################################################################
