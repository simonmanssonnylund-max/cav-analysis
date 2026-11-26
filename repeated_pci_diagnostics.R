# Repeated PCI: descriptive diagnostics for the FYND observation process
# Requires fit_adj, msm_dat9, and pci_pairs_max from the preceding analyses.
# Time is measured in years since transplant. Same-time PCI records are
# treated as one occasion. PCI at the current examination is excluded from
# its prior history; PCI at the preceding examination is included.

library(dplyr)
library(purrr)

stopifnot(
  all(c("lopnr", "years", "included_pci", "state") %in% names(msm_dat9)),
  all(c("lopnr", "years", "prev_years", "prev_state", "prev_max_state",
        "prior_pci", "years_since_pci", "pci_group", "can_decrease",
        "downward_consecutive", "persistent_decrease", "downward_max") %in%
        names(pci_pairs_max))
)

# ==============================================================================
# 1. IDENTIFY ALL PCI OCCASIONS ----
# ==============================================================================

cohort_ids <- unique(as.character(model.frame(fit_adj)[["(subject)"]]))

pci_events <- msm_dat9 %>%
  ungroup() %>%
  filter(included_pci %in% TRUE, is.finite(years)) %>%
  transmute(lopnr = as.character(lopnr), pci_year = years) %>%
  filter(lopnr %in% cohort_ids) %>%
  distinct(lopnr, pci_year) %>%
  group_by(lopnr) %>%
  arrange(pci_year, .by_group = TRUE) %>%
  mutate(
    pci_number = row_number(),
    interval_since_previous_pci = pci_year - lag(pci_year)
  ) %>%
  ungroup()

# ==============================================================================
# 2. DISTRIBUTION OF NUMBER OF PCIs PER PATIENT ----
# ==============================================================================

pci_number_summary <- pci_events %>%
  count(lopnr, name = "n_pci") %>%
  count(n_pci, name = "n_patients") %>%
  mutate(pct_patients = 100 * n_patients / sum(n_patients))

print(pci_number_summary, n = Inf)

# ==============================================================================
# 3. INTERVALS BETWEEN CONSECUTIVE PCIs ----
# ==============================================================================

pci_interval_summary <- pci_events %>%
  filter(pci_number > 1) %>%
  group_by(pci_number) %>%
  summarise(
    n_intervals = n(),
    min_years = min(interval_since_previous_pci),
    q25 = quantile(interval_since_previous_pci, 0.25),
    median_years = median(interval_since_previous_pci),
    q75 = quantile(interval_since_previous_pci, 0.75),
    max_years = max(interval_since_previous_pci),
    .groups = "drop"
  )

print(pci_interval_summary, n = Inf)

# ==============================================================================
# 4. NUMBER OF PREVIOUS PCIs AT EACH OBSERVATION ----
# ==============================================================================

pci_histories <- pci_events %>%
  group_by(lopnr) %>%
  summarise(pci_times = list(sort(pci_year)), .groups = "drop")

pci_pairs_repeat <- pci_pairs_max %>%
  ungroup() %>%
  mutate(lopnr = as.character(lopnr)) %>%
  left_join(pci_histories, by = "lopnr") %>%
  mutate(
    # Ensure that patients with no PCI have an empty numeric history.
    pci_times = map(pci_times, ~ sort(as.numeric(.x[!is.na(.x)]))),
    n_prior_pci = map2_int(years, pci_times, ~ sum(.y < .x)),
    n_pci_between = n_prior_pci -
      map2_int(prev_years, pci_times, ~ sum(.y < .x)),
    last_pci_year = map2_dbl(pci_times, n_prior_pci,
                             ~ if (.y > 0) .x[.y] else NA_real_),
    pci_history = case_when(
      n_prior_pci == 0 ~ "No previous PCI",
      n_prior_pci == 1 ~ "One previous PCI",
      n_prior_pci >= 2 ~ "Multiple previous PCIs"
    ),
    pci_history = factor(pci_history, levels = c(
      "No previous PCI", "One previous PCI", "Multiple previous PCIs"
    ))
  )

# ==============================================================================
# 5. VERIFY CONSISTENCY WITH EXISTING PCI VARIABLES ----
# ==============================================================================

pci_check <- pci_pairs_repeat %>%
  summarise(
    n_disagreement_pci = sum((n_prior_pci > 0) != (prior_pci == 1),
                            na.rm = TRUE),
    n_disagreement_time = sum(n_prior_pci > 0 &
      (is.na(years_since_pci) |
       abs(years - last_pci_year - years_since_pci) > 1e-8),
      na.rm = TRUE),
    n_missing_prior_pci = sum(is.na(prior_pci))
  )

print(pci_check)

# Do not silently mix incompatible definitions of PCI history.
if (any(unlist(pci_check) > 0)) {
  warning("PCI histories disagree with the existing model covariates. Check pci_history_disagreements before interpreting the summaries.")
}

pci_history_disagreements <- pci_pairs_repeat %>%
  filter(is.na(prior_pci) |
    (n_prior_pci > 0) != (prior_pci == 1) |
    (n_prior_pci > 0 & (is.na(years_since_pci) |
      abs(years - last_pci_year - years_since_pci) > 1e-8)))

# ==============================================================================
# 6. STRATIFY BY NUMBER OF PREVIOUS PCIs ----
# ==============================================================================

pci_summary_repeat <- pci_pairs_repeat %>%
  group_by(pci_history, pci_group) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_possible_decrease = sum(can_decrease),
    n_eligible_patients = n_distinct(lopnr[can_decrease]),
    n_downward_consecutive = sum(downward_consecutive),
    n_downward_persistent = sum(persistent_decrease),
    n_downward_max = sum(downward_max),
    n_patients_downward = n_distinct(lopnr[downward_max]),
    pct_downward = if_else(n_possible_decrease > 0,
      100 * n_downward_max / n_possible_decrease, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(pci_history, pci_group)

print(pci_summary_repeat, n = Inf, width = Inf)

# ==============================================================================
# 7. ADDITIONALLY STRATIFY BY HISTORICAL MAXIMUM ----
# ==============================================================================

pci_summary_repeat_strat <- pci_pairs_repeat %>%
  filter(can_decrease) %>%
  group_by(pci_history, pci_group, prev_max_state) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_downward = sum(downward_max),
    pct_downward = 100 * n_downward / n_pairs,
    .groups = "drop"
  ) %>%
  arrange(pci_history, prev_max_state, pci_group)

print(pci_summary_repeat_strat, n = Inf, width = Inf)

# ==============================================================================
# 8. CONSECUTIVE DECREASES FOLLOWING FIRST VS REPEATED PCI ----
# ==============================================================================
# One row per eligible consecutive assessment pair with intervening PCI.
# If several PCIs occurred between exams, the pair is classified according
# to the most recent intervention; n_multiple_intervening_pci flags these pairs.

pci_first_followup_summary <- pci_pairs_repeat %>%
  filter(n_pci_between > 0, prev_state %in% 2:3) %>%
  mutate(pci_type = if_else(n_prior_pci == 1, "First PCI", "Repeated PCI")) %>%
  group_by(pci_type, prev_state) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_multiple_intervening_pci = sum(n_pci_between > 1),
    n_downward = sum(downward_consecutive),
    pct_downward = 100 * n_downward / n_pairs,
    median_time_since_pci = median(years_since_pci),
    median_exam_interval = median(years - prev_years),
    .groups = "drop"
  )

print(pci_first_followup_summary, n = Inf, width = Inf)

# ==============================================================================
# 9. LONGITUDINAL HISTORIES OF REPEATED PCI PATIENTS ----
# ==============================================================================

repeat_pci_ids <- pci_events %>%
  count(lopnr) %>%
  filter(n > 1) %>%
  pull(lopnr)

repeat_pci_trajectories <- msm_dat9 %>%
  ungroup() %>%
  mutate(lopnr = as.character(lopnr)) %>%
  filter(lopnr %in% repeat_pci_ids,
         state %in% 1:3 | included_pci %in% TRUE) %>%
  select(lopnr, years, state, included_pci) %>%
  arrange(lopnr, years)

if (interactive()) View(repeat_pci_trajectories)

# Primary outputs: pci_number_summary, pci_check, pci_summary_repeat,
# pci_first_followup_summary. Review discrepancies before interpreting.
# These descriptive analyses do not identify causal PCI effects or latent
# misclassification probabilities, and assessment pairs are not independent.
