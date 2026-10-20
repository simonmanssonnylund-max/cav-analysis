# SCAAR segment data audit - updated 2026-10-20
# Compares test2.6 with the segment-level SCAAR data.
# All patients in test2.6 are retained for cohort-level checks, even without dates.
# Dated patient visits are used only when an examination date is required.
# Run after loading test2.6 and the segment dataset into your R session.

library(dplyr)
library(tidyr)

# 1. Set the input dataset and check variable names

segment_data <- scaar_segments  # Change to your segment dataset object name

required_visits <- c("lopnr", "angio_date", "from_strax", "event", "angio_fynd", "cav_proxy")
required_segments <- c("LopNr", "pseudo_MCEID", "pseudo_SID", "INTERDAT", "SEGMENT", "FINDING")

stopifnot(all(required_visits %in% names(test2.6)))
stopifnot(all(required_segments %in% names(segment_data)))

# The angiography flag has appeared under two names in earlier code.
angio_flag <- if("included_angiography" %in% names(test2.6)) {
  "included_angiography"
} else if("included_angio" %in% names(test2.6)) {
  "included_angio"
} else {
  stop("Neither included_angiography nor included_angio exists in test2.6.")
}
stopifnot("included_pci" %in% names(test2.6))

# Keep every source row, including those without angiography dates.
# Dates must be Date objects or ISO date strings.
visits_raw <- test2.6 %>%
  ungroup() %>%
  mutate(
    lopnr = as.character(lopnr),
    angio_date = as.Date(angio_date),
    included_angiography = .data[[angio_flag]] %in% TRUE,
    included_pci = included_pci %in% TRUE,
    angio_fynd = suppressWarnings(as.integer(as.character(angio_fynd))),
    cav_proxy = suppressWarnings(as.numeric(as.character(cav_proxy)))
  )

seg_raw <- segment_data %>%
  transmute(
    lopnr = as.character(LopNr),
    pseudo_MCEID, pseudo_SID,
    angio_date = as.Date(INTERDAT),
    SEGMENT = suppressWarnings(as.integer(as.character(SEGMENT))),
    FINDING = suppressWarnings(as.integer(as.character(FINDING)))
  )

# Report missing identifiers/dates before excluding them from date comparisons.
missing_keys <- bind_rows(
  visits_raw %>% summarise(
    source = "test2.6", rows = n(),
    missing_id = sum(is.na(lopnr) | lopnr == ""),
    missing_date = sum(is.na(angio_date))
  ),
  seg_raw %>% summarise(
    source = "Segments", rows = n(),
    missing_id = sum(is.na(lopnr) | lopnr == ""),
    missing_date = sum(is.na(angio_date))
  )
)

# Do not filter on angio_date when defining the patient cohort.
# The raw objects above still contain every original row, including missing IDs.
visits <- visits_raw %>% filter(!is.na(lopnr), lopnr != "")
cohort_ids <- visits %>% distinct(lopnr)

# Explicitly identify patients with no dated registration in test2.6.
patient_registration_status <- visits %>%
  group_by(lopnr) %>%
  summarise(
    n_rows = n(),
    n_dated_rows = sum(!is.na(angio_date)),
    n_unique_dates = n_distinct(angio_date[!is.na(angio_date)]),
    has_dated_registration = n_dated_rows > 0,
    has_included_angiography = any(included_angiography),
    has_included_pci = any(included_pci),
    .groups = "drop"
  )

cohort_summary <- patient_registration_status %>%
  summarise(
    n_patients = n(),
    n_without_dated_registration = sum(!has_dated_registration),
    n_with_dated_registration = sum(has_dated_registration),
    n_with_included_angiography = sum(has_included_angiography),
    n_with_included_pci = sum(has_included_pci)
  )

cohort_registration_distribution <- patient_registration_status %>%
  count(has_dated_registration, has_included_angiography, has_included_pci,
        name = "n_patients")

# These separate datasets are only for analyses requiring an examination date.
visits_dated <- visits %>% filter(!is.na(angio_date))
visits_undated <- visits %>% filter(is.na(angio_date))
seg_ids <- seg_raw %>%
  filter(!is.na(lopnr), lopnr != "") %>%
  distinct(lopnr)
seg <- seg_raw %>% filter(!is.na(lopnr), lopnr != "", !is.na(angio_date))
keys <- c("lopnr", "angio_date")

# 2. Check identifiers, stenosis codes, and duplicate registrations

invalid_codes <- seg %>%
  filter(
    (!is.na(SEGMENT) & !SEGMENT %in% 1:20) |
    (!is.na(FINDING) & !FINDING %in% 0:6)
  )

# There may be several segment records for a single examination.
segment_duplicates <- seg %>%
  count(lopnr, angio_date, SEGMENT, name = "n_rows") %>%
  filter(n_rows > 1) %>%
  arrange(desc(n_rows), lopnr, angio_date)

segment_conflicts <- seg %>%
  group_by(lopnr, angio_date, SEGMENT) %>%
  summarise(
    n_rows = n(),
    n_distinct_findings = n_distinct(FINDING, na.rm = TRUE),
    findings = paste(sort(unique(FINDING[!is.na(FINDING)])), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_distinct_findings > 1)

# Does an examination date have several MCEID or SID values?
exam_identifiers <- seg %>%
  group_by(lopnr, angio_date) %>%
  summarise(
    n_MCEID = n_distinct(pseudo_MCEID, na.rm = TRUE),
    n_SID = n_distinct(pseudo_SID, na.rm = TRUE),
    .groups = "drop"
  )

identifier_distribution <- exam_identifiers %>% count(n_MCEID, n_SID, sort = TRUE)

# 3. Aggregate the segment records to examination-level summaries

seg_visits <- seg %>%
  group_by(lopnr, angio_date) %>%
  summarise(
    n_segment_rows = n(),
    n_segments = n_distinct(SEGMENT, na.rm = TRUE),
    n_with_finding = n_distinct(SEGMENT[!is.na(FINDING)], na.rm = TRUE),
    n_coded = n_distinct(SEGMENT[FINDING %in% 1:6], na.rm = TRUE),
    n_ge50 = n_distinct(SEGMENT[FINDING %in% 3:6], na.rm = TRUE),
    n_ge70 = n_distinct(SEGMENT[FINDING %in% 4:6], na.rm = TRUE),
    n_code0 = sum(FINDING %in% 0),
    .groups = "drop"
  )

segment_count_distribution <- seg_visits %>% count(n_segments, name = "n_exams")

segment_completeness <- seg_visits %>%
  summarise(
    n_exams = n(),
    n_all20 = sum(n_segments == 20),
    pct_all20 = 100 * mean(n_segments == 20),
    n_any_coded = sum(n_coded > 0),
    pct_any_coded = 100 * mean(n_coded > 0),
    median_segments = median(n_segments),
    q1_segments = unname(quantile(n_segments, 0.25)),
    q3_segments = unname(quantile(n_segments, 0.75))
  )

# Number of examinations where each specific segment is documented.
# A missing segment record does not necessarily mean the vessel was not imaged.
segment_by_visit <- seg %>%
  filter(!is.na(SEGMENT), SEGMENT %in% 1:20) %>%
  group_by(lopnr, angio_date, SEGMENT) %>%
  summarise(
    n_rows = n(),
    n_codes = n_distinct(FINDING, na.rm = TRUE),
    FINDING = if(all(is.na(FINDING))) NA_integer_ else max(FINDING, na.rm = TRUE),
    .groups = "drop"
  )
# Maximum is used only for descriptive checks if a segment has multiple codes.
# Review segment_conflicts before using these values for formal CAV classification.

segment_frequency <- segment_by_visit %>%
  group_by(SEGMENT) %>%
  summarise(
    n_exams = n(),
    n_coded = sum(FINDING %in% 1:6),
    n_ge50 = sum(FINDING %in% 3:6),
    n_ge70 = sum(FINDING %in% 4:6),
    .groups = "drop"
  ) %>%
  mutate(pct_exams = 100 * n_exams / nrow(seg_visits)) %>%
  arrange(SEGMENT)

finding_distribution <- seg %>% count(FINDING, sort = TRUE)

# 4. Compare patient and examination dates between datasets

# Keep all patients here; do not filter from_strax in the main comparison.
# source-only subsets below are additional diagnostics, not exclusions.
only_one_grade <- function(x) {
  u <- unique(x[!is.na(x)])
  if(length(u) == 1) u else NA_integer_
}

visits_by_date <- visits_dated %>%
  group_by(lopnr, angio_date) %>%
  summarise(
    n_source_rows = n(),
    has_strax = any(from_strax %in% 1),
    has_scaar = any(from_strax %in% 0),
    included_angiography = any(included_angiography),
    included_pci = any(included_pci),
    n_fynd_codes = n_distinct(angio_fynd, na.rm = TRUE),
    angio_fynd = only_one_grade(angio_fynd),
    .groups = "drop"
  )

# Also inspect date duplication and conflicting FYND codes in test2.6.
duplicate_visit_dates <- visits_by_date %>% filter(n_source_rows > 1)
conflicting_visit_fynd <- visits_by_date %>% filter(n_fynd_codes > 1)

scaar_dates <- visits_dated %>%
  filter(from_strax %in% 0) %>%
  distinct(lopnr, angio_date)
all_dates <- visits_by_date %>% select(all_of(keys))
segment_dates <- seg_visits %>% select(all_of(keys))
# cohort_ids was created from ALL patients above, including those with no dates.

# Compare both all test2.6 dates and the SCAAR-source dates separately.
match_summary <- bind_rows(
  tibble(
    comparison = "All test2.6 dates",
    n_reference = nrow(all_dates),
    n_matched = nrow(inner_join(all_dates, segment_dates, by = keys))
  ),
  tibble(
    comparison = "SCAAR-source dates only",
    n_reference = nrow(scaar_dates),
    n_matched = nrow(inner_join(scaar_dates, segment_dates, by = keys))
  )
) %>%
  mutate(pct_matched = if_else(n_reference > 0, 100 * n_matched / n_reference, NA_real_))

only_test26 <- anti_join(all_dates, segment_dates, by = keys)
only_scaar_dates <- anti_join(scaar_dates, segment_dates, by = keys)

# Restrict segment-only dates to patients represented in test2.6.
# Some can still lie outside the follow-up period, which needs separate checking.
only_segment_dates <- segment_dates %>%
  semi_join(cohort_ids, by = "lopnr") %>%
  anti_join(all_dates, by = keys)

# Patient matching is independent of the availability of examination dates.
patient_overlap <- tibble(
  patients_test26 = nrow(cohort_ids),
  patients_segments = nrow(seg_ids),
  patients_in_both = nrow(inner_join(cohort_ids, seg_ids, by = "lopnr")),
  patients_only_test26 = nrow(anti_join(cohort_ids, seg_ids, by = "lopnr")),
  patients_only_segments = nrow(anti_join(seg_ids, cohort_ids, by = "lopnr"))
)

# Distinguish patients without any dated examination from patients whose
# dated examinations fail to match segment-level registrations.
patient_segment_coverage <- patient_registration_status %>%
  mutate(has_segment_record = lopnr %in% seg_ids$lopnr) %>%
  count(has_dated_registration, has_segment_record, name = "n_patients")

# Explore whether unmatched SCAAR-source dates have nearby segment dates.
# A near match is not automatically the same examination.
nearby_dates <- only_scaar_dates %>%
  inner_join(
    segment_dates %>% rename(segment_date = angio_date),
    by = "lopnr", relationship = "many-to-many"
  ) %>%
  mutate(days_apart = as.integer(segment_date - angio_date)) %>%
  filter(abs(days_apart) <= 30) %>%
  group_by(lopnr, angio_date) %>%
  slice_min(abs(days_apart), n = 1, with_ties = FALSE) %>%
  ungroup()

nearby_day_distribution <- nearby_dates %>% count(days_apart, sort = TRUE)

# 5. Does segment coverage depend on source, procedure, or FYND?

comparison <- visits_by_date %>%
  left_join(seg_visits, by = keys, relationship = "one-to-one") %>%
  mutate(
    has_segments = !is.na(n_segments),
    source = case_when(
      has_scaar & has_strax ~ "Both",
      has_scaar ~ "SCAAR only",
      has_strax ~ "STRAX only",
      TRUE ~ "Other"
    ),
    procedure = case_when(
      included_angiography & included_pci ~ "Angiography + PCI",
      included_angiography ~ "Angiography only",
      included_pci ~ "PCI only",
      TRUE ~ "Neither"
    )
  )

coverage_summary <- comparison %>%
  group_by(source, procedure) %>%
  summarise(
    n_exams = n(),
    n_with_segments = sum(has_segments),
    pct_with_segments = 100 * mean(has_segments),
    .groups = "drop"
  )

coverage_by_fynd <- comparison %>%
  group_by(angio_fynd) %>%
  summarise(
    n_exams = n(),
    n_with_segments = sum(has_segments),
    pct_with_segments = 100 * mean(has_segments),
    median_n_segments = if(all(is.na(n_segments))) NA_real_ else median(n_segments, na.rm = TRUE),
    n_with_ge50 = sum(n_ge50 > 0, na.rm = TRUE),
    n_with_ge70 = sum(n_ge70 > 0, na.rm = TRUE),
    .groups = "drop"
  )

# Restrict to original SCAAR registrations when examining consistency with FYND.
# This avoids mistaking a STRAX-only finding for a SCAAR-coded examination.
scaar_fynd <- visits_dated %>%
  filter(from_strax %in% 0) %>%
  group_by(lopnr, angio_date) %>%
  summarise(
    n_fynd_codes = n_distinct(angio_fynd, na.rm = TRUE),
    angio_fynd = only_one_grade(angio_fynd),
    .groups = "drop"
  )

scaar_fynd_comparison <- scaar_fynd %>%
  left_join(seg_visits, by = keys, relationship = "one-to-one")

# Review these, but do not automatically label them registration errors.
# A threshold/segment mismatch or procedure-related coding can be relevant.
potential_discrepancies <- scaar_fynd_comparison %>%
  filter(angio_fynd == 1, n_ge50 > 0)

# 6. Explore segment records around decreases in cav_proxy

# Only non-missing CAV-proxy observations define transitions.
# All source rows remain available in visits; only dated rows enter this analysis.
reversals <- visits_dated %>%
  arrange(lopnr, angio_date, desc(from_strax)) %>%
  group_by(lopnr) %>%
  filter(!is.na(cav_proxy)) %>%
  mutate(
    previous_date = lag(angio_date),
    previous_cav = lag(cav_proxy),
    previous_source = lag(from_strax)
  ) %>%
  filter(!is.na(previous_cav), cav_proxy < previous_cav) %>%
  ungroup()

reversal_summary <- reversals %>%
  group_by(previous_cav, cav_proxy) %>%
  summarise(
    n_decreases = n(),
    n_patients = n_distinct(lopnr),
    .groups = "drop"
  ) %>%
  arrange(desc(previous_cav), cav_proxy)

reversal_sources <- reversals %>% count(previous_source, from_strax, sort = TRUE)

# Compare how many segments were recorded before and after each decrease.
reversal_coverage <- reversals %>%
  left_join(
    seg_visits %>%
      transmute(lopnr, previous_date = angio_date,
                previous_n_segments = n_segments,
                previous_n_ge50 = n_ge50,
                previous_n_ge70 = n_ge70),
    by = c("lopnr", "previous_date"), relationship = "many-to-one"
  ) %>%
  left_join(
    seg_visits %>%
      transmute(lopnr, angio_date,
                current_n_segments = n_segments,
                current_n_ge50 = n_ge50,
                current_n_ge70 = n_ge70),
    by = keys, relationship = "many-to-one"
  )

reversal_coverage_summary <- reversal_coverage %>%
  summarise(
    n_decreases = n(),
    n_with_both_segment_visits = sum(
      !is.na(previous_n_segments) & !is.na(current_n_segments)
    ),
    n_missing_previous_segments = sum(is.na(previous_n_segments)),
    n_missing_current_segments = sum(is.na(current_n_segments))
  )

# Compare each segment on both examinations for every reversal.
# NA means no recorded grade for that segment; it is not proof of no stenosis.
reversal_pairs <- reversals %>% distinct(lopnr, previous_date, angio_date)

previous_segments <- reversal_pairs %>%
  inner_join(
    segment_by_visit %>%
      transmute(lopnr, previous_date = angio_date, SEGMENT,
                previous_finding = FINDING),
    by = c("lopnr", "previous_date"), relationship = "many-to-many"
  )

current_segments <- reversal_pairs %>%
  inner_join(
    segment_by_visit %>%
      transmute(lopnr, angio_date, SEGMENT,
                current_finding = FINDING),
    by = keys, relationship = "many-to-many"
  )

reversal_segments <- full_join(
  previous_segments, current_segments,
  by = c("lopnr", "previous_date", "angio_date", "SEGMENT")
) %>%
  mutate(
    previously_ge50 = previous_finding %in% 3:6,
    currently_ge50 = current_finding %in% 3:6,
    previously_ge70 = previous_finding %in% 4:6,
    currently_ge70 = current_finding %in% 4:6,
    missing_on_current = !is.na(previous_finding) & is.na(current_finding),
    lower_recorded_category = !is.na(previous_finding) &
      !is.na(current_finding) & current_finding < previous_finding
  ) %>%
  arrange(lopnr, previous_date, angio_date, SEGMENT)

# Potentially informative segment changes for closer inspection.
reversal_segment_changes <- reversal_segments %>%
  filter(
    missing_on_current | lower_recorded_category |
    (previously_ge50 & !currently_ge50)
  )

# 7. Print concise results (the detailed objects remain in your environment)

cat("\nMissing patient IDs and dates:\n")
print(missing_keys)
cat("\nCohort counts (including patients without examination dates):\n")
print(cohort_summary)
cat("\nCohort registration distribution:\n")
print(cohort_registration_distribution, n = Inf)
cat("\nPatient segment coverage, with undated patients included:\n")
print(patient_segment_coverage, n = Inf)
cat("\nPatient overlap:\n")
print(patient_overlap)
cat("\nDate overlap:\n")
print(match_summary)
cat("\nNearby unmatched dates (within 30 days):\n")
print(nearby_day_distribution)
cat("\nSegment completeness:\n")
print(segment_completeness)
cat("\nNumber of segments per examination:\n")
print(segment_count_distribution, n = Inf)
cat("\nCoverage by source and procedure:\n")
print(coverage_summary, n = Inf)
cat("\nCoverage by overall angiographic finding:\n")
print(coverage_by_fynd, n = Inf)
cat("\nFrequency of individual segments:\n")
print(segment_frequency, n = Inf)
cat("\nFINDING coding:\n")
print(finding_distribution, n = Inf)
cat("\nDuplicate segment groups and conflicting findings:\n")
print(tibble(duplicate_groups = nrow(segment_duplicates), conflicting_groups = nrow(segment_conflicts)))
cat("\nDecreases by type:\n")
print(reversal_summary, n = Inf)
cat("\nSegment coverage around decreases:\n")
print(reversal_coverage_summary)

# Optional exports; uncomment if useful.
# write.csv(patient_registration_status, "scaar_patient_registration_status.csv", row.names = FALSE)
# write.csv(cohort_summary, "scaar_cohort_summary.csv", row.names = FALSE)
# write.csv(visits_undated, "scaar_undated_visit_rows.csv", row.names = FALSE)
# write.csv(match_summary, "scaar_segment_date_overlap.csv", row.names = FALSE)
# write.csv(only_scaar_dates, "scaar_dates_without_segments.csv", row.names = FALSE)
# write.csv(only_segment_dates, "segment_dates_without_test26.csv", row.names = FALSE)
# write.csv(comparison, "scaar_segment_coverage_per_visit.csv", row.names = FALSE)
# write.csv(reversal_coverage, "scaar_reversal_segment_coverage.csv", row.names = FALSE)
# write.csv(reversal_segments, "scaar_reversal_segment_details.csv", row.names = FALSE)
