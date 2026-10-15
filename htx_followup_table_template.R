# ============================================================================
# HTx FOLLOW-UP TABLE: R CALCULATIONS -> SAS PRESENTATION
# ============================================================================
# Purpose
#   Create all rows shown in the follow-up / angiography / PCI table.
#   Perform the analyses in R; export a formatted table for presentation in SAS.
#
# What YOU need to edit
#   1. Load/prepare the relevant patient-, angiography-, and PCI-level datasets.
#   2. Replace NA_real_ in the THREE "CALCULATIONS" sections below.
#   3. Review denominator definitions and which columns are applicable.
#   4. Change output_dir to the desired directory.
#
# The script runs as-is, producing a COMPLETE BLANK TEMPLATE. As calculations
# are added, the relevant cells populate; uncalculated cells remain blank.
#
# Important
#   * For patient counts, count unique PATIENTS, not procedure records.
#   * For procedure counts, count procedures using the relevant procedures data.
#   * Pick the correct denominator separately for every percentage.
#   * "Entire follow-up" must be calculated independently where appropriate.
#     In particular, do not add distinct-patient counts from the two periods.
#   * All numbers are ILLUSTRATIVE PLACEHOLDERS (NA), not estimated results.
#
# File encoding: UTF-8 (labels include ≥ and †).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

# ============================================================================
# 1. LOAD AND PREPARE YOUR ANALYSIS DATA
# ============================================================================
# Example (replace with your own data-loading code):
# patients <- readRDS("path/to/patients.rds")
# angiographies <- readRDS("path/to/angiographies.rds")
# pci <- readRDS("path/to/pci.rds")
#
# Prepare any period-specific person-time, patient status, or procedure datasets
# here. The four sections below do not require a single common unit of analysis.
#
# NOTE: The period 'p0_7' means 0 <= follow-up < 7 years, and 'ge7' means
# follow-up >= 7 years. Define the boundary and censoring consistently with
# the Statistical Analysis Plan.


# ============================================================================
# 2. DEFINE TABLE LAYOUT / ROW METADATA -- NO ANALYSIS IN THIS SECTION
# ============================================================================
# tribble(): columns are marked with ~, then each subsequent line is ONE row.
# row_id    : stable link between the definition and the calculated statistic.
# section   : section code; full title is defined in section_def.
# label     : text for SAS to print in the first column.
# indent    : 1 for indented detail rows; 0 otherwise.
# row_type  : 'subheading' is a text-only row, not a statistic.
# footnote  : optional footnote marker to format in SAS.

section_def <- tribble(
  ~section,    ~section_title,
  "patients",  "Patients and follow-up",
  "angio",     "Investigation: coronary angiography",
  "diagnosis", "Diagnosis: highest angiographic finding during period, n (%)†",
  "pci",       "Treatment: PCI"
)

row_def <- tribble(
  ~row_id,         ~section,     ~label,                                            ~indent, ~row_type,   ~footnote,

  # Patients and follow-up
  "n_start",       "patients",    "Patients at start of period, N",                   0,       "data",      "",
  "age",           "patients",    "Age at start of period, years, median (IQR)",      0,       "data",      "",
  "male",          "patients",    "Male sex (%)",                                    0,       "data",      "",
  "fu",            "patients",    "Follow-up years, median (IQR)",                   0,       "data",      "",
  "py",            "patients",    "Total patient years",                             0,       "data",      "",
  "end_reason",    "patients",    "Reason follow-up ended, n (%)",                   0,       "subheading", "",
  "death",         "patients",    "Death during period, n (%)",                      1,       "data",      "",
  "retx",          "patients",    "Re-HTx",                                          1,       "data",      "",
  "reached7",      "patients",    "Reached 7 years alive without re-HTx",             1,       "data",      "",
  "alive_cutoff",  "patients",    "Alive at data cutoff",                            1,       "data",      "",

  # Investigation: coronary angiography
  "angio_total",   "angio",       "Total angiographies, n",                          0,       "data",      "",
  "angio_concl",   "angio",       "Conclusive, n",                                   0,       "data",      "",
  "angio_pat",     "angio",       "Number of patients with ≥1 angiography",           0,       "data",      "",
  "angio_rate",    "angio",       "Angiographies per patient year",                   0,       "data",      "",

  # Diagnosis: highest angiographic finding during period
  "fynd1",         "diagnosis",   "Fynd 1: normal or non-obstructive atheroma",        0,       "data",      "",
  "fynd_ge2",      "diagnosis",   "Fynd ≥2: single, multivessel or left main",        0,       "data",      "",
  "fynd2",         "diagnosis",   "Fynd 2: one-vessel disease",                      0,       "data",      "",
  "fynd_ge3",      "diagnosis",   "Fynd ≥3: multivessel or left main disease",        0,       "data",      "",

  # Treatment: PCI
  "pci_total",     "pci",         "Total PCI procedures, n",                         0,       "data",      "",
  "pci_success",   "pci",         "Successful, n (%)",                              0,       "data",      "",
  "pci_pat",       "pci",         "Patients with ≥1 PCI, n (%)",                     0,       "data",      "",
  "pci_rate",      "pci",         "PCI per 100 patient-years",                      0,       "data",      "",
  "pci_fynd",      "pci",         "Patients with Fynd ≥2 treated with PCI, n (%)",   0,       "data",      "c"
) %>%
  mutate(row_order = row_number()) %>%
  left_join(section_def, by = "section") %>%
  relocate(row_order, .before = row_id)

stopifnot(!anyDuplicated(row_def$row_id))
stopifnot(!anyNA(row_def$section_title))

# Map machine-readable period names to the intended printed SAS headers.
period_def <- tribble(
  ~period, ~column_heading,
  "p0_7",   "Follow-up 0–<7 years",
  "ge7",    "Follow-up ≥7 years",
  "entire", "Entire follow-up"
)


# ============================================================================
# 3. HELPER FUNCTIONS TO STORE NUMERIC RESULTS AND FORMAT DISPLAY STRINGS
# ============================================================================
# Your calculated results retain n, denominator, median, q1, q3, and value.
# 'display' is what gets printed in SAS. Missing components -> blank display.
# All helpers return a one-row tibble; bind_rows() stacks these vertically.

fmt_n <- function(n) {
  if (length(n) != 1L || is.na(n) || !is.finite(n)) return(NA_character_)
  sprintf("%.0f", n)
}

fmt_n_pct <- function(n, denominator, digits = 1) {
  if (length(n) != 1L || length(denominator) != 1L ||
      is.na(n) || is.na(denominator) ||
      !is.finite(n) || !is.finite(denominator) || denominator <= 0) {
    return(NA_character_)
  }
  sprintf(paste0("%s (%.", digits, "f%%)"), fmt_n(n), 100 * n / denominator)
}

fmt_median_iqr <- function(median, q1, q3, digits = 1) {
  x <- c(median, q1, q3)
  if (length(x) != 3L || anyNA(x) || any(!is.finite(x))) return(NA_character_)
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"), median, q1, q3)
}

fmt_value <- function(value, digits = 2) {
  if (length(value) != 1L || is.na(value) || !is.finite(value)) return(NA_character_)
  sprintf(paste0("%.", digits, "f"), value)
}

stat_row <- function(row_id, n = NA_real_, denominator = NA_real_,
                     median = NA_real_, q1 = NA_real_, q3 = NA_real_,
                     value = NA_real_, display = NA_character_) {
  tibble(row_id, n, denominator, median, q1, q3, value, display)
}

stat_count <- function(row_id, n = NA_real_) {
  stat_row(row_id, n = n, display = fmt_n(n))
}

stat_n_pct <- function(row_id, n = NA_real_, denominator = NA_real_, digits = 1) {
  stat_row(row_id, n = n, denominator = denominator,
           display = fmt_n_pct(n, denominator, digits = digits))
}

stat_median_iqr <- function(row_id, median = NA_real_, q1 = NA_real_,
                            q3 = NA_real_, digits = 1) {
  stat_row(row_id, median = median, q1 = q1, q3 = q3,
           display = fmt_median_iqr(median, q1, q3, digits = digits))
}

stat_value <- function(row_id, value = NA_real_, digits = 2) {
  stat_row(row_id, value = value, display = fmt_value(value, digits = digits))
}

stat_blank <- function(row_id) {
  stat_row(row_id)
}

# EXAMPLES of how to calculate and then enter values (commented out):
# n_male <- sum(patients_0_7$sex == "Male", na.rm = TRUE)
# denom_male <- sum(!is.na(patients_0_7$sex))
# age_q <- quantile(patients_0_7$age, c(.25, .50, .75), na.rm = TRUE)
# Then use stat_n_pct("male", n_male, denom_male)
# and stat_median_iqr("age", age_q[2], age_q[1], age_q[3]) below.
# You may also write your R calculations directly INSIDE a stat_*() call.


# ============================================================================
# 4. CALCULATIONS: FOLLOW-UP 0 TO <7 YEARS
# ============================================================================
# TODO: Add your period-specific computation code HERE, above the tibble.
# Replace each NA_real_ below with the corresponding calculated value,
# or with a named scalar created in this calculation section.
# Keep NA_real_ for not-yet-calculated / inapplicable cells.

results_p0_7 <- bind_rows(
  # --- PATIENTS AND FOLLOW-UP ---
  stat_count("n_start", n = NA_real_),  # Unique patients eligible at period start
  stat_median_iqr("age", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Age at start of period, not necessarily age at first HTx
  stat_n_pct("male", n = NA_real_, denominator = NA_real_),  # Male patients / patients with reported sex
  stat_median_iqr("fu", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Follow-up duration in THIS period
  stat_value("py", value = NA_real_, digits = 1),  # Sum of patient-years contributed in this period
  stat_blank("end_reason"),  # A subheading row only; no numeric result
  stat_n_pct("death", n = NA_real_, denominator = NA_real_),  # Deaths / appropriate at-risk patients
  stat_n_pct("retx", n = NA_real_, denominator = NA_real_),  # Re-HTx / appropriate at-risk patients
  stat_n_pct("reached7", n = NA_real_, denominator = NA_real_),  # Reached 7y alive w/o re-HTx; NA if not applicable
  stat_n_pct("alive_cutoff", n = NA_real_, denominator = NA_real_),  # Alive at data cutoff; define eligibility

  # --- INVESTIGATION: CORONARY ANGIOGRAPHY ---
  stat_count("angio_total", n = NA_real_),  # Number of angiography PROCEDURES
  stat_count("angio_concl", n = NA_real_),  # Number of conclusive angiography PROCEDURES
  stat_count("angio_pat", n = NA_real_),  # DISTINCT patients with >=1 angiography
  stat_value("angio_rate", value = NA_real_, digits = 2),  # Example: n_angiographies / patient_years

  # --- DIAGNOSIS: HIGHEST ANGIOGRAPHIC FINDING ---
  stat_n_pct("fynd1", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 1 / chosen denominator
  stat_n_pct("fynd_ge2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=2 / chosen denominator
  stat_n_pct("fynd2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 2 / chosen denominator
  stat_n_pct("fynd_ge3", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=3 / chosen denominator

  # --- TREATMENT: PCI ---
  stat_count("pci_total", n = NA_real_),  # Number of PCI PROCEDURES
  stat_n_pct("pci_success", n = NA_real_, denominator = NA_real_),  # Successful PCI / all eligible PCI procedures
  stat_n_pct("pci_pat", n = NA_real_, denominator = NA_real_),  # Distinct patients with >=1 PCI / eligible patients
  stat_value("pci_rate", value = NA_real_, digits = 2),  # Example: 100 * n_PCI_procedures / patient_years
  stat_n_pct("pci_fynd", n = NA_real_, denominator = NA_real_)  # Fynd >=2 patients receiving PCI / eligible Fynd >=2 patients
) %>%
  mutate(period = "p0_7")


# ============================================================================
# 5. CALCULATIONS: FOLLOW-UP >=7 YEARS
# ============================================================================
# TODO: Add your period-specific computation code HERE, above the tibble.
# Replace each NA_real_ below with the corresponding calculated value,
# or with a named scalar created in this calculation section.
# Keep NA_real_ for not-yet-calculated / inapplicable cells.
# Only include patients who actually enter the >=7y period.

results_ge7 <- bind_rows(
  # --- PATIENTS AND FOLLOW-UP ---
  stat_count("n_start", n = NA_real_),  # Unique patients eligible at period start
  stat_median_iqr("age", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Age at start of period, not necessarily age at first HTx
  stat_n_pct("male", n = NA_real_, denominator = NA_real_),  # Male patients / patients with reported sex
  stat_median_iqr("fu", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Follow-up duration in THIS period
  stat_value("py", value = NA_real_, digits = 1),  # Sum of patient-years contributed in this period
  stat_blank("end_reason"),  # A subheading row only; no numeric result
  stat_n_pct("death", n = NA_real_, denominator = NA_real_),  # Deaths / appropriate at-risk patients
  stat_n_pct("retx", n = NA_real_, denominator = NA_real_),  # Re-HTx / appropriate at-risk patients
  stat_n_pct("reached7", n = NA_real_, denominator = NA_real_),  # Reached 7y alive w/o re-HTx; NA if not applicable
  stat_n_pct("alive_cutoff", n = NA_real_, denominator = NA_real_),  # Alive at data cutoff; define eligibility

  # --- INVESTIGATION: CORONARY ANGIOGRAPHY ---
  stat_count("angio_total", n = NA_real_),  # Number of angiography PROCEDURES
  stat_count("angio_concl", n = NA_real_),  # Number of conclusive angiography PROCEDURES
  stat_count("angio_pat", n = NA_real_),  # DISTINCT patients with >=1 angiography
  stat_value("angio_rate", value = NA_real_, digits = 2),  # Example: n_angiographies / patient_years

  # --- DIAGNOSIS: HIGHEST ANGIOGRAPHIC FINDING ---
  stat_n_pct("fynd1", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 1 / chosen denominator
  stat_n_pct("fynd_ge2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=2 / chosen denominator
  stat_n_pct("fynd2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 2 / chosen denominator
  stat_n_pct("fynd_ge3", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=3 / chosen denominator

  # --- TREATMENT: PCI ---
  stat_count("pci_total", n = NA_real_),  # Number of PCI PROCEDURES
  stat_n_pct("pci_success", n = NA_real_, denominator = NA_real_),  # Successful PCI / all eligible PCI procedures
  stat_n_pct("pci_pat", n = NA_real_, denominator = NA_real_),  # Distinct patients with >=1 PCI / eligible patients
  stat_value("pci_rate", value = NA_real_, digits = 2),  # Example: 100 * n_PCI_procedures / patient_years
  stat_n_pct("pci_fynd", n = NA_real_, denominator = NA_real_)  # Fynd >=2 patients receiving PCI / eligible Fynd >=2 patients
) %>%
  mutate(period = "ge7")


# ============================================================================
# 6. CALCULATIONS: ENTIRE FOLLOW-UP
# ============================================================================
# TODO: Add your period-specific computation code HERE, above the tibble.
# Replace each NA_real_ below with the corresponding calculated value,
# or with a named scalar created in this calculation section.
# Keep NA_real_ for not-yet-calculated / inapplicable cells.
# Recalculate distinct-patient counts across the ENTIRE follow-up.
# DO NOT sum period-specific patient counts for distinct patients.

results_entire <- bind_rows(
  # --- PATIENTS AND FOLLOW-UP ---
  stat_count("n_start", n = NA_real_),  # Unique patients eligible at period start
  stat_median_iqr("age", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Age at start of period, not necessarily age at first HTx
  stat_n_pct("male", n = NA_real_, denominator = NA_real_),  # Male patients / patients with reported sex
  stat_median_iqr("fu", median = NA_real_, q1 = NA_real_, q3 = NA_real_),  # Follow-up duration in THIS period
  stat_value("py", value = NA_real_, digits = 1),  # Sum of patient-years contributed in this period
  stat_blank("end_reason"),  # A subheading row only; no numeric result
  stat_n_pct("death", n = NA_real_, denominator = NA_real_),  # Deaths / appropriate at-risk patients
  stat_n_pct("retx", n = NA_real_, denominator = NA_real_),  # Re-HTx / appropriate at-risk patients
  stat_n_pct("reached7", n = NA_real_, denominator = NA_real_),  # Reached 7y alive w/o re-HTx; NA if not applicable
  stat_n_pct("alive_cutoff", n = NA_real_, denominator = NA_real_),  # Alive at data cutoff; define eligibility

  # --- INVESTIGATION: CORONARY ANGIOGRAPHY ---
  stat_count("angio_total", n = NA_real_),  # Number of angiography PROCEDURES
  stat_count("angio_concl", n = NA_real_),  # Number of conclusive angiography PROCEDURES
  stat_count("angio_pat", n = NA_real_),  # DISTINCT patients with >=1 angiography
  stat_value("angio_rate", value = NA_real_, digits = 2),  # Example: n_angiographies / patient_years

  # --- DIAGNOSIS: HIGHEST ANGIOGRAPHIC FINDING ---
  stat_n_pct("fynd1", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 1 / chosen denominator
  stat_n_pct("fynd_ge2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=2 / chosen denominator
  stat_n_pct("fynd2", n = NA_real_, denominator = NA_real_),  # Patients in Fynd 2 / chosen denominator
  stat_n_pct("fynd_ge3", n = NA_real_, denominator = NA_real_),  # Patients in Fynd >=3 / chosen denominator

  # --- TREATMENT: PCI ---
  stat_count("pci_total", n = NA_real_),  # Number of PCI PROCEDURES
  stat_n_pct("pci_success", n = NA_real_, denominator = NA_real_),  # Successful PCI / all eligible PCI procedures
  stat_n_pct("pci_pat", n = NA_real_, denominator = NA_real_),  # Distinct patients with >=1 PCI / eligible patients
  stat_value("pci_rate", value = NA_real_, digits = 2),  # Example: 100 * n_PCI_procedures / patient_years
  stat_n_pct("pci_fynd", n = NA_real_, denominator = NA_real_)  # Fynd >=2 patients receiving PCI / eligible Fynd >=2 patients
) %>%
  mutate(period = "entire")


# ============================================================================
# 7. COMBINE RESULTS AND CHECK THAT NO TABLE ROWS HAVE BEEN OMITTED
# ============================================================================
results_long <- bind_rows(results_p0_7, results_ge7, results_entire)

# The unique key for the analytical results is (row_id, period).
if (anyDuplicated(results_long[c("row_id", "period")])) {
  stop("Duplicate (row_id, period) combinations: check the three results blocks.")
}

expected_keys <- tidyr::expand_grid(
  row_id = row_def$row_id,
  period = period_def$period
)

missing_keys <- anti_join(expected_keys, results_long, by = c("row_id", "period"))
extra_keys <- anti_join(results_long, expected_keys, by = c("row_id", "period"))

if (nrow(missing_keys) > 0L || nrow(extra_keys) > 0L) {
  print(missing_keys)
  print(extra_keys)
  stop("Rows are missing or unrecognized. Compare row_id values with row_def.")
}

# Optional audit view: all numeric components + the formatted text.
# View(results_long)


# ============================================================================
# 8. CREATE THE WIDE TABLE FOR SAS
# ============================================================================
results_wide <- results_long %>%
  select(row_id, period, display) %>%
  pivot_wider(
    names_from = period,
    values_from = display
  )

# All 23 detail rows, ordered and accompanied by metadata for SAS.
table_detail <- row_def %>%
  left_join(results_wide, by = "row_id") %>%
  arrange(row_order) %>%
  mutate(across(all_of(period_def$period), ~coalesce(.x, "")))

# Add 4 physical section-heading rows to simplify SAS presentation.
# SAS can use row_type == "section" for grey shading and spanning labels.
section_rows <- row_def %>%
  group_by(section, section_title) %>%
  summarise(row_order = min(row_order) - 0.5, .groups = "drop") %>%
  transmute(
    row_order = as.numeric(row_order),
    row_id = paste0("section_", section),
    section,
    label = section_title,
    indent = 0L,
    row_type = "section",
    footnote = "",
    section_title,
    p0_7 = "",
    ge7 = "",
    entire = ""
  )

table_sas <- bind_rows(table_detail, section_rows) %>%
  arrange(row_order) %>%
  select(
    row_order, row_id, section, row_type, label, indent, footnote,
    p0_7, ge7, entire
  )

stopifnot(nrow(table_detail) == 23L)
stopifnot(nrow(table_sas) == 27L)  # 23 characteristics + 4 section headers


# ============================================================================
# 9. EXPORT
# ============================================================================
# Set this to an existing folder you can access from RStudio Server.
# getwd() uses your current R working directory by default.
output_dir <- getwd()

# One CSV ready for presentation in SAS, and one audit CSV retaining raw values.
# readr::write_csv() saves UTF-8. Import using UTF-8 in SAS if necessary.
readr::write_csv(table_sas, file.path(output_dir, "htx_followup_table_for_sas.csv"), na = "")
readr::write_csv(results_long, file.path(output_dir, "htx_followup_table_numeric_audit.csv"), na = "")

message("Saved SAS-ready table to: ", file.path(output_dir, "htx_followup_table_for_sas.csv"))
message("Saved numeric audit file to: ", file.path(output_dir, "htx_followup_table_numeric_audit.csv"))

# View(table_sas)
# The three SAS result columns are p0_7, ge7, entire; customize their
# printed headers using period_def above.
# row_type: 'section' = shaded section header, 'subheading' = plain subheading,
#           'data' = ordinary statistic.
# indent: 1 = indent label, 0 = no indent.
# footnote: 'c' for the final PCI row, to be superscripted by SAS if desired.
