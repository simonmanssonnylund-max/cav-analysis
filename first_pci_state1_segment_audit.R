# SCAAR segment findings at first PCI with current and previous state 1
# Requires first_pci and scaar_segments in the R session.
# FINDING thresholds follow the earlier SCAAR segment audit coding.

library(dplyr)

# ==============================================================================
# 1. IDENTIFY RELEVANT FIRST PCI OCCASIONS ----
# ==============================================================================

pci_state1 <- first_pci %>%
  ungroup() %>%
  filter(state == 1, previous_FYND == 1) %>%
  mutate(
    lopnr = as.character(lopnr),
    angio_date = as.Date(angio_date)
  )

# ==============================================================================
# 2. RETRIEVE SEGMENT DATA ON THE PCI DATE ----
# ==============================================================================

pci_state1_segments <- pci_state1 %>%
  select(lopnr, angio_date, years, state, previous_FYND) %>%
  left_join(
    scaar_segments %>%
      transmute(
        lopnr = as.character(LopNr),
        angio_date = as.Date(INTERDAT),
        pseudo_MCEID,
        pseudo_SID,
        SEGMENT = as.integer(as.character(SEGMENT)),
        FINDING = as.integer(as.character(FINDING))
      ),
    by = c("lopnr", "angio_date"),
    na_matches = "never"
  ) %>%
  arrange(lopnr, SEGMENT)

View(pci_state1_segments)

# ==============================================================================
# 3. SUMMARISE FINDINGS PER PATIENT AND PCI DATE ----
# ==============================================================================

pci_state1_summary <- pci_state1_segments %>%
  group_by(lopnr, angio_date) %>%
  summarise(
    n_segments = n_distinct(SEGMENT, na.rm = TRUE),
    n_with_finding = n_distinct(
      SEGMENT[!is.na(FINDING)], na.rm = TRUE
    ),
    n_ge50 = n_distinct(
      SEGMENT[FINDING %in% 3:6], na.rm = TRUE
    ),
    n_ge70 = n_distinct(
      SEGMENT[FINDING %in% 4:6], na.rm = TRUE
    ),
    .groups = "drop"
  )

View(pci_state1_summary)

# Note: zero affected segments is not evidence of no stenosis if segment
# records or conclusive FINDING codes are missing. Same-date matches may
# include multiple procedures and do not establish pre/post-PCI timing.
