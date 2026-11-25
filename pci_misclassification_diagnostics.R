library(ggplot2)
library(dplyr)

# Requires the fitted msm model, fit_adj, in the current R session.

# ==============================================================================
# 1. EXTRACT PCI COEFFICIENTS AND COVARIANCE MATRIX ----
# ==============================================================================

fit <- fit_adj
p <- fit$paramdata
pci_vars <- c("prior_pci", "pci_logtime")

# Downward misclassification e21; shared with e31 and e32.
beta <- vapply(
  pci_vars,
  function(v) fit$Ematrices[[v]][2, 1],
  numeric(1)
)

idx_free <- intersect(which(p$plabs == "hcov"), p$optpars)

if (!isTRUE(fit$foundse) || length(idx_free) < 2L) {
  stop("Check the Hessian and the free HMM covariate parameter indices.")
}

idx <- vapply(
  beta,
  function(b) idx_free[which.min(abs(fit$estimates[idx_free] - b))],
  integer(1)
)

stopifnot(
  length(unique(idx)) == 2L,
  isTRUE(all.equal(unname(fit$estimates[idx]), unname(beta)))
)

V <- fit$covmat[idx, idx, drop = FALSE]
dimnames(V) <- list(pci_vars, pci_vars)

print(beta)
print(V)
print(cov2cor(V))

# ==============================================================================
# 2. ESTIMATE OR CURVE AND POINTWISE 95% CIs ----
# ==============================================================================

pci_or_df <- data.frame(
  years_since_pci = seq(0, 7, length.out = 500)
)

X <- cbind(1, log1p(pci_or_df$years_since_pci))

pci_or_df <- pci_or_df %>%
  mutate(
    log_OR = as.numeric(X %*% beta),
    SE = sqrt(rowSums((X %*% V) * X)),
    OR = exp(log_OR),
    LCL = exp(log_OR - qnorm(0.975) * SE),
    UCL = exp(log_OR + qnorm(0.975) * SE)
  )

# ==============================================================================
# 3. PLOT ----
# ==============================================================================

p_pci <- ggplot(pci_or_df, aes(x = years_since_pci, y = OR)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL),
              fill = "#009fc3", alpha = 0.18) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_line(colour = "#009fc3", linewidth = 1) +
  scale_x_continuous(breaks = 0:7) +
  scale_y_log10(labels = scales::label_number()) +
  labs(
    x = "Years since most recent PCI",
    y = "Odds ratio",
    title = "Downward misclassification following PCI",
    subtitle = "Relative to no previous PCI",
    caption = "Shaded area represents the pointwise 95% confidence interval."
  ) +
  theme_bw(base_size = 12, base_family = "Aptos") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(colour = "grey40")
  )

print(p_pci)

# ==============================================================================
# 4. CONSTRUCT CONSECUTIVE OBSERVED STATE PAIRS ----
# ==============================================================================

# These are observed downward transitions, not known misclassifications.
# Exclude the assumed baseline, unknown states, and the absorbing state.
mf <- model.frame(fit_adj)

pci_pairs <- mf %>%
  transmute(
    lopnr = `(subject)`,
    years = `(time)`,
    state = as.integer(as.character(`(state)`)),
    prior_pci = prior_pci,
    pci_logtime = pci_logtime
  ) %>%
  filter(state %in% 1:3, years > 0) %>%
  group_by(lopnr) %>%
  arrange(years, .by_group = TRUE) %>%
  mutate(
    prev_state = lag(state),
    prev_years = lag(years),
    interval_years = years - prev_years,
    downward = state < prev_state,
    can_decrease = prev_state %in% 2:3,
    years_since_pci = if_else(
      prior_pci == 1, expm1(pci_logtime), NA_real_
    )
  ) %>%
  filter(!is.na(prev_state)) %>%
  ungroup()

# ==============================================================================
# 5. CATEGORISE TIME SINCE PCI ----
# ==============================================================================

pci_groups <- c(
  "No previous PCI",
  "0–<0.5 years",
  "0.5–<1 year",
  "1–<2 years",
  "2–<3 years",
  "3+ years"
)

pci_pairs <- pci_pairs %>%
  mutate(
    pci_group = case_when(
      prior_pci == 0        ~ "No previous PCI",
      years_since_pci < 0.5 ~ "0–<0.5 years",
      years_since_pci < 1   ~ "0.5–<1 year",
      years_since_pci < 2   ~ "1–<2 years",
      years_since_pci < 3   ~ "2–<3 years",
      years_since_pci >= 3  ~ "3+ years"
    ),
    pci_group = factor(pci_group, levels = pci_groups)
  )

# ==============================================================================
# 6. SUMMARISE DOWNWARD TRANSITIONS ----
# ==============================================================================

pci_summary <- pci_pairs %>%
  group_by(pci_group, .drop = FALSE) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_possible_decrease = sum(can_decrease),
    n_downward = sum(downward),
    n_patients_downward = n_distinct(lopnr[downward]),
    pct_downward_all = if (n_pairs > 0) 100 * n_downward / n_pairs else NA_real_,
    pct_downward_eligible = if (n_possible_decrease > 0)
      100 * n_downward / n_possible_decrease else NA_real_,
    .groups = "drop"
  )

print(pci_summary, width = Inf)

# ==============================================================================
# 7. STRATIFY BY PRECEDING OBSERVED STATE ----
# ==============================================================================

pci_summary_strat <- pci_pairs %>%
  filter(can_decrease) %>%
  group_by(pci_group, prev_state) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_downward = sum(downward),
    pct_downward = 100 * n_downward / n_pairs,
    median_interval = median(interval_years),
    .groups = "drop"
  ) %>%
  arrange(prev_state, pci_group)

print(pci_summary_strat, n = Inf, width = Inf)

# ==============================================================================
# 8. PCI BETWEEN CONSECUTIVE ASSESSMENTS ----
# ==============================================================================

# Assumes PCI at the preceding angiography occurred after its FYND assessment.
# Uses a small numerical tolerance for reconstructing PCI times from log1p(u).
pci_pairs <- pci_pairs %>%
  mutate(
    last_pci_year = years - years_since_pci,
    pci_between = case_when(
      prior_pci == 0 ~ FALSE,
      prior_pci == 1 ~ last_pci_year >= prev_years - 1e-9,
      TRUE ~ NA
    )
  )

pci_between_summary <- pci_pairs %>%
  filter(can_decrease) %>%
  group_by(pci_between, prev_state) %>%
  summarise(
    n_patients = n_distinct(lopnr),
    n_pairs = n(),
    n_downward = sum(downward),
    pct_downward = 100 * n_downward / n_pairs,
    median_interval = median(interval_years),
    .groups = "drop"
  )

print(pci_between_summary, n = Inf, width = Inf)
