# ======================================================================
# FYND | Sahlgrenska-only informal goodness-of-fit diagnostic
# Uses the existing pooled fit_constant (no refitting).
# Local patches to msm 1.8.2, not to the installed package.
# Required objects: fit_constant, msm_dat_final4, cens_time
# ======================================================================

library(msm)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

stopifnot(packageVersion("msm") == "1.8.2")
stopifnot(inherits(fit_constant, "msm"))

# 1. Settings ----------------------------------------------------------
cut_time <- 45 / 365.25
max_years <- 7
times <- sort(unique(c(seq(0, max_years, by = 0.05), cut_time)))
q_names <- c("sex", "age", "tx_era", "post_cut")
e_names <- "prior_pci"
col_obs <- "#9EB8B7"
col_pred <- "#651220"
col_rec <- "#405C5B"
state_labels <- c("FYND = 1", "FYND = 2", "FYND >= 3", "Death/Re-HTx")

# 2. Identify Sahlgrenska subjects in the fitted model's order -----------
mf <- fit_constant$data$mf
stopifnot(all(c(q_names, e_names) %in% attr(mf, "covnames")),
          all(c(q_names, e_names) %in% names(mf)),
          fit_constant$qmodel$nstates == 4,
          fit_constant$emodel$misc,
          is.null(fit_constant$pci),
          all(c("Early", "After") %in% unique(as.character(mf$post_cut))))

mf_dat <- data.frame(lopnr = mf[["(subject)"]],
                     years = mf[["(time)"]], state = mf[["(state)"]])
ids_all <- unique(mf_dat$lopnr)

centres <- msm_dat_final4 %>%
  filter(lopnr %in% ids_all) %>%
  distinct(lopnr, unit)
stopifnot(!anyDuplicated(centres$lopnr), !anyNA(centres$unit),
          all(ids_all %in% centres$lopnr))
ids_sahl <- ids_all[centres$unit[match(ids_all, centres$lopnr)] == "Sahlgrenska"]
stopifnot(length(ids_sahl) > 0L)

stopifnot(all(c("lopnr", "obs_yrs") %in% names(cens_time)),
          !anyDuplicated(cens_time$lopnr))
censtime_vec <- cens_time$obs_yrs[match(ids_sahl, cens_time$lopnr)]
stopifnot(is.numeric(censtime_vec), !anyNA(censtime_vec), all(censtime_vec >= 0))

patient_dat <- mf_dat %>%
  filter(lopnr %in% ids_sahl) %>%
  group_by(lopnr) %>%
  summarise(first_time = min(years), last_time = max(years),
            death_time = if (any(state == 4)) min(years[state == 4]) else Inf,
            .groups = "drop") %>%
  slice(match(ids_sahl, lopnr)) %>%
  mutate(obs_yrs = censtime_vec)
stopifnot(identical(as.character(patient_dat$lopnr), as.character(ids_sahl)),
          all(abs(patient_dat$first_time) < 1e-10))
if (any(patient_dat$obs_yrs + 1e-8 < patient_dat$last_time)) {
  stop("Some obs_yrs values precede actual fitted observations. Check cens_time.")
}
message("Sahlgrenska patients in diagnostic: ", length(ids_sahl),
        " (of ", length(ids_all), " in fitted model)")

# Check that the selected cutoff agrees with actual recorded coding.
# At exactly the cutoff allow either coding (floating-point/calendar ties).
away_from_cut <- abs(mf[["(time)"]] - cut_time) > 1e-8
expected_cut <- ifelse(mf[["(time)"]] < cut_time, "Early", "After")
if (any(as.character(mf$post_cut)[away_from_cut] != expected_cut[away_from_cut])) {
  stop("Recorded post_cut does not agree with cut_time. Check the original coding before running.")
}

# 3. Local copies of original prevalence functions ----------------------
ns <- asNamespace("msm")
patch_env <- new.env(parent = ns)
expected_fixed <- get("expected.msm", envir = ns)
observed_fixed <- get("observed.msm", envir = ns)
covhist_fixed <- get("get.covhist", envir = ns)
prevalence_fixed <- msm::prevalence.msm

# Patch complete R statements, not individual deparsed lines.
# Only traverse statement positions, not arbitrary function arguments.
# Every change must match exactly once; otherwise stop with an error.
patch_statement <- function(fun, old, new) {
  n <- 0L
  walk <- function(node) {
    if (identical(node, old, ignore.srcref = TRUE)) {
      n <<- n + 1L
      return(new)
    }
    if (!is.call(node)) return(node)
    h <- node[[1L]]
    inds <- if (identical(h, quote(`{`))) {
      seq_along(node)[-1L]
    } else if (identical(h, quote(`if`))) {
      seq.int(3L, length(node))
    } else if (identical(h, quote(`for`))) {
      4L
    } else if (identical(h, quote(`while`))) {
      3L
    } else integer(0)
    for (k in inds) node[[k]] <- walk(node[[k]])
    node
  }
  body(fun) <- walk(body(fun))
  if (n != 1L) stop("Patch matched ", n, " statements: ",
                    paste(deparse(old), collapse = " "))
  fun
}

# 4. Correct grouping of distinct covariate histories -------------------
# Keep a list entry per subject, including subjects with no changes.
# Include the final observation time in the key: otherwise subjects with
# different terminal times may be grouped although their future age updates
# would differ. The original grouping/representative-row logic is retained.
covhist_fixed <- patch_statement(
  covhist_fixed,
  old = quote(change.times <- split(change.times[keep & (!(first | last))],
                                    subj[keep & (!(first | last))])),
  new = quote(change.times <- split(mf[["(time)"]][keep & !first],
                                    factor(subj[keep & !first], levels = unique(subj)),
                                    drop = FALSE))
)

# 5. Complete deterministic covariate histories -------------------------
# Preserve actual covariate updates; add the day-45 switch and annual age
# updates beyond final observation, even for patients who died early.
# prior_pci is held at its last observed value after follow-up ends.
complete_covhist <- function(dat, horizon, cutoff) {
  dat <- dat[order(dat$time), , drop = FALSE]
  # At duplicate times use the last recorded covariate profile.
  dat <- dat[!duplicated(dat$time, fromLast = TRUE), , drop = FALSE]
  original_times <- dat$time
  stopifnot(length(original_times) >= 1L, abs(original_times[1]) < 1e-10,
            !anyNA(dat[, c(q_names, e_names), drop = FALSE]))
  last_time <- max(original_times)
  age_at_tx <- dat$age[1] - original_times[1]
  annual <- seq(0, ceiling(horizon), by = 1)
  future_age <- annual[annual > last_time & annual <= horizon]
  new_times <- sort(unique(c(original_times, future_age,
                             if (cutoff <= horizon) cutoff)))
  ind <- findInterval(new_times, original_times)
  out <- dat[ind, , drop = FALSE]
  out$time <- new_times
  inserted <- !new_times %in% original_times
  out$age[inserted] <- age_at_tx + new_times[inserted]
  early_late <- ifelse(new_times < cutoff, "Early", "After")
  out$post_cut <- if (is.factor(dat$post_cut)) {
    factor(early_late, levels = levels(dat$post_cut))
  } else early_late
  stopifnot(!anyNA(out$post_cut), !anyDuplicated(out$time))
  out
}

# 6. Correct expected.msm() --------------------------------------------
# Retain the package's population weighting and piecewise transition
# probability computation. Replace only covariate-history preparation,
# the arguments of the matrix computation, and the emission profile.
expected_fixed <- patch_statement(
  expected_fixed,
  old = quote(ctimes <- covmat$time[covmat$subject == subji][-c(1, ni)]),
  new = quote({
    hist_i <- complete_covhist(covmat[covmat$subject == subji, , drop = FALSE],
                               horizon = max(times), cutoff = cut_time)
    ctimes <- hist_i$time[-1L]
  })
)
expected_fixed <- patch_statement(
  expected_fixed,
  old = quote(covs <- covmat[covmat$subject == subji,
                            attr(x$data$mf, "covnames"), drop = FALSE][-ni, , drop = FALSE]),
  new = quote(covs <- hist_i[, attr(x$data$mf, "covnames"), drop = FALSE])
)
expected_fixed <- patch_statement(
  expected_fixed,
  old = quote(pmat <- pmatrix.piecewise.msm(x, t1 = timezero, t2 = times[j],
                                           times = ctimes, covariates = ccovs)),
  new = quote(pmat <- pmatrix.piecewise.msm(x, t1 = timezero, t2 = times[j],
                                           times = ctimes,
                                           covariates = lapply(ccovs, function(z) z[q_names])))
)
expected_fixed <- patch_statement(
  expected_fixed,
  old = quote(ecovs <- if (length(ctimes) == 0) ccovs else ccovs[[length(ccovs)]]),
  new = quote(ecovs <- ccovs[[findInterval(times[j], ctimes) + 1L]][e_names])
)
# The original convention of applying E(t) even at t = 0 is retained.

# 7. Correct absorbing-state carry-forward in observed.msm() ------------
# Prevent death/re-HTx from being carried beyond potential observation time.
# The original function otherwise retains the absorbing state at every
# subsequent grid point after checking censtime at only the first one.
observed_fixed <- patch_statement(
  observed_fixed,
  old = quote(states.expand[i, j:(length(times))] <- state.i[length(time.i)]),
  new = quote(states.expand[i, j:(length(times))] <-
                ifelse(times[j:length(times)] <= censtime[i],
                       state.i[length(time.i)], NA))
)

# 8. Link local functions and check patched history construction --------
patch_env$complete_covhist <- complete_covhist
patch_env$cut_time <- cut_time
patch_env$q_names <- q_names
patch_env$e_names <- e_names

environment(covhist_fixed) <- patch_env
environment(observed_fixed) <- patch_env
environment(expected_fixed) <- patch_env
environment(prevalence_fixed) <- patch_env
patch_env$get.covhist <- covhist_fixed
patch_env$observed.msm <- observed_fixed
patch_env$expected.msm <- expected_fixed

ch <- covhist_fixed(fit_constant, subset = ids_sahl)
stopifnot(length(ch$hist) == length(ids_sahl),
          length(unique(ch$example$subject)) == length(unique(ch$hist)))
# Verify the extrapolation on a truncated history (baseline only).
h0 <- ch$example[ch$example$subject == ch$example$subject[1], , drop = FALSE]
h0 <- h0[h0$time == 0, , drop = FALSE]
h0 <- h0[1, , drop = FALSE]
h_check <- complete_covhist(h0, horizon = max_years, cutoff = cut_time)
stopifnot(all(as.character(h_check$post_cut[h_check$time >= cut_time]) == "After"),
          all(as.character(h_check$post_cut[h_check$time < cut_time]) == "Early"),
          all(1:max_years %in% h_check$time))

# 9. Original prevalence.msm() with local patches -----------------------
prev_constant <- prevalence_fixed(fit_constant,
                                  times = times, timezero = 0,
                                  covariates = "population", interp = "start",
                                  subset = ids_sahl, censtime = censtime_vec,
                                  ci = "none", plot = FALSE)

obs_n <- prev_constant$Observed[, "Total"]
exp_n <- prev_constant$Expected[, "Total"]
stopifnot(isTRUE(all.equal(as.numeric(obs_n), as.numeric(exp_n), tolerance = 1e-8)))
nonempty <- obs_n > 0
stopifnot(all(abs(rowSums(prev_constant$`Expected percentages`)[nonempty] - 100) < 1e-6),
          all(abs(rowSums(prev_constant$`Observed percentages`)[nonempty] - 100) < 1e-6))
# Independent check of denominators against fitted observation times and
# potential administrative follow-up (including patients who died).
expected_n <- vapply(times, function(t) {
  sum(t <= patient_dat$last_time |
        (is.finite(patient_dat$death_time) & t <= patient_dat$obs_yrs))
}, integer(1))
stopifnot(all(as.numeric(obs_n) == expected_n))
print(head(prev_constant$`Observed percentages`))
print(head(prev_constant$`Expected percentages`))

# 10. Identify genuine conclusive angiographies -------------------------
# Use only retained model observations; look up the artificial flag in
# the original data when available. Baseline and state 99 never count.
rec_source <- msm_dat_final4 %>% filter(lopnr %in% ids_sahl)
if ("artificial" %in% names(rec_source)) {
  assessment_dat <- rec_source %>%
    filter(years > 0, state %in% 1:3, !coalesce(as.logical(artificial), FALSE))
} else {
  warning("No artificial flag in msm_dat_final4: assuming post-baseline states 1-3 are genuine assessments.")
  assessment_dat <- rec_source %>% filter(years > 0, state %in% 1:3)
}
assessment_dat <- assessment_dat %>%
  distinct(lopnr, years, state) %>%
  semi_join(mf_dat, by = c("lopnr", "years", "state")) %>%
  distinct(lopnr, years)
assessment_list <- split(assessment_dat$years, as.character(assessment_dat$lopnr))

# 11. Descriptive recency and absence, independent of Viterbi ------------
# At each t: alive, re-HTx-free, and under observation.
# Recency percentiles: only those with at least one prior assessment.
# Absence percentage: denominator includes everyone eligible/alive.
recency_patient_times <- bind_rows(lapply(seq_len(nrow(patient_dat)), function(i) {
  pat <- patient_dat[i, ]
  eval_times <- times[times <= min(pat$obs_yrs, pat$last_time) & times < pat$death_time]
  if (!length(eval_times)) return(NULL)
  a_times <- sort(unique(as.numeric(assessment_list[[as.character(pat$lopnr)]])))
  idx <- findInterval(eval_times, a_times)
  assessed <- idx > 0L
  recency <- rep(NA_real_, length(eval_times))
  recency[assessed] <- eval_times[assessed] - a_times[idx[assessed]]
  data.frame(lopnr = pat$lopnr, years = eval_times, assessed = assessed, recency = recency)
}))

safe_quantile <- function(z, p) {
  if (all(is.na(z))) return(NA_real_)
  as.numeric(quantile(z, probs = p, na.rm = TRUE))
}
recency_summary <- recency_patient_times %>%
  group_by(years) %>%
  summarise(n_alive = n(), n_assessed = sum(assessed), n_unassessed = sum(!assessed),
            pct_unassessed = 100 * n_unassessed / n_alive,
            median = safe_quantile(recency, 0.50),
            q10 = safe_quantile(recency, 0.10), q25 = safe_quantile(recency, 0.25),
            q75 = safe_quantile(recency, 0.75), q90 = safe_quantile(recency, 0.90),
            .groups = "drop")
# Validate the descriptive denominator against observed nonabsorbing states.
rec_check <- match(recency_summary$years, times)
stopifnot(all(recency_summary$n_alive == rowSums(prev_constant$Observed[, 1:3, drop = FALSE])[rec_check]))

# 12. Prepare four-panel prevalence data -------------------------------
obs <- as.data.frame(prev_constant$`Observed percentages`)
pred <- as.data.frame(prev_constant$`Expected percentages`)
stopifnot(ncol(obs) == 4, ncol(pred) == 4)
names(obs) <- names(pred) <- paste("State", 1:4)
prev_df <- bind_rows(mutate(obs, years = times, type = "Observed"),
                     mutate(pred, years = times, type = "Model-predicted")) %>%
  pivot_longer(cols = starts_with("State"), names_to = "state", values_to = "prevalence") %>%
  mutate(state = factor(state, levels = paste("State", 1:4), labels = state_labels),
         type = factor(type, levels = c("Observed", "Model-predicted")))

# 13. Common plotting theme --------------------------------------------
theme_common <- theme_minimal(base_size = 11, base_family = "Aptos") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#E8E8E8", linewidth = 0.35),
        panel.border = element_rect(colour = "#BDBDBD", fill = NA, linewidth = 0.5),
        strip.text = element_text(face = "bold", size = 11),
        axis.text = element_text(colour = "#333333"),
        plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = grid::unit(1.3, "cm"), plot.margin = margin(6, 8, 6, 8))
common_x <- scale_x_continuous(breaks = 0:max_years, limits = c(0, max_years),
                               expand = expansion(mult = c(0.01, 0.01)))

p_prev <- ggplot(prev_df, aes(years, prevalence, colour = type, linetype = type)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  facet_wrap(~ state, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("Observed" = col_obs, "Model-predicted" = col_pred)) +
  scale_linetype_manual(values = c("Observed" = "solid", "Model-predicted" = "dashed")) +
  common_x + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Prevalence (%)") + theme_common

# 14. Recency: median, IQR and 10th-90th percentiles --------------------
p_recency <- ggplot(recency_summary, aes(x = years)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = col_obs, alpha = 0.20, na.rm = TRUE) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = col_obs, alpha = 0.55, na.rm = TRUE) +
  geom_line(aes(y = median), colour = col_rec, linewidth = 0.9, na.rm = TRUE) +
  common_x + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Recency of CAV ascertainment", subtitle = "Median, IQR and 10th-90th percentiles",
       x = NULL, y = "Time since last\nassessment (years)") + theme_common

# 15. Absence of ascertainment -----------------------------------------
p_unassessed <- ggplot(recency_summary, aes(years, pct_unassessed)) +
  geom_area(fill = col_obs, alpha = 0.25, na.rm = TRUE) +
  geom_line(colour = col_rec, linewidth = 0.9, na.rm = TRUE) + common_x +
  scale_y_continuous(breaks = seq(0, 100, by = 25), limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Absence of CAV ascertainment", x = "Years since HTx",
       y = "Patients without prior\nassessment (%)") + theme_common

# 16. Combine and export -----------------------------------------------
# wrap_plots avoids the previous overloaded '/' operator problem.
p_gof <- patchwork::wrap_plots(p_prev, p_recency, p_unassessed,
                               ncol = 1, heights = c(3.2, 1.3, 1))
print(p_gof)
ggsave("FYND_prevalence_GoF_Sahlgrenska.png", plot = p_gof,
       width = 10, height = 10, units = "in", dpi = 400, bg = "white")
