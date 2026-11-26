# ============================================================================== 
# PCI MISCLASSIFICATION: CONSTANT VS PERSISTENT + EXPONENTIAL DECAY
# msm 1.8.2; requires fit_adj and msm_dat_final4 in the R session.
# For downward misclassification, conditional on the latent state:
#   Constant:           log(OR) = beta * Z
#   Persistent + decay: log(OR) = beta * Z + gamma * Z * 2^(-u / h)
# Z = prior PCI indicator; u = years since most recent PCI; h = fixed half-life.
# Both coefficients are unconstrained in sign.
# ============================================================================== 

library(msm)
library(dplyr)
library(ggplot2)

# 1. CONSTRUCT EXPONENTIAL PCI COVARIATE ----

half_life <- 1  # years; fixed before fitting
stopifnot(is.numeric(half_life), length(half_life) == 1L,
          is.finite(half_life), half_life > 0)

msm_dat_compare <- msm_dat_final4 %>%
  ungroup() %>%
  mutate(
    pci_decay = if_else(
      prior_pci == 0, 0,
      2^(-expm1(pci_logtime) / half_life)
    )
  )

# Check covariates rather than silently dropping observations.
stopifnot(
  all(msm_dat_compare$prior_pci %in% c(0, 1)),
  all(is.finite(msm_dat_compare$pci_logtime)),
  all(msm_dat_compare$pci_logtime[msm_dat_compare$prior_pci == 1] >= 0),
  all(is.finite(msm_dat_compare$pci_decay)),
  all(msm_dat_compare$pci_decay >= 0 & msm_dat_compare$pci_decay <= 1)
)

# 2. FIT CONSTANT PCI ASSOCIATION ----
# Preserve Q, E structure, covariate constraints and observation scheme in fit_adj.
# The fixedpars indices assume the same parameter ordering as the original fit.

fit_constant <- update(
  fit_adj,
  data = msm_dat_compare,
  misccovariates = ~ prior_pci,
  miscconstraint = list(
    prior_pci = c(1, 2, 1, 2, 2)
  ),
  misccovinits = list(
    prior_pci = rep(0, 5)
  ),
  fixedpars = 18,
  method = "BFGS",
  control = list(trace = 1, REPORT = 1, maxit = 1000),
  hessian = TRUE
)

# 3. FIT PERSISTENT + EXPONENTIALLY DECAYING ASSOCIATION ----

fit_decay_persistent <- update(
  fit_adj,
  data = msm_dat_compare,
  misccovariates = ~ prior_pci + pci_decay,
  miscconstraint = list(
    prior_pci = c(1, 2, 1, 2, 2),
    pci_decay = c(1, 2, 1, 2, 2)
  ),
  misccovinits = list(
    prior_pci = rep(0, 5),
    pci_decay = rep(0, 5)
  ),
  fixedpars = c(18, 19),
  method = "BFGS",
  control = list(trace = 1, REPORT = 1, maxit = 1000),
  hessian = TRUE
)

# 4. VERIFY PARAMETER CONSTRAINTS ----
# Forward log-odds coefficients fixed at zero; all downward coefficients shared.
# This does not imply that absolute forward misclassification probabilities
# remain unchanged (the multinomial denominator may change).

check_constraints <- function(fit, variable, tol = 1e-10) {
  E <- fit$Ematrices[[variable]]
  stopifnot(is.matrix(E), all(dim(E) == c(4, 4)))
  print(E)
  stopifnot(
    abs(E[1, 2]) < tol,
    abs(E[2, 3]) < tol,
    abs(E[2, 1] - E[3, 1]) < tol,
    abs(E[2, 1] - E[3, 2]) < tol
  )
}

check_constraints(fit_constant, "prior_pci")
check_constraints(fit_decay_persistent, "prior_pci")
check_constraints(fit_decay_persistent, "pci_decay")

# 5. VERIFY IDENTICAL ESTIMATION SAMPLES AND PARAMETER COUNTS ----

vars <- c("(subject)", "(time)", "(state)")

stopifnot(
  isTRUE(all.equal(
    model.frame(fit_constant)[, vars, drop = FALSE],
    model.frame(fit_decay_persistent)[, vars, drop = FALSE],
    check.attributes = FALSE
  )),
  fit_decay_persistent$paramdata$nopt == fit_constant$paramdata$nopt + 1L
)

# 6. COMPARE MODEL FIT ----

models <- list(
  Constant = fit_constant,
  "Persistent + exponential decay" = fit_decay_persistent
)

model_comparison <- bind_rows(lapply(names(models), function(nm) {
  fit <- models[[nm]]
  data.frame(
    model = nm,
    n_parameters = fit$paramdata$nopt,
    loglik = -fit$minus2loglik / 2,
    minus2loglik = fit$minus2loglik,
    AIC = AIC(fit),
    convergence = fit$opt$convergence,
    foundse = isTRUE(fit$foundse)
  )
})) %>%
  mutate(delta_AIC = AIC - min(AIC))

print(model_comparison, width = Inf)

# Nested models, with gamma = 0 under H0. The half-life is fixed and
# gamma can be positive or negative, so there is no sign boundary under H0.
# The usual chi-square approximation still requires regularity conditions.

if (any(model_comparison$convergence != 0) ||
    any(!model_comparison$foundse)) {
  warning("Convergence or Hessian problem: interpret comparison with caution.")
}

LR <- fit_constant$minus2loglik - fit_decay_persistent$minus2loglik

if (LR < -1e-6) {
  warning("The larger model has a worse likelihood. Check optimization before using the LRT.")
} else {
  lrt_result <- data.frame(
    LR = max(0, LR),
    df = 1L,
    p_value = pchisq(max(0, LR), df = 1L, lower.tail = FALSE)
  )
  print(lrt_result)
}

# Individual covariate effects (not the combined OR at time u).
print(odds.msm(fit_constant))
print(odds.msm(fit_decay_persistent))

# 7. EXTRACT AND PLOT ESTIMATED PCI ASSOCIATIONS ----
# e21 represents the common coefficient for e21, e31 and e32.

beta_constant <- fit_constant$Ematrices$prior_pci[2, 1]
beta_persistent <- fit_decay_persistent$Ematrices$prior_pci[2, 1]
gamma <- fit_decay_persistent$Ematrices$pci_decay[2, 1]

print(data.frame(
  beta_constant = beta_constant,
  beta_persistent = beta_persistent,
  gamma = gamma,
  OR_constant = exp(beta_constant),
  OR_immediately_after_PCI = exp(beta_persistent + gamma),
  OR_long_run = exp(beta_persistent)
))

u <- seq(0, 7, length.out = 500)

pci_comparison_df <- bind_rows(
  data.frame(
    years_since_pci = u,
    log_OR = beta_constant,
    model = "Constant"
  ),
  data.frame(
    years_since_pci = u,
    log_OR = beta_persistent + gamma * 2^(-u / half_life),
    model = "Persistent + exponential decay"
  )
) %>%
  mutate(OR = exp(log_OR))

p_pci_comparison <- ggplot(
  pci_comparison_df,
  aes(x = years_since_pci, y = OR, colour = model, linetype = model)
) +
  geom_hline(yintercept = 1, colour = "grey45", linetype = "dashed") +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 0:7) +
  scale_y_log10() +
  scale_colour_manual(values = c(
    "Constant" = "#009fc3",
    "Persistent + exponential decay" = "#651220"
  )) +
  labs(
    x = "Years since most recent PCI",
    y = "Odds ratio",
    colour = NULL,
    linetype = NULL,
    title = "PCI-associated downward misclassification",
    subtitle = paste0("Exponential half-life: ", half_life, " year(s)"),
    caption = "Point estimates; relative to no previous PCI."
  ) +
  theme_bw(base_size = 12, base_family = "Aptos") +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey85"),
    plot.title = element_text(face = "bold")
  )

print(p_pci_comparison)

# Positive gamma implies a decreasing OR over time; negative gamma implies
# an increasing OR over time. The long-run OR is exp(beta_persistent).
# The plot displays point estimates, not confidence intervals.
# Ultimately compare population-standardized latent state occupancy curves.
