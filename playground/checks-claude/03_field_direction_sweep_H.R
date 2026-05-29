# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXPLORATORY: Field-direction sweep (Camino-H)
# 03_field_direction_sweep_H.R
#
# The main analysis (script 07) follows Camino-T: fix IUL (~external field H),
# sweep MSP (~temperature T) across MSP_c, extract susceptibility exponent gamma.
#
# This exploratory script does the COMPLEMENTARY Camino-H: fix MSP, sweep IUL
# (the field), and probe the response. Two physics targets:
#   (a) Critical isotherm exponent delta:  Phi ~ |IUL - IUL_c|^(1/delta)
#       at MSP = MSP_c.  Mean-field predicts delta = 3.
#   (b) Field-direction susceptibility:    chi_H = Var(Phi) along the IUL axis.
#
# Widom scaling relation (mean-field check):  gamma = beta * (delta - 1).
#   With beta = 1/2, delta = 3  =>  gamma = 1. Consistent with the Camino-T result.
#
# NOTE: In the social model "negative field" (IUL < 0) is meaningless, so the
# classic first-order field reversal cannot be traced. We instead sit on the
# critical isotherm and sweep IUL > 0. This is exploratory.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

PLOTS_DIR <- "playground/checks-claude/plots/"
DATA_OUT  <- "playground/checks-claude/data/"
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT,  showWarnings = FALSE, recursive = TRUE)

MEAN_KEY <- "mean_0.40"   # matches main-text Figure 1 (tau_mean = 0.4)
FILES <- list(
  GSS = "output/04_GSS_diffusion_sims/results_m03-06_sd0.12_random.rds",
  DP  = "output/04_GSS_ER_diffusion_sims/results_m03-06_sd0.12_random.rds"
)

# ---------------------------------------------------------------------------
# Helper: load and aggregate to (IUL, MSP) cells
# ---------------------------------------------------------------------------
load_agg <- function(path) {
  if (!file.exists(path)) { cat("MISSING:", path, "\n"); return(NULL) }
  raw <- readRDS(path)
  if (is.list(raw) && !is.data.frame(raw)) {
    raw <- if (!is.null(raw[[MEAN_KEY]])) raw[[MEAN_KEY]] else raw[[1]]
  }
  raw %>%
    mutate(phi = num_adopters / N_nodes_actual) %>%
    group_by(IUL = innovation_iul_Gamma, MSP = social_distance_h) %>%
    summarise(avg_phi = mean(phi, na.rm = TRUE),
              var_phi = var(phi, na.rm = TRUE),
              .groups = "drop")
}

# ---------------------------------------------------------------------------
# Helper: extract delta on the critical isotherm (fixed MSP, sweep IUL)
#   Phi - Phi_c ~ (IUL - IUL_c)^(1/delta)   for IUL > IUL_c
# IUL_c is taken as the INFLECTION point (max slope d phi / d IUL), which is the
# field-direction analogue of the critical point on this isotherm.
# ---------------------------------------------------------------------------
fit_isotherm <- function(slice) {
  slice <- slice %>% arrange(IUL)
  d_phi <- c(NA, diff(slice$avg_phi))
  infl <- which.max(d_phi)
  iul_c <- slice$IUL[infl]
  phi_c <- slice$avg_phi[infl]
  up <- slice %>% filter(IUL > iul_c, avg_phi > phi_c) %>%
    mutate(dx = IUL - iul_c, dy = pmax(avg_phi - phi_c, 1e-6))
  delta <- NA; r2 <- NA
  if (nrow(up) >= 3) {
    fit <- lm(log10(dy) ~ log10(dx), data = up)
    slope <- coef(fit)[2]              # slope = 1/delta
    delta <- if (is.finite(slope) && slope > 1e-6) 1 / slope else NA
    r2 <- summary(fit)$r.squared
  }
  list(iul_c = iul_c, phi_c = phi_c, delta = delta, r2 = r2, n = nrow(up))
}

# ---------------------------------------------------------------------------
# Helper: field-direction susceptibility exponent
#   chi_H = Var(Phi) ~ |IUL - IUL_c|^(-gamma_H), super-critical side (IUL > IUL_c)
# ---------------------------------------------------------------------------
fit_chi_H <- function(slice) {
  slice <- slice %>% arrange(IUL)
  iul_c <- slice$IUL[which.max(slice$var_phi)]
  up <- slice %>% filter(IUL > iul_c) %>%
    mutate(dx = IUL - iul_c, chi = pmax(var_phi, 1e-7)) %>%
    filter(dx > 0)
  gamma_H <- NA; r2 <- NA
  if (nrow(up) >= 3) {
    fit <- lm(log10(chi) ~ log10(dx), data = up)
    gamma_H <- -coef(fit)[2]
    r2 <- summary(fit)$r.squared
  }
  list(iul_c = iul_c, gamma_H = gamma_H, r2 = r2, n = nrow(up))
}

# ---------------------------------------------------------------------------
# Critical MSP from the Camino-T analysis (script 07 / abstract Table 2).
# These are the temperatures at which the T-path located the continuous
# transition; the critical isotherm of the H-path must sit at the same MSP_c.
# ---------------------------------------------------------------------------
MSP_C_TPATH <- c(GSS = 0.25, DP = 0.4167)

# ---------------------------------------------------------------------------
# Main loop over topologies
# ---------------------------------------------------------------------------
results <- list()
iso_family_plots <- list()
iso_crit_plots <- list()

for (topo in names(FILES)) {
  agg <- load_agg(FILES[[topo]])
  if (is.null(agg)) next

  msp_c <- MSP_C_TPATH[[topo]]
  msp_sorted <- sort(unique(agg$MSP))
  msp_c <- msp_sorted[which.min(abs(msp_sorted - msp_c))]  # snap to grid

  # ---- (1) Family-of-isotherms plot: Phi vs IUL, one line per MSP ----
  fam <- agg %>% filter(IUL <= 0.6)  # focus on the active region
  p_fam <- ggplot(fam, aes(IUL, avg_phi, color = factor(round(MSP, 2)), group = MSP)) +
    geom_line(linewidth = 0.5) +
    scale_color_viridis_d(option = "plasma", name = "MSP (T)") +
    labs(title = paste0(topo, ": isotherm family (Camino-H)"),
         x = "IUL  (external field H)", y = expression(Phi)) +
    theme_minimal(base_size = 10)
  iso_family_plots[[topo]] <- p_fam

  # ---- (2) Exponents at sub / critical / super isotherms ----
  msp_sub   <- msp_sorted[max(1, which(msp_sorted == msp_c) - 2)]
  msp_super <- msp_sorted[min(length(msp_sorted), which(msp_sorted == msp_c) + 2)]
  slices <- c(sub = msp_sub, critical = msp_c, super = msp_super)

  for (lab in names(slices)) {
    mspv <- slices[[lab]]
    slice <- agg %>% filter(abs(MSP - mspv) < 1e-6)
    iso <- fit_isotherm(slice)
    chi <- fit_chi_H(slice)
    results[[paste(topo, lab)]] <- data.frame(
      topology = topo, slice = lab, MSP = round(mspv, 3),
      IUL_c = iso$iul_c, delta = iso$delta, delta_r2 = iso$r2,
      gamma_H = chi$gamma_H, gamma_H_r2 = chi$r2
    )

    if (lab == "critical") {
      p <- ggplot(slice, aes(IUL, avg_phi)) +
        geom_line(color = "grey50") + geom_point(size = 1.3) +
        geom_vline(xintercept = iso$iul_c, linetype = "dashed", color = "red") +
        labs(title = paste0(topo, ": critical isotherm (MSP=", round(mspv,3), ")"),
             subtitle = paste0("IUL_c=", round(iso$iul_c,3),
                               "  delta=", ifelse(is.na(iso$delta),"NA",round(iso$delta,2)),
                               "  (R2=", ifelse(is.na(iso$r2),"NA",round(iso$r2,2)), ")"),
             x = "IUL  (external field H)", y = expression(Phi)) +
        theme_minimal(base_size = 10)
      iso_crit_plots[[topo]] <- p
    }
  }
}

summary_df <- bind_rows(results)
cat("\n================ CAMINO-H EXPLORATORY RESULTS ================\n")
print(summary_df, row.names = FALSE)
cat("\nMean-field targets: delta = 3.0 ; Widom: gamma = beta(delta-1) = 0.5*2 = 1.0\n")
cat("==============================================================\n")
write.csv(summary_df, file.path(DATA_OUT, "camino_H_exponents.csv"), row.names = FALSE)

if (length(iso_family_plots) > 0) {
  combined_fam <- wrap_plots(iso_family_plots, ncol = length(iso_family_plots))
  ggsave(file.path(PLOTS_DIR, "camino_H_isotherm_family.pdf"), combined_fam,
         width = 6 * length(iso_family_plots), height = 4)
  cat("Saved isotherm-family plot.\n")
}
if (length(iso_crit_plots) > 0) {
  combined_crit <- wrap_plots(iso_crit_plots, ncol = length(iso_crit_plots))
  ggsave(file.path(PLOTS_DIR, "camino_H_critical_isotherm.pdf"), combined_crit,
         width = 5 * length(iso_crit_plots), height = 4)
  cat("Saved critical-isotherm plot.\n")
}
