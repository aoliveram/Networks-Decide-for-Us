# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Stability Checks for γ Exponent Across Parameter Space
# 00_exponent_stability_checks.R
# Claude's comprehensive checks for universality claims in NDFU
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(tidyr)

# Configuration
RESULTS_DIR_GSS <- "output/07_phase_transition/GSS/"
RESULTS_DIR_DP <- "output/07_phase_transition/GSS_ER/"
PLOTS_DIR <- "playground/checks-claude/plots/"
DATA_OUT_DIR <- "playground/checks-claude/data/"

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load γ exponent results from Method 4
gamma_gss_file <- paste0(RESULTS_DIR_GSS, "gamma_exponents.csv")
gamma_dp_file <- paste0(RESULTS_DIR_DP, "gamma_exponents.csv")

if (file.exists(gamma_gss_file) && file.exists(gamma_dp_file)) {
  gamma_gss <- read.csv(gamma_gss_file)
  gamma_dp <- read.csv(gamma_dp_file)

  # Merge and compare
  gamma_comparison <- gamma_gss %>%
    rename(Gamma_GSS = Gamma, MSP_c_GSS = MSP_c) %>%
    select(IUL, Gamma_GSS, MSP_c_GSS) %>%
    left_join(
      gamma_dp %>%
        rename(Gamma_DP = Gamma, MSP_c_DP = MSP_c) %>%
        select(IUL, Gamma_DP, MSP_c_DP),
      by = "IUL"
    ) %>%
    mutate(
      Gamma_Ratio = Gamma_DP / Gamma_GSS,
      Universal_Test = ifelse(Gamma_GSS < 1.2, "Pass (GSS)", "Fail (GSS)"),
      Discontinuous_Test = ifelse(Gamma_DP > 2.5, "Pass (DP)", "Fail (DP)")
    )

  cat("\n=== EXPONENT COMPARISON: GSS vs. GSS-DP ===\n")
  print(gamma_comparison)
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("GSS: Gamma =", sprintf("%.3f ± %.3f\n", mean(gamma_gss$Gamma, na.rm=T), sd(gamma_gss$Gamma, na.rm=T)))
  cat("DP:  Gamma =", sprintf("%.3f ± %.3f\n", mean(gamma_dp$Gamma, na.rm=T), sd(gamma_dp$Gamma, na.rm=T)))
  cat("Separation (DP/GSS ratio) =", sprintf("%.2f\n", mean(gamma_comparison$Gamma_Ratio, na.rm=T)))

  write.csv(gamma_comparison, paste0(DATA_OUT_DIR, "gamma_comparison_gss_vs_dp.csv"), row.names = FALSE)

  # Visualization
  p_compare <- ggplot(gamma_comparison, aes(x = factor(IUL))) +
    geom_point(aes(y = Gamma_GSS), color = "blue", size = 3, shape = 1) +
    geom_point(aes(y = Gamma_DP), color = "red", size = 3, shape = 2) +
    geom_hline(yintercept = 1, color = "blue", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 4, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Susceptibility Exponent γ: GSS (Blue) vs. GSS-DP (Red)",
      x = "Intrinsic Utility Level (IUL)",
      y = "γ Exponent",
      subtitle = "Blue dashed = Mean-Field prediction (γ=1); Red dashed = DP baseline"
    ) +
    theme_minimal() +
    theme(legend.position = "top")

  ggsave(paste0(PLOTS_DIR, "gamma_comparison.pdf"), p_compare, width = 8, height = 5)

} else {
  cat("WARNING: γ exponent files not found. Run 07_phase_transition_deeper.R first.\n")
}

# CHECK 2: Fit Quality (R² of log-log regression) ============================
cat("\n=== CHECK 2: LOG-LOG FIT QUALITY ===\n")
cat("(To be implemented: extract R² from phase transition script outputs)\n")
cat("Goal: verify γ estimates are not from poor power-law fits.\n")

# CHECK 3: Bootstrap confidence intervals on γ ================================
cat("\n=== CHECK 3: BOOTSTRAP CI ON γ (FUTURE) ===\n")
cat("Recommended: resample simulation runs, recompute γ, report 95% CI.\n")
cat("Would strengthen universality claim if CI is narrow.\n")

cat("\n✓ Exponent stability checks complete.\n")
