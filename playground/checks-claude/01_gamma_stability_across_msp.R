# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Check: Is γ stable as we vary MSP around the critical point?
# 01_gamma_stability_across_msp.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(patchwork)

# This script extracts γ for a FIXED IUL but sliding window of MSP values
# Goal: if true universality, γ should be ~constant in the critical regime

RESULTS_DIR_GSS <- "output/04_GSS_diffusion_sims/"
RESULTS_DIR_DP <- "output/04_GSS_ER_diffusion_sims/"  # Note: actual folder named with "ER" suffix
PLOTS_DIR <- "playground/checks-claude/plots/"
DATA_OUT_DIR <- "playground/checks-claude/data/"

TARGET_IUL <- 0.250  # Pick one IUL and study in detail
TAU_MEAN <- 0.4
TAU_SD <- 0.12
SEEDING <- "random"

# ============================================================================
# Function: Compute γ for a sliding window around MSP
# ============================================================================
compute_gamma_sliding <- function(df_agg, target_iul, window_size = 0.15) {
  df_filtered <- df_agg %>%
    filter(abs(innovation_iul_Gamma - target_iul) < 1e-4) %>%
    arrange(social_distance_h)

  if (nrow(df_filtered) < 5) return(data.frame())

  # Find critical point (max variance)
  msp_c_idx <- which.max(df_filtered$sd_adoption)
  msp_c <- df_filtered$social_distance_h[msp_c_idx]

  # Sliding windows centered on points around MSP_c
  msp_values <- unique(df_filtered$social_distance_h)
  results_list <- list()

  for (msp_center in msp_values) {
    window_data <- df_filtered %>%
      filter(
        social_distance_h >= msp_center - window_size / 2,
        social_distance_h <= msp_center + window_size / 2,
        social_distance_h > msp_c  # Only super-critical regime (above critical point)
      )

    if (nrow(window_data) >= 3) {
      fit_data <- window_data %>%
        mutate(
          abs_dist = abs(social_distance_h - msp_c),
          sd_adoption_safe = pmax(sd_adoption, 1e-5)
        ) %>%
        filter(abs_dist > 0)

      if (nrow(fit_data) >= 2) {
        lm_fit <- lm(log10(sd_adoption_safe) ~ log10(abs_dist), data = fit_data)
        gamma <- -coef(lm_fit)[2]
        r_squared <- summary(lm_fit)$r.squared

        results_list[[length(results_list) + 1]] <- data.frame(
          MSP_Window_Center = msp_center,
          Gamma = gamma,
          R_squared = r_squared,
          N_Points = nrow(fit_data),
          MSP_c = msp_c
        )
      }
    }
  }

  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    return(data.frame())
  }
}

# ============================================================================
# Load and process data
# ============================================================================

# GSS data
file_gss <- paste0(RESULTS_DIR_GSS, "results_m03-06_all_random.rds")  # Try all_random first
if (!file.exists(file_gss)) {
  file_gss <- paste0(RESULTS_DIR_GSS, "results_m03-06_sd0.12_random.rds")  # Fallback
}
if (file.exists(file_gss)) {
  cat("Loading GSS data from:", file_gss, "\n")
  raw_gss <- readRDS(file_gss)

  # Extract the correct dataframe
  if (inherits(raw_gss, "list")) {
    mean_key <- "mean_0.40"
    if (!is.null(raw_gss[[mean_key]])) {
      raw_gss <- raw_gss[[mean_key]]
    } else {
      # If mean_0.40 doesn't exist, try to grab the first available
      raw_gss <- raw_gss[[1]]
    }
  }

  # Ensure it's a dataframe
  if (!is.data.frame(raw_gss)) {
    cat("ERROR: raw_gss is not a dataframe. Structure:\n")
    str(raw_gss)
    raw_gss <- NULL
  }

  if (!is.null(raw_gss)) {
    agg_gss <- raw_gss %>%
      mutate(adopters_prop = num_adopters / N_nodes_actual) %>%
      group_by(innovation_iul_Gamma, social_distance_h) %>%
      summarise(
        sd_adoption = sd(adopters_prop, na.rm = TRUE),
        avg_adoption = mean(adopters_prop, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    agg_gss <- NULL
  }

  if (!is.null(agg_gss)) {
    gamma_window_gss <- compute_gamma_sliding(agg_gss, TARGET_IUL)
  } else {
    gamma_window_gss <- data.frame()
  }

  cat("\n=== GSS: γ Stability Across MSP (IUL = ", TARGET_IUL, ") ===\n")
  if (nrow(gamma_window_gss) > 0) {
    print(gamma_window_gss)
    cat("\nGamma Mean:", sprintf("%.3f ± %.3f\n",
        mean(gamma_window_gss$Gamma, na.rm=T),
        sd(gamma_window_gss$Gamma, na.rm=T)))

    write.csv(gamma_window_gss,
              paste0(DATA_OUT_DIR, "gamma_stability_gss_iul", sprintf("%.3f", TARGET_IUL), ".csv"),
              row.names = FALSE)

    p_gss <- ggplot(gamma_window_gss, aes(x = MSP_Window_Center, y = Gamma)) +
      geom_point(color = "blue", size = 2) +
      geom_line(color = "blue", alpha = 0.5) +
      geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
      geom_ribbon(aes(ymin = Gamma - sd(Gamma, na.rm=T),
                      ymax = Gamma + sd(Gamma, na.rm=T)),
                  alpha = 0.2, fill = "blue") +
      labs(
        title = paste0("GSS: γ Stability (IUL = ", sprintf("%.3f", TARGET_IUL), ")"),
        x = "MSP Window Center",
        y = "γ Exponent",
        subtitle = "Red dashed = Mean-Field prediction"
      ) +
      theme_minimal()

    print(p_gss)
    ggsave(paste0(PLOTS_DIR, "gamma_stability_gss.pdf"), p_gss, width = 8, height = 5)
  } else {
    cat("ERROR: No data generated for GSS stability check.\n")
  }
} else {
  cat("WARNING: GSS results file not found at ", file_gss, "\n")
}

# DP data
file_dp <- paste0(RESULTS_DIR_DP, "results_m03-06_all_random.rds")  # Try all_random first
if (!file.exists(file_dp)) {
  file_dp <- paste0(RESULTS_DIR_DP, "results_m03-06_sd0.12_random.rds")  # Fallback
}
if (file.exists(file_dp)) {
  cat("Loading DP data from:", file_dp, "\n")
  raw_dp <- readRDS(file_dp)

  # Extract the correct dataframe
  if (inherits(raw_dp, "list")) {
    mean_key <- "mean_0.40"
    if (!is.null(raw_dp[[mean_key]])) {
      raw_dp <- raw_dp[[mean_key]]
    } else {
      # If mean_0.40 doesn't exist, try to grab the first available
      raw_dp <- raw_dp[[1]]
    }
  }

  # Ensure it's a dataframe
  if (!is.data.frame(raw_dp)) {
    cat("ERROR: raw_dp is not a dataframe. Structure:\n")
    str(raw_dp)
    raw_dp <- NULL
  }

  if (!is.null(raw_dp)) {
    agg_dp <- raw_dp %>%
      mutate(adopters_prop = num_adopters / N_nodes_actual) %>%
      group_by(innovation_iul_Gamma, social_distance_h) %>%
      summarise(
        sd_adoption = sd(adopters_prop, na.rm = TRUE),
        avg_adoption = mean(adopters_prop, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    agg_dp <- NULL
  }

  if (!is.null(agg_dp)) {
    gamma_window_dp <- compute_gamma_sliding(agg_dp, TARGET_IUL)
  } else {
    gamma_window_dp <- data.frame()
  }

  cat("\n=== DP: γ Stability Across MSP (IUL = ", TARGET_IUL, ") ===\n")
  if (nrow(gamma_window_dp) > 0) {
    print(gamma_window_dp)
    cat("\nGamma Mean:", sprintf("%.3f ± %.3f\n",
        mean(gamma_window_dp$Gamma, na.rm=T),
        sd(gamma_window_dp$Gamma, na.rm=T)))

    write.csv(gamma_window_dp,
              paste0(DATA_OUT_DIR, "gamma_stability_dp_iul", sprintf("%.3f", TARGET_IUL), ".csv"),
              row.names = FALSE)

    p_dp <- ggplot(gamma_window_dp, aes(x = MSP_Window_Center, y = Gamma)) +
      geom_point(color = "red", size = 2) +
      geom_line(color = "red", alpha = 0.5) +
      geom_hline(yintercept = 1, color = "blue", linetype = "dashed") +
      geom_ribbon(aes(ymin = Gamma - sd(Gamma, na.rm=T),
                      ymax = Gamma + sd(Gamma, na.rm=T)),
                  alpha = 0.2, fill = "red") +
      labs(
        title = paste0("DP: γ Stability (IUL = ", sprintf("%.3f", TARGET_IUL), ")"),
        x = "MSP Window Center",
        y = "γ Exponent",
        subtitle = "Blue dashed = Mean-Field prediction"
      ) +
      theme_minimal()

    print(p_dp)
    ggsave(paste0(PLOTS_DIR, "gamma_stability_dp.pdf"), p_dp, width = 8, height = 5)
  } else {
    cat("ERROR: No data generated for DP stability check.\n")
  }
} else {
  cat("WARNING: DP results file not found at ", file_dp, "\n")
}

cat("\n✓ γ Stability checks complete. See playground/checks-claude/plots/ and data/\n")
