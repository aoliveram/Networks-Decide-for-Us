# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to generate Variance/SD Heatmap vs ER Baseline for CPIN2
# REDUCED VERSION: Only Empirical Networks to save height
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# --- Configuration ---
RESULTS_DIR_EMPIRICAL <- "output/04_GSS_diffusion_sims/"
OUTPUT_PDF <- "paper-NDFU/abstract-cpin2/figure.pdf"
SEEDING_STRATEGY <- "random"
THRESHOLD_MEAN_TARGET <- 0.4
TAU_SD_TARGETS <- c(0.12)
PHASE_TRANSITION_THRESHOLD_JUMP <- 1 / 3

# Parameters for plotting
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0 / 12, 12 / 12, by = 1 / 12)

# Function to load and process data
process_data <- function(results_dir, label) {
  file_path <- paste0(results_dir, "results_m03-06_sd", sprintf("%.2f", TAU_SD_TARGETS[1]), "_", SEEDING_STRATEGY, ".rds")
  if (!file.exists(file_path)) stop("File not found: ", file_path)

  raw_df <- readRDS(file_path)

  mean_label <- paste0("mean_", sprintf("%.2f", THRESHOLD_MEAN_TARGET))
  if (!is.null(raw_df[[mean_label]])) raw_df <- raw_df[[mean_label]] # Handle nested list structure if present

  if (inherits(raw_df, "list")) {
    # If it's a list with mean_labels, extract target
    raw_df <- raw_df[[mean_label]]
  }

  base_run_summary <- raw_df %>%
    group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
    summarise(
      adopters_prop_at_cell = first(num_adopters) / first(N_nodes_actual),
      prop_rational_of_pop_at_cell = first(num_adopted_rational) / first(N_nodes_actual),
      prop_social_of_pop_at_cell = first(num_adopted_social) / first(N_nodes_actual),
      .groups = "drop"
    ) %>%
    arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
    group_by(run_id, social_distance_h) %>%
    mutate(
      jump_at_step = adopters_prop_at_cell - lag(adopters_prop_at_cell),
      is_transition_vs_prev_gamma = ifelse(!is.na(jump_at_step) & jump_at_step >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
    ) %>%
    ungroup()

  # Calculate SD for Phase Transition and Total Adoption
  sd_df <- base_run_summary %>%
    group_by(innovation_iul_Gamma, social_distance_h) %>%
    summarise(
      sd_adoption = sd(adopters_prop_at_cell, na.rm = TRUE),
      sd_transition = sd(is_transition_vs_prev_gamma, na.rm = TRUE),
      avg_adoption = mean(adopters_prop_at_cell, na.rm = T),
      prop_transition = mean(is_transition_vs_prev_gamma, na.rm = T),
      avg_rational = mean(prop_rational_of_pop_at_cell, na.rm = T),
      avg_social = mean(prop_social_of_pop_at_cell, na.rm = T),
      .groups = "drop"
    ) %>%
    mutate(network = label)

  return(sd_df)
}

df_emp <- process_data(RESULTS_DIR_EMPIRICAL, "Empirical")

# --- Plotting ---
create_heatmap <- function(df, fill_col, title, limit_max, show_y = TRUE, show_legend = FALSE, color_opt = "viridis") {
  df <- df %>% filter(innovation_iul_Gamma <= 0.75)
  h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
  df$social_distance_h_factor <- factor(sprintf("%.2f", df$social_distance_h), levels = h_levels_sorted)

  y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)])
  x_breaks <- seq(0, 0.75, 0.25)

  p <- ggplot(df, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_col]])) +
    geom_tile(color = "white", lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0, limit_max), option = color_opt, na.value = "grey90", name = "") +
    labs(x = expression(paste("IUL (", Gamma, ")")), y = if (show_y) "MSP (h)" else NULL, title = title) +
    scale_x_continuous(breaks = x_breaks, labels = sprintf("%.2f", x_breaks), expand = c(0, 0)) +
    scale_y_discrete(drop = FALSE, breaks = y_breaks, labels = if (show_y) y_breaks else NULL) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(size = 6, color = if (show_y) "black" else "transparent"),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold", angle = 90),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.text = element_text(size = 5),
      plot.margin = margin(t = 2, r = 2, b = 2, l = if (show_y) 2 else 0)
    )
  return(p)
}

max_unified <- 1.0
max_sd <- max(df_emp$sd_adoption, na.rm = TRUE)

# Plots Empirical
p_mean_emp <- create_heatmap(df_emp, "avg_adoption", "Avg Total Adopt.", max_unified, show_y = TRUE, show_legend = FALSE, color_opt = "viridis")
p_rat_emp <- create_heatmap(df_emp, "avg_rational", "Rational Choice", max_unified, show_y = FALSE, show_legend = FALSE, color_opt = "viridis")
p_soc_emp <- create_heatmap(df_emp, "avg_social", "Social Influence", max_unified, show_y = FALSE, show_legend = TRUE, color_opt = "viridis")
p_sd_emp <- create_heatmap(df_emp, "sd_adoption", "Susceptibility", max_sd, show_y = FALSE, show_legend = TRUE, color_opt = "magma")

design <- "ABCD"
final_plot <- (p_mean_emp + p_rat_emp + p_soc_emp + p_sd_emp) +
  plot_layout(design = design, widths = c(1, 1, 1.15, 1.15))

# Extremely condensed layout
ggsave(OUTPUT_PDF, final_plot, width = 7.0, height = 1.7)
