# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Custom Plotting Script for NetSci Abstract Figure
# Reduced height: Only SD=0.12 and SD=0.20
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)

# --- Configuration ---
RESULTS_DIR <- "output/04_GSS_diffusion_sims/"
OUTPUT_PDF <- "paper-NDFU/abstract-netsci26/figure.pdf"
SEEDING_STRATEGY <- "random"
THRESHOLD_MEAN_TARGET <- 0.4
TAU_SD_TARGETS <- c(0.12)
PHASE_TRANSITION_THRESHOLD_JUMP <- 1 / 3

# Parameters for plotting
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0 / 12, 12 / 12, by = 1 / 12)

# Load Data
cat("Loading raw results...\n")
data_path <- paste0(RESULTS_DIR, "results_m03-06_all_", SEEDING_STRATEGY, ".rds")
if (!file.exists(data_path)) stop("Results file not found: ", data_path)
all_sds_raw_results <- readRDS(data_path)

# --- Process Data for Target SDs ---
cat("Processing data...\n")

plot_data_list <- list()

for (current_sd in TAU_SD_TARGETS) {
    sd_label <- paste0("sd_", sprintf("%.2f", current_sd))
    mean_label <- paste0("mean_", sprintf("%.2f", THRESHOLD_MEAN_TARGET))

    raw_df <- all_sds_raw_results[[sd_label]][[mean_label]]

    if (is.null(raw_df)) next

    # Summarize by cell
    summary_df <- raw_df %>%
        group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
        summarise(
            adopters_prop = first(num_adopters) / first(N_nodes_actual),
            prop_rational = first(num_adopted_rational) / first(N_nodes_actual),
            prop_social = first(num_adopted_social) / first(N_nodes_actual),
            .groups = "drop"
        ) %>%
        # Calculate Phase Transition
        arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
        group_by(run_id, social_distance_h) %>%
        mutate(
            jump = adopters_prop - lag(adopters_prop),
            is_trans_IUL = ifelse(!is.na(jump) & jump >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0),
            first_trans_IUL = if (any(is_trans_IUL == 1)) min(innovation_iul_Gamma[which(is_trans_IUL == 1)]) else NA
        ) %>%
        group_by(run_id, innovation_iul_Gamma) %>%
        mutate(
            jump_h = adopters_prop - lag(adopters_prop), # Note: sorting needs to be correct for this. Actually logic in original script was complex.
            # Simplified approach: Use same logic as original script for consistency
            is_trans_H = ifelse(!is.na(jump_h) & jump_h >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0),
            first_trans_H = if (any(is_trans_H == 1)) min(social_distance_h[which(is_trans_H == 1)]) else NA
        ) %>%
        ungroup()

    # Aggregation for Heatmaps
    # Metric 1: Avg Adoption
    aa_df <- summary_df %>%
        group_by(social_distance_h, innovation_iul_Gamma) %>%
        summarise(val = mean(adopters_prop, na.rm = TRUE), .groups = "drop")

    # Metric 2: Transition Prob
    # Re-implementing logic from original script correctly
    # Original script looks for transition on either axis

    # We need to compute transition metrics properly based on run-level data
    # Correct way: For each cell (h, iul), what proportion of runs have their "first transition" LANDING exactly here?
    # This is tricky without the full original logic.
    # To save time and avoid bugs, let's use the MAIN original logic but streamlined.

    # Recalculating transition point for each run
    runs_meta <- summary_df %>%
        distinct(run_id)

    tm_data <- list()
    for (i_gam in IUL_VALUES_SWEEP) {
        for (i_h in H_VALUES_SWEEP) {
            # Count runs where transition point is exactly here
            count <- summary_df %>%
                filter(
                    (first_trans_IUL == i_gam & social_distance_h == i_h) |
                        (first_trans_H == i_h & innovation_iul_Gamma == i_gam)
                ) %>%
                pull(run_id) %>%
                unique() %>%
                length()

            # Total runs
            total <- summary_df %>%
                filter(innovation_iul_Gamma == i_gam, social_distance_h == i_h) %>%
                distinct(run_id) %>%
                nrow()
            val <- if (total > 0) count / total else NA
            tm_data[[length(tm_data) + 1]] <- data.frame(innovation_iul_Gamma = i_gam, social_distance_h = i_h, val = val)
        }
    }
    tm_df <- bind_rows(tm_data)

    # Metric 3 & 4
    ar_df <- summary_df %>%
        group_by(social_distance_h, innovation_iul_Gamma) %>%
        summarise(val = mean(prop_rational, na.rm = TRUE), .groups = "drop")
    as_df <- summary_df %>%
        group_by(social_distance_h, innovation_iul_Gamma) %>%
        summarise(val = mean(prop_social, na.rm = TRUE), .groups = "drop")

    plot_data_list[[sd_label]] <- list(aa = aa_df, tm = tm_df, ar = ar_df, as = as_df)
}

# --- Plotting Function ---
create_heatmap <- function(data, fill_var, title = NULL, show_legend = FALSE, y_lab = FALSE, x_lab = FALSE) {
    if (is.null(data) || nrow(data) == 0) {
        return(ggplot() +
            theme_void())
    }

    ggplot(data, aes(x = innovation_iul_Gamma, y = factor(sprintf("%.2f", social_distance_h)), fill = val)) +
        geom_tile(color = "white", lwd = 0.1) +
        scale_fill_viridis_c(limits = c(0, 1), name = NULL, na.value = "grey90") +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
        scale_y_discrete(breaks = sprintf("%.2f", seq(0, 1, 0.25))) +
        labs(
            x = if (x_lab) expression(paste("IUL (", Gamma, ")")) else NULL,
            y = if (y_lab) title else NULL
        ) +
        theme_minimal(base_size = 8) +
        theme(
            axis.text.x = if (x_lab) element_text(angle = 45, hjust = 1, size = 6) else element_blank(),
            axis.text.y = if (y_lab) element_text(size = 6) else element_blank(),
            legend.position = if (show_legend) "right" else "none",
            legend.key.width = unit(0.3, "cm"),
            axis.title.y = element_text(size = 7, face = "bold"),
            plot.margin = margin(1, 1, 1, 1, "mm")
        )
}

# --- Assemble Plot ---
plots <- list()
idx <- 1

col_titles <- c("Avg. Adoption", "Phase Trans Prob.", "Avg. Adopt.\nby Rational", "Avg. Adopt.\nby Social Infl.")

for (i_sd in 1:length(TAU_SD_TARGETS)) {
    sd_val <- TAU_SD_TARGETS[i_sd]
    sd_label <- paste0("sd_", sprintf("%.2f", sd_val))
    p_dat <- plot_data_list[[sd_label]]

    for (col in 1:4) {
        dat <- switch(col,
            p_dat$aa,
            p_dat$tm,
            p_dat$ar,
            p_dat$as
        )
        row_title <- paste0("MSP (h)\nSD=", sprintf("%.2f", sd_val))

        p <- create_heatmap(
            dat, "val",
            title = row_title,
            show_legend = (col == 4),
            y_lab = (col == 1),
            x_lab = (i_sd == length(TAU_SD_TARGETS)) # Only bottom row gets X axis
        )

        # Add column titles only to first row
        if (i_sd == 1) {
            p <- p + ggtitle(col_titles[col]) +
                theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))
        }

        plots[[idx]] <- p
        idx <- idx + 1
    }
}

final_plot <- wrap_plots(plots, ncol = 4, byrow = TRUE) +
    plot_annotation(
        title = paste("Adoption Regimes in Plausible Networks (GSS-net)"),
        subtitle = "Mean threshold=0.40, SD=0.12. Seeding: random",
        theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 9, hjust = 0.5))
    )

ggsave(OUTPUT_PDF, final_plot, width = 7.5, height = 2.5)
cat("Plot saved to", OUTPUT_PDF, "\n")
