# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to Analyze Raw Results and Generate PDF with Heatmaps (GSS)
#
# Helper script that exposes 'generate_diffusion_plots' (GSS version).
# Expects a main driver script to call this function.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)
library(doParallel)

generate_diffusion_plots <- function(
    RESULTS_DIR,
    PLOTS_DIR,
    SEEDING_STRATEGY_FIXED,
    NETWORK_LABEL = "GSS-net",
    NUM_CORES = 8) {
  # --- Visualization Parameters ---
  IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
  H_VALUES_SWEEP <- seq(0 / 12, 12 / 12, by = 1 / 12)
  THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
  TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20)
  PHASE_TRANSITION_THRESHOLD_JUMP <- 1 / 3

  dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

  cat(paste0("\n====================================================================\n"))
  cat(paste0("Processing Plots for Strategy: ", SEEDING_STRATEGY_FIXED, "\n"))
  cat(paste0("Results Source: ", RESULTS_DIR, "\n"))
  cat(paste0("====================================================================\n"))

  cat("Loading raw results...\n")
  # Note: GSS filenames have 'GSS_phase_transition_GRAND_COMBINED_raw_results_' prefix
  # updated to match the commit for GSS diffusion sims (checked in Git history/previous context)
  # Actually, looking at 04_GSS_diffusion_sims.R commit:
  # "GSS_phase_transition_GRAND_COMBINED_raw_results_" + SEEDING_STRATEGY_FIXED + ".rds"

  grand_raw_results_path <- paste0(RESULTS_DIR, "results_m03-06_all_", SEEDING_STRATEGY_FIXED, ".rds")
  # Wait, in the GSS script 04_GSS_diffusion_sims main content I wrote:
  # saveRDS(all_sds_raw_results_list, paste0(RESULTS_DIR, "results_m03-06_all_", SEEDING_STRATEGY_FIXED, ".rds"))
  # So the name matches ATP structure now!

  if (!file.exists(grand_raw_results_path)) {
    # Fallback to check if it's the old name just in case? No, I updated the script.
    cat(paste0("ERROR: Raw results file not found: ", grand_raw_results_path, "\n"))
    return(NULL)
  }
  all_sds_raw_results_list_from_file <- readRDS(grand_raw_results_path)
  cat("Results loaded.\n")

  # --- STEP 0: Pre-process raw data (PARALLEL) ---
  cat(paste0("\nStarting parallel processing (", NUM_CORES, " cores)...\n"))

  # Using FORK cluster for efficiency on Unix-like systems
  cl <- makeCluster(NUM_CORES, type = "FORK")
  registerDoParallel(cl)

  processed_results_by_sd <- foreach(
    sd_idx = 1:length(TAU_NORMAL_SD_SWEEP_LIST),
    .packages = c("dplyr", "readr"),
    .export = c("PHASE_TRANSITION_THRESHOLD_JUMP", "IUL_VALUES_SWEEP", "H_VALUES_SWEEP", "THRESHOLD_MEAN_SWEEP_LIST")
  ) %dopar% {
    current_tau_sd_proc <- TAU_NORMAL_SD_SWEEP_LIST[sd_idx]
    sd_label_proc <- paste0("sd_", sprintf("%.2f", current_tau_sd_proc))
    current_sd_all_means_raw_results_proc <- all_sds_raw_results_list_from_file[[sd_label_proc]]

    if (is.null(current_sd_all_means_raw_results_proc)) {
      return(NULL)
    }

    heatmap_data_tm_this_sd_list_proc <- list()
    heatmap_data_aa_this_sd_list_proc <- list()
    heatmap_data_ar_pop_this_sd_list_proc <- list()
    heatmap_data_as_pop_this_sd_list_proc <- list()

    for (current_threshold_mean_proc in THRESHOLD_MEAN_SWEEP_LIST) {
      mean_label_proc <- paste0("mean_", sprintf("%.2f", current_threshold_mean_proc))
      raw_df_this_mean_sd_proc <- current_sd_all_means_raw_results_proc[[mean_label_proc]]
      if (is.null(raw_df_this_mean_sd_proc) || nrow(raw_df_this_mean_sd_proc) == 0) {
        next
      }

      NUM_RUNS_THIS_COMBO_ACTUAL_PROC <- length(unique(raw_df_this_mean_sd_proc$run_id))
      if (NUM_RUNS_THIS_COMBO_ACTUAL_PROC == 0) next

      # General Pre-calc
      base_run_summary_proc <- raw_df_this_mean_sd_proc %>%
        group_by(run_id, social_distance_h, innovation_iul_Gamma) %>%
        summarise(
          adopters_prop_at_cell = first(num_adopters) / first(N_nodes_actual),
          prop_rational_of_pop_at_cell = first(num_adopted_rational) / first(N_nodes_actual),
          prop_social_of_pop_at_cell = first(num_adopted_social) / first(N_nodes_actual),
          num_adopters_at_cell = first(num_adopters),
          n_nodes_at_cell = first(N_nodes_actual),
          .groups = "drop"
        ) %>%
        arrange(run_id, social_distance_h, innovation_iul_Gamma) %>%
        group_by(run_id, social_distance_h) %>%
        mutate(
          jump_at_step = adopters_prop_at_cell - lag(adopters_prop_at_cell),
          is_transition_vs_prev_gamma = ifelse(!is.na(jump_at_step) & jump_at_step >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
        ) %>%
        group_by(run_id, social_distance_h) %>%
        mutate(first_transition_IUL_for_series = if (any(is_transition_vs_prev_gamma == 1, na.rm = TRUE)) {
          min(innovation_iul_Gamma[which(is_transition_vs_prev_gamma == 1)])
        } else {
          NA_real_
        }) %>%
        ungroup() %>%
        # Transition by H
        arrange(run_id, innovation_iul_Gamma, social_distance_h) %>%
        group_by(run_id, innovation_iul_Gamma) %>%
        mutate(
          jump_vs_prev_h = adopters_prop_at_cell - lag(adopters_prop_at_cell),
          is_transition_vs_prev_h = ifelse(!is.na(jump_vs_prev_h) & jump_vs_prev_h >= PHASE_TRANSITION_THRESHOLD_JUMP, 1, 0)
        ) %>%
        group_by(run_id, innovation_iul_Gamma) %>%
        mutate(first_transition_H_for_series = if (any(is_transition_vs_prev_h == 1, na.rm = TRUE)) {
          min(social_distance_h[which(is_transition_vs_prev_h == 1)])
        } else {
          NA_real_
        }) %>%
        ungroup()


      panel_data_tm_list_proc <- list()
      panel_data_aa_list_proc <- list()
      panel_data_ar_pop_list_proc <- list()
      panel_data_as_pop_list_proc <- list()

      for (iul_val_proc in IUL_VALUES_SWEEP) {
        for (h_val_proc in H_VALUES_SWEEP) {
          # Metric 1: Transition Metric
          runs_transitioned_to_this_cell_by_gamma <- base_run_summary_proc %>%
            filter(
              social_distance_h == h_val_proc,
              !is.na(first_transition_IUL_for_series),
              first_transition_IUL_for_series == iul_val_proc
            ) %>%
            pull(run_id)
          runs_transitioned_to_this_cell_by_h <- base_run_summary_proc %>%
            filter(
              innovation_iul_Gamma == iul_val_proc,
              !is.na(first_transition_H_for_series),
              first_transition_H_for_series == h_val_proc
            ) %>%
            pull(run_id)
          unique_runs_transitioned_to_cell <- unique(c(runs_transitioned_to_this_cell_by_gamma, runs_transitioned_to_this_cell_by_h))
          prop_tm_cell_combined <- length(unique_runs_transitioned_to_cell) / NUM_RUNS_THIS_COMBO_ACTUAL_PROC
          panel_data_tm_list_proc[[length(panel_data_tm_list_proc) + 1]] <- data.frame(iul = iul_val_proc, h = h_val_proc, val = prop_tm_cell_combined)

          # Metrics 2, 3, 4
          current_cell_data_proc <- base_run_summary_proc %>%
            filter(innovation_iul_Gamma == iul_val_proc, social_distance_h == h_val_proc)

          # Metric 2: Average Adoption (total)
          prop_aa_cell <- NA_real_
          if (nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {
            prop_aa_cell <- mean(current_cell_data_proc$adopters_prop_at_cell, na.rm = TRUE)
            if (is.nan(prop_aa_cell)) prop_aa_cell <- NA_real_
          }
          panel_data_aa_list_proc[[length(panel_data_aa_list_proc) + 1]] <- data.frame(iul = iul_val_proc, h = h_val_proc, val = prop_aa_cell)

          # Metric 3: Avg. Adoption by Rational Choice (prop. of population)
          prop_ar_pop_cell <- NA_real_
          if (nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {
            prop_ar_pop_cell <- mean(current_cell_data_proc$prop_rational_of_pop_at_cell, na.rm = TRUE)
            if (is.nan(prop_ar_pop_cell)) prop_ar_pop_cell <- NA_real_
          }
          panel_data_ar_pop_list_proc[[length(panel_data_ar_pop_list_proc) + 1]] <- data.frame(iul = iul_val_proc, h = h_val_proc, val = prop_ar_pop_cell)

          # Metric 4: Avg. Adoption by Social Influence (prop. of population)
          prop_as_pop_cell <- NA_real_
          if (nrow(current_cell_data_proc) == NUM_RUNS_THIS_COMBO_ACTUAL_PROC) {
            prop_as_pop_cell <- mean(current_cell_data_proc$prop_social_of_pop_at_cell, na.rm = TRUE)
            if (is.nan(prop_as_pop_cell)) prop_as_pop_cell <- NA_real_
          }
          panel_data_as_pop_list_proc[[length(panel_data_as_pop_list_proc) + 1]] <- data.frame(iul = iul_val_proc, h = h_val_proc, val = prop_as_pop_cell)
        }
      }
      heatmap_data_tm_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_tm_list_proc) %>%
        mutate(tau_mean_param = current_threshold_mean_proc, proportion_value_to_plot = val, tau_sd_param = current_tau_sd_proc) %>%
        select(-val)
      heatmap_data_aa_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_aa_list_proc) %>%
        mutate(tau_mean_param = current_threshold_mean_proc, mean_adopters_prop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>%
        select(-val)
      heatmap_data_ar_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_ar_pop_list_proc) %>%
        mutate(tau_mean_param = current_threshold_mean_proc, avg_rational_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>%
        select(-val)
      heatmap_data_as_pop_this_sd_list_proc[[mean_label_proc]] <- bind_rows(panel_data_as_pop_list_proc) %>%
        mutate(tau_mean_param = current_threshold_mean_proc, avg_social_adopt_pop_to_plot = val, tau_sd_param = current_tau_sd_proc) %>%
        select(-val)
    }

    # Return combined list for this SD
    list(
      tm = bind_rows(heatmap_data_tm_this_sd_list_proc[!sapply(heatmap_data_tm_this_sd_list_proc, is.null)]),
      aa = bind_rows(heatmap_data_aa_this_sd_list_proc[!sapply(heatmap_data_aa_this_sd_list_proc, is.null)]),
      ar = bind_rows(heatmap_data_ar_pop_this_sd_list_proc[!sapply(heatmap_data_ar_pop_this_sd_list_proc, is.null)]),
      as = bind_rows(heatmap_data_as_pop_this_sd_list_proc[!sapply(heatmap_data_as_pop_this_sd_list_proc, is.null)])
    )
  }

  stopCluster(cl)
  cat("Parallel processing completed. Aggregating results...\n")

  all_sds_transition_metric_heatmap_df_list <- list()
  all_sds_avg_adoption_heatmap_df_list <- list()
  all_sds_avg_rational_adopt_pop_heatmap_df_list <- list()
  all_sds_avg_social_adopt_pop_heatmap_df_list <- list()

  # Unpack results
  for (i in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
    current_tau_sd <- TAU_NORMAL_SD_SWEEP_LIST[i]
    sd_label <- paste0("sd_", sprintf("%.2f", current_tau_sd))
    res <- processed_results_by_sd[[i]]

    if (!is.null(res)) {
      all_sds_transition_metric_heatmap_df_list[[sd_label]] <- res$tm
      all_sds_avg_adoption_heatmap_df_list[[sd_label]] <- res$aa
      all_sds_avg_rational_adopt_pop_heatmap_df_list[[sd_label]] <- res$ar
      all_sds_avg_social_adopt_pop_heatmap_df_list[[sd_label]] <- res$as
    }
  }
  cat("Data preparation complete.\n")

  # --- Single Heatmap Function ---
  create_single_heatmap_v3 <- function(df_plot_data, fill_col_name, legend_title_text, viridis_option,
                                       show_legend = TRUE, y_axis_label_on = TRUE, x_axis_label_on = TRUE,
                                       panel_row_title = "") {
    if (is.null(df_plot_data) || nrow(df_plot_data) == 0 || all(is.na(df_plot_data[[fill_col_name]]))) {
      return(ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No plottable data") +
        theme_void() +
        labs(title = NULL, y = if (y_axis_label_on) panel_row_title else NULL) +
        theme(axis.title.y = element_text(size = 8, face = "bold", angle = 90, margin = margin(r = 5))))
    }

    if ("iul" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, innovation_iul_Gamma = iul)
    if ("h" %in% names(df_plot_data)) df_plot_data <- rename(df_plot_data, social_distance_h = h)

    # Ensure factor levels
    h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
    df_plot_data$social_distance_h_factor <- factor(sprintf("%.2f", df_plot_data$social_distance_h), levels = h_levels_sorted)

    # Axes definitions
    y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)])
    y_labels <- y_breaks
    x_breaks <- seq(0, 1, 0.25)
    x_labels <- sprintf("%.2f", x_breaks)

    plot_obj <- ggplot(df_plot_data, aes(x = innovation_iul_Gamma, y = social_distance_h_factor, fill = .data[[fill_col_name]])) +
      geom_tile(color = "white", lwd = 0.1) +
      scale_fill_viridis_c(
        name = if (show_legend) legend_title_text else NULL,
        limits = c(0, 1), option = viridis_option, n.breaks = 4, na.value = "grey90"
      ) +
      labs(
        x = if (x_axis_label_on) expression(paste("IUL (", Gamma, ")")) else NULL,
        y = if (y_axis_label_on) panel_row_title else NULL,
        title = NULL
      ) +
      scale_x_continuous(breaks = x_breaks, labels = if (x_axis_label_on) x_labels else NULL, expand = c(0, 0)) +
      scale_y_discrete(drop = FALSE, breaks = y_breaks, labels = if (y_axis_label_on) y_labels else NULL) +
      theme_minimal(base_size = 7) +
      theme(
        axis.text.x = element_text(
          angle = 45, hjust = 1, size = 6,
          color = if (x_axis_label_on) "black" else "transparent"
        ),
        axis.text.y = element_text(
          size = 6,
          color = if (y_axis_label_on) "black" else "transparent"
        ),
        axis.title.x = element_text(size = 8, face = "bold", margin = margin(t = 2, unit = "mm")),
        axis.title.y = element_text(size = 8, face = "bold", angle = 90, margin = margin(r = 2, unit = "mm")),
        legend.position = if (show_legend) "right" else "none",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.6, "lines"),
        panel.grid = element_line(linewidth = 0.15),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm")
      )
    return(plot_obj)
  }


  # --- Main Plotting Loop ---
  for (current_threshold_mean_plot in THRESHOLD_MEAN_SWEEP_LIST) {
    cat(paste0("Generating consolidated PDF for Mean Tau = ", current_threshold_mean_plot, "\n"))

    plot_list_for_this_mean_pdf_ordered <- list()

    metric_titles_for_cols <- c(
      "Avg. Adoption",
      "Phase Trans Prob.",
      "Avg. Adopt.\nby Rational",
      "Avg. Adopt.\nby Social Infl."
    )

    metric_fill_vars_for_cols <- c(
      "mean_adopters_prop_to_plot",
      "proportion_value_to_plot",
      "avg_rational_adopt_pop_to_plot",
      "avg_social_adopt_pop_to_plot"
    )

    metric_color_options_for_cols <- rep("viridis", 4)

    data_sources_for_cols <- list(
      all_sds_avg_adoption_heatmap_df_list,
      all_sds_transition_metric_heatmap_df_list,
      all_sds_avg_rational_adopt_pop_heatmap_df_list,
      all_sds_avg_social_adopt_pop_heatmap_df_list
    )

    for (row_idx in 1:length(TAU_NORMAL_SD_SWEEP_LIST)) {
      current_tau_sd_plot <- TAU_NORMAL_SD_SWEEP_LIST[row_idx]
      sd_label_plot <- paste0("sd_", sprintf("%.2f", current_tau_sd_plot))

      for (col_idx in 1:4) {
        df_all_means_for_sd_metric <- data_sources_for_cols[[col_idx]][[sd_label_plot]]

        current_df_for_panel <- if (!is.null(df_all_means_for_sd_metric) && nrow(df_all_means_for_sd_metric) > 0) {
          df_all_means_for_sd_metric %>% filter(tau_mean_param == current_threshold_mean_plot)
        } else {
          data.frame(innovation_iul_Gamma = numeric(0), social_distance_h = numeric(0)) %>%
            mutate(!!metric_fill_vars_for_cols[col_idx] := numeric(0))
        }

        row_title_str <- if (col_idx == 1) paste0("MSP (h) - SD=", sprintf("%.2f", current_tau_sd_plot)) else ""

        y_label_visible <- (col_idx == 1)
        x_label_visible <- (row_idx == length(TAU_NORMAL_SD_SWEEP_LIST))
        legend_visible <- (col_idx == 4)

        plot_index_in_list <- ((row_idx - 1) * 4) + col_idx

        plot_list_for_this_mean_pdf_ordered[[plot_index_in_list]] <-
          create_single_heatmap_v3(
            df_plot_data = current_df_for_panel,
            fill_col_name = metric_fill_vars_for_cols[col_idx],
            legend_title_text = NULL,
            viridis_option = metric_color_options_for_cols[col_idx],
            show_legend = legend_visible,
            y_axis_label_on = y_label_visible,
            x_axis_label_on = x_label_visible,
            panel_row_title = row_title_str
          )
      }
    }

    # --- Assembly with Patchwork ---
    if (length(plot_list_for_this_mean_pdf_ordered) == (length(TAU_NORMAL_SD_SWEEP_LIST) * 4)) {
      col_titles_plots <- lapply(metric_titles_for_cols, function(title) {
        ggplot() +
          labs(title = title) +
          theme_void() +
          theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold", margin = margin(b = 0, t = 2, unit = "mm")))
      })

      column_titles_row_layout <- Reduce(`+`, col_titles_plots) + plot_layout(ncol = 4)
      heatmaps_grid_layout <- wrap_plots(plot_list_for_this_mean_pdf_ordered, ncol = 4, byrow = TRUE)

      final_combined_layout <- column_titles_row_layout / heatmaps_grid_layout +
        plot_layout(heights = c(0.05, 1))

      # Try to get number of runs for subtitle
      num_runs_val_subtitle <- NA
      first_valid_sd_label <- names(all_sds_raw_results_list_from_file)[which(!sapply(all_sds_raw_results_list_from_file, is.null))[1]]
      if (!is.na(first_valid_sd_label) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]])) {
        first_valid_mean_label_for_runs <- names(all_sds_raw_results_list_from_file[[first_valid_sd_label]])[which(!sapply(all_sds_raw_results_list_from_file[[first_valid_sd_label]], is.null))[1]]
        if (!is.na(first_valid_mean_label_for_runs) && !is.null(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label_for_runs]])) {
          num_runs_val_subtitle <- length(unique(all_sds_raw_results_list_from_file[[first_valid_sd_label]][[first_valid_mean_label_for_runs]]$run_id))
        }
      }
      if (is.na(num_runs_val_subtitle)) num_runs_val_subtitle <- "N/A"

      final_plot_with_annotation <- final_combined_layout +
        plot_annotation(
          title = paste("Consolidated Heatmaps for", NETWORK_LABEL, "- Mean threshold =", sprintf("%.2f", current_threshold_mean_plot)),
          subtitle = paste("Thresholds ~ N(mu=", sprintf("%.2f", current_threshold_mean_plot), ", SD=var). ",
            num_runs_val_subtitle, " runs per cell. Seeding strategy: ", SEEDING_STRATEGY_FIXED,
            sep = ""
          ),
          theme = theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = element_text(hjust = 0.5, size = 9.5)
          )
        )

      # Filename format: plots_mXX_strategy.pdf
      mean_formatted <- sprintf("%02d", as.integer(current_threshold_mean_plot * 10)) # 0.3 -> 03
      plot_filename <- paste0(PLOTS_DIR, "plots_m", mean_formatted, "_", SEEDING_STRATEGY_FIXED, ".pdf")

      ggsave(plot_filename, final_plot_with_annotation, width = 7.5, height = 7.0, limitsize = FALSE)
      cat(paste0("  Saved: ", plot_filename, "\n"))
    } else {
      cat(paste0("  Warning: Not enough plots generated for Mean Tau = ", current_threshold_mean_plot, "\n"))
    }
  }
  cat(paste0("\nAll PDF generations finished for strategy: ", SEEDING_STRATEGY_FIXED, "\n"))
}
