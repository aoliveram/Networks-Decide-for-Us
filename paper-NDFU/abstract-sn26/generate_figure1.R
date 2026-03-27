# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate 5-panel Figure 1 for abstract_SN26
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# --- Configuration ---
RESULTS_DIR <- "output/04_GSS_ER_diffusion_sims/"
OUTPUT_PDF <- "paper-NDFU/abstract-sn26/figure1_5panels_GSS_ER.pdf"
SEEDING_STRATEGY <- "random"
THRESHOLD_MEAN_TARGET <- 0.4
TAU_SD_TARGET <- 0.12

# --- Data Loading ---
file_path <- paste0(
  RESULTS_DIR, "results_m03-06_sd",
  sprintf("%.2f", TAU_SD_TARGET), "_", SEEDING_STRATEGY, ".rds"
)
if (!file.exists(file_path)) stop("File not found: ", file_path)

raw_df <- readRDS(file_path)
mean_label <- paste0("mean_", sprintf("%.2f", THRESHOLD_MEAN_TARGET))
if (inherits(raw_df, "list") && !is.null(raw_df[[mean_label]])) {
  raw_df <- raw_df[[mean_label]]
}

base_df <- raw_df %>%
  mutate(
    adopters_prop = num_adopters / N_nodes_actual,
    prop_rational = num_adopted_rational / N_nodes_actual,
    prop_social   = num_adopted_social / N_nodes_actual
  )

# --- Aggregation ---
agg_df <- base_df %>%
  group_by(innovation_iul_Gamma, social_distance_h) %>%
  summarise(
    avg_adoption = mean(adopters_prop, na.rm = TRUE),
    avg_rational = mean(prop_rational, na.rm = TRUE),
    avg_social = mean(prop_social, na.rm = TRUE),
    susceptibility = var(adopters_prop, na.rm = TRUE), # chi via FDT
    sd_steps = sd(num_steps, na.rm = TRUE), # Method 2 CSD
    .groups = "drop"
  )

# --- Plotting Helpers ---
h_levels <- sprintf("%.2f", sort(unique(agg_df$social_distance_h)))
y_breaks <- h_levels[seq(1, length(h_levels), by = 3)]

make_heatmap <- function(df, fill_col, title, legend_label,
                         show_y = FALSE, palette = "viridis", limits = NULL) {
  df2 <- df %>%
    mutate(h_factor = factor(sprintf("%.2f", social_distance_h), levels = h_levels))

  p <- ggplot(df2, aes(x = innovation_iul_Gamma, y = h_factor, fill = .data[[fill_col]])) +
    geom_tile(color = "white", lwd = 0.05) +
    scale_x_continuous(
      expand = c(0, 0),
      limits = c(0, 0.75),
      breaks = c(0, 0.25, 0.5, 0.75),
      labels = c("0", "0.25", "0.5", "0.75")
    ) +
    scale_y_discrete(expand = c(0, 0), breaks = y_breaks) +
    labs(
      title = title,
      x     = expression(paste("IUL (", Gamma, ")")),
      y     = if (show_y) "MSP (h)" else NULL,
      fill  = legend_label
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title      = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 20)),
      axis.text       = element_text(size = 6),
      axis.title      = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.title    = element_text(size = 6),
      legend.text     = element_text(size = 5),
      panel.grid      = element_blank(),
      aspect.ratio    = 0.66,
      plot.margin     = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    )

  if (!is.null(limits)) {
    p <- p + scale_fill_viridis_c(option = palette, limits = limits, name = legend_label)
  } else {
    p <- p + scale_fill_viridis_c(option = palette, name = legend_label)
  }
  p
}

# --- Generate Panels ---
p1 <- make_heatmap(agg_df, "avg_adoption", "Total Adoption", "Prop.", show_y = TRUE, limits = c(0, 1))
p2 <- make_heatmap(agg_df, "avg_rational", "Rational Adoption", "Prop.", limits = c(0, 1))
p3 <- make_heatmap(agg_df, "avg_social", "Social Adoption", "Prop.", limits = c(0, 1))

p4 <- make_heatmap(agg_df, "susceptibility", expression(chi), expression(chi), show_y = TRUE, palette = "magma")
p5 <- make_heatmap(agg_df, "sd_steps", "SD Steps", "SD", palette = "magma")

# --- Layout Assembly ---
# Patchwork row wrappers to manage guides properly
row1 <- (p1 | p2 | p3) +
  plot_layout(guides = "collect")

row2 <- (plot_spacer() | p4 | plot_spacer() | p5 | plot_spacer()) +
  plot_layout(widths = c(0.1, 1, 0.01, 1, 0.1))

final_fig <- (row1 / row2) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Emergence of Adoption and Phase Transition Signatures",
    theme = theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 15))
    )
  )

# --- Save Output ---
ggsave(OUTPUT_PDF, final_fig, width = 8, height = 4.5)
cat("Figure saved to:", OUTPUT_PDF, "\n")
