# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to explore Phase Transition metrics deeply
# 07_phase_transition_deeper.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# --- Configuration ---
RESULTS_DIR <- "output/04_GSS_ER_diffusion_sims/"
PLOTS_DIR <- "plots/07_phase_transition/GSS_ER/"
DATA_OUT_DIR <- "output/07_phase_transition/GSS_ER/"
SEEDING_STRATEGY <- "random"
THRESHOLD_MEAN_TARGET <- 0.4
TAU_SD_TARGETS <- c(0.12)

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0, 1, by = 1 / 12)

# Load data
file_path <- paste0(RESULTS_DIR, "results_m03-06_sd", sprintf("%.2f", TAU_SD_TARGETS[1]), "_", SEEDING_STRATEGY, ".rds")
if (!file.exists(file_path)) stop("File not found")

raw_df <- readRDS(file_path)
mean_label <- paste0("mean_", sprintf("%.2f", THRESHOLD_MEAN_TARGET))
if (inherits(raw_df, "list") && !is.null(raw_df[[mean_label]])) {
  raw_df <- raw_df[[mean_label]]
}

base_df <- raw_df %>%
  mutate(
    adopters_prop = num_adopters / N_nodes_actual,
    prop_rational = num_adopted_rational / N_nodes_actual,
    prop_social = num_adopted_social / N_nodes_actual
  )

# Calculate aggregations
agg_df <- base_df %>%
  group_by(innovation_iul_Gamma, social_distance_h) %>%
  summarise(
    avg_steps = mean(num_steps, na.rm = TRUE),
    sd_steps = sd(num_steps, na.rm = TRUE),
    sd_adoption = sd(adopters_prop, na.rm = TRUE),
    avg_adoption = mean(adopters_prop, na.rm = TRUE),
    .groups = "drop"
  )

# --- Correlation Table for User ---
cat("\n=== CORRELATION BETWEEN AVG TOTAL ADOPTION AND AVG STEPS ===\n")
cor_table <- agg_df %>%
  group_by(social_distance_h) %>%
  summarise(
    Correlation = cor(avg_adoption, avg_steps, use = "complete.obs"),
    Avg_Steps_Overall = mean(avg_steps),
    .groups = "drop"
  )
print(as.data.frame(cor_table))
cat("============================================================\n\n")

write.csv(cor_table, paste0(DATA_OUT_DIR, "correlation_adoption_steps.csv"), row.names = FALSE)

# --- Method 2: Critical Slowing Down ---
h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)])

agg_df_plot <- agg_df %>%
  mutate(h_factor = factor(sprintf("%.2f", social_distance_h), levels = h_levels_sorted)) %>%
  filter(innovation_iul_Gamma <= 0.75)

p_steps1 <- ggplot(agg_df_plot, aes(x = innovation_iul_Gamma, y = h_factor, fill = avg_steps)) +
  geom_tile(color = "white", lwd = 0.1) +
  scale_fill_viridis_c(option = "magma", name = "Avg Steps") +
  labs(title = "Method 2: Critical Slowing Down (Avg Steps)", x = "IUL", y = "MSP") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.0, 0.25, 0.5, 0.75), labels = c("0.00", "0.25", "0.50", "0.75")) +
  scale_y_discrete(breaks = y_breaks) +
  theme_minimal()

p_steps2 <- ggplot(agg_df_plot, aes(x = innovation_iul_Gamma, y = h_factor, fill = sd_steps)) +
  geom_tile(color = "white", lwd = 0.1) +
  scale_fill_viridis_c(option = "magma", name = "SD Steps") +
  labs(title = "Method 2: Critical Slowing Down (SD of Steps)", x = "IUL", y = "MSP") +
  scale_x_continuous(expand = c(0, 0), breaks = c(0.0, 0.25, 0.5, 0.75), labels = c("0.00", "0.25", "0.50", "0.75")) +
  scale_y_discrete(breaks = y_breaks) +
  theme_minimal()

ggsave(paste0(PLOTS_DIR, "method2_critical_slowing_down.pdf"), p_steps1 + p_steps2, width = 10, height = 3)

# Find the Global Maximum of Susceptibility (original Method 3 anchoring)
global_max_idx <- which.max(agg_df$sd_adoption)
IUL_max_suscep <- agg_df$innovation_iul_Gamma[global_max_idx]

# Find the Global Maximum Jump of Adoption (Method 4 H=0 effective baseline)
jump_df <- agg_df %>%
  group_by(innovation_iul_Gamma) %>%
  arrange(social_distance_h) %>%
  summarise(max_jump = max(diff(avg_adoption), na.rm = TRUE), .groups = "drop")

IUL_jump <- jump_df$innovation_iul_Gamma[which.max(jump_df$max_jump)]
delta_iul <- 0.025

# --- Method 4: Order Parameter & Susceptibility Exponents for MULTIPLE IULs ---
# Effective Zero-Field Analogue swept around IUL = 0.200
selected_iuls <- c(0.200, 0.225, 0.250, 0.275)

# Snap to available grid values
selected_iuls <- IUL_VALUES_SWEEP[sapply(selected_iuls, function(x) which.min(abs(IUL_VALUES_SWEEP - x)))]
selected_iuls <- unique(selected_iuls)

plot_list <- list()
exponent_results <- data.frame(IUL = numeric(), MSP_c = numeric(), Gamma = numeric())

for (target_iul in selected_iuls) {
  m4_df <- agg_df %>% filter(abs(innovation_iul_Gamma - target_iul) < 1e-4)
  if (nrow(m4_df) == 0) next

  # Find critical point MSP_c where susceptibility (var) peaks
  msp_c <- m4_df$social_distance_h[which.max(m4_df$sd_adoption)]

  m4_plot_df <- m4_df %>%
    mutate(
      dist_to_c = social_distance_h - msp_c,
      abs_dist = abs(dist_to_c)
    )

  p_beta <- ggplot(m4_plot_df, aes(x = dist_to_c, y = avg_adoption)) +
    geom_point(color = "black", size = 1) +
    geom_line(color = "grey50") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = paste0("IUL=", sprintf("%.3f", target_iul), " | MSP_c=", sprintf("%.2f", msp_c)),
      x = "MSP - MSP_c", y = "Avg Adoption (Phi)"
    ) +
    theme_minimal(base_size = 10)

  m4_gamma_df <- m4_plot_df %>% filter(abs_dist > 0, dist_to_c > 0, sd_adoption > 0)

  if (nrow(m4_gamma_df) >= 2) {
    fit <- lm(log10(sd_adoption) ~ log10(abs_dist), data = m4_gamma_df)
    gamma_val <- -coef(fit)[2]
  } else {
    gamma_val <- NA
  }

  exponent_results <- rbind(exponent_results, data.frame(IUL = target_iul, MSP_c = msp_c, Gamma = gamma_val))

  p_gamma <- ggplot(m4_gamma_df, aes(x = abs_dist, y = sd_adoption)) +
    geom_point(color = "purple", size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = paste0("Gamma (", sprintf("%.3f", gamma_val), ") | IUL=", sprintf("%.3f", target_iul)),
      x = "log(|MSP - MSP_c|)", y = "log(Susceptibility)"
    ) +
    theme_minimal(base_size = 10)

  plot_list[[length(plot_list) + 1]] <- p_beta
  plot_list[[length(plot_list) + 1]] <- p_gamma
}

write.csv(exponent_results, paste0(DATA_OUT_DIR, "gamma_exponents.csv"), row.names = FALSE)

final_m4 <- wrap_plots(plot_list, ncol = 2)
ggsave(paste0(PLOTS_DIR, "method4_exponents.pdf"), final_m4, width = 10, height = 2.5 * (length(plot_list) / 2))

# --- Method 3: Avalanche Power Laws 3x3 Array ---
IUL_opt <- IUL_max_suscep
m3_slice <- agg_df %>% filter(abs(innovation_iul_Gamma - IUL_opt) < 1e-4)

# 1. Critical
msp_c_exact <- m3_slice$social_distance_h[which.max(m3_slice$sd_adoption)]

# 2. Super-Critical (Disordered Phase / High MSP ~0.80)
msp_super <- H_VALUES_SWEEP[which.min(abs(H_VALUES_SWEEP - 0.80))]

# 3. Sub-Critical (Ordered Phase / Low MSP ~0.20)
msp_sub <- H_VALUES_SWEEP[which.min(abs(H_VALUES_SWEEP - 0.20))]

selected_msps_m3 <- c(msp_c_exact, msp_super, msp_sub)
labels_m3 <- c("Critical", "Super-Critical", "Sub-Critical")

# Function to generate a row of 3 plots
build_m3_row <- function(msp_val, label_text) {
  df_runs <- base_df %>% filter(
    abs(innovation_iul_Gamma - IUL_opt) < 1e-4,
    abs(social_distance_h - msp_val) < 1e-4
  )

  make_hist <- function(col_name, title_txt, fill_col) {
    ggplot(df_runs, aes(x = .data[[col_name]])) +
      geom_histogram(bins = 100, fill = fill_col, color = NA, alpha = 0.9) +
      scale_x_continuous(limits = c(-5, 1005)) +
      labs(
        title = paste0(title_txt, "\n(", label_text, ")"),
        subtitle = paste0("IUL=", sprintf("%.3f", IUL_opt), ", MSP=", sprintf("%.2f", msp_val)),
        x = "Size", y = "Frequency"
      ) +
      theme_minimal(base_size = 10)
  }

  list(
    make_hist("num_adopters", "Total Adopters", "coral"),
    make_hist("num_adopted_rational", "Rational Adopters", "mediumseagreen"),
    make_hist("num_adopted_social", "Social Adopters", "orchid")
  )
}

all_m3_plots <- c(
  build_m3_row(selected_msps_m3[1], labels_m3[1]),
  build_m3_row(selected_msps_m3[2], labels_m3[2]),
  build_m3_row(selected_msps_m3[3], labels_m3[3])
)

final_m3 <- wrap_plots(all_m3_plots, ncol = 3)
ggsave(paste0(PLOTS_DIR, "method3_avalanches.pdf"), final_m3, width = 12, height = 9)

print("All deeper phase transition analysis plots generated successfully in plots/07_phase_transition/")
