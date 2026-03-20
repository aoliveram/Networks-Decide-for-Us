# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to explore Phase Transition metrics deeply
# 07_phase_transition_deeper.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# --- Configuration ---
RESULTS_DIR <- "output/04_GSS_diffusion_sims/"
PLOTS_DIR <- "plots/07_phase_transition/"
SEEDING_STRATEGY <- "random"
THRESHOLD_MEAN_TARGET <- 0.4
TAU_SD_TARGETS <- c(0.12)

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP <- seq(0, 1, by = 1 / 12)

# Load data
file_path <- paste0(RESULTS_DIR, "results_m03-06_sd", sprintf("%.2f", TAU_SD_TARGETS[1]), "_", SEEDING_STRATEGY, ".rds")
if(!file.exists(file_path)) stop("File not found")

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
    avg_steps = mean(num_steps, na.rm=TRUE),
    max_steps = max(num_steps, na.rm=TRUE),
    sd_adoption = sd(adopters_prop, na.rm=TRUE),
    avg_adoption = mean(adopters_prop, na.rm=TRUE),
    .groups = "drop"
  )

# --- Method 2: Critical Slowing Down ---
h_levels_sorted <- sprintf("%.2f", sort(unique(H_VALUES_SWEEP)))
y_breaks <- sprintf("%.2f", H_VALUES_SWEEP[seq(1, length(H_VALUES_SWEEP), by = 2)])

agg_df_plot <- agg_df %>%
  mutate(h_factor = factor(sprintf("%.2f", social_distance_h), levels = h_levels_sorted)) %>%
  filter(innovation_iul_Gamma <= 0.75)

p_steps1 <- ggplot(agg_df_plot, aes(x = innovation_iul_Gamma, y = h_factor, fill = avg_steps)) +
  geom_tile(color="white", lwd=0.1) +
  scale_fill_viridis_c(option="magma", name="Avg Steps") +
  labs(title="Method 2: Critical Slowing Down (Avg Steps)", x="IUL", y="MSP") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(breaks=y_breaks) + theme_minimal()

p_steps2 <- ggplot(agg_df_plot, aes(x = innovation_iul_Gamma, y = h_factor, fill = max_steps)) +
  geom_tile(color="white", lwd=0.1) +
  scale_fill_viridis_c(option="magma", name="Max Steps") +
  labs(title="Method 2: Critical Slowing Down (Max Steps)", x="IUL", y="MSP") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(breaks=y_breaks) + theme_minimal()

ggsave(paste0(PLOTS_DIR, "method2_critical_slowing_down.pdf"), p_steps1 + p_steps2, width=11, height=4)

# --- Method 4: Order Parameter & Susceptibility Exponents ---
# Pick a fixed MSP, target MSP = 0.5
target_msp <- H_VALUES_SWEEP[which.min(abs(H_VALUES_SWEEP - 0.5))]

m4_df <- agg_df %>% filter(abs(social_distance_h - target_msp) < 1e-4)

# Find critical point IUL_c where susceptibility (var) peaks
iul_c <- m4_df$innovation_iul_Gamma[which.max(m4_df$sd_adoption)]

m4_plot_df <- m4_df %>%
  mutate(
    dist_to_c = innovation_iul_Gamma - iul_c,
    abs_dist = abs(dist_to_c)
  )

# Plot beta (Adoption vs IUL)
p_beta <- ggplot(m4_plot_df, aes(x=innovation_iul_Gamma, y=avg_adoption)) +
  geom_point(color="black") + geom_line(color="grey50") +
  geom_vline(xintercept=iul_c, color="red", linetype="dashed") +
  labs(title=paste0("Method 4: Order Parameter (Phi) [MSP=", sprintf("%.2f", target_msp), "]"),
       subtitle=paste0("IUL_c = ", iul_c, " (red dashed line)"),
       x="IUL", y="Average Total Adoption (Phi)") + theme_minimal()

# Plot gamma (Susceptibility vs |IUL-IULc|) in log-log
m4_gamma_df <- m4_plot_df %>% filter(abs_dist > 0, dist_to_c > 0) # Approached from RHS
p_gamma <- ggplot(m4_gamma_df, aes(x=abs_dist, y=sd_adoption)) +
  geom_point(color="purple") + 
  geom_smooth(method="lm", se=FALSE, color="blue", linewidth=0.5) +
  scale_x_log10() + scale_y_log10() +
  labs(title="Method 4: Susceptibility Exponent", subtitle="Log-Log Scale for IUL > IUL_c",
       x="log(|IUL - IUL_c|)", y="log(Susceptibility)") + theme_minimal()

ggsave(paste0(PLOTS_DIR, "method4_exponents.pdf"), p_beta + p_gamma, width=11, height=4)

# --- Method 3: Avalanche Power Laws ---
# Extract the specific runs at (IUL_c, MSP target)
m3_runs <- base_df %>% 
  filter(abs(social_distance_h - target_msp) < 1e-4, 
         abs(innovation_iul_Gamma - iul_c) < 1e-4)

# Also get a non-critical point for comparison (e.g. IUL = 0.75)
iul_nc <- m4_df$innovation_iul_Gamma[which.min(abs(m4_df$innovation_iul_Gamma - 0.75))]
m3_runs_nc <- base_df %>%
  filter(abs(social_distance_h - target_msp) < 1e-4, 
         abs(innovation_iul_Gamma - iul_nc) < 1e-4)

# Create a combined data frame to plot histograms side by side via patchwork
p_avl1 <- ggplot(m3_runs, aes(x=num_adopted_social)) +
  geom_histogram(bins=15, fill="coral", color="black", alpha=0.7) +
  scale_y_continuous(transform="log1p", breaks=c(0, 1, 5, 10, 20)) +
  labs(title=paste0("Method 3: Avalanche Dist. (Critical)"), 
       subtitle=paste0("IUL=", iul_c, ", MSP=", sprintf("%.2f", target_msp)), 
       x="Social Adopters (Avalanche Size)", y="Log(Frequency)") + theme_minimal()

p_avl2 <- ggplot(m3_runs_nc, aes(x=num_adopted_social)) +
  geom_histogram(bins=15, fill="steelblue", color="black", alpha=0.7) +
  scale_y_continuous(transform="log1p", breaks=c(0, 1, 5, 10, 20)) +
  labs(title=paste0("Avalanche Dist. (Non-Critical)"), 
       subtitle=paste0("IUL=", iul_nc, ", MSP=", sprintf("%.2f", target_msp)), 
       x="Social Adopters (Avalanche Size)", y="Log(Frequency)") + theme_minimal()

ggsave(paste0(PLOTS_DIR, "method3_avalanches.pdf"), p_avl1 + p_avl2, width=11, height=4)

print("All deeper phase transition analysis plots generated successfully in plots/07_phase_transition/")
