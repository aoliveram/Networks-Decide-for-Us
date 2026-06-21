# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Execution Script for GSS Diffusion Plots
#
# Generates heatmaps for both GSS-net and GSS-net-ER topologies.
# Iterates through all seeding strategies.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load the plotting function
source("scripts/05_GSS_diffusion_plot.R")

# Configuration
strategies_to_plot <- c("random", "central", "marginal", "eigen", "closeness")
graph_types <- c("GSS", "ER", "SH")

cat("Starting massive plot generation for GSS datasets...\n")

for (graph_type in graph_types) {
    # Determine paths based on graph type
    if (graph_type == "GSS") {
        RESULTS_DIR <- "output/04_GSS_diffusion_sims/"
        PLOTS_DIR <- "plots/04_GSS_diffusion_sims/"
        NET_LABEL <- "GSS-net"
    } else if (graph_type == "ER") {
        RESULTS_DIR <- "output/04_GSS_ER_diffusion_sims/"
        PLOTS_DIR <- "plots/04_GSS_ER_diffusion_sims/"
        NET_LABEL <- "GSS-net-ER"
    } else if (graph_type == "SH") {
        RESULTS_DIR <- "output/04_GSS_SH_diffusion_sims/"
        PLOTS_DIR <- "plots/04_GSS_SH_diffusion_sims/"
        NET_LABEL <- "GSS-net-SH"
    }

    cat(paste0("\n---------------------------------------------------------------\n"))
    cat(paste0("Processing Graph Type: ", graph_type, " (Label: ", NET_LABEL, ")\n"))
    cat(paste0("---------------------------------------------------------------\n"))

    for (strategy in strategies_to_plot) {
        generate_diffusion_plots(
            RESULTS_DIR = RESULTS_DIR,
            PLOTS_DIR = PLOTS_DIR,
            SEEDING_STRATEGY_FIXED = strategy,
            NETWORK_LABEL = NET_LABEL,
            NUM_CORES = 8 # Adjust if needed
        )
    }
}

cat("\nAll GSS plots generated successfully.\n")
