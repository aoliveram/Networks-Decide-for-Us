# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Script to Run Diffusion Simulations (GSS - ER Topology)
#
# This script configures parameters and runs the simulation for ALL
# defined seeding strategies.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load simulation function
source("scripts/04_GSS_diffusion_sims.R")

# -----------------------------------------------------------------------------
# Parameter Configuration
# -----------------------------------------------------------------------------

# 1. Social Influence Thresholds (Normal Distribution)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6) # Means
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08, 0.12, 0.16, 0.20) # Standard Deviations

# 2. Graph Type
CURRENT_GRAPH_TYPE_LABEL <- "ER" # Options: "GSS", "ER"

# 3. Innovation and Social Flexibility Parameters (Internal Sweep)
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025) # Intrinsic Utility
H_VALUES_SWEEP <- seq(0 / 12, 12 / 12, by = 1 / 12) # Social Distance

# 4. Seeding Strategies to Run
strategies_to_run <- c("random", "central", "marginal", "eigen", "closeness")

# 5. Output Directory Configuration
if (CURRENT_GRAPH_TYPE_LABEL == "GSS") {
    RESULTS_DIR <- "output/04_GSS_diffusion_sims/"
} else if (CURRENT_GRAPH_TYPE_LABEL == "ER") {
    RESULTS_DIR <- "output/04_GSS_ER_diffusion_sims/"
} else {
    RESULTS_DIR <- paste0("output/04_", CURRENT_GRAPH_TYPE_LABEL, "_diffusion_sims/")
}

# -----------------------------------------------------------------------------
# Main Loop Execution
# -----------------------------------------------------------------------------

cat("Starting massive simulation execution (GSS - ER)...\n")

for (strategy in strategies_to_run) {
    cat(paste0("\n\n>>> Starting simulations for strategy: ", strategy, " <<<\n"))

    run_diffusion_simulation(
        SEEDING_STRATEGY_FIXED = strategy,
        THRESHOLD_MEAN_SWEEP_LIST = THRESHOLD_MEAN_SWEEP_LIST,
        TAU_NORMAL_SD_SWEEP_LIST = TAU_NORMAL_SD_SWEEP_LIST,
        CURRENT_GRAPH_TYPE_LABEL = CURRENT_GRAPH_TYPE_LABEL,
        IUL_VALUES_SWEEP = IUL_VALUES_SWEEP,
        H_VALUES_SWEEP = H_VALUES_SWEEP,
        NUM_SEED_RUNS_TOTAL = 24, # Adapted to 24
        NUM_CORES_TO_USE = 8, # Adjust according to machine capacity
        RESULTS_DIR = RESULTS_DIR, # Pass the dynamic directory
        NETWORKS_DIR = "data/02_GSS_network_ergm/"
    )
}

cat("\n\nAll simulations have finished.\n")
