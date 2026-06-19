# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main Script — Diffusion Simulations on GSS-SH (attribute-SHuffle null model)
#
# GSS-SH = the plausible GSS networks with their TOPOLOGY KEPT BYTE-IDENTICAL
# but node attributes randomly PERMUTED, so ONLY homophily is destroyed
# (clustering / triangles / modularity / degree all intact). This is the clean
# counterfactual the GSS-DP (degree-preserving rewiring) baseline cannot give:
# if the structural premium survives DP but vanishes here, the cause is homophily,
# not clustering.
#
# Runs the SAME sweep as 04_GSS_diffusion_sims_main.R, but with the "SH" topology,
# writing results to a SEPARATE directory so nothing collides with GSS or GSS-DP.
#
#   Run it with a single command:
#       Rscript scripts/04_GSS_SH_diffusion_sims_main.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load simulation function (now supports CURRENT_GRAPH_TYPE_LABEL == "SH")
source("scripts/04_GSS_diffusion_sims.R")

# -----------------------------------------------------------------------------
# Parameter Configuration — identical sweep to the GSS / GSS-DP runs
# -----------------------------------------------------------------------------
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)
TAU_NORMAL_SD_SWEEP_LIST  <- c(0.08, 0.12, 0.16, 0.20)

CURRENT_GRAPH_TYPE_LABEL  <- "SH"                 # <-- attribute-shuffle topology

IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)
H_VALUES_SWEEP   <- seq(0 / 12, 12 / 12, by = 1 / 12)

strategies_to_run <- c("random", "central", "marginal", "eigen", "closeness")

# SEPARATE output directory — does not touch GSS or GSS-DP
RESULTS_DIR <- "output/04_GSS_SH_diffusion_sims/"

# Source networks: the ORIGINAL plausible GSS networks. The "SH" branch in the
# simulation permutes their attributes per run (topology untouched), so we read
# from the same folder; nothing is overwritten.
NETWORKS_DIR <- "data/02_GSS_network_ergm/"

# -----------------------------------------------------------------------------
# Main Loop Execution
# -----------------------------------------------------------------------------
cat("Starting massive simulation execution (GSS-SH / attribute shuffle)...\n")
cat("Results -> ", RESULTS_DIR, "\n")

for (strategy in strategies_to_run) {
  cat(paste0("\n\n>>> Starting simulations for strategy: ", strategy, " <<<\n"))
  run_diffusion_simulation(
    SEEDING_STRATEGY_FIXED    = strategy,
    THRESHOLD_MEAN_SWEEP_LIST = THRESHOLD_MEAN_SWEEP_LIST,
    TAU_NORMAL_SD_SWEEP_LIST  = TAU_NORMAL_SD_SWEEP_LIST,
    CURRENT_GRAPH_TYPE_LABEL  = CURRENT_GRAPH_TYPE_LABEL,
    IUL_VALUES_SWEEP          = IUL_VALUES_SWEEP,
    H_VALUES_SWEEP            = H_VALUES_SWEEP,
    NUM_SEED_RUNS_TOTAL       = 24,
    NUM_CORES_TO_USE          = 8,        # adjust to the machine
    RESULTS_DIR               = RESULTS_DIR,
    NETWORKS_DIR              = NETWORKS_DIR
  )
}

cat("\n\nAll GSS-SH simulations have finished. Results in ", RESULTS_DIR, "\n")
