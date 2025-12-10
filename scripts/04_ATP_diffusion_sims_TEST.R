# TEST VERSION of 04_ATP_diffusion_sims.R

library(igraph)
library(doParallel)
library(dplyr)
library(cluster)
library(netdiffuseR)

# -----------------------------------------------------------------------------
# 1. Core Parameters & Setup (TEST CONFIG)
# -----------------------------------------------------------------------------
cat("Setting up core parameters for TEST...\n")

NUM_SEED_RUNS_TOTAL <- 2 
N_NODES_GLOBAL <- 50 # Small graph for testing

# Reduced Sweeps
THRESHOLD_MEAN_SWEEP_LIST <- c(0.4)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.12)
IUL_VALUES_SWEEP <- c(0.5)
H_VALUES_SWEEP <- c(0.5)

NUM_CORES_TO_USE <- 1 

RESULTS_DIR <- "output/test_results/"
if(!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Helper to generate dummy graph
# -----------------------------------------------------------------------------
generate_dummy_graph <- function(n) {
  g <- erdos.renyi.game(n, p = 0.1)
  V(g)$name <- 1:n
  V(g)$age <- sample(18:80, n, replace = TRUE)
  V(g)$educ_num <- sample(1:5, n, replace = TRUE)
  V(g)$race <- sample(c("White", "Black", "Hispanic", "Other"), n, replace = TRUE)
  V(g)$relig <- sample(c("Prot", "Cath", "None", "Other"), n, replace = TRUE)
  V(g)$sex <- sample(c("Male", "Female"), n, replace = TRUE)
  V(g)$q_i <- runif(n, 0, 1) # Intrinsic utility
  return(g)
}

# -----------------------------------------------------------------------------
# 3. Main Simulation Loop
# -----------------------------------------------------------------------------
cat("Starting TEST simulation...\n")

SEEDING_STRATEGY_FIXED <- "random"

for (current_tau_sd in TAU_NORMAL_SD_SWEEP_LIST) {
  for (current_threshold_mean in THRESHOLD_MEAN_SWEEP_LIST) {
    
    cat(paste0("Running for Mean=", current_threshold_mean, ", SD=", current_tau_sd, "\n"))
    
    # Sequential loop for testing
    results_list_all_runs <- list()
    
    for (run_idx in 1:NUM_SEED_RUNS_TOTAL) {
      cat(paste0("  Run ", run_idx, "/", NUM_SEED_RUNS_TOTAL, "\n"))
      
      # Generate dummy graph instead of loading
      graph_for_this_run <- generate_dummy_graph(N_NODES_GLOBAL)
      
      N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
      node_mur_q_specific <- V(graph_for_this_run)$q_i
      node_degrees_specific <- igraph::degree(graph_for_this_run)
      
      attributes_for_distance_specific <- data.frame(
        age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
        race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
        sex = as.factor(V(graph_for_this_run)$sex)
      )
      d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
      
      # Thresholds
      set.seed(run_idx * 1000) 
      node_thresholds_tau_frac_specific <- rnorm(
        n = N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd
      )
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
      node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
      
      # Seeding
      primary_seed_for_this_run <- sample(1:N_NODES_SPECIFIC_GRAPH, 1)
      initial_infectors_for_this_sim_run <- c(primary_seed_for_this_run)
      
      # --- REFACTORING: Use netdiffuseR ---
      adj_mat <- as_adjacency_matrix(graph_for_this_run, sparse = FALSE)
      
      for (current_Gamma in IUL_VALUES_SWEEP) {
        # Calculate effective thresholds
        effective_thresholds <- node_thresholds_tau_frac_specific
        rational_indices <- which(node_mur_q_specific <= current_Gamma)
        effective_thresholds[rational_indices] <- 1e-6 # Epsilon
        
        for (current_h in H_VALUES_SWEEP) {
           # Calculate probability matrix W
           W <- 1 / (1 + exp((d_ij_matrix - current_h) / 0.02))
           W <- W * adj_mat 
           
           # Run rdiffnet
           # Run rdiffnet
           diff_model <- rdiffnet(
             n = N_NODES_SPECIFIC_GRAPH,
             seed.nodes = initial_infectors_for_this_sim_run,
             threshold.dist = effective_thresholds,
             seed.graph = W, 
             exposure.mode = "stochastic",
             t = 20 
           )
           
           adopters <- which(!is.na(diff_model$toa))
           non_seed_adopters <- setdiff(adopters, initial_infectors_for_this_sim_run)
           
           num_adopted_rational <- sum(non_seed_adopters %in% rational_indices)
           num_adopted_social <- sum(!(non_seed_adopters %in% rational_indices))
           
           res_row <- data.frame(
             run_id = run_idx,
             innovation_iul_Gamma = current_Gamma,
             social_distance_h = current_h,
             num_adopters = length(adopters),
             num_adopted_rational = num_adopted_rational,
             num_adopted_social = num_adopted_social
           )
           results_list_all_runs[[length(results_list_all_runs) + 1]] <- res_row
        }
      }
    }
    
    final_df <- bind_rows(results_list_all_runs)
    print(final_df)
    saveRDS(final_df, paste0(RESULTS_DIR, "test_results.rds"))
  }
}

cat("Test script completed successfully.\n")
