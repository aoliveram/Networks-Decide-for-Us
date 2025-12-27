# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT DE EFICIENCIA: SIMULACIÓN DE DIFUSIÓN (ATP)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Objetivo: Correr simulaciones de difusión de manera eficiente usando FORK y matrices dispersas.
# Parámetros reducidos para prueba: 50 redes, 50 corridas, 8 cores.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(igraph)
library(doParallel)
library(dplyr)
library(readr)
library(intergraph)
library(cluster)
library(netdiffuseR)
library(Matrix)

# Evitar conflictos de OpenMP con FORK
Sys.setenv(OMP_NUM_THREADS="1")

# -----------------------------------------------------------------------------
# 1. Parámetros
# -----------------------------------------------------------------------------
NETWORKS_DIR <- "data/02_ATP_network_ergm/"
RESULTS_DIR <- "output/04_ATP_diffusion_sims_efficiency/"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# Configuración de la prueba
NUM_NETWORK_INSTANCES_AVAILABLE <- 50
NUM_SEED_RUNS_TOTAL <- 50
NUM_CORES_TO_USE <- 8
N_NODES_GLOBAL <- 1000
MAX_TIME_STEPS <- 200 

# Parámetros del modelo (Solo 1 configuración de Mean/SD para la prueba)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3)
TAU_NORMAL_SD_SWEEP_LIST <- c(0.08)

# Barrido completo interno
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025) # 41 valores
H_VALUES_SWEEP <- seq(0/12, 12/12, by = 1/12) # 13 valores

# Estrategia de semilla
SEEDING_STRATEGY_FIXED <- "central"
CURRENT_GRAPH_TYPE_LABEL <- "ATP"

# -----------------------------------------------------------------------------
# 2. Funciones Auxiliares
# -----------------------------------------------------------------------------
# (Necesarias para exportar a los workers si se usara PSOCK, pero útiles en general)

graph_density_target <- function(g) {
  igraph::ecount(g) / (igraph::vcount(g) * (igraph::vcount(g) - 1) / 2)
}

random_er_same_density <- function(g_base) {
  n <- igraph::vcount(g_base)
  p <- graph_density_target(g_base)
  gr <- igraph::erdos.renyi.game(n = n, p.or.m = p, type = "gnp", directed = FALSE, loops = FALSE)
  if (!is.null(igraph::V(g_base)$name)) igraph::V(gr)$name <- igraph::V(g_base)$name
  for (attr in igraph::vertex_attr_names(g_base)) {
    igraph::vertex_attr(gr, attr) <- igraph::vertex_attr(g_base, attr)
  }
  gr
}

random_degree_preserving <- function(g_base, niter_factor = 20) {
  m <- igraph::ecount(g_base)
  gr <- igraph::rewire(g_base, with = igraph::keeping_degseq(niter = niter_factor * m))
  gr
}

# Función Core de Simulación (para usar tanto en Benchmark como en Producción)
run_single_simulation_task <- function(run_idx, network_file_idx, current_threshold_mean, current_tau_sd) {
  
  # 1. Cargar Red
  current_network_path <- paste0(NETWORKS_DIR, "ATP_net_sim_1000_", sprintf("%03d", network_file_idx), ".rds")
  if (!file.exists(current_network_path)) return(NULL)
  
  graph_for_this_run_ergm <- readRDS(current_network_path)
  graph_for_this_run <- asIgraph(graph_for_this_run_ergm)
  
  # (Omitimos lógica de cambio de topología para simplificar la prueba ATP, pero se puede añadir)
  N_NODES_SPECIFIC_GRAPH <- vcount(graph_for_this_run)
  
  # 2. Configurar Atributos y Semillas
  set.seed(run_idx * 3000 + round(current_threshold_mean * 100) + round(current_tau_sd * 100))
  node_mur_q_specific <- V(graph_for_this_run)$mur_score
  node_degrees_specific <- igraph::degree(graph_for_this_run)
  
  attributes_for_distance_specific <- data.frame(
    age = V(graph_for_this_run)$age, educ_num = V(graph_for_this_run)$educ_num,
    race = as.factor(V(graph_for_this_run)$race), relig = as.factor(V(graph_for_this_run)$relig),
    sex = as.factor(V(graph_for_this_run)$sex)
  )
  d_ij_matrix <- as.matrix(daisy(attributes_for_distance_specific, metric = "gower"))
  
  # Thresholds
  set.seed(run_idx * 1000 + round(current_threshold_mean * 100) + round(current_tau_sd * 1000)) 
  node_thresholds_tau_frac_specific <- rnorm(N_NODES_SPECIFIC_GRAPH, mean = current_threshold_mean, sd = current_tau_sd)
  node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific < 0] <- 0
  node_thresholds_tau_frac_specific[node_thresholds_tau_frac_specific > 1] <- 1
  
  # Seeding Logic (Simplified for 'central')
  primary_seed_for_this_run <- which.max(igraph::degree(graph_for_this_run))
  if(length(primary_seed_for_this_run) > 1) primary_seed_for_this_run <- primary_seed_for_this_run[1]
  
  # Cluster seeding logic (simplified)
  node_thresholds_count <- round(node_thresholds_tau_frac_specific * node_degrees_specific)
  node_thresholds_count[node_thresholds_count <= 0 & node_thresholds_tau_frac_specific > 0] <- 1
  node_thresholds_count[node_thresholds_count <= 0 & node_thresholds_tau_frac_specific == 0] <- 0
  
  num_seeds <- node_thresholds_count[primary_seed_for_this_run]
  num_seeds <- min(num_seeds, N_NODES_SPECIFIC_GRAPH, node_degrees_specific[primary_seed_for_this_run] + 1)
  if (num_seeds < 1) num_seeds <- 1
  
  initial_infectors <- c(primary_seed_for_this_run)
  if (num_seeds > 1) {
    neighbors_of_primary <- as.numeric(neighbors(graph_for_this_run, primary_seed_for_this_run, mode="total"))
    needed <- num_seeds - 1
    if (length(neighbors_of_primary) >= needed) {
      initial_infectors <- c(initial_infectors, sample(neighbors_of_primary, needed))
    } else {
      initial_infectors <- c(initial_infectors, neighbors_of_primary)
    }
  }
  initial_infectors <- unique(initial_infectors)
  
  # 3. Simulación (Bucles H y IUL)
  results_list <- list()
  
  # Pre-cálculo disperso
  el <- as_edgelist(graph_for_this_run, names = FALSE)
  if (nrow(el) > 0) {
    edge_dists <- d_ij_matrix[el]
  } else {
    edge_dists <- numeric(0)
  }
  
  total_h <- length(H_VALUES_SWEEP)
  num_iul <- length(IUL_VALUES_SWEEP)
  
  # Preparar listas para rdiffnet multi-behavior (vectorización de IUL)
  # seed.p.adopt debe ser una lista para activar el modo multi-behavior
  seed_p_adopt_list <- as.list(rep(0, num_iul)) 
  # seed.nodes debe ser una lista de vectores numéricos (mismos seeds para todos los IUL)
  seed_nodes_list <- replicate(num_iul, initial_infectors, simplify = FALSE)
  
  for (h_idx in 1:total_h) {
    current_h <- H_VALUES_SWEEP[h_idx]
    
    # Feedback de progreso (solo visible si no se captura stdout)
    # cat(sprintf("\r    [Run %d] Procesando H=%.2f (%d/%d)...", run_idx, current_h, h_idx, total_h))
    
    # Cálculo W disperso
    if (length(edge_dists) > 0) {
      edge_weights <- 1 / (1 + exp((edge_dists - current_h) / 0.02))
      W_sparse <- Matrix::sparseMatrix(
        i = c(el[,1], el[,2]), 
        j = c(el[,2], el[,1]), 
        x = c(edge_weights, edge_weights),
        dims = c(N_NODES_SPECIFIC_GRAPH, N_NODES_SPECIFIC_GRAPH)
      )
    } else {
      W_sparse <- Matrix::sparseMatrix(i=integer(0), j=integer(0), dims=c(N_NODES_SPECIFIC_GRAPH, N_NODES_SPECIFIC_GRAPH))
    }
    
    # Construir Matriz de Thresholds para todos los valores de IUL simultáneamente
    # Filas: Nodos, Columnas: Comportamientos (Valores de IUL)
    threshold_matrix <- matrix(0, nrow = N_NODES_SPECIFIC_GRAPH, ncol = num_iul)
    
    for (i in 1:num_iul) {
      current_Gamma <- IUL_VALUES_SWEEP[i]
      effective_thresholds <- node_thresholds_tau_frac_specific
      rational_indices <- which(node_mur_q_specific <= current_Gamma)
      effective_thresholds[rational_indices] <- 1e-6 
      threshold_matrix[, i] <- effective_thresholds
    }
    
    # Ejecutar rdiffnet una sola vez para todos los IUL (Multi-behavior)
    # Esto reduce drásticamente el overhead de llamadas a función
    diff_model <- rdiffnet(
      n = N_NODES_SPECIFIC_GRAPH,
      seed.nodes = seed_nodes_list,
      seed.p.adopt = seed_p_adopt_list, # Activa multi-behavior
      threshold.dist = threshold_matrix,
      seed.graph = W_sparse, 
      exposure.mode = "stochastic",
      exposure.args = list(valued = TRUE),
      t = MAX_TIME_STEPS, # Usar límite razonable en lugar de N_NODES
      stop.no.diff = FALSE # No detener si un IUL no difunde
    )
    
    # Extraer resultados para cada IUL
    # diff_model$toa es ahora una matriz (N x num_iul)
    
    for (i in 1:num_iul) {
      current_Gamma <- IUL_VALUES_SWEEP[i]
      
      # Obtener TOA para este comportamiento (columna i)
      toa_col <- diff_model$toa[, i]
      adopters <- which(!is.na(toa_col))
      
      # Recalcular conteos racional/social (necesitamos saber quiénes eran racionales para este Gamma)
      rational_indices <- which(node_mur_q_specific <= current_Gamma)
      non_seed_adopters <- setdiff(adopters, initial_infectors)
      
      # Manejo de Inf/-Inf en max() si no hay adoptantes
      num_steps_val <- if(length(adopters) > 0) max(toa_col, na.rm = TRUE) else 0
      
      res_row <- data.frame(
        innovation_iul_Gamma = current_Gamma,
        social_distance_h = current_h,
        seed = primary_seed_for_this_run,
        num_adopters = length(adopters),
        num_steps = num_steps_val,
        num_adopted_rational = sum(non_seed_adopters %in% rational_indices),
        num_adopted_social = sum(!(non_seed_adopters %in% rational_indices)),
        initial_cluster_size = length(initial_infectors)
      )
      results_list[[length(results_list) + 1]] <- res_row
    }
  }
  
  return(bind_rows(results_list))
}

# -----------------------------------------------------------------------------
# 3. Estimación de Tiempo (Benchmark)
# -----------------------------------------------------------------------------
cat("\n====================================================================\n")
cat("  FASE 1: ESTIMACIÓN DE TIEMPO (Benchmark)\n")
cat("====================================================================\n")
cat("Ejecutando 1 corrida completa (todas las H, todas las IUL) secuencialmente...\n")

t_bench_start <- Sys.time()
# Corremos la simulación 1 (red 1)
bench_res <- run_single_simulation_task(
  run_idx = 1, 
  network_file_idx = 1, 
  current_threshold_mean = THRESHOLD_MEAN_SWEEP_LIST[1], 
  current_tau_sd = TAU_NORMAL_SD_SWEEP_LIST[1]
)
t_bench_end <- Sys.time()

duration_bench_secs <- as.numeric(difftime(t_bench_end, t_bench_start, units = "secs"))
cat(paste0("Tiempo de 1 corrida completa: ", round(duration_bench_secs, 2), " segundos.\n"))

# Cálculo de estimación
# Total corridas = NUM_SEED_RUNS_TOTAL (50)
# Cores = NUM_CORES_TO_USE (8)
# Lotes paralelos = ceiling(50 / 8) = 7 lotes
# Tiempo total estimado = duración_1_corrida * lotes

total_batches <- ceiling(NUM_SEED_RUNS_TOTAL / NUM_CORES_TO_USE)
estimated_total_secs <- duration_bench_secs * total_batches
estimated_total_mins <- estimated_total_secs / 60

cat(paste0("\n--- ESTIMACIÓN PARA ", NUM_SEED_RUNS_TOTAL, " CORRIDAS EN ", NUM_CORES_TO_USE, " CORES ---\n"))
cat(paste0("Lotes secuenciales necesarios: ", total_batches, "\n"))
cat(paste0("TIEMPO TOTAL ESTIMADO: ", round(estimated_total_mins, 2), " MINUTOS.\n"))
cat("====================================================================\n\n")

# -----------------------------------------------------------------------------
# 4. Ejecución Principal
# -----------------------------------------------------------------------------
cat("  FASE 2: EJECUCIÓN PARALELA\n")
cat("====================================================================\n")

# Usamos FORK como solicitado
cl <- makeCluster(NUM_CORES_TO_USE, type = "FORK")
registerDoParallel(cl)
cat(paste0("Cluster FORK iniciado con ", NUM_CORES_TO_USE, " cores.\n"))

time_start_parallel <- Sys.time()

# Barra de progreso "manual" imprimiendo desde el master loop
# Nota: foreach con %dopar% no devuelve control hasta terminar, 
# pero podemos usar .verbose=TRUE para debug o simplemente confiar en la estimación.
# Para tener una barra real se requerirían paquetes extra o callbacks complejos.
# Imprimiremos al inicio y final.

cat("Iniciando simulaciones...\n")

final_results_list <- foreach(
  run_idx = 1:NUM_SEED_RUNS_TOTAL,
  .combine = 'list',
  .multicombine = TRUE,
  .packages = c('igraph', 'dplyr', 'readr', 'intergraph', 'cluster', 'netdiffuseR', 'Matrix'),
  .errorhandling = 'pass'
) %dopar% {
  
  # Cálculo del índice de red
  network_file_idx <- ((run_idx - 1) %% NUM_NETWORK_INSTANCES_AVAILABLE) + 1
  
  # Ejecutar tarea
  res <- run_single_simulation_task(
    run_idx = run_idx,
    network_file_idx = network_file_idx,
    current_threshold_mean = THRESHOLD_MEAN_SWEEP_LIST[1],
    current_tau_sd = TAU_NORMAL_SD_SWEEP_LIST[1]
  )
  
  # Añadir metadatos
  if (!is.null(res)) {
    res$run_id <- run_idx
    res$network_instance_file_idx <- network_file_idx
    res$threshold_mean_param <- THRESHOLD_MEAN_SWEEP_LIST[1]
    res$threshold_sd_param <- TAU_NORMAL_SD_SWEEP_LIST[1]
  }
  
  return(res)
}

stopCluster(cl)

time_end_parallel <- Sys.time()
actual_duration_mins <- as.numeric(difftime(time_end_parallel, time_start_parallel, units = "mins"))

cat(paste0("\nSimulaciones terminadas.\n"))
cat(paste0("Tiempo real: ", round(actual_duration_mins, 2), " minutos.\n"))

# -----------------------------------------------------------------------------
# 5. Guardado
# -----------------------------------------------------------------------------
cat("Guardando resultados...\n")
valid_results <- final_results_list[!sapply(final_results_list, is.null)]
if (length(valid_results) > 0) {
  combined_df <- bind_rows(valid_results)
  filename <- paste0(RESULTS_DIR, "efficiency_test_results_sd", TAU_NORMAL_SD_SWEEP_LIST[1], ".rds")
  saveRDS(combined_df, filename)
  cat(paste0("Archivo guardado: ", filename, "\n"))
} else {
  cat("ADVERTENCIA: No se generaron resultados válidos.\n")
}

cat("Script finalizado.\n")
