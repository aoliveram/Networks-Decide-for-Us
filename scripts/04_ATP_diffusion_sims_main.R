# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script Principal para Ejecutar Simulaciones de Difusión (ATP)
#
# Este script configura los parámetros y ejecuta la simulación para TODAS
# las estrategias de seeding definidas.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cargar la función de simulación
source("scripts/04_ATP_diffusion_sims.R")

# -----------------------------------------------------------------------------
# Configuración de Parámetros
# -----------------------------------------------------------------------------

# 1. Umbrales de Influencia Social (Distribución Normal)
THRESHOLD_MEAN_SWEEP_LIST <- c(0.3, 0.4, 0.5, 0.6)     # Medias
TAU_NORMAL_SD_SWEEP_LIST  <- c(0.08, 0.12, 0.16, 0.20) # Desviaciones Estándar

# 2. Tipo de Grafo
CURRENT_GRAPH_TYPE_LABEL <- "ATP" # Opciones: "ATP", "ER"

# 3. Parámetros de Innovación y Flexibilidad Social (Barrido Interno)
IUL_VALUES_SWEEP <- seq(0.0, 1.0, by = 0.025)   # Utilidad Intrínseca
H_VALUES_SWEEP   <- seq(0/12, 12/12, by = 1/12) # Distancia Social

# 4. Estrategias de Seeding a Ejecutar
# Nota: "random" y "central" ya se ejecutaron.
strategies_to_run <- c("marginal", "eigen", "closeness")

# -----------------------------------------------------------------------------
# Ejecución del Bucle Principal
# -----------------------------------------------------------------------------

cat("Iniciando ejecución masiva de simulaciones...\n")

for (strategy in strategies_to_run) {
  
  cat(paste0("\n\n>>> Iniciando simulaciones para estrategia: ", strategy, " <<<\n"))
  
  run_diffusion_simulation(
    SEEDING_STRATEGY_FIXED = strategy,
    THRESHOLD_MEAN_SWEEP_LIST = THRESHOLD_MEAN_SWEEP_LIST,
    TAU_NORMAL_SD_SWEEP_LIST = TAU_NORMAL_SD_SWEEP_LIST,
    CURRENT_GRAPH_TYPE_LABEL = CURRENT_GRAPH_TYPE_LABEL,
    IUL_VALUES_SWEEP = IUL_VALUES_SWEEP,
    H_VALUES_SWEEP = H_VALUES_SWEEP,
    NUM_SEED_RUNS_TOTAL = 24, # Ajustado a 24 según solicitud
    NUM_CORES_TO_USE = 8      # Ajustar según capacidad de la máquina
  )
  
}

cat("\n\nTodas las simulaciones han finalizado.\n")
