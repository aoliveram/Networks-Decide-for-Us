# ------------------------------------------------------------------------------
# Calcula y asigna, a los nodos de las redes simuladas 'GSS', el índice de
# PROPENSIÓN A LA ACCIÓN COLECTIVA y el requisito mínimo de utilidad (MUR) que
# usa el modelo.
#
# Teoría: la regla de adopción es  Gamma + lambda * E_i >= q_i,  donde q_i es un
# REQUISITO MÍNIMO de utilidad: q alto = más difícil de convencer (aversión).
# Los 9 ítems GSS miden lo contrario (propensión: alto = más dispuesto), así que
#   propensity_score = suma de los 9 ítems recodificados (4 - x) / 27   [0, 1]
#   mur_score        = 1 - propensity_score                             [0, 1]
# El motor (scripts/05) lee SOLO 'mur_score'. 'propensity_score' se guarda tal
# como se construye, para que el score siga siendo transparente respecto a los
# datos (corr(mur_score, propensity_score) = -1 exactamente).
#
# Objetivo:
#   1. Leer las redes simuladas existentes en 'data/02_GSS_network_ergm/'.
#   2. Calcular 'propensity_score' (9 ítems) y 'mur_score' = 1 - propensity_score.
#   3. Actualizar las redes existentes añadiendo ambos atributos a los nodos.
#   4. Generar gráficos de diagnóstico.
#
# Entradas:
#   - data/02_GSS_network_ergm/GSS_net_sim_1000_XXX.rds
# Salidas:
#   - data/02_GSS_network_ergm/GSS_net_sim_1000_XXX.rds (Sobreescrito con los atributos)
#   - plots/03_MUR_calculation/GSS_propensity_vars_distribution.pdf
#   (el gráfico de la distribución MUR de GSS y ATP lo hace 03_MUR_distribution_plot.R)
# ------------------------------------------------------------------------------

library(network)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(psych)

# --- Configuración ---
networks_dir <- "data/02_GSS_network_ergm/"
plots_dir    <- "plots/03_MUR_calculation/"
N_networks   <- 100

# Crear directorio de plots si no existe
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Las 9 variables "ingrediente" para nuestro score
propensity_ingredient_vars <- c("signdpet", "avoidbuy", "joindem", "attrally", 
                                "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")

# ==============================================================================
# 1. Diagnóstico y Visualización (Usando la primera red como muestra)
# ==============================================================================

sample_net_path <- file.path(networks_dir, "GSS_net_sim_1000_001.rds")
sample_net <- readRDS(sample_net_path)

# Extraer atributos a un dataframe
df_attr <- data.frame(vertex_id = 1:network.size(sample_net))
for (var in propensity_ingredient_vars) {
  if (var %in% list.vertex.attributes(sample_net)) {
    df_attr[[var]] <- get.vertex.attribute(sample_net, var)
  } else {
    warning(paste("Atributo", var, "no encontrado en la red. Se llenará con NA."))
    df_attr[[var]] <- NA
  }
}

# --- A) Distribución de variables originales ---
plots_propensity <- list()
for (p_var in propensity_ingredient_vars) {
  # Asegurar que la variable es un factor para el gráfico de barras
  # Los valores originales van de 1 a 4
  valid_data <- na.omit(df_attr[[p_var]])
  
  if (length(valid_data) > 0) {
    freq_table <- as.data.frame(table(factor(valid_data, levels = 1:4)))
    colnames(freq_table) <- c("Respuesta", "Frecuencia")
    
    p <- ggplot(freq_table, aes(x = Respuesta, y = Frecuencia, fill = Respuesta)) +
      geom_bar(stat = "identity") +
      scale_x_discrete(drop = FALSE) +
      labs(title = p_var, x = "Response (1-4)", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    
    plots_propensity[[p_var]] <- p
  }
}

pdf(file.path(plots_dir, "GSS_propensity_vars_distribution.pdf"), width = 12, height = 9)
do.call(grid.arrange, c(plots_propensity, ncol = 3))
invisible(dev.off())

# --- B) Cálculo y distribución de propensity_score y mur_score (Muestra) ---
# Codificación original GSS: 1 = "lo hice el último año" ... 4 = "nunca lo haría".
# Recodificación a propensión (alto = más dispuesto): 4 - x
#   1 -> 3, 2 -> 2, 3 -> 1, 4 -> 0
# propensity_score = suma / 27 (max posible = 9 * 3);  mur_score = 1 - propensity_score

df_attr <- df_attr %>%
  mutate(
    across(all_of(propensity_ingredient_vars), 
           ~ 4 - .,
           .names = "recod_{.col}")
  ) %>%
  mutate(
    raw_sum = rowSums(select(., starts_with("recod_"))),
    propensity_score = raw_sum / 27,       # Propensión [0, 1]
    mur_score = 1 - propensity_score       # Requisito mínimo de utilidad (aversión) [0, 1]
  )


# ==============================================================================
# Cronbach α: Internal Consistency of the propensity construct (mur_score = 1 - it)
# ==============================================================================

# Extract recoded values for all 9 items
recoded_items <- df_attr %>%
  select(starts_with("recod_")) %>%
  rename_with(~ gsub("recod_", "", .))

# Remove any rows with missing values for alpha calculation
recoded_items_complete <- recoded_items[complete.cases(recoded_items), ]

# Calculate Cronbach's alpha
cronbach_result <- psych::alpha(recoded_items_complete, warnings = FALSE)$total$raw_alpha
cat("\n========== CRONBACH'S ALPHA INTERNAL CONSISTENCY ==========\n")
cat("Construct: GSS Collective Action Propensity (propensity_score)\n")
cat("Items: signdpet, avoidbuy, joindem, attrally, cntctgov, polfunds, usemedia, interpol, actlaw (9 items)\n")
cat("Cronbach's α =", sprintf("%.4f\n", cronbach_result))
cat("Sample size (complete cases) = ", nrow(recoded_items_complete), "\n")
cat("Interpretation: ")
if (cronbach_result >= 0.70) {
  cat("PASS - Sufficient internal consistency (α ≥ 0.70)\n")
} else {
  cat("WARNING - Low internal consistency (α < 0.70)\n")
}
cat("===========================================================\n\n")

# ==============================================================================
# 2. Procesamiento Masivo: Actualizar Redes
# ==============================================================================

for (i in 1:N_networks) {
  filename <- sprintf("GSS_net_sim_1000_%03d.rds", i)
  full_path <- file.path(networks_dir, filename)
  
  if (!file.exists(full_path)) {
    stop("Archivo no encontrado: ", full_path)
  }
  
  # Cargar red
  net <- readRDS(full_path)
  
  # Extraer matriz de valores originales (1-4), una columna por ítem
  vals_matrix <- matrix(NA, nrow = network.size(net), ncol = length(propensity_ingredient_vars))
  for (k in seq_along(propensity_ingredient_vars)) {
    vals_matrix[, k] <- get.vertex.attribute(net, propensity_ingredient_vars[k])
  }
  
  # Recodificar a propensión: 4 - x  (1 -> 3, 2 -> 2, 3 -> 1, 4 -> 0)
  propensity_vals <- rowSums(4 - vals_matrix) / 27
  
  # Asignar atributos: propensión (construcción) y MUR = 1 - propensión (modelo)
  set.vertex.attribute(net, "propensity_score", propensity_vals)
  set.vertex.attribute(net, "mur_score", 1 - propensity_vals)
  
  # Guardar (Sobreescribir)
  saveRDS(net, full_path)
  
  if (i %% 10 == 0) cat(sprintf("  Procesada red %d/%d\n", i, N_networks))
}

cat("\nProceso completado exitosamente.\n")
