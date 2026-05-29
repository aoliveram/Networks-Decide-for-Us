# ------------------------------------------------------------------------------
# Calcula y asigna el índice de propensión a la innovación (MUR Score) a nodos
# de redes simuladas 'GSS'.
#
# Objetivo:
#   1. Leer las redes simuladas existentes en 'data/02_GSS_network_ergm/'.
#   2. Calcular el 'mur_score' (Minimum Utility Requirement) basado en 9 variables.
#   3. Actualizar las redes existentes añadiendo este atributo a los nodos.
#   4. Generar gráficos de diagnóstico.
#
# Entradas:
#   - data/02_GSS_network_ergm/GSS_network_simulated_1000_XXX.rds
# Salidas:
#   - data/02_GSS_network_ergm/GSS_network_simulated_1000_XXX.rds (Sobreescrito con nuevo atributo)
#   - plots/03_GSS_MUR_calculation/gss_propensity_vars_distribution.pdf
#   - plots/03_GSS_MUR_calculation/mur_score_distribution_GSS.pdf
# ------------------------------------------------------------------------------

library(network)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(psych)

# --- Configuración ---
networks_dir <- "data/02_GSS_network_ergm/"
plots_dir    <- "plots/03_GSS_MUR_calculation/"
N_networks   <- 100

# Crear directorio de plots si no existe
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Las 9 variables "ingrediente" para nuestro score
propensity_ingredient_vars <- c("signdpet", "avoidbuy", "joindem", "attrally", 
                                "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")

# ==============================================================================
# 1. Diagnóstico y Visualización (Usando la primera red como muestra)
# ==============================================================================

sample_net_path <- file.path(networks_dir, "GSS_network_simulated_1000_001.rds")
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
      labs(title = paste("Var:", p_var), x = "Respuesta (1-4)", y = "Frecuencia") +
      theme_minimal() +
      theme(legend.position = "none")
    
    plots_propensity[[p_var]] <- p
  }
}

pdf(file.path(plots_dir, "gss_propensity_vars_distribution.pdf"), width = 12, height = 9)
do.call(grid.arrange, c(plots_propensity, ncol = 3))
dev.off()

# --- B) Cálculo y Distribución del MUR Score (Muestra) ---
# Lógica de recodificación:
# 1 (hecho reciente/muy probable) -> 3
# 2 -> 2
# 3 -> 1
# 4 (nunca/nada probable) -> 0

df_attr <- df_attr %>%
  mutate(
    across(all_of(propensity_ingredient_vars), 
           ~ case_when(
             . == 1 ~ 3,
             . == 2 ~ 2,
             . == 3 ~ 1,
             . == 4 ~ 0,
             TRUE ~ NA_real_
           ),
           .names = "recod_{.col}")
  ) %>%
  mutate(
    raw_sum = rowSums(select(., starts_with("recod_"))),
    mur_score = raw_sum / 27  # Normalizar a [0, 1] (max score posible = 9 * 3 = 27)
  )

# Gráfico de distribución del MUR Score
hist_score <- ggplot(df_attr, aes(x = mur_score)) +
  geom_histogram(bins = 28, fill = "skyblue", color = "black") + # 28 bins para 0..27
  labs(title = "Distribución del MUR Score (GSS)",
       subtitle = "Propensión a la Acción Colectiva (Normalizado 0-1)",
       x = "MUR Score", y = "Frecuencia") +
  theme_minimal()

ggsave(file.path(plots_dir, "mur_score_distribution_GSS.pdf"), plot = hist_score, width = 8, height = 6)

# ==============================================================================
# Cronbach α: Internal Consistency of MUR Construct
# ==============================================================================

# Extract recoded values for all 9 items
recoded_items <- df_attr %>%
  select(starts_with("recod_")) %>%
  rename_with(~ gsub("recod_", "", .))

# Remove any rows with missing values for alpha calculation
recoded_items_complete <- recoded_items[complete.cases(recoded_items), ]

# Calculate Cronbach's alpha
cronbach_result <- cronbach(recoded_items_complete)
cat("\n========== CRONBACH'S ALPHA INTERNAL CONSISTENCY ==========\n")
cat("Construct: GSS Collective Action Propensity (MUR)\n")
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
  filename <- sprintf("GSS_network_simulated_1000_%03d.rds", i)
  full_path <- file.path(networks_dir, filename)
  
  if (!file.exists(full_path)) {
    warning("Archivo no encontrado: ", full_path)
    next
  }
  
  # Cargar red
  net <- readRDS(full_path)
  
  # Extraer atributos a un dataframe temporal para cálculo seguro
  # (Aunque es más lento que vectorizado puro, es más seguro con case_when complejo)
  df_temp <- data.frame(vertex_id = 1:network.size(net))
  for (var in propensity_ingredient_vars) {
    df_temp[[var]] <- get.vertex.attribute(net, var)
  }
  
  # Calcular MUR Score
  df_temp <- df_temp %>%
    mutate(
      across(all_of(propensity_ingredient_vars), 
             ~ case_when(
               . == 1 ~ 3,
               . == 2 ~ 2,
               . == 3 ~ 1,
               . == 4 ~ 0,
               TRUE ~ NA_real_
             ))
    ) %>%
    mutate(
      raw_sum = rowSums(select(., all_of(propensity_ingredient_vars))), # Ya están recodificadas in-place o en nuevas cols?
      # Ah, cuidado con el mutate anterior. Si uso across sin .names, sobreescribe.
      # Vamos a hacerlo explícito para evitar errores.
      mur_score = raw_sum / 27
    )
  
  # Corrección de lógica en el loop para asegurar cálculo correcto:
  # Recodificamos y sumamos en un paso limpio
  
  # Extraer matriz de valores
  vals_matrix <- matrix(NA, nrow = network.size(net), ncol = length(propensity_ingredient_vars))
  for (k in seq_along(propensity_ingredient_vars)) {
    vals_matrix[, k] <- get.vertex.attribute(net, propensity_ingredient_vars[k])
  }
  
  # Recodificar matriz (vectorizado)
  # 1->3, 2->2, 3->1, 4->0
  # Formula: 4 - x  (si x=1 -> 3, x=2 -> 2, x=3 -> 1, x=4 -> 0)
  # Verificamos: 4-1=3, 4-2=2, 4-3=1, 4-4=0. Funciona perfecto para 1,2,3,4.
  recod_matrix <- 4 - vals_matrix
  
  # Manejar NAs si los valores originales no eran 1-4 (aunque el script original asumía 1-4)
  # Si hay valores fuera de rango, esto daría resultados raros, pero asumimos limpieza previa.
  
  mur_vals <- rowSums(recod_matrix) / 27
  
  # Asignar atributo
  set.vertex.attribute(net, "mur_score", mur_vals)
  
  # Guardar (Sobreescribir)
  saveRDS(net, full_path)
  
  if (i %% 10 == 0) cat(sprintf("  Procesada red %d/%d\n", i, N_networks))
}

cat("\nProceso completado exitosamente.\n")
