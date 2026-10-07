# ------------------------------------------------------------------------------
# Calcula y asigna, a los nodos de las redes simuladas 'ATP', el índice de
# PROPENSIÓN A LA INNOVACIÓN y el requisito mínimo de utilidad (MUR) que usa el
# modelo.
#
# Teoría: la regla de adopción es  Gamma + lambda * E_i >= q_i,  donde q_i es un
# REQUISITO MÍNIMO de utilidad: q alto = más difícil de convencer (aversión).
# Los 6 ítems METECH miden lo contrario (propensión: alto = más dispuesto), así
# que
#   propensity_score = promedio de los 6 ítems (anti-innovación invertidos) [0, 1]
#   mur_score        = 1 - propensity_score                                 [0, 1]
# El motor (scripts/05) lee SOLO 'mur_score'. 'propensity_score' se guarda tal
# como se construye, para que el score siga siendo transparente respecto a los
# datos (corr(mur_score, propensity_score) = -1 exactamente).
#
# Objetivo:
#   1. Leer las redes simuladas existentes en 'data/02_ATP_network_ergm/'.
#   2. Calcular 'propensity_score' (6 variables METECH) y 'mur_score' = 1 - propensity_score.
#   3. Actualizar las redes existentes añadiendo ambos atributos a los nodos.
#   4. Generar gráficos de diagnóstico.
#
# Entradas:
#   - data/02_ATP_network_ergm/ATP_net_sim_1000_XXX.rds
# Salidas:
#   - data/02_ATP_network_ergm/ATP_net_sim_1000_XXX.rds (Sobreescrito con los atributos)
#   - plots/03_MUR_calculation/ATP_metech_vars_distribution.pdf
#   (el gráfico de la distribución MUR de GSS y ATP lo hace 03_MUR_distribution_plot.R)
# ------------------------------------------------------------------------------

library(network)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(psych)

# --- Configuración ---
networks_dir <- "data/02_ATP_network_ergm/"
plots_dir    <- "plots/03_MUR_calculation/"
N_networks   <- 100

# Crear directorio de plots si no existe
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

metech_vars <- c("metech_a", "metech_b", "metech_c", 
                 "metech_d", "metech_e", "metech_f")

# ==============================================================================
# 1. Diagnóstico y Visualización (Usando la primera red como muestra)
# ==============================================================================

sample_net_path <- file.path(networks_dir, "ATP_net_sim_1000_001.rds")
sample_net <- readRDS(sample_net_path)

# Extraer atributos a un dataframe
df_attr <- data.frame(vertex_id = 1:network.size(sample_net))
for (var in metech_vars) {
  if (var %in% list.vertex.attributes(sample_net)) {
    df_attr[[var]] <- get.vertex.attribute(sample_net, var)
  } else {
    warning(paste("Atributo", var, "no encontrado en la red. Se llenará con NA."))
    df_attr[[var]] <- NA
  }
}

# --- A) Distribución de variables originales METECH ---
plots_metech <- list()
for (col_name in metech_vars) {
  valid_data <- na.omit(df_attr[[col_name]])
  
  if (length(valid_data) > 0) {
    freq_table <- as.data.frame(table(factor(valid_data, levels = c(0, 1))))
    colnames(freq_table) <- c("Respuesta", "Frecuencia")
    
    p <- ggplot(freq_table, aes(x = Respuesta, y = Frecuencia, fill = Respuesta)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("0" = "skyblue", "1" = "coral")) +
      labs(title = col_name, x = "(0 = No, 1 = Yes)", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
    
    plots_metech[[col_name]] <- p
  }
}

pdf(file.path(plots_dir, "ATP_metech_vars_distribution.pdf"), width = 10, height = 7)
do.call(grid.arrange, c(plots_metech, ncol = 3))
invisible(dev.off())

# --- B) Cálculo y distribución de propensity_score y mur_score (Muestra) ---
# Lógica de cálculo ("marque todas las que apliquen", 1 = aplica):
# a. Usually try new products before others do (METECH_A): 1 = pro-innovación
# b. Prefer my tried and trusted brands (METECH_B): 1 = ANTI-innovación (invertir)
# c. Like being able to tell others about new brands (METECH_C): 1 = pro-innovación
# d. Like the variety of trying new products (METECH_D): 1 = pro-innovación
# e. Feel more comfortable using familiar brands (METECH_E): 1 = ANTI-innovación (invertir)
# f. Wait until I hear about others' experiences (METECH_F): 1 = ANTI-innovación (invertir)
# propensity_score = promedio de los 6 ítems;  mur_score = 1 - propensity_score

df_attr <- df_attr %>%
  mutate(
    score_A = ifelse(is.na(metech_a), NA, metech_a),       # Pro
    score_B = ifelse(is.na(metech_b), NA, 1 - metech_b),   # Anti -> Invertir
    score_C = ifelse(is.na(metech_c), NA, metech_c),       # Pro
    score_D = ifelse(is.na(metech_d), NA, metech_d),       # Pro
    score_E = ifelse(is.na(metech_e), NA, 1 - metech_e),   # Anti -> Invertir
    score_F = ifelse(is.na(metech_f), NA, 1 - metech_f)    # Anti -> Invertir
  ) %>%
  mutate(
    raw_sum = rowSums(select(., starts_with("score_"))),
    propensity_score = raw_sum / 6,        # Propensión [0, 1]
    mur_score = 1 - propensity_score       # Requisito mínimo de utilidad (aversión) [0, 1]
  )


# ==============================================================================
# Cronbach α: Internal Consistency of the propensity construct (mur_score = 1 - it)
# ==============================================================================

# Items recodificados en la misma dirección (anti-innovación invertidos); sobre
# los binarios crudos, con ítems en sentidos opuestos, el alpha no tiene sentido.
binary_items <- df_attr %>%
  select(starts_with("score_"))

# Remove any rows with missing values for alpha calculation
binary_items_complete <- binary_items[complete.cases(binary_items), ]

# Calculate Cronbach's alpha
# Note: For binary items, alpha is still valid and useful
cronbach_result_atp <- psych::alpha(binary_items_complete, warnings = FALSE)$total$raw_alpha
cat("\n========== CRONBACH'S ALPHA INTERNAL CONSISTENCY ==========\n")
cat("Construct: ATP Innovation Propensity (propensity_score)\n")
cat("Items: metech_a, metech_b, metech_c, metech_d, metech_e, metech_f (6 binary items)\n")
cat("Cronbach's α =", sprintf("%.4f\n", cronbach_result_atp))
cat("Sample size (complete cases) = ", nrow(binary_items_complete), "\n")
cat("Interpretation: ")
if (cronbach_result_atp >= 0.70) {
  cat("PASS - Sufficient internal consistency (α ≥ 0.70)\n")
} else if (cronbach_result_atp >= 0.60) {
  cat("ACCEPTABLE - Marginal internal consistency (0.60 ≤ α < 0.70)\n")
} else {
  cat("WARNING - Low internal consistency (α < 0.60)\n")
}
cat("===========================================================\n\n")

# ==============================================================================
# 2. Procesamiento Masivo: Actualizar Redes
# ==============================================================================

for (i in 1:N_networks) {
  filename <- sprintf("ATP_net_sim_1000_%03d.rds", i)
  full_path <- file.path(networks_dir, filename)
  
  if (!file.exists(full_path)) {
    stop("Archivo no encontrado: ", full_path)
  }
  
  # Cargar red
  net <- readRDS(full_path)
  
  # Extraer atributos necesarios
  # Nota: Podríamos usar el df_attr calculado arriba si los nodos son idénticos en orden y atributos
  # en todas las redes simuladas (lo cual es usual en ERGM si no cambian los nodos).
  # Sin embargo, para ser seguros, recalculamos por red o extraemos de la red actual.
  
  # Extracción rápida vectorizada
  vals_a <- get.vertex.attribute(net, "metech_a")
  vals_b <- get.vertex.attribute(net, "metech_b")
  vals_c <- get.vertex.attribute(net, "metech_c")
  vals_d <- get.vertex.attribute(net, "metech_d")
  vals_e <- get.vertex.attribute(net, "metech_e")
  vals_f <- get.vertex.attribute(net, "metech_f")
  
  # Cálculo vectorizado (NA se propaga: NA + 1 = NA)
  s_a <- vals_a
  s_b <- 1 - vals_b
  s_c <- vals_c
  s_d <- vals_d
  s_e <- 1 - vals_e
  s_f <- 1 - vals_f
  
  propensity_vals <- (s_a + s_b + s_c + s_d + s_e + s_f) / 6
  
  # Asignar atributos: propensión (construcción) y MUR = 1 - propensión (modelo)
  set.vertex.attribute(net, "propensity_score", propensity_vals)
  set.vertex.attribute(net, "mur_score", 1 - propensity_vals)
  
  # Guardar (Sobreescribir)
  saveRDS(net, full_path)
  
  if (i %% 10 == 0) cat(sprintf("  Procesada red %d/%d\n", i, N_networks))
}

cat("\nProceso completado exitosamente.\n")
