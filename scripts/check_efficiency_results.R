
library(dplyr)
library(readr)

results_file <- "output/04_ATP_diffusion_sims_efficiency/efficiency_test_results_sd0.08.rds"

if (!file.exists(results_file)) {
  stop("El archivo de resultados no existe.")
}

df <- readRDS(results_file)

cat("=== Estructura de los datos ===\n")
glimpse(df)

cat("\n=== Resumen de num_adopters ===\n")
summary(df$num_adopters)

cat("\n=== Adopción promedio por H (Social Distance) ===\n")
df %>%
  group_by(social_distance_h) %>%
  summarise(
    mean_adopters = mean(num_adopters),
    min_adopters = min(num_adopters),
    max_adopters = max(num_adopters),
    n = n()
  ) %>%
  print(n = 20)

cat("\n=== Adopción promedio por Gamma (IUL) - Muestra ===\n")
# Muestreamos algunos valores de Gamma para no imprimir los 41
gammas_sample <- unique(df$innovation_iul_Gamma)[seq(1, length(unique(df$innovation_iul_Gamma)), length.out = 10)]

df %>%
  filter(innovation_iul_Gamma %in% gammas_sample) %>%
  group_by(innovation_iul_Gamma) %>%
  summarise(
    mean_adopters = mean(num_adopters),
    min_adopters = min(num_adopters),
    max_adopters = max(num_adopters),
    n = n()
  ) %>%
  print()

cat("\n=== Distribución de pasos (num_steps) ===\n")
summary(df$num_steps)

cat("\n=== Verificación de Racional vs Social ===\n")
df %>%
  summarise(
    total_rational = sum(num_adopted_rational),
    total_social = sum(num_adopted_social),
    ratio_social_rational = sum(num_adopted_social) / sum(num_adopted_rational)
  ) %>%
  print()
