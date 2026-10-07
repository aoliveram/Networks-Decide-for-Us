# ------------------------------------------------------------------------------
# Distribución del MUR score (q = 1 - propensity_score) en GSS y ATP, en un solo
# gráfico. Requiere haber corrido antes 03_GSS_MUR_calculation.R y
# 03_ATP_MUR_calculation.R (lee la red 001 de cada carpeta).
#
# Salida: plots/03_MUR_calculation/MUR_distribution.pdf
# ------------------------------------------------------------------------------

library(network)
library(ggplot2)

plots_dir <- "plots/03_MUR_calculation/"
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

mur <- function(path) get.vertex.attribute(readRDS(path), "mur_score")

# The scores are discrete (GSS: 28 levels = k/27; ATP: 7 levels = k/6), so we draw
# one bar per level, with width proportional to the spacing between levels
# (a histogram with fixed bins gives very thin bars for ATP).
level_counts <- function(x, survey, n_levels) {
  step <- 1 / (n_levels - 1)
  tab  <- table(round(x / step))
  data.frame(survey = survey, mur_score = as.numeric(names(tab)) * step,
             n = as.vector(tab), width = 0.9 * step)
}
df <- rbind(
  level_counts(mur("data/02_GSS_network_ergm/GSS_net_sim_1000_001.rds"), "GSS", 28),
  level_counts(mur("data/02_ATP_network_ergm/ATP_net_sim_1000_001.rds"), "ATP", 7)
)
df$survey <- factor(df$survey, levels = c("GSS", "ATP"))

p <- ggplot(df, aes(x = mur_score, y = n, fill = survey, width = width)) +
  geom_col(color = "black") +
  facet_wrap(~ survey, scales = "free_y") +
  scale_fill_manual(values = c(GSS = "skyblue", ATP = "coral")) +
  labs(x = "MUR score (q)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, "MUR_distribution.pdf"), plot = p, width = 10, height = 4.5)
