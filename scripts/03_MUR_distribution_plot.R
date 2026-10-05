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
df <- rbind(
  data.frame(survey = "GSS 2004 (collective action)",
             mur_score = mur("data/02_GSS_network_ergm/GSS_net_sim_1000_001.rds")),
  data.frame(survey = "ATP 2014 (innovation)",
             mur_score = mur("data/02_ATP_network_ergm/ATP_net_sim_1000_001.rds"))
)

p <- ggplot(df, aes(x = mur_score)) +
  geom_histogram(bins = 28, fill = "coral", color = "black") +
  facet_wrap(~ survey, scales = "free_y") +
  labs(x = "MUR score (q)", y = "Count") +
  theme_minimal()

ggsave(file.path(plots_dir, "MUR_distribution.pdf"), plot = p, width = 10, height = 4.5)
