# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Criticality with the SOCIAL adoption as order parameter
# 04_social_order_parameter.R
#
# Motivation: total adoption Phi = Phi_rational + Phi_social + seeds.
#   - Phi_rational (q_i <= Gamma) is monotone in IUL, FLAT in MSP, topology-blind:
#     a non-critical "background" that dilutes the transition.
#   - Phi_social carries ALL the network physics (selective, MSP-dependent).
# So using Phi_social as the order parameter should sharpen the transition and
# give a cleaner susceptibility exponent.
#
# We replicate EXACTLY the main-text Method-4 extraction (script 07):
#   MSP_c = MSP just before the largest jump in the order parameter;
#   fit log10(SD(order)) ~ log10(|MSP - MSP_c|) on the super-critical side;
#   gamma = -slope.  (Done on SD to match the published gamma_total; we also
#   report the Var-based exponent for completeness, since FDT chi = Var.)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(dplyr); library(ggplot2); library(patchwork); library(viridis)
})

PLOTS_DIR <- "playground/checks-claude/plots/"
DATA_OUT  <- "playground/checks-claude/data/"
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT,  showWarnings = FALSE, recursive = TRUE)

MEAN_KEY <- "mean_0.40"
SELECTED_IULS <- c(0.200, 0.225, 0.250, 0.275)
FILES <- list(
  GSS = "output/04_GSS_diffusion_sims/results_m03-06_sd0.12_random.rds",
  DP  = "output/04_GSS_ER_diffusion_sims/results_m03-06_sd0.12_random.rds"
)

load_agg <- function(path) {
  r <- readRDS(path); if (is.list(r) && !is.data.frame(r)) r <- r[[MEAN_KEY]]
  r %>% mutate(
      p_total = num_adopters / N_nodes_actual,
      p_social = num_adopted_social / N_nodes_actual,
      p_rational = num_adopted_rational / N_nodes_actual
    ) %>%
    group_by(IUL = innovation_iul_Gamma, MSP = social_distance_h) %>%
    summarise(
      avg_total = mean(p_total), sd_total = sd(p_total), var_total = var(p_total),
      avg_social = mean(p_social), sd_social = sd(p_social), var_social = var(p_social),
      avg_rational = mean(p_rational),
      .groups = "drop")
}

# Method-4 gamma for a chosen order parameter -------------------------------
# order = "total" or "social"; disp = "sd" or "var"
extract_gamma <- function(agg, target_iul, order = "social", disp = "sd") {
  avg_col <- paste0("avg_", order)
  d_col   <- paste0(disp, "_", order)
  sl <- agg %>% filter(abs(IUL - target_iul) < 1e-4) %>% arrange(MSP)
  if (nrow(sl) < 4) return(list(msp_c = NA, gamma = NA, r2 = NA, n = 0))
  jump_idx <- which.max(diff(sl[[avg_col]]))
  msp_c <- sl$MSP[jump_idx]
  up <- sl %>% mutate(dist = MSP - msp_c, disp_v = pmax(.data[[d_col]], 1e-6)) %>%
    filter(dist > 0)
  if (nrow(up) < 2) return(list(msp_c = msp_c, gamma = NA, r2 = NA, n = nrow(up)))
  fit <- lm(log10(disp_v) ~ log10(dist), data = up)
  list(msp_c = msp_c, gamma = -coef(fit)[2], r2 = summary(fit)$r.squared, n = nrow(up))
}

rows <- list(); aggs <- list()
for (topo in names(FILES)) {
  agg <- load_agg(FILES[[topo]]); aggs[[topo]] <- agg
  for (iul in SELECTED_IULS) {
    gs  <- extract_gamma(agg, iul, "social", "sd")
    gsv <- extract_gamma(agg, iul, "social", "var")
    gt  <- extract_gamma(agg, iul, "total",  "sd")
    rows[[paste(topo, iul)]] <- data.frame(
      topology = topo, IUL = iul,
      MSP_c_social = gs$msp_c,
      gamma_social_SD = gs$gamma, r2_social_SD = gs$r2,
      gamma_social_Var = gsv$gamma,
      gamma_total_SD = gt$gamma)
  }
}
res <- bind_rows(rows)

cat("\n================ SOCIAL vs TOTAL ORDER PARAMETER (gamma) ================\n")
print(res, row.names = FALSE, digits = 3)
cat("\nMean-field prediction: gamma = 1 (SD-based, as in the published Method 4).\n")
summ <- res %>% group_by(topology) %>%
  summarise(gamma_social_SD = mean(gamma_social_SD, na.rm = TRUE),
            gamma_social_Var = mean(gamma_social_Var, na.rm = TRUE),
            gamma_total_SD = mean(gamma_total_SD, na.rm = TRUE), .groups = "drop")
cat("\n--- means across IUL ---\n"); print(as.data.frame(summ), digits = 3)
cat("========================================================================\n")
write.csv(res, file.path(DATA_OUT, "gamma_social_vs_total.csv"), row.names = FALSE)

# Heatmaps: social order parameter + its susceptibility ---------------------
make_heat <- function(agg, topo) {
  d <- agg %>% filter(IUL <= 0.75)
  h_lv <- sprintf("%.2f", sort(unique(agg$MSP)))
  d <- d %>% mutate(hf = factor(sprintf("%.2f", MSP), levels = h_lv))
  ybr <- h_lv[seq(1, length(h_lv), by = 2)]
  p1 <- ggplot(d, aes(IUL, hf, fill = avg_social)) + geom_tile() +
    scale_fill_viridis_c(name = expression(Phi[soc])) +
    scale_y_discrete(breaks = ybr) +
    labs(title = paste0(topo, ": Social adoption"), x = "IUL", y = "MSP") +
    theme_minimal(base_size = 9)
  p2 <- ggplot(d, aes(IUL, hf, fill = sd_social)) + geom_tile() +
    scale_fill_viridis_c(option = "magma", name = "SD") +
    scale_y_discrete(breaks = ybr) +
    labs(title = paste0(topo, ": Social susceptibility"), x = "IUL", y = "MSP") +
    theme_minimal(base_size = 9)
  p1 + p2
}
heat <- (make_heat(aggs$GSS, "GSS")) / (make_heat(aggs$DP, "GSS-DP"))
ggsave(file.path(PLOTS_DIR, "social_order_heatmaps.pdf"), heat, width = 10, height = 6)

# Exponent log-log panels (social), GSS vs DP ------------------------------
panels <- list()
for (topo in names(FILES)) {
  agg <- aggs[[topo]]
  for (iul in SELECTED_IULS) {
    g <- extract_gamma(agg, iul, "social", "sd")
    sl <- agg %>% filter(abs(IUL - iul) < 1e-4) %>%
      mutate(dist = MSP - g$msp_c, disp_v = pmax(sd_social, 1e-6)) %>% filter(dist > 0)
    panels[[length(panels) + 1]] <- ggplot(sl, aes(dist, disp_v)) +
      geom_point(color = "orchid", size = 1) +
      geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 0.4) +
      scale_x_log10() + scale_y_log10() +
      labs(title = paste0(topo, " IUL=", sprintf("%.3f", iul),
                          " | g=", ifelse(is.na(g$gamma), "NA", sprintf("%.2f", g$gamma))),
           x = "log|MSP-MSP_c|", y = "log SD(Phi_soc)") +
      theme_minimal(base_size = 8)
  }
}
ggsave(file.path(PLOTS_DIR, "gamma_social_exponents.pdf"),
       wrap_plots(panels, ncol = 4), width = 12, height = 5)

cat("\nSaved: gamma_social_vs_total.csv, social_order_heatmaps.pdf, gamma_social_exponents.pdf\n")
