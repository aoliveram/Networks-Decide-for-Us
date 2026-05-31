# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regenerate network diagnostic plots (English, LaTeX math labels, professional)
# 02b_network_plots_EN.R
#
# Self-contained replacement for the plotting blocks of
# 02_GSS_network_ergm_tests.R / 02_ATP_network_ergm_tests.R, which are stale
# (they expect 'GSS_network_simulated_1000_*.rds'; the actual files are
# 'GSS_net_sim_1000_*.rds') and carry heavy unrelated dependencies.
#
# Produces, for GSS and ATP:
#   - <PREFIX>_degree_dist.pdf : degree distribution (linear + log-log),
#     mean over 100 imputed nets with +/-1 SD band, vs ER and BA theory.
#     Legend: network identifiers (BA / ER / GSS|ATP) first, the SD band below.
#   - <PREFIX>_Ck_vs_k.pdf : local clustering C(k) vs k (median, [p25,p75]).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(network)
  library(igraph)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# colours (consistent with the original deck)
col_ER  <- "#2ECC71"   # green
col_BA  <- "#6A1B9A"   # purple
col_EMP <- "#E74C3C"   # red (empirical: GSS / ATP)

theme_pro <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
        legend.key.height = unit(11, "pt"),
        plot.title = element_text(face = "bold"))

# ---------------------------------------------------------------------------
# load one folder of statnet 'network' objects -> list of igraph graphs
# ---------------------------------------------------------------------------
load_graphs <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  lapply(files, function(f) {
    net <- readRDS(f)
    igraph::graph_from_adjacency_matrix(as.matrix.network(net, matrix.type = "adjacency"),
                                        mode = "undirected", diag = FALSE)
  })
}

# ---------------------------------------------------------------------------
# degree-distribution statistics across nets
# ---------------------------------------------------------------------------
degree_stats <- function(graphs) {
  pk_list <- lapply(seq_along(graphs), function(i) {
    d <- igraph::degree(graphs[[i]])
    tibble(id = i, k = d) %>% count(k) %>% mutate(pk = n / sum(n))
  })
  bind_rows(pk_list) %>% group_by(k) %>%
    summarise(mean_pk = mean(pk), sd_pk = sd(pk), .groups = "drop")
}

clustering_stats <- function(graphs) {
  ck_list <- lapply(graphs, function(g) {
    tibble(k = igraph::degree(g), Ck = igraph::transitivity(g, type = "local"))
  })
  bind_rows(ck_list) %>% filter(!is.na(Ck), k >= 1) %>%
    group_by(k) %>%
    summarise(Ck_med = median(Ck), Ck_p25 = quantile(Ck, .25),
              Ck_p75 = quantile(Ck, .75), n = n(), .groups = "drop") %>%
    filter(n >= 5)
}

# ---------------------------------------------------------------------------
# build the two plots for one population
# ---------------------------------------------------------------------------
make_plots <- function(graphs, label, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  emp_lab <- paste0(label, " (mean)")
  sd_lab  <- paste0(label, " (±1 SD)")

  N0     <- igraph::vcount(graphs[[1]])
  k_mean <- mean(sapply(graphs, function(g) mean(igraph::degree(g))))
  p_er   <- k_mean / (N0 - 1)

  deg <- degree_stats(graphs)

  # ER theory: Binomial(N-1, p)
  k_er <- 0:(N0 - 1)
  er_t <- tibble(k = k_er, pk = dbinom(k_er, size = N0 - 1, prob = p_er)) %>%
    filter(pk > 1e-12)
  # BA theory: P(k) = 2 m (m+1) / [k (k+1)(k+2)]
  m_ba <- max(1L, round(k_mean / 2))
  k_ba <- m_ba:(N0 - 1)
  ba_t <- tibble(k = k_ba, pk = (2 * m_ba * (m_ba + 1)) / (k_ba * (k_ba + 1) * (k_ba + 2)))

  col_map  <- setNames(c(col_BA, col_ER, col_EMP), c("BA (theory)", "ER (theory)", emp_lab))
  fill_map <- setNames(col_EMP, sd_lab)
  brk_col  <- c("BA (theory)", "ER (theory)", emp_lab)

  base_layers <- function() list(
    geom_ribbon(data = deg, aes(x = k, ymin = pmax(mean_pk - sd_pk, 0),
                                ymax = mean_pk + sd_pk, fill = sd_lab),
                alpha = 0.45, inherit.aes = FALSE),
    geom_line(data = ba_t, aes(x = k, y = pk, color = "BA (theory)"), linewidth = 0.9),
    geom_line(data = er_t, aes(x = k, y = pk, color = "ER (theory)"), linewidth = 0.9),
    geom_line(data = deg, aes(x = k, y = mean_pk, color = emp_lab), linewidth = 1.1),
    scale_color_manual(values = col_map, breaks = brk_col, name = NULL),
    scale_fill_manual(values = fill_map, name = NULL),
    guides(color = guide_legend(order = 1), fill = guide_legend(order = 2))
  )

  p_lin <- ggplot() + base_layers() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, max(deg$mean_pk) * 1.25)) +
    labs(title = "(a) Linear scale", x = expression(italic(k)),
         y = expression(italic(p)[italic(k)])) + theme_pro

  p_log <- ggplot() + base_layers() +
    scale_x_log10() + scale_y_log10() +
    coord_cartesian(xlim = c(5, 150), ylim = c(1e-5, 1)) +
    labs(title = "(b) Log-log scale",
         x = expression(log[10]~italic(k)),
         y = expression(log[10]~italic(p)[italic(k)])) + theme_pro

  ggsave(file.path(out_dir, paste0(label, "_degree_dist.pdf")),
         p_lin | p_log, width = 11, height = 4.8, dpi = 300)

  # ---- C(k) vs k ----
  ck <- clustering_stats(graphs)
  C_ER <- k_mean / N0
  p_ck <- ggplot(ck, aes(x = k, y = Ck_med)) +
    geom_ribbon(aes(ymin = Ck_p25, ymax = Ck_p75), fill = col_EMP, alpha = 0.4) +
    geom_line(color = col_EMP, linewidth = 1.1) +
    geom_hline(yintercept = C_ER, linetype = "dashed", color = col_ER, linewidth = 1) +
    annotate("text", x = max(ck$k) * 0.62, y = C_ER, vjust = -0.6, color = col_ER,
             label = deparse(bquote(italic(C)[ER] %~~% .(round(C_ER, 3)))), parse = TRUE) +
    labs(title = paste0("Local clustering: ", label, " (median, [p25, p75])"),
         x = expression(italic(k)), y = expression(italic(C)(italic(k)))) + theme_pro
  ggsave(file.path(out_dir, paste0(label, "_Ck_vs_k.pdf")),
         p_ck, width = 6.6, height = 4.2, dpi = 300)

  cat(sprintf("[%s] N=%d  <k>=%.2f  p=%.4f  m_BA=%d  -> 2 PDFs written\n",
              label, N0, k_mean, p_er, m_ba))
}

# ---------------------------------------------------------------------------
cat("Loading GSS networks...\n")
gss <- load_graphs("data/02_GSS_network_ergm", "^GSS_net_sim_1000_\\d{3}\\.rds$")
make_plots(gss, "GSS", "plots/02_GSS_network_ergm_tests")

cat("Loading ATP networks...\n")
atp <- load_graphs("data/02_ATP_network_ergm", "^ATP_net_sim_1000_\\d{3}\\.rds$")
make_plots(atp, "ATP", "plots/02_ATP_network_ergm_tests")

cat("\nDone.\n")
