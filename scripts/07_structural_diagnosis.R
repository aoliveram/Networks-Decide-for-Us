# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structural diagnosis: what does "breaking homophily" (GSS -> GSS-DP) destroy?
# 08_structural_diagnosis.R   (GSS vs DP: what does randomization destroy?)
#
# For the 100 imputed GSS plausible networks and 100 degree-preserving
# randomized (DP) counterparts (one DP per plausible, via the SAME double-edge
# swap used in the diffusion sims: igraph::rewire(keeping_degseq)), compute a
# panel of network statistics in three families:
#   A. Controls the DP preserves by construction  (<k>, density)  -> sanity check
#   B. Higher-order topology (the candidate topological "culprit")
#        clustering(global), clustering(avg local), #triangles, degree
#        assortativity, mean path length, modularity (Louvain), giant component
#   C. Homophily (the sociological "culprit")
#        attribute assortativity per Blau axis (race, sex, relig, age, educ)
#        + MEAN GOWER DISTANCE OVER EDGES (the model's own d_ij)
#
# Significance: 100 vs 100 -> Mann-Whitney U (two-sided) + Cohen's d effect size.
# Outputs:
#   - data/structural_diagnosis_summary.csv  (per-stat: means, p, d, % change)
#   - plots/structural_diagnosis_B.pdf  (3x3 grid, family B)
#   - plots/structural_diagnosis_C.pdf  (1x2: per-axis assortativity + Gower dist)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(network); library(igraph); library(cluster)   # cluster::daisy for Gower
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

set.seed(20260619)
NET_DIR  <- "data/02_GSS_network_ergm"
PAT      <- "^GSS_net_sim_1000_\\d{3}\\.rds$"
N_NETS   <- 100
SWAP_FAC <- 20                 # same niter_factor as the diffusion sims
PLOTS    <- "plots/07_structural_diagnosis"
DATA_OUT <- "output/07_structural_diagnosis"
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)

BLAU_CAT <- c("race", "sex", "relig")    # categorical Blau axes
BLAU_NUM <- c("age", "educ_num")          # continuous Blau axes
GOWER_VARS <- c("age", "educ_num", "race", "relig", "sex")  # = model's d_ij inputs

# ---------------------------------------------------------------------------
# load a statnet network -> igraph, carrying the demographic attributes
# ---------------------------------------------------------------------------
load_igraph <- function(path) {
  net <- readRDS(path)
  A <- network::as.matrix.network(net, matrix.type = "adjacency")
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", diag = FALSE)
  for (a in c(BLAU_CAT, BLAU_NUM, "mur_score")) {
    if (a %in% network::list.vertex.attributes(net))
      g <- igraph::set_vertex_attr(g, a, value = network::get.vertex.attribute(net, a))
  }
  g
}

# degree-preserving randomization (identical to 04_GSS_diffusion_sims.R)
make_dp <- function(g) {
  igraph::rewire(g, with = igraph::keeping_degseq(niter = SWAP_FAC * igraph::ecount(g)))
}

# mixed-type assortativity per attribute: nominal (Newman) for categorical,
# numeric (Pearson over edges) for continuous
assort_attr <- function(g, attr, categorical) {
  v <- igraph::vertex_attr(g, attr)
  if (categorical) {
    igraph::assortativity_nominal(g, types = as.integer(factor(v)) , directed = FALSE)
  } else {
    igraph::assortativity(g, values = as.numeric(v), directed = FALSE)
  }
}

# mean Gower distance over the existing edges (the model's d_ij on ties)
mean_edge_gower <- function(g) {
  df <- data.frame(
    age      = as.numeric(igraph::vertex_attr(g, "age")),
    educ_num = as.numeric(igraph::vertex_attr(g, "educ_num")),
    race     = factor(igraph::vertex_attr(g, "race")),
    relig    = factor(igraph::vertex_attr(g, "relig")),
    sex      = factor(igraph::vertex_attr(g, "sex"))
  )
  D <- as.matrix(cluster::daisy(df, metric = "gower"))
  el <- igraph::as_edgelist(g, names = FALSE)
  mean(D[el], na.rm = TRUE)
}

# all statistics for one graph -> named numeric vector
stats_one <- function(g) {
  comps <- igraph::components(g)
  out <- c(
    # A. controls
    mean_degree        = mean(igraph::degree(g)),
    density            = igraph::edge_density(g),
    # B. higher-order topology
    clustering_global  = igraph::transitivity(g, type = "global"),
    clustering_local   = igraph::transitivity(g, type = "localaverage", isolates = "zero"),
    triangles          = sum(igraph::count_triangles(g)) / 3,
    assort_degree      = igraph::assortativity_degree(g, directed = FALSE),
    mean_path_length   = igraph::mean_distance(g, directed = FALSE),
    modularity         = igraph::modularity(igraph::cluster_louvain(g)),
    giant_frac         = max(comps$csize) / igraph::vcount(g),
    # C. homophily
    assort_race        = assort_attr(g, "race",  TRUE),
    assort_sex         = assort_attr(g, "sex",   TRUE),
    assort_relig       = assort_attr(g, "relig", TRUE),
    assort_age         = assort_attr(g, "age",      FALSE),
    assort_educ        = assort_attr(g, "educ_num", FALSE),
    mean_edge_gower    = mean_edge_gower(g)
  )
  out
}

# ---------------------------------------------------------------------------
# run over the 100 networks
# ---------------------------------------------------------------------------
files <- list.files(NET_DIR, pattern = PAT, full.names = TRUE)[1:N_NETS]
cat("Computing statistics on", length(files), "GSS networks and their DP counterparts...\n")

# GSS-SH networks (topology byte-identical, attributes permuted) — already built
# and saved by 04b_make_SH_networks.R; load them so the panel can place an
# SH column BETWEEN GSS and GSS-DP.
SH_DIR <- "data/02_GSS_SH_network"
SH_PAT <- "^GSS_sh_sim_1000_(\\d{3})\\.rds$"
sh_files <- list.files(SH_DIR, pattern = SH_PAT, full.names = TRUE)

rows <- list()
for (i in seq_along(files)) {
  g  <- load_igraph(files[i])
  gd <- make_dp(g)
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = "GSS",    t(stats_one(g)))
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = "GSS-SH", t(stats_one(load_igraph(sh_files[i]))))
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = "GSS-DP", t(stats_one(gd)))
  if (i %% 10 == 0) cat("  done", i, "/", length(files), "\n")
}
allstats <- bind_rows(rows)
# fixed display order: GSS (real) -> GSS-SH (homophily gone) -> GSS-DP (topology gone)
allstats$topo <- factor(allstats$topo, levels = c("GSS", "GSS-SH", "GSS-DP"))
write.csv(allstats, file.path(DATA_OUT, "structural_diagnosis_raw.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# significance: Mann-Whitney U + Cohen's d, per statistic
# ---------------------------------------------------------------------------
stat_cols <- setdiff(names(allstats), c("net", "topo"))
cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx-1)*var(x) + (ny-1)*var(y)) / (nx+ny-2))
  if (sp == 0 || is.na(sp)) return(NA_real_)
  (mean(x) - mean(y)) / sp
}
summ <- lapply(stat_cols, function(s) {
  x  <- allstats[[s]][allstats$topo == "GSS"]
  sh <- allstats[[s]][allstats$topo == "GSS-SH"]
  y  <- allstats[[s]][allstats$topo == "GSS-DP"]
  p     <- tryCatch(wilcox.test(x, y)$p.value,  error = function(e) NA_real_)  # GSS vs DP
  p_sh  <- tryCatch(wilcox.test(x, sh)$p.value, error = function(e) NA_real_)  # GSS vs SH
  data.frame(statistic = s,
             mean_GSS = mean(x, na.rm = TRUE), mean_SH = mean(sh, na.rm = TRUE),
             mean_DP = mean(y, na.rm = TRUE),
             pct_change_DP = 100 * (mean(y, na.rm = TRUE)  - mean(x, na.rm = TRUE)) / abs(mean(x, na.rm = TRUE)),
             pct_change_SH = 100 * (mean(sh, na.rm = TRUE) - mean(x, na.rm = TRUE)) / abs(mean(x, na.rm = TRUE)),
             cohens_d = cohens_d(x, y), p_value = p, p_value_SH = p_sh)
}) %>% bind_rows()
write.csv(summ, file.path(DATA_OUT, "structural_diagnosis_summary.csv"), row.names = FALSE)
cat("\n=============== STRUCTURAL DIAGNOSIS SUMMARY ===============\n")
print(summ, row.names = FALSE, digits = 3)
cat("===========================================================\n")

# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------
# GSS-DP now uses the deck's dark green (ndgreen = RGB 0,118,54) to match the
# "Two ways to break the network" concept map.
COLS <- c("GSS" = "#E74C3C", "GSS-SH" = "#2980B9", "GSS-DP" = "#007636")
TOPO_LV <- c("GSS", "GSS-SH", "GSS-DP")

# single boxplot, 3 topologies, labels above each box (no legend, no p-value)
box_one <- function(stat, title) {
  d <- allstats[, c("topo", stat)]; names(d)[2] <- "y"
  ytop <- max(d$y, na.rm = TRUE); yr <- diff(range(d$y, na.rm = TRUE))
  labs_df <- data.frame(topo = factor(TOPO_LV, levels = TOPO_LV), y = ytop + 0.12 * yr)
  ggplot(d, aes(topo, y, fill = topo)) +
    geom_boxplot(width = 0.6, outlier.size = 0.35, alpha = 0.85) +
    geom_text(data = labs_df, aes(label = topo, color = topo), vjust = 0, size = 2.5,
              fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = COLS, guide = "none") +
    scale_color_manual(values = COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

# --- Plot 1: family B, 3x2 grid (drop the control row: giant comp / degree / density) ---
B_stats  <- c("clustering_global","clustering_local","triangles",
              "assort_degree","mean_path_length","modularity")
B_titles <- c("Global clustering","Avg local clustering","Triangles",
              "Degree assortativity","Mean path length","Modularity")
pB <- wrap_plots(Map(box_one, B_stats, B_titles), ncol = 3)
ggsave(file.path(PLOTS, "structural_diagnosis_B.pdf"), pB, width = 9, height = 5.6)

# --- Plot 2: family C, 1x2 (3 topologies). Left panel: dashed separators between
#     Blau axes, per-box GSS/GSS-SH/GSS-DP labels (no legend), compact. ---
AXIS_ORDER <- c("Age", "Education", "Race", "Religion", "Sex")
assort_long <- allstats %>%
  select(topo, assort_race, assort_sex, assort_relig, assort_age, assort_educ) %>%
  pivot_longer(-topo, names_to = "axis", values_to = "val") %>%
  mutate(axis = factor(recode(axis, assort_race="Race", assort_sex="Sex", assort_relig="Religion",
                       assort_age="Age", assort_educ="Education"), levels = AXIS_ORDER))
DODGE <- 0.7
ytopC <- max(assort_long$val, na.rm = TRUE)
# per-box vertical labels just above each of the 3 boxes in every axis group
labC <- do.call(rbind, lapply(seq_along(AXIS_ORDER), function(k) data.frame(
  axis = factor(AXIS_ORDER[k], levels = AXIS_ORDER),
  topo = factor(TOPO_LV, levels = TOPO_LV),
  x    = k + c(-DODGE/3, 0, DODGE/3),
  val  = ytopC + 0.03)))
# dashed vertical separators between axis groups (at the half-integer positions)
seps <- data.frame(x = seq_len(length(AXIS_ORDER) - 1) + 0.5)
pC1 <- ggplot(assort_long, aes(axis, val, fill = topo)) +
  geom_vline(data = seps, aes(xintercept = x), linetype = "dashed",
             color = "grey70", linewidth = 0.3) +
  geom_boxplot(width = 0.62, outlier.size = 0.25, alpha = 0.85,
               position = position_dodge(DODGE)) +
  geom_text(data = labC, aes(x = x, label = topo, color = topo), vjust = 0, hjust = 0,
            angle = 90, size = 1.7, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_color_manual(values = COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
  labs(title = "Attribute assortativity by Blau axis", x = NULL, y = "Assortativity") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
pC2 <- box_one("mean_edge_gower", "Mean social distance")
pC <- pC1 + pC2 + plot_layout(widths = c(4, 1))
# smaller overall canvas so the fonts read larger relative to the plot
ggsave(file.path(PLOTS, "structural_diagnosis_C.pdf"), pC, width = 6.6, height = 2.7)

cat("\nSaved: structural_diagnosis_summary.csv, _B.pdf (3 topos, 3x2), _C.pdf (3 topos, 1x2)\n")
