# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structural diagnosis: what does "breaking the network" destroy?
# 07_structural_diagnosis.R
#
# For the 24 networks of a survey (02_network_calibrated.R) compares three
# versions of each network:
#   x2  the main networks (calibrated + rewired to twice the triangles)
#   x1  the same nodes, degrees and mixing, before the rewiring: the extra
#       closure is gone, the homophily stays
#   DP  the degree-preserving twin of x2, built exactly as in 05
#       (set.seed(500000 + k), rewire(keeping_degseq, niter = 20 m)): the
#       homophily is gone too
# (Until 2026-11 the middle column was the SH network, attributes shuffled on
# the same topology; SH is no longer used.)
# Panel of statistics in three families:
#   A. Controls the DP preserves by construction  (<k>, density)  -> sanity check
#   B. Higher-order topology: clustering (global, average local), triangles,
#      degree assortativity, mean path length, modularity (Louvain), giant component
#   C. Homophily: attribute assortativity per Blau axis (race, sex, relig, age,
#      educ; ELSOC: sex, relig, ideol, age, educ; GSS14 adds ideol, which is not in
#      its tie model) + the MEAN GOWER DISTANCE OVER EDGES (the model's own d_ij)
#
# Significance: 24 vs 24 -> Mann-Whitney U (two-sided) + Cohen's d (x2 vs DP),
# and Mann-Whitney x2 vs x1.
# Outputs (<S> = GSS, ATP, GSS14 or ELSOC):
#   output/07_structural_diagnosis/<S>/structural_diagnosis_{raw,summary}.csv
#   plots/07_structural_diagnosis/<S>/structural_diagnosis_B.pdf  (3x2 grid, family B)
#   plots/07_structural_diagnosis/<S>/structural_diagnosis_C.pdf  (assortativity + Gower distance)
# Usage: Rscript scripts/07_structural_diagnosis.R --survey=GSS     (~1 min)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressWarnings(suppressMessages({
  library(network); library(igraph); library(cluster)   # cluster::daisy for Gower
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
}))

ARGS <- commandArgs(trailingOnly = TRUE)
SURVEY <- toupper(sub("^--survey=", "", grep("^--survey=", ARGS, value = TRUE)))
if (length(SURVEY) == 0) SURVEY <- "GSS"
stopifnot(SURVEY %in% c("GSS", "ATP", "GSS14", "ELSOC"))
NET_DIR  <- "data/02_calibrated_networks"
SWAP_FAC <- 20                 # same niter factor as the DP twins of 05
PLOTS    <- file.path("plots/07_structural_diagnosis", SURVEY)
DATA_OUT <- file.path("output/07_structural_diagnosis", SURVEY)
dir.create(PLOTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)
TOPO_LV <- paste0(SURVEY, c(" x2", " x1", "-DP"))      # display order: main -> no extra closure -> no homophily

BLAU_CAT <- switch(SURVEY, ELSOC = c("sex", "relig", "ideol"), GSS14 = c("race", "sex", "relig", "ideol"),
                   c("race", "sex", "relig"))    # categorical Blau axes
BLAU_NUM <- c("age", "educ_num")          # continuous Blau axes
AXIS_LAB <- c(race = "Race", sex = "Sex", relig = "Religion", ideol = "Ideology", age = "Age", educ = "Education")

# DP twin of network k, identical to 05's (make_pair in 05_unified_engine.R)
make_dp <- function(g, k) {
  set.seed(500000 + k)
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
  if (is.null(g$blau)) {        # GSS, ATP: age, educ_num, race, relig, sex
  df <- data.frame(
    age      = as.numeric(igraph::vertex_attr(g, "age")),
    educ_num = as.numeric(igraph::vertex_attr(g, "educ_num")),
    race     = factor(igraph::vertex_attr(g, "race")),
    relig    = factor(igraph::vertex_attr(g, "relig")),
    sex      = factor(igraph::vertex_attr(g, "sex"))
  )
  } else {                      # the network's own Blau space (02 writes g$blau; as gower_of() of 05)
    df <- as.data.frame(lapply(setNames(g$blau, g$blau), function(a) {
      x <- igraph::vertex_attr(g, a); if (a %in% BLAU_NUM) as.numeric(x) else factor(x) }))
  }
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
    modularity         = { set.seed(1); igraph::modularity(igraph::cluster_louvain(g)) },
    giant_frac         = max(comps$csize) / igraph::vcount(g),
    # C. homophily (GSS, ATP: assort_race, _sex, _relig, _age, _educ)
    setNames(vapply(BLAU_CAT, function(a) assort_attr(g, a, TRUE), 0), paste0("assort_", BLAU_CAT)),
    assort_age         = assort_attr(g, "age",      FALSE),
    assort_educ        = assort_attr(g, "educ_num", FALSE),
    mean_edge_gower    = mean_edge_gower(g)
  )
  out
}

# ---------------------------------------------------------------------------
# run over the 24 networks of each level
# ---------------------------------------------------------------------------
g2s <- readRDS(file.path(NET_DIR, paste0(SURVEY, "_x2.rds")))
g1s <- readRDS(file.path(NET_DIR, paste0(SURVEY, "_x1.rds")))
cat("Computing statistics on", length(g2s), SURVEY, "networks (x2, x1) and the DP twins of x2...\n")
rows <- list()
for (i in seq_along(g2s)) {
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = TOPO_LV[1], t(stats_one(g2s[[i]])))
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = TOPO_LV[2], t(stats_one(g1s[[i]])))
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = TOPO_LV[3], t(stats_one(make_dp(g2s[[i]], i))))
  if (i %% 6 == 0) cat("  done", i, "/", length(g2s), "\n")
}
allstats <- bind_rows(rows)
allstats$topo <- factor(allstats$topo, levels = TOPO_LV)
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
  x  <- allstats[[s]][allstats$topo == TOPO_LV[1]]
  x1 <- allstats[[s]][allstats$topo == TOPO_LV[2]]
  y  <- allstats[[s]][allstats$topo == TOPO_LV[3]]
  p    <- tryCatch(wilcox.test(x, y)$p.value,  error = function(e) NA_real_)  # x2 vs DP
  p_x1 <- tryCatch(wilcox.test(x, x1)$p.value, error = function(e) NA_real_)  # x2 vs x1
  data.frame(statistic = s,
             mean_x2 = mean(x, na.rm = TRUE), mean_x1 = mean(x1, na.rm = TRUE),
             mean_DP = mean(y, na.rm = TRUE),
             pct_change_DP = 100 * (mean(y, na.rm = TRUE)  - mean(x, na.rm = TRUE)) / abs(mean(x, na.rm = TRUE)),
             pct_change_x1 = 100 * (mean(x1, na.rm = TRUE) - mean(x, na.rm = TRUE)) / abs(mean(x, na.rm = TRUE)),
             cohens_d = cohens_d(x, y), p_value = p, p_value_x1 = p_x1)
}) %>% bind_rows()
write.csv(summ, file.path(DATA_OUT, "structural_diagnosis_summary.csv"), row.names = FALSE)
cat("\n=============== STRUCTURAL DIAGNOSIS SUMMARY ===============\n")
print(summ, row.names = FALSE, digits = 3)
cat("===========================================================\n")

# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------
# DP uses the deck's dark green (ndgreen = RGB 0,118,54) to match the
# "Two ways to break the network" concept map.
COLS <- setNames(c("#E74C3C", "#2980B9", "#007636"), TOPO_LV)

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
#     Blau axes, per-box x2/x1/DP labels (no legend), compact. ---
AX <- c(BLAU_CAT, "age", "educ")
AXIS_ORDER <- sort(unname(AXIS_LAB[AX]))          # GSS: Age, Education, Race, Religion, Sex
assort_long <- allstats %>%
  select(topo, all_of(paste0("assort_", AX))) %>%
  pivot_longer(-topo, names_to = "axis", values_to = "val") %>%
  mutate(axis = factor(unname(AXIS_LAB[sub("^assort_", "", axis)]), levels = AXIS_ORDER))
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
