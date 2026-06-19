# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structural diagnosis: what does "breaking homophily" (GSS -> GSS-DP) destroy?
# 05_structural_diagnosis.R   (Paper 4.1.1)
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
PLOTS    <- "playground/checks-claude/plots"
DATA_OUT <- "playground/checks-claude/data"
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

rows <- list()
for (i in seq_along(files)) {
  g  <- load_igraph(files[i])
  gd <- make_dp(g)
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = "GSS",    t(stats_one(g)))
  rows[[length(rows) + 1]] <- data.frame(net = i, topo = "GSS-DP", t(stats_one(gd)))
  if (i %% 10 == 0) cat("  done", i, "/", length(files), "\n")
}
allstats <- bind_rows(rows)
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
  x <- allstats[[s]][allstats$topo == "GSS"]
  y <- allstats[[s]][allstats$topo == "GSS-DP"]
  p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  data.frame(statistic = s,
             mean_GSS = mean(x, na.rm = TRUE), mean_DP = mean(y, na.rm = TRUE),
             pct_change = 100 * (mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE)) / abs(mean(x, na.rm = TRUE)),
             cohens_d = cohens_d(x, y), p_value = p)
}) %>% bind_rows()
write.csv(summ, file.path(DATA_OUT, "structural_diagnosis_summary.csv"), row.names = FALSE)
cat("\n=============== STRUCTURAL DIAGNOSIS SUMMARY ===============\n")
print(summ, row.names = FALSE, digits = 3)
cat("===========================================================\n")

# ---------------------------------------------------------------------------
# plotting helpers
# ---------------------------------------------------------------------------
COLS <- c("GSS" = "#E74C3C", "GSS-DP" = "#7F8C8D")

# single boxplot with the GSS / GSS-DP labels placed ABOVE each box (no legend, no p-value)
box_one <- function(stat, title) {
  d <- allstats[, c("topo", stat)]; names(d)[2] <- "y"
  ytop <- max(d$y, na.rm = TRUE); yr <- diff(range(d$y, na.rm = TRUE))
  labs_df <- data.frame(topo = c("GSS", "GSS-DP"), y = ytop + 0.10 * yr)
  ggplot(d, aes(topo, y, fill = topo)) +
    geom_boxplot(width = 0.55, outlier.size = 0.4, alpha = 0.85) +
    geom_text(data = labs_df, aes(label = topo, color = topo), vjust = 0, size = 3,
              fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = COLS, guide = "none") +
    scale_color_manual(values = COLS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
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

# --- Plot 2: family C, 1x2 (assortativity panel wider, social-distance panel thinner) ---
# left: per-axis attribute assortativity; GSS/GSS-DP labels above EACH axis pair (no legend)
AXIS_ORDER <- c("Age", "Education", "Race", "Religion", "Sex")
assort_long <- allstats %>%
  select(topo, assort_race, assort_sex, assort_relig, assort_age, assort_educ) %>%
  pivot_longer(-topo, names_to = "axis", values_to = "val") %>%
  mutate(axis = factor(recode(axis, assort_race="Race", assort_sex="Sex", assort_relig="Religion",
                       assort_age="Age", assort_educ="Education"), levels = AXIS_ORDER))
ytopC <- max(assort_long$val, na.rm = TRUE)
# one GSS + one GSS-DP label above every axis pair (dodge half-width ~0.2)
labC <- do.call(rbind, lapply(seq_along(AXIS_ORDER), function(k) data.frame(
  axis = factor(AXIS_ORDER[k], levels = AXIS_ORDER),
  topo = c("GSS", "GSS-DP"),
  x    = k + c(-0.2, 0.2),
  val  = ytopC + 0.085)))
pC1 <- ggplot(assort_long, aes(axis, val, fill = topo)) +
  geom_boxplot(width = 0.55, outlier.size = 0.3, alpha = 0.85, position = position_dodge(0.55)) +
  geom_text(data = labC, aes(x = x, label = topo, color = topo), vjust = 0,
            size = 2.1, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_color_manual(values = COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(title = "Attribute assortativity by Blau axis", x = NULL, y = "Assortativity") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 9, face = "bold"))
pC2 <- box_one("mean_edge_gower", "Mean social distance")
pC <- pC1 + pC2 + plot_layout(widths = c(2.6, 1))
ggsave(file.path(PLOTS, "structural_diagnosis_C.pdf"), pC, width = 8.2, height = 3.3)

cat("\nSaved: structural_diagnosis_summary.csv, _B.pdf (3x2), _C.pdf (1x2)\n")
