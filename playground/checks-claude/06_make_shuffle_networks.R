# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create the GSS-SH (attribute-SHuffle) networks and verify them
# 06_make_shuffle_networks.R   (Paper 4.1.2 — null model #3)
#
# GSS-SH = take each plausible GSS network, keep its topology BYTE-IDENTICAL,
# and randomly PERMUTE the node attributes (demographics + propensity) across
# vertices. This destroys ONLY homophily (the attribute<->position alignment),
# leaving clustering / triangles / modularity / degree exactly intact — the
# clean counterfactual that the GSS-DP (rewiring) baseline cannot give.
#
# This script ONLY builds + saves the networks and checks the statistics; it
# does NOT run any diffusion. (Building a shuffle net ~ 0.7 ms.)
#
# Output:
#   - data/02_GSS_SH_network/GSS_sh_sim_1000_XXX.rds   (statnet 'network' objects,
#       same class/format as the GSS originals, so the sim pipeline can load them)
#   - playground/checks-claude/data/shuffle_verification.csv
#   - playground/checks-claude/plots/shuffle_verification.pdf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(network); library(igraph); library(intergraph); library(cluster)
  library(dplyr); library(tidyr); library(ggplot2)
})

set.seed(20260620)
SRC_DIR  <- "data/02_GSS_network_ergm"
SRC_PAT  <- "^GSS_net_sim_1000_(\\d{3})\\.rds$"
OUT_DIR  <- "data/02_GSS_SH_network"
DATA_OUT <- "playground/checks-claude/data"
PLOTS    <- "playground/checks-claude/plots"
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS,    showWarnings = FALSE, recursive = TRUE)

# attributes to permute together (one common permutation keeps each person's
# full profile intact, only moving the WHOLE person to another node slot)
SHUFFLE_ATTRS <- c("age", "educ_num", "race", "relig", "sex", "mur_score",
                   "signdpet","avoidbuy","joindem","attrally","cntctgov",
                   "polfunds","usemedia","interpol","actlaw")
BLAU_CAT <- c("race","sex","relig"); BLAU_NUM <- c("age","educ_num")
GOWER_VARS <- c("age","educ_num","race","relig","sex")

# ---------------------------------------------------------------------------
# make ONE GSS-SH network from a statnet 'network' object: permute node attrs
# (topology untouched). Returns a statnet 'network' object (same format in/out).
# ---------------------------------------------------------------------------
make_shuffle <- function(net) {
  n <- network::network.size(net)
  perm <- sample(n)
  present <- intersect(SHUFFLE_ATTRS, network::list.vertex.attributes(net))
  for (a in present) {
    v <- network::get.vertex.attribute(net, a)
    network::set.vertex.attribute(net, a, v[perm])
  }
  net
}

# ---------------------------------------------------------------------------
# stats to verify: a few B (must be ~identical) and C (must collapse)
# ---------------------------------------------------------------------------
assort_attr <- function(g, attr, categorical) {
  v <- igraph::vertex_attr(g, attr)
  if (categorical) igraph::assortativity_nominal(g, as.integer(factor(v)), directed = FALSE)
  else igraph::assortativity(g, as.numeric(v), directed = FALSE)
}
mean_edge_gower <- function(g) {
  df <- data.frame(age=as.numeric(igraph::vertex_attr(g,"age")),
                   educ_num=as.numeric(igraph::vertex_attr(g,"educ_num")),
                   race=factor(igraph::vertex_attr(g,"race")),
                   relig=factor(igraph::vertex_attr(g,"relig")),
                   sex=factor(igraph::vertex_attr(g,"sex")))
  D <- as.matrix(cluster::daisy(df, metric="gower"))
  el <- igraph::as_edgelist(g, names = FALSE)
  mean(D[el], na.rm = TRUE)
}
stats_one <- function(g) c(
  # B: topology — MUST be identical to GSS (shuffle doesn't touch edges)
  clustering_global = igraph::transitivity(g, "global"),
  triangles         = sum(igraph::count_triangles(g))/3,
  modularity        = igraph::modularity(igraph::cluster_louvain(g)),
  assort_degree     = igraph::assortativity_degree(g, directed = FALSE),
  mean_degree       = mean(igraph::degree(g)),
  # C: homophily — MUST collapse toward 0 / rise
  assort_race  = assort_attr(g,"race",TRUE),
  assort_age   = assort_attr(g,"age",FALSE),
  assort_relig = assort_attr(g,"relig",TRUE),
  mean_edge_gower = mean_edge_gower(g)
)

# ---------------------------------------------------------------------------
# build, save, and collect verification stats over all source networks
# ---------------------------------------------------------------------------
files <- list.files(SRC_DIR, pattern = SRC_PAT, full.names = TRUE)
cat("Building", length(files), "GSS-SH networks (topology fixed, attributes permuted)...\n")

t0 <- Sys.time()
rows <- list()
for (f in files) {
  idx <- sub(SRC_PAT, "\\1", basename(f))
  net  <- readRDS(f)
  gss_g <- intergraph::asIgraph(net)                 # original (for the GSS row)
  sh_net <- make_shuffle(net)                        # permute attrs in place
  sh_g   <- intergraph::asIgraph(sh_net)
  saveRDS(sh_net, file.path(OUT_DIR, sprintf("GSS_sh_sim_1000_%s.rds", idx)))
  rows[[length(rows)+1]] <- data.frame(net = idx, topo = "GSS",    t(stats_one(gss_g)))
  rows[[length(rows)+1]] <- data.frame(net = idx, topo = "GSS-SH", t(stats_one(sh_g)))
  if (as.integer(idx) %% 20 == 0) cat("  ", idx, "\n")
}
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
ver <- bind_rows(rows)
write.csv(ver, file.path(DATA_OUT, "shuffle_verification.csv"), row.names = FALSE)

cat(sprintf("\nBuilt + saved %d GSS-SH networks in %.1f s (%.1f ms each incl. stats).\n",
            length(files), elapsed, 1000*elapsed/length(files)))

# ---------------------------------------------------------------------------
# verification summary: B should match GSS, C should collapse
# ---------------------------------------------------------------------------
summ <- ver %>% group_by(topo) %>%
  summarise(across(-net, ~mean(.x, na.rm = TRUE)), .groups = "drop")
cat("\n=============== GSS  vs  GSS-SH  (means) ===============\n")
print(as.data.frame(summ), digits = 4)
cat("\nExpected: B (clustering/triangles/modularity/assort_degree/mean_degree) IDENTICAL;\n")
cat("          C (attribute assortativities -> ~0, mean_edge_gower up).\n")
cat("=======================================================\n")

# quick visual: topology controls (should overlap) + homophily (should split)
vlong <- ver %>%
  select(topo, clustering_global, modularity, assort_degree,
         assort_race, assort_age, mean_edge_gower) %>%
  pivot_longer(-topo, names_to = "stat", values_to = "val") %>%
  mutate(family = ifelse(stat %in% c("clustering_global","modularity","assort_degree"),
                         "B: topology (expect identical)", "C: homophily (expect collapse)"))
pV <- ggplot(vlong, aes(topo, val, fill = topo)) +
  geom_boxplot(width = 0.6, outlier.size = 0.4, alpha = 0.85) +
  facet_wrap(~ stat, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("GSS"="#E74C3C","GSS-SH"="#2980B9"), guide = "none") +
  labs(title = "GSS vs GSS-SH: topology preserved, homophily destroyed",
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face="bold", size=10))
ggsave(file.path(PLOTS, "shuffle_verification.pdf"), pV, width = 9, height = 5)

cat("\nSaved: data/02_GSS_SH_network/ (100 nets), shuffle_verification.csv + .pdf\n")
