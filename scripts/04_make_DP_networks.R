# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Persist the GSS-DP (degree-preserving) twin networks to disk
# 04_make_DP_networks.R
#
# WHY: 05_unified_diffusion_sweep.R builds the DP twins in memory with the recipe
#
#     set.seed(500000 + k); rewire(g, with = keeping_degseq(niter = 20 * m))
#
# which is deterministic (verified: identical edge-list digests across separate
# R processes). This script reproduces the SAME twins with the SAME recipe and
# saves them, so that
#   (a) later runs (other seeding strategies, gamma-form, criticality) reuse the
#       EXACT same twins even if igraph's RNG changes across versions;
#   (b) descriptive statistics can be computed on the twins actually used.
#
# Format mirrors data/02_GSS_SH_network/ (statnet 'network' objects), so the
# legacy pipeline can load them directly.
#
# Outputs:
#   data/02_GSS_DP_network/GSS_dp_sim_1000_XXX.rds        (96 networks)
#   output/04_null_models/dp_verification.csv     (GSS vs DP stats)
#
# Runtime: ~2-4 min single-threaded (safe to run alongside the sweep).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(igraph); library(network); library(intergraph); library(cluster)
  library(dplyr)
})

SRC_DIR  <- "data/02_GSS_network_ergm"
OUT_DIR  <- "data/02_GSS_DP_network"
DATA_OUT <- "output/04_null_models"
N_NETS   <- 96
NITER_F  <- 20        # niter = NITER_F * |E|   (identical to script 12)
SEED_OFF <- 500000    # set.seed(SEED_OFF + k)  (identical to script 12)

dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(DATA_OUT, showWarnings = FALSE, recursive = TRUE)

GOWER_VARS <- c("age", "educ_num", "race", "relig", "sex")

net_stats <- function(g, d_mat) {
  el <- as_edgelist(g, names = FALSE)
  ed <- d_mat[cbind(c(el[, 1], el[, 2]), c(el[, 2], el[, 1]))]
  c(mean_degree   = mean(degree(g)),
    density       = edge_density(g),
    clustering    = transitivity(g, type = "global"),
    triangles     = sum(count_triangles(g)) / 3,
    assort_degree = assortativity_degree(g),
    assort_race   = assortativity_nominal(g, as.integer(as.factor(V(g)$race))),
    assort_relig  = assortativity_nominal(g, as.integer(as.factor(V(g)$relig))),
    assort_age    = assortativity(g, V(g)$age),
    mean_gower_on_edges = mean(ed))
}

message(format(Sys.time(), "%H:%M:%S"), "  Building and saving ", N_NETS, " DP twins...")
rows <- vector("list", N_NETS)
for (k in seq_len(N_NETS)) {
  path <- file.path(SRC_DIR, sprintf("GSS_net_sim_1000_%03d.rds", k))
  if (!file.exists(path)) stop("Missing network file: ", path)
  g <- asIgraph(readRDS(path))

  attrs <- data.frame(
    age = V(g)$age, educ_num = V(g)$educ_num,
    race = as.factor(V(g)$race), relig = as.factor(V(g)$relig),
    sex = as.factor(V(g)$sex)
  )
  d_mat <- as.matrix(daisy(attrs, metric = "gower"))

  # --- the exact recipe used by 05_unified_diffusion_sweep.R ---
  set.seed(SEED_OFF + k)
  g_dp <- rewire(g, with = keeping_degseq(niter = NITER_F * ecount(g)))

  # save as statnet 'network' object (same class/format as the GSS and SH files)
  saveRDS(asNetwork(g_dp),
          file.path(OUT_DIR, sprintf("GSS_dp_sim_1000_%03d.rds", k)))

  rows[[k]] <- bind_rows(
    data.frame(net = k, topology = "GSS", t(net_stats(g,    d_mat))),
    data.frame(net = k, topology = "DP",  t(net_stats(g_dp, d_mat)))
  )
  if (k %% 12 == 0) message("  ...", k, "/", N_NETS)
}

verif <- bind_rows(rows)
write.csv(verif, file.path(DATA_OUT, "dp_verification.csv"), row.names = FALSE)

summ <- verif |>
  group_by(topology) |>
  summarise(across(-net, mean), .groups = "drop") |>
  mutate(topology = factor(topology, levels = c("GSS", "DP"))) |>
  arrange(topology)

message("\n=== GSS vs GSS-DP (means over ", N_NETS, " network pairs) ===")
print(as.data.frame(summ), digits = 3)
message("\nExpected signature: degree/density IDENTICAL; clustering, triangles and\n",
        "attribute assortativities DOWN; mean Gower distance on ties UP (~0.23 -> ~0.33).")
message("Saved ", OUT_DIR, "/ (", N_NETS, " files) and ",
        file.path(DATA_OUT, "dp_verification.csv"))
