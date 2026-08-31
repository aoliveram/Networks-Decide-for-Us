# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE STRUCTURAL PREMIUM under the unified (bottom-up) linear rule:
#
#   adopt_i(t)  <=>  Gamma + lambda * E~_i(t) >= q_i    ( & at least 1 contact )
#
# GSS plausible networks vs. their DEGREE-PRESERVING twins (GSS-DP), built here
# with igraph::rewire(keeping_degseq, niter = 20*m) — the true DP null,
# constructed fresh per network (this also sidesteps the legacy "ER" naming
# ambiguity of the old engine's outputs: the new-engine premium is genuinely DP).
#
# Design (per MAKING-IT-BOTTOM-UP.tex §8 and chat decisions, 2026-08-31):
#   - IUL: 41 levels (0..1 by 0.025)          — legacy grid
#   - MSP: 25 levels (0..1 by 1/24)           — NEW denser grid; nests the
#                                               legacy 13-level (1/12) grid
#   - lambda: {0.5, ..., 1.0}                 — the identified high-fit band plus
#                                               lambda = 1, the immunity boundary
#                                               (q <= 1 => no immune agents)
#   - Networks: GSS 001..096, one run each; DP twin per network
#   - Seeding: one strategy per invocation (--seeding=), 1 primary node +
#     a neighbour cluster of ~ 0.40 * degree(primary):
#       random    a uniformly drawn node
#       central   the highest-degree node
#       closeness the highest-closeness node
#       eigen     the highest-eigenvector-centrality node
#       marginal  a node drawn from the bottom 10% by degree
#   - Gate: sigma((h - d_ij)/0.02), stochastic per step; denominator = degree
#
# Total per seeding: 2 topologies x 6 lambda x 41 x 25 x 96 = 1,180,800 sims.
# Runtime: ~50 min on 8 cores. Checkpoints per (topology, lambda) — safe to
# kill and relaunch; completed combos are skipped.
#
# Outputs (<S> = the seeding strategy):
#   output/05_unified_diffusion/checkpoints/<S>/<topo>_lambda_X.XX.rds
#   output/05_unified_diffusion/unified_premium_results_<S>.csv
#   output/05_unified_diffusion/unified_premium_bam_<S>.csv   (beta_DP per lambda)
#   plots/05_unified_diffusion/unified_premium_heatmaps_<S>.pdf
#
# Usage:
#   Rscript scripts/05_unified_diffusion_sweep.R                     # random
#   Rscript scripts/05_unified_diffusion_sweep.R --seeding=central
#   Rscript scripts/05_unified_diffusion_sweep.R --test              # smoke
#   Rscript scripts/05_unified_diffusion_sweep_main.R                # all five
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(igraph); library(intergraph); library(cluster)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(parallel); library(mgcv)
})
Sys.setenv(OMP_NUM_THREADS = "1")

ARGS      <- commandArgs(trailingOnly = TRUE)
TEST_MODE <- "--test" %in% ARGS

# Seeding strategy: --seeding=<random|central|marginal|closeness|eigen>
SEEDING <- sub("^--seeding=", "", grep("^--seeding=", ARGS, value = TRUE))
if (length(SEEDING) == 0) SEEDING <- "random"
stopifnot(SEEDING %in% c("random", "central", "marginal", "closeness", "eigen"))

# ------------------------------- configuration -------------------------------
NETWORKS_DIR <- "data/02_GSS_network_ergm/"
# --test writes everything (checkpoints, CSVs, plots) into its own sandbox, so a
# smoke run can never overwrite or get skipped by the real sweep's results
DATA_OUT     <- if (TEST_MODE) "output/05_unified_diffusion/_test" else "output/05_unified_diffusion"
PLOTS_OUT    <- if (TEST_MODE) "output/05_unified_diffusion/_test" else "plots/05_unified_diffusion"
# one checkpoint sub-directory and one set of result files per seeding strategy,
# so strategies never overwrite each other and each can be resumed on its own
CKPT_DIR     <- file.path(DATA_OUT, "checkpoints", SEEDING)
SUFFIX       <- paste0("_", SEEDING)

N_RUNS       <- if (TEST_MODE) 2 else 96
N_CORES      <- if (TEST_MODE) 2 else 8
IUL_SWEEP    <- if (TEST_MODE) seq(0, 1, by = 0.25) else seq(0, 1, by = 0.025)
H_SWEEP      <- if (TEST_MODE) seq(0, 1, by = 0.25) else seq(0, 1, by = 1/24)
LAMBDAS      <- if (TEST_MODE) c(0.7) else c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
TOPOLOGIES   <- c("GSS", "DP")
GATE_T       <- 0.02
TAU_REF_SEED <- 0.40
MAX_STEPS    <- 200
PATIENCE     <- 25

dir.create(CKPT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_OUT, showWarnings = FALSE, recursive = TRUE)

sigmoid_gate <- function(d, h, T = GATE_T) 1 / (1 + exp((d - h) / T))

edge_arrays <- function(g, d_mat) {
  el <- as_edgelist(g, names = FALSE)
  src <- c(el[, 1], el[, 2]); dst <- c(el[, 2], el[, 1])
  deg <- degree(g)
  out <- list(g = g, src = src, dst = dst, edge_d = d_mat[cbind(src, dst)],
              deg = pmax(deg, 1L), deg_raw = deg)
  # position-based primary seeds: computed once per network and topology
  # (only the one the current strategy needs — closeness/eigen are not free)
  if (SEEDING == "central")   out$primary_central   <- which.max(deg)[1]
  if (SEEDING == "closeness") out$primary_closeness <- which.max(closeness(g))[1]
  if (SEEDING == "eigen")     out$primary_eigen     <- which.max(eigen_centrality(g)$vector)[1]
  out
}

# ------------------- load GSS networks & build DP twins -----------------------
message(format(Sys.time(), "%H:%M:%S"), "  Loading ", N_RUNS,
        " GSS networks, computing Gower distances, building DP twins...")
nets <- vector("list", N_RUNS)
for (k in seq_len(N_RUNS)) {
  path <- file.path(NETWORKS_DIR, sprintf("GSS_net_sim_1000_%03d.rds", k))
  if (!file.exists(path)) stop("Missing network file: ", path)
  g <- asIgraph(readRDS(path))
  attrs <- data.frame(
    age = V(g)$age, educ_num = V(g)$educ_num,
    race = as.factor(V(g)$race), relig = as.factor(V(g)$relig),
    sex = as.factor(V(g)$sex)
  )
  d_mat <- as.matrix(daisy(attrs, metric = "gower"))
  # DP twin: degree-preserving double-edge swaps (the engine's own null recipe)
  set.seed(500000 + k)
  g_dp <- rewire(g, with = keeping_degseq(niter = 20 * ecount(g)))
  nets[[k]] <- list(n = vcount(g), q = V(g)$mur_score,
                    GSS = edge_arrays(g, d_mat), DP = edge_arrays(g_dp, d_mat))
  if (k %% 12 == 0) message("  ...", k, "/", N_RUNS)
}

# quick sanity: DP preserves degrees, destroys homophily (edge distance rises)
dg <- sapply(nets, function(x) all(sort(x$GSS$deg_raw) == sort(x$DP$deg_raw)))
md <- sapply(nets, function(x) c(mean(x$GSS$edge_d), mean(x$DP$edge_d)))
message(sprintf("Sanity — degree sequences preserved in %d/%d twins; mean tie distance %.3f (GSS) vs %.3f (DP).",
                sum(dg), length(dg), mean(md[1, ]), mean(md[2, ])))

draw_seeds <- function(top, lam, run, topo) {
  set.seed(run * 2000 + round(lam * 100) + ifelse(topo == "DP", 777L, 0L))
  n <- length(top$deg_raw)
  primary <- switch(SEEDING,
    # deterministic, position-based strategies (precomputed at load time)
    central   = top$primary_central,
    closeness = top$primary_closeness,
    eigen     = top$primary_eigen,
    # stochastic strategies (seeded above, so reproducible per run/lambda/topology)
    random    = sample.int(n, 1),
    marginal  = as.integer(sample(
                  order(top$deg_raw)[seq_len(ceiling(n * 0.1))], 1)))
  n_seeds <- max(1, min(round(TAU_REF_SEED * top$deg_raw[primary]),
                        top$deg_raw[primary] + 1))
  seeds <- primary
  if (n_seeds > 1) {
    nb <- as.integer(neighbors(top$g, primary))
    if (length(nb) > 0)
      seeds <- c(primary, sample(nb, min(length(nb), n_seeds - 1)))
  }
  unique(seeds)
}

run_unified <- function(top, q, n, seeds, h, Gamma, lambda, sim_seed) {
  tau_star <- pmax(0, (q - Gamma) / lambda)
  w_edge   <- sigmoid_gate(top$edge_d, h)
  a <- rep(FALSE, n); a[seeds] <- TRUE
  set.seed(sim_seed)
  stall <- 0L; steps <- 0L
  for (t in seq_len(MAX_STEPS)) {
    idx <- which(a[top$src] & !a[top$dst])
    if (length(idx) == 0L) break
    pass <- runif(length(idx)) <= w_edge[idx]
    counts <- tabulate(top$dst[idx][pass], nbins = n)
    E <- counts / top$deg
    new <- !a & counts > 0L & E >= tau_star
    if (any(new)) { a[new] <- TRUE; stall <- 0L; steps <- t } else stall <- stall + 1L
    if (all(a) || stall >= PATIENCE) break
  }
  non_seed <- setdiff(which(a), seeds)
  c(num_adopters = sum(a),
    num_adopted_rational = sum(q[non_seed] <= Gamma),
    num_adopted_social   = sum(q[non_seed] >  Gamma),
    num_steps = steps, initial_cluster_size = length(seeds))
}

# --------------------------------- the sweep ----------------------------------
sims_per_combo <- N_RUNS * length(H_SWEEP) * length(IUL_SWEEP)
n_combos <- length(TOPOLOGIES) * length(LAMBDAS)
message(format(Sys.time(), "%H:%M:%S"), "  Seeding: ", toupper(SEEDING), " | Sweep: ",
        n_combos, " (topology, lambda) combos x ", sims_per_combo, " sims = ",
        n_combos * sims_per_combo, " total simulations.")
t_all <- Sys.time()

combos <- expand.grid(topo = TOPOLOGIES, lam = LAMBDAS, stringsAsFactors = FALSE)
for (ci in seq_len(nrow(combos))) {
  topo <- combos$topo[ci]; lam <- combos$lam[ci]
  ck <- file.path(CKPT_DIR, sprintf("%s_lambda_%.2f.rds", topo, lam))
  if (file.exists(ck)) {
    message(format(Sys.time(), "%H:%M:%S"), "  ", topo, " lambda=", sprintf("%.2f", lam),
            " checkpoint exists — skipping.")
    next
  }
  t_c <- Sys.time()
  per_run <- mclapply(seq_len(N_RUNS), function(run) {
    nt <- nets[[run]]; top <- nt[[topo]]
    seeds <- draw_seeds(top, lam, run, topo)
    rows <- vector("list", length(H_SWEEP) * length(IUL_SWEEP)); ri <- 0L
    for (h in H_SWEEP) {
      for (Gamma in IUL_SWEEP) {
        ri <- ri + 1L
        sim_seed <- run * 1e6 + round(h * 24) * 1e4 + round(Gamma * 1000) +
          round(lam * 100) * 7L + ifelse(topo == "DP", 3e8, 0)
        out <- run_unified(top, nt$q, nt$n, seeds, h, Gamma, lam, sim_seed)
        rows[[ri]] <- data.frame(topology = topo, lambda = lam, run_id = run,
                                 social_distance_h = h, innovation_iul_Gamma = Gamma,
                                 t(out), N_nodes = nt$n)
      }
    }
    do.call(rbind, rows)
  }, mc.cores = N_CORES)
  bad <- vapply(per_run, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("Worker failure at ", topo, " lambda=", lam,
                     " (runs: ", paste(which(bad), collapse = ","), ")")
  saveRDS(do.call(rbind, per_run), ck)
  el_min <- as.numeric(difftime(Sys.time(), t_c, units = "mins"))
  done <- sum(file.exists(file.path(
    CKPT_DIR, sprintf("%s_lambda_%.2f.rds", combos$topo, combos$lam))))
  message(format(Sys.time(), "%H:%M:%S"), "  ", topo, " lambda=", sprintf("%.2f", lam),
          " done in ", round(el_min, 1), " min  (", done, "/", n_combos,
          "; ~", round(el_min * (n_combos - done)), " min remaining)")
}

# ------------------------------ combine & save --------------------------------
results <- bind_rows(lapply(
  file.path(CKPT_DIR, sprintf("%s_lambda_%.2f.rds", combos$topo, combos$lam)), readRDS))
write.csv(results, file.path(DATA_OUT, paste0("unified_premium_results", SUFFIX, ".csv")), row.names = FALSE)
message(format(Sys.time(), "%H:%M:%S"), "  Combined results saved (", nrow(results),
        " rows; total ", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1),
        " min).")

# ---------------- THE PREMIUM: BAM beta_DP per lambda + pooled ----------------
results <- results |>
  mutate(topology = factor(topology, levels = c("GSS", "DP")),
         fail = N_nodes - num_adopters)

fit_one <- function(df, kk = if (TEST_MODE) 4 else 10) {
  m <- bam(cbind(num_adopters, fail) ~ s(innovation_iul_Gamma, social_distance_h, k = kk)
           + topology,
           family = binomial(), data = df, discrete = TRUE, nthreads = N_CORES)
  co <- summary(m)$p.table["topologyDP", ]
  data.frame(beta_DP = co[1], se = co[2], OR = exp(co[1]),
             odds_drop_pct = 100 * (exp(co[1]) - 1),
             dev_expl = summary(m)$dev.expl)
}

message("\n=== STRUCTURAL PREMIUM under the unified rule (ref. = GSS) ===")
per_lambda <- results |>
  group_by(lambda) |> group_modify(~ fit_one(.x)) |> ungroup()
print(as.data.frame(per_lambda), digits = 3)

pooled <- fit_one(results |> mutate(lambda_f = factor(lambda)))
message("Pooled (all lambdas): beta_DP = ", round(pooled$beta_DP, 3),
        "  OR = ", round(pooled$OR, 3),
        "  odds drop = ", round(pooled$odds_drop_pct, 1), "%")
write.csv(bind_rows(per_lambda |> mutate(lambda = as.character(lambda)),
                    pooled |> mutate(lambda = "pooled")),
          file.path(DATA_OUT, paste0("unified_premium_bam", SUFFIX, ".csv")), row.names = FALSE)
message("(legacy engine reference: beta_DP = -0.850, OR = 0.43, -57% odds)")

# ---------------------------------- plots -------------------------------------
agg <- results |>
  group_by(topology, lambda, Gamma = innovation_iul_Gamma, h = social_distance_h) |>
  summarise(phi = mean(num_adopters / N_nodes), .groups = "drop")

heat <- function(df, title) {
  ggplot(df, aes(Gamma, h, fill = phi)) +
    geom_tile(color = "white", lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis",
                         n.breaks = 4, na.value = "grey90", name = "Adoption") +
    labs(title = title, x = "IUL (Γ)", y = "MSP (h)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))
}

cairo_pdf(file.path(PLOTS_OUT, paste0("unified_premium_heatmaps", SUFFIX, ".pdf")),
          width = 11, height = 4.2, onefile = TRUE)
for (l in LAMBDAS) {
  print(
    heat(filter(agg, lambda == l, topology == "GSS"),
         sprintf("GSS — Φ_total (unified rule, λ = %.1f)", l)) +
    heat(filter(agg, lambda == l, topology == "DP"),
         sprintf("GSS-DP (rewired) — Φ_total (λ = %.1f)", l)) +
    plot_layout(nrow = 1, guides = "collect")
  )
}
dev.off()
message(format(Sys.time(), "%H:%M:%S"), "  Saved ",
        file.path(PLOTS_OUT, paste0("unified_premium_heatmaps", SUFFIX, ".pdf")), " — DONE.")
