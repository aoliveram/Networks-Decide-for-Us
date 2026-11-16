# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE STRUCTURAL PREMIUM under the unified (bottom-up) linear rule:
#
#   adopt_i  <=>  Gamma + lambda * E_i >= q_i    ( & >= 1 adopting contact behind an open tie )
#
# Networks with empirical homophily (NET) vs their DEGREE-PRESERVING twins (DP),
# built here with igraph::rewire(keeping_degseq, niter = 20 m), set.seed(500000 + k).
# The twin is rewired from the x2 network also in the x1 sweep (same nodes and
# degrees): both levels share ONE twin, so x2 vs x1 is paired.
#
# Design (since 2026-11, PR #3 / #4):
#   - Networks: data/02_calibrated_networks/<SURVEY>_<LEVEL>.rds (02_network_calibrated.R),
#     24 per set, N = 1000, mean degree 30, nodes resampled per network.
#     SURVEY: GSS (2004, collective action), ATP (2014, new products), GSS14 (2014,
#     collective action, GSS 2004 targets) or ELSOC (2018 willingness, Chilean targets; PR #5).
#     LEVEL:  x2 (main: calibrated + rewired to twice the triangles) or x1
#             (calibrated, before the rewiring; same nodes, degrees and mixing).
#   - Engine: Rcpp, FROZEN gate (one coin per tie per simulation, w = sigma((h - d_ij)/0.02));
#     each cascade stops at its fixed point (no step limit, no patience).
#   - Exposure: the base model E_i = a_i / d_i, and the robustness axis (one formula,
#     two dials; EXPOSURE_CONDITIONS in 05_unified_engine.R): independent
#     reinforcement theta (race x religion circles) and cohesive reinforcement kappa
#     (pairs of contacts who know each other), at c = -1, -0.5, +0.5, +1, each with
#     its permuted control. Networks with 'ideol' (ELSOC, GSS14) add the theta arm with
#     ideology circles (c = -1, -0.5 and permuted controls): 13 conditions.
#   - Grid: IUL Gamma 0..1 by 0.025 (41), MSP h 0..1 by 1/24 (25), lambda 0.5..1.0 by 0.1.
#   - Seeding: one strategy per invocation (--seeding=) picks the PRIMARY node:
#       random    a uniformly drawn node
#       central   the highest-degree node
#       closeness the highest-closeness node
#       eigen     the highest-eigenvector-centrality node
#       marginal  a node drawn from the bottom 10% by degree
#     and the seed cluster is the primary plus its nearest nodes, 0.40 x the
#     network's maximum degree in total (~28 seeds), the same in the network and its twin.
#
# Per seeding and set: 9 conditions x 2 topologies x 6 lambda x 1025 cells x 24 networks
# = 2,656,800 simulations (~0.15-0.3 ms each): ~2-3 min on 8 cores.
# Checkpoints per (condition, network): resumable; existing ones are skipped unless
# they were computed on other networks or another grid (fingerprint check).
#
# Outputs (<SET> = <SURVEY>_<LEVEL>, <S> = seeding):
#   output/05_unified_diffusion/<SET>/checkpoints/<S>/<condition>/run_XX.rds
#   output/05_unified_diffusion/<SET>/unified_premium_results_<S>.csv   base exposure, one row
#       per simulation; read by 06 and 06b. topology = "NET" is the network with empirical
#       homophily (whatever the survey; "GSS" until 2026-11), "DP" its twin.
#   output/05_unified_diffusion/<SET>/exposure_axis_<S>.rds   all conditions, compact
#       (integer array cells x topology x lambda x run x output); read by 06c
#   output/05_unified_diffusion/<SET>/unified_premium_bam_<S>.csv   beta_DP per lambda (s(k=10))
#   plots/05_unified_diffusion/<SET>/unified_premium_heatmaps_<S>.pdf
#
# Usage (repo root):
#   Rscript scripts/05_unified_diffusion_sweep.R                          # GSS x2, random, all conditions
#   Rscript scripts/05_unified_diffusion_sweep.R --survey=ATP --level=x1 --seeding=central
#   Rscript scripts/05_unified_diffusion_sweep.R --conditions=base        # base exposure only
#   Rscript scripts/05_unified_diffusion_sweep.R --test                   # smoke test, sandbox
#   Rscript scripts/05_unified_diffusion_sweep_main.R                     # everything, then 06-07
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressWarnings(suppressMessages({
  library(igraph); library(dplyr); library(ggplot2); library(patchwork)
  library(parallel); library(mgcv)
}))
Sys.setenv(OMP_NUM_THREADS = "1")
source("scripts/05_unified_engine.R")

ARGS <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) sub(".*=", "", a[1]) else d }
TEST_MODE <- "--test" %in% ARGS
SURVEY  <- toupper(arg("survey", "GSS")); stopifnot(SURVEY %in% c("GSS", "ATP", "GSS14", "ELSOC"))
LEVEL   <- tolower(arg("level", "x2"))    # x2 (main), x1; other rewired levels of 02 (e.g. x2n) if their file exists
SEEDING <- arg("seeding", "random")
stopifnot(SEEDING %in% c("random", "central", "marginal", "closeness", "eigen"))
COND_ARG <- arg("conditions", "all")
SET <- paste0(SURVEY, "_", LEVEL)

# ------------------------------- configuration -------------------------------
NET_FILE  <- file.path("data/02_calibrated_networks", paste0(SET, ".rds"))
DATA_OUT  <- file.path(if (TEST_MODE) "output/05_unified_diffusion/_test" else "output/05_unified_diffusion", SET)
PLOTS_OUT <- if (TEST_MODE) DATA_OUT else file.path("plots/05_unified_diffusion", SET)
CKPT_DIR  <- file.path(DATA_OUT, "checkpoints", SEEDING)
SUFFIX    <- paste0("_", SEEDING)

N_RUNS    <- if (TEST_MODE) 2L else 24L
N_CORES   <- as.integer(arg("cores", if (TEST_MODE) 2 else 8))
IUL_SWEEP <- if (TEST_MODE) seq(0, 1, by = 0.25) else seq(0, 1, by = 0.025)
H_SWEEP   <- if (TEST_MODE) seq(0, 1, by = 0.25) else seq(0, 1, by = 1/24)
LAMBDAS   <- if (TEST_MODE) c(0.7) else c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
TOPOLOGIES <- c("NET", "DP")
OUTS <- c("num_adopters", "num_adopted_rational", "t50", "num_steps")
for (d in c(CKPT_DIR, PLOTS_OUT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(NET_FILE)) stop("Missing ", NET_FILE, ": run scripts/02_network_calibrated.R --survey=", SURVEY,
                                 if (!LEVEL %in% c("x1", "x2")) paste0(" with the options of level ", LEVEL))
nets <- readRDS(NET_FILE)
ALL_CONDS <- exposure_conditions(nets[[1]])          # + the ideology theta arm if the nodes carry 'ideol'
CONDS <- if (COND_ARG == "all") ALL_CONDS else ALL_CONDS[strsplit(COND_ARG, ",")[[1]]]
if (any(vapply(CONDS, is.null, logical(1))) || !"base" %in% names(CONDS))
  stop("--conditions must include 'base' and only names of the network's conditions: ", paste(names(ALL_CONDS), collapse = ", "))
for (d in file.path(CKPT_DIR, names(CONDS))) dir.create(d, showWarnings = FALSE, recursive = TRUE)
DP_FILE <- file.path("data/02_calibrated_networks", paste0(SURVEY, "_x2.rds"))   # the twins come from x2
dp_nets <- if (LEVEL == "x2") nets else readRDS(DP_FILE)
stopifnot(length(nets) >= N_RUNS, length(dp_nets) >= N_RUNS)
# fingerprint of what a checkpoint depends on: network, twin source, grid, lambdas, seeding
edges_md5 <- function(g) digest::digest(as_edgelist(g, names = FALSE), algo = "md5")
FP <- vapply(seq_len(N_RUNS), function(r) digest::digest(list(edges_md5(nets[[r]]), edges_md5(dp_nets[[r]]),
                                                              H_SWEEP, IUL_SWEEP, LAMBDAS, SEEDING), algo = "md5"), "")

# --------------------------------- the sweep ----------------------------------
ck_path <- function(cn, run) file.path(CKPT_DIR, cn, sprintf("run_%02d.rds", run))
ck_done <- function(cn, run) { f <- ck_path(cn, run); file.exists(f) && identical(readRDS(f)$fp, FP[[run]]) }
n_cells <- length(H_SWEEP) * length(IUL_SWEEP)
message(format(Sys.time(), "%H:%M:%S"), "  ", SET, " | seeding ", toupper(SEEDING), " | ", length(CONDS),
        " conditions x 2 topologies x ", length(LAMBDAS), " lambda x ", n_cells, " cells x ", N_RUNS, " networks = ",
        format(length(CONDS) * 2 * length(LAMBDAS) * n_cells * N_RUNS, big.mark = ","), " simulations")
t_all <- Sys.time()

one_run <- function(run) {
  todo <- names(CONDS)[!vapply(names(CONDS), ck_done, logical(1), run = run)]
  if (!length(todo)) return(NULL)
  nt <- make_pair(nets[[run]], run, SEEDING, dp_from = dp_nets[[run]])
  for (topo in TOPOLOGIES) {
    top <- nt[[topo]]
    seeds <- draw_seeds(top, run, topo, SEEDING)
    W <- gate_matrix(top, H_SWEEP)
    ctx <- build_exposure_ctx(top, run, as.integer(topo == "DP"), CONDS[todo])
    nt[[topo]] <- c(top, list(seeds = seeds, W = W, ctx = ctx))
  }
  for (cn in todo) {
    arr <- array(0L, dim = c(n_cells, 2L, length(LAMBDAS), length(OUTS)),
                 dimnames = list(NULL, TOPOLOGIES, sprintf("%.2f", LAMBDAS), OUTS))
    for (topo in TOPOLOGIES) for (li in seq_along(LAMBDAS)) {
      tp <- nt[[topo]]
      arr[, topo, li, ] <- run_condition(tp, tp$W, tp$ctx, CONDS[[cn]], nt$q, tp$seeds, run, LAMBDAS[li], topo,
                                         H_SWEEP, IUL_SWEEP)[, OUTS]
    }
    saveRDS(list(run = run, condition = cn, N = nt$n, fp = FP[[run]],
                 seeds = c(NET = length(nt$NET$seeds), DP = length(nt$DP$seeds)), res = arr), ck_path(cn, run))
  }
  c(sanity_degrees = all(sort(nt$NET$deg_raw) == sort(nt$DP$deg_raw)),
    tie_d_net = mean(nt$NET$edge_d), tie_d_dp = mean(nt$DP$edge_d))
}
out <- mclapply(seq_len(N_RUNS), function(r) tryCatch(one_run(r), error = function(e) paste("ERROR", conditionMessage(e))),
                mc.cores = N_CORES, mc.preschedule = FALSE)
err <- Filter(function(x) inherits(x, "try-error") || (is.character(x) && grepl("^ERROR", x[1])), out)
if (length(err)) stop("worker failure: ", err[[1]])
san <- do.call(rbind, Filter(Negate(is.null), out))
if (!is.null(san))
  message(sprintf("Sanity: degree sequences preserved in %d/%d twins; mean tie distance %.3f (network) vs %.3f (DP)",
                  sum(san[, "sanity_degrees"]), nrow(san), mean(san[, "tie_d_net"]), mean(san[, "tie_d_dp"])))
message(format(Sys.time(), "%H:%M:%S"), "  sweep done in ",
        round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min")

# ------------------------------ combine & save --------------------------------
ck <- lapply(names(CONDS), function(cn) lapply(seq_len(N_RUNS), function(run) readRDS(ck_path(cn, run))))
names(ck) <- names(CONDS)
stopifnot(all(vapply(unlist(ck, recursive = FALSE), function(x) identical(x$fp, FP[[x$run]]), logical(1))))
N_NODES <- ck$base[[1]]$N
axis <- list(set = SET, survey = SURVEY, level = LEVEL, seeding = SEEDING, H = H_SWEEP, G = IUL_SWEEP, lambdas = LAMBDAS,
             N = N_NODES, seeds = t(vapply(ck$base, `[[`, integer(2), "seeds")),
             conditions = do.call(rbind, lapply(names(CONDS), function(cn) data.frame(
               condition = cn, theta = CONDS[[cn]]$theta, kappa = CONDS[[cn]]$kappa_p, c = CONDS[[cn]]$c,
               circles = CONDS[[cn]]$circles, pairs = CONDS[[cn]]$pairs, control_of = CONDS[[cn]]$control))),
             res = lapply(ck, function(L) simplify2array(lapply(L, `[[`, "res"))))   # cells x topo x lambda x out x run
saveRDS(axis, file.path(DATA_OUT, paste0("exposure_axis", SUFFIX, ".rds")))

# base condition in 05's long format (one row per simulation)
cell_h <- rep(H_SWEEP, each = length(IUL_SWEEP)); cell_g <- rep(IUL_SWEEP, times = length(H_SWEEP))
results <- bind_rows(lapply(ck$base, function(x) bind_rows(lapply(TOPOLOGIES, function(topo) bind_rows(lapply(seq_along(LAMBDAS), function(li) {
  m <- x$res[, topo, li, ]
  data.frame(topology = topo, lambda = LAMBDAS[li], run_id = x$run,
             social_distance_h = cell_h, innovation_iul_Gamma = cell_g,
             num_adopters = m[, "num_adopters"], num_adopted_rational = m[, "num_adopted_rational"],
             num_adopted_social = m[, "num_adopters"] - x$seeds[[topo]] - m[, "num_adopted_rational"],
             num_steps = m[, "num_steps"], initial_cluster_size = x$seeds[[topo]], t50 = m[, "t50"],
             N_nodes = x$N, survey = SURVEY, level = LEVEL)
}))))))
write.csv(results, file.path(DATA_OUT, paste0("unified_premium_results", SUFFIX, ".csv")), row.names = FALSE)
message(format(Sys.time(), "%H:%M:%S"), "  saved ", nrow(results), " rows (base) and the exposure axis (",
        length(CONDS), " conditions)")

# ---------------- THE PREMIUM: BAM beta_DP per lambda + pooled ----------------
# Continuity with 05/06: s(k = 10) surface + topology shift. 06b gives the
# te(12, 12) estimates with run-bootstrap intervals.
results <- results |>
  mutate(topology = factor(topology, levels = c("NET", "DP")), fail = N_nodes - num_adopters)
fit_one <- function(df, pooled = FALSE, kk = if (TEST_MODE) 4 else 10) {
  form <- if (pooled && length(unique(df$lambda)) > 1) {
    cbind(num_adopters, fail) ~ s(innovation_iul_Gamma, social_distance_h, k = kk) + topology + lambda_f
  } else {
    cbind(num_adopters, fail) ~ s(innovation_iul_Gamma, social_distance_h, k = kk) + topology
  }
  m <- bam(form, family = binomial(), data = df, discrete = TRUE)
  co <- summary(m)$p.table["topologyDP", ]
  data.frame(beta_DP = co[1], se = co[2], OR = exp(co[1]), odds_drop_pct = 100 * (exp(co[1]) - 1),
             dev_expl = summary(m)$dev.expl)
}
message("\n=== STRUCTURAL PREMIUM, ", SET, ", base exposure (ref. = network) ===")
per_lambda <- results |> group_by(lambda) |> group_modify(~ fit_one(.x)) |> ungroup()
print(as.data.frame(per_lambda), digits = 3)
gap <- results |> group_by(lambda, run_id, social_distance_h, innovation_iul_Gamma) |>
  summarise(gap = (num_adopters[topology == "NET"] - num_adopters[topology == "DP"]) / N_nodes[1], .groups = "drop") |>
  group_by(lambda) |> summarise(gap_pp = 100 * mean(gap), .groups = "drop")
print(as.data.frame(gap), digits = 3)
pooled <- fit_one(results |> mutate(lambda_f = factor(lambda)), pooled = TRUE)
write.csv(bind_rows(per_lambda |> mutate(lambda = as.character(lambda)), pooled |> mutate(lambda = "pooled")),
          file.path(DATA_OUT, paste0("unified_premium_bam", SUFFIX, ".csv")), row.names = FALSE)

# ---------------------------------- plots -------------------------------------
agg <- results |>
  group_by(topology, lambda, Gamma = innovation_iul_Gamma, h = social_distance_h) |>
  summarise(phi = mean(num_adopters / N_nodes), .groups = "drop")
heat <- function(df, title) {
  ggplot(df, aes(Gamma, h, fill = phi)) +
    geom_tile(color = "white", lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis", n.breaks = 4, na.value = "grey90", name = "Adoption") +
    labs(title = title, x = "IUL (Γ)", y = "MSP (h)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title = element_text(size = 10, face = "bold"))
}
cairo_pdf(file.path(PLOTS_OUT, paste0("unified_premium_heatmaps", SUFFIX, ".pdf")), width = 11, height = 4.2, onefile = TRUE)
for (l in LAMBDAS) {
  print(heat(filter(agg, lambda == l, topology == "NET"), sprintf("%s %s — Φ_total (λ = %.1f)", SURVEY, LEVEL, l)) +
        heat(filter(agg, lambda == l, topology == "DP"), sprintf("%s-DP — Φ_total (λ = %.1f)", SURVEY, l)) +
        plot_layout(nrow = 1, guides = "collect"))
}
invisible(dev.off())
message(format(Sys.time(), "%H:%M:%S"), "  DONE (", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min)")
