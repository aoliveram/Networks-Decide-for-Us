# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTENDED sweep of the unified (bottom-up) adoption rule — Option 1 of
# MAKING-IT-BOTTOM-UP.tex — sized for a ~1–1.5 h background run, single core.
#
#   adopt_i(t)  ⟺  Γ + λ·Ẽ_i(t) ≥ q_i   ∧   Ẽ_i(t) > 0
#
# Goal: replicate the engine's heatmap outputs (total / rational / social
# adoption over the IUL × MSP plane) under the unified rule, with λ as the
# new axis playing the role the (τ_mean, τ_sd) configurations play today.
#
# Design (mirrors scripts/04_GSS_diffusion_sims.R conventions):
#   - Topology: GSS plausible networks 001..096 (one run per network instance;
#     96 replicas per cell = 4x the engine's 24), parallel over runs on 8 cores
#   - Seeding: RANDOM only — 1 random primary + neighbor cluster of size
#     ≈ round(0.40 · degree(primary)) (reference τ_mean, constant across Γ)
#   - Gate: P(contact counts) = 1/(1+exp((d_ij − h)/0.02)), redrawn per step
#   - Exposure denominator = degree (paper Eq. 2)
#   - MAX_STEPS = 200 (engine value); patience cutoff 25 stalled steps
#
# Grid: 41 IUL × 13 MSP × 13 λ × 96 runs = 665,184 simulations.
#
# Checkpointing: results are saved per λ under data/unified_sweep/ — the script
# can be killed and relaunched and will resume, skipping completed λ values.
#
# Outputs:
#   output/lambda_calibration/unified_sweep/lambda_X.XX.rds   (checkpoints)
#   output/lambda_calibration/unified_sweep_results.csv       (combined)
#   output/lambda_calibration/unified_sweep_vs_current.csv    (fit table)
#   plots/lambda_calibration/unified_sweep_heatmaps.pdf     (1 page per λ:
#       total / rational / social triptych, engine-style; final comparison page)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(igraph); library(intergraph); library(cluster)
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(parallel)
})

# Avoid OpenMP conflicts with FORK (same guard as the engine)
Sys.setenv(OMP_NUM_THREADS = "1")

# ------------------------------- configuration -------------------------------
NETWORKS_DIR <- "data/02_GSS_network_ergm/"
DATA_OUT     <- "output/lambda_calibration"
CKPT_DIR     <- file.path(DATA_OUT, "checkpoints")
PLOTS_OUT    <- "plots/lambda_calibration"
CURRENT_RAW  <- "legacy/output/04_GSS_diffusion_sims/results_m03-06_sd0.12_random.rds"

N_RUNS       <- 96                        # network instances 001..096, 1 run each
N_CORES      <- 8                         # 8 P-cores, FORK workers (as the engine)
IUL_SWEEP    <- seq(0, 1, by = 0.025)     # 41 levels — identical to engine
H_SWEEP      <- seq(0, 1, by = 1/12)      # 13 levels — identical to engine
LAMBDAS      <- seq(0.3, 1.5, by = 0.1)   # 13 values of the social-proof coupling
GATE_T       <- 0.02
TAU_REF_SEED <- 0.40
MAX_STEPS    <- 200
PATIENCE     <- 25

dir.create(CKPT_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_OUT, showWarnings = FALSE, recursive = TRUE)

sigmoid_gate <- function(d, h, T = GATE_T) 1 / (1 + exp((d - h) / T))

# ------------------------- load networks & precompute -------------------------
message(format(Sys.time(), "%H:%M:%S"), "  Loading ", N_RUNS,
        " GSS networks and precomputing Gower distances...")
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
  el <- as_edgelist(g, names = FALSE)
  src <- c(el[, 1], el[, 2]); dst <- c(el[, 2], el[, 1])
  deg <- degree(g)
  nets[[k]] <- list(g = g, n = vcount(g), src = src, dst = dst,
                    edge_d = d_mat[cbind(src, dst)],
                    q = V(g)$mur_score, deg = pmax(deg, 1L), deg_raw = deg)
  if (k %% 12 == 0) message("  ...", k, "/", N_RUNS)
}

# ---- implied derived-threshold diagnostic (ported from the removed minimal
# ---- test script): tau* = max(0,(q - Gamma)/lambda) across all loaded nets
message("\nImplied derived-threshold tau* = max(0,(q-Gamma)/lambda):")
q_all <- unlist(lapply(nets, `[[`, "q"))
diag_tab <- expand.grid(lambda = c(0.5, 0.7, 1.0), Gamma = c(0.20, 0.30, 0.50)) |>
  rowwise() |>
  mutate(ts = list(pmax(0, (q_all - Gamma) / lambda)),
         tau_star_mean = mean(unlist(ts)), tau_star_sd = sd(unlist(ts)),
         pct_instigators = mean(unlist(ts) == 0) * 100,
         pct_immune = mean(unlist(ts) > 1) * 100) |>
  select(-ts) |> ungroup()
print(as.data.frame(diag_tab), digits = 3)
message("(reference: the legacy engine's main-text config is tau ~ Normal(0.40, 0.12))\n")

draw_seeds <- function(net, lam, run) {
  # mirror engine: random primary, cluster sized by reference τ · degree;
  # redrawn per (λ, run) as the engine redraws per (τ-config, run)
  set.seed(run * 2000 + round(lam * 100))
  primary <- sample.int(net$n, 1)
  n_seeds <- max(1, min(round(TAU_REF_SEED * net$deg_raw[primary]),
                        net$deg_raw[primary] + 1))
  seeds <- primary
  if (n_seeds > 1) {
    nb <- as.integer(neighbors(net$g, primary))
    if (length(nb) > 0)
      seeds <- c(primary, sample(nb, min(length(nb), n_seeds - 1)))
  }
  unique(seeds)
}

run_unified <- function(net, seeds, h, Gamma, lambda, sim_seed) {
  tau_star <- pmax(0, (net$q - Gamma) / lambda)
  w_edge   <- sigmoid_gate(net$edge_d, h)
  a <- rep(FALSE, net$n); a[seeds] <- TRUE
  set.seed(sim_seed)
  stall <- 0L; steps <- 0L
  for (t in seq_len(MAX_STEPS)) {
    idx <- which(a[net$src] & !a[net$dst])
    if (length(idx) == 0L) break
    pass <- runif(length(idx)) <= w_edge[idx]
    counts <- tabulate(net$dst[idx][pass], nbins = net$n)
    E <- counts / net$deg
    new <- !a & counts > 0L & E >= tau_star
    if (any(new)) { a[new] <- TRUE; stall <- 0L; steps <- t } else stall <- stall + 1L
    if (all(a) || stall >= PATIENCE) break
  }
  non_seed <- setdiff(which(a), seeds)
  c(num_adopters = sum(a),
    num_adopted_rational = sum(net$q[non_seed] <= Gamma),
    num_adopted_social   = sum(net$q[non_seed] >  Gamma),
    num_steps = steps, initial_cluster_size = length(seeds))
}

# --------------------------------- the sweep ----------------------------------
sims_per_lambda <- N_RUNS * length(H_SWEEP) * length(IUL_SWEEP)
message(format(Sys.time(), "%H:%M:%S"), "  Sweep: ", length(LAMBDAS),
        " lambdas x ", sims_per_lambda, " sims = ",
        length(LAMBDAS) * sims_per_lambda, " total simulations.")
t_all <- Sys.time()

for (lam in LAMBDAS) {
  ck <- file.path(CKPT_DIR, sprintf("lambda_%.2f.rds", lam))
  if (file.exists(ck)) {
    message(format(Sys.time(), "%H:%M:%S"), "  lambda=", sprintf("%.2f", lam),
            " checkpoint exists — skipping.")
    next
  }
  t_lam <- Sys.time()
  # Parallel over runs (FORK): nets[] is shared copy-on-write; every simulation
  # sets its own deterministic seed, so results are identical to the sequential
  # version regardless of worker scheduling.
  per_run <- mclapply(seq_len(N_RUNS), function(run) {
    net <- nets[[run]]
    seeds <- draw_seeds(net, lam, run)
    rows <- vector("list", length(H_SWEEP) * length(IUL_SWEEP)); ri <- 0L
    for (h in H_SWEEP) {
      for (Gamma in IUL_SWEEP) {
        ri <- ri + 1L
        sim_seed <- run * 1e6 + round(h * 12) * 1e4 +
          round(Gamma * 1000) + round(lam * 100) * 7L
        out <- run_unified(net, seeds, h, Gamma, lam, sim_seed)
        rows[[ri]] <- c(lambda = lam, run_id = run, social_distance_h = h,
                        innovation_iul_Gamma = Gamma, out, N_nodes = net$n)
      }
    }
    do.call(rbind, rows)
  }, mc.cores = N_CORES)
  bad <- vapply(per_run, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("Worker failure at lambda=", lam, " (runs: ",
                     paste(which(bad), collapse = ","), ")")
  df <- as.data.frame(do.call(rbind, per_run))
  saveRDS(df, ck)
  el_min <- as.numeric(difftime(Sys.time(), t_lam, units = "mins"))
  done <- sum(file.exists(file.path(CKPT_DIR, sprintf("lambda_%.2f.rds", LAMBDAS))))
  eta <- el_min * (length(LAMBDAS) - done)
  message(format(Sys.time(), "%H:%M:%S"), "  lambda=", sprintf("%.2f", lam),
          " done in ", round(el_min, 1), " min  (", done, "/", length(LAMBDAS),
          "; ~", round(eta), " min remaining)")
}

# ------------------------------ combine & save --------------------------------
results <- bind_rows(lapply(
  file.path(CKPT_DIR, sprintf("lambda_%.2f.rds", LAMBDAS)), readRDS))
write.csv(results, file.path(DATA_OUT, "unified_sweep_results.csv"), row.names = FALSE)
message(format(Sys.time(), "%H:%M:%S"), "  Combined results saved (",
        nrow(results), " rows; total ",
        round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min).")

agg <- results |>
  mutate(phi_total    = num_adopters / N_nodes,
         phi_rational = num_adopted_rational / N_nodes,
         phi_social   = num_adopted_social / N_nodes) |>
  group_by(lambda, Gamma = innovation_iul_Gamma, h = social_distance_h) |>
  summarise(across(c(phi_total, phi_rational, phi_social), mean),
            mean_steps = mean(num_steps), .groups = "drop")

# ------------------- comparison table vs the current model --------------------
fit_tab <- NULL
if (file.exists(CURRENT_RAW)) {
  raw <- readRDS(CURRENT_RAW)
  cur <- raw[["mean_0.40"]] |>
    mutate(G = round(innovation_iul_Gamma, 3), H = round(social_distance_h, 3)) |>
    group_by(G, H) |>
    summarise(phi_cur = mean(num_adopters / N_nodes_actual), .groups = "drop")
  fit_tab <- agg |>
    mutate(G = round(Gamma, 3), H = round(h, 3)) |>
    inner_join(cur, by = c("G", "H")) |>
    group_by(lambda) |>
    summarise(r = cor(phi_total, phi_cur),
              MAD = mean(abs(phi_total - phi_cur)), n = n(), .groups = "drop")
  write.csv(fit_tab, file.path(DATA_OUT, "unified_sweep_vs_current.csv"), row.names = FALSE)
  message("\nFit vs CURRENT model (tau~N(0.40,0.12), random seeding), by lambda:")
  print(as.data.frame(fit_tab), digits = 3)
}

# ---------------------------------- plots -------------------------------------
heat <- function(df, fill, title) {
  ggplot(df, aes(x = Gamma, y = h, fill = .data[[fill]])) +
    geom_tile(color = "white", lwd = 0.1) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis",
                         n.breaks = 4, na.value = "grey90", name = "Adoption") +
    labs(title = title, x = "IUL (Γ)", y = "MSP (h)") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))
}

cairo_pdf(file.path(PLOTS_OUT, "unified_sweep_heatmaps.pdf"),
          width = 11, height = 4.2, onefile = TRUE)
for (lam in LAMBDAS) {
  d <- filter(agg, lambda == lam)
  print(
    heat(d, "phi_total",    sprintf("Unified — Φ_total (λ = %.2f)", lam)) +
    heat(d, "phi_rational", "Φ_rational (q ≤ Γ, non-seed)") +
    heat(d, "phi_social",   "Φ_social (q > Γ, non-seed)") +
    plot_layout(nrow = 1, guides = "collect")
  )
}
if (!is.null(fit_tab)) {
  best <- fit_tab$lambda[which.min(fit_tab$MAD)]
  cur_plot <- raw[["mean_0.40"]] |>
    group_by(Gamma = innovation_iul_Gamma, h = social_distance_h) |>
    summarise(phi_total = mean(num_adopters / N_nodes_actual), .groups = "drop")
  print(
    heat(cur_plot, "phi_total", "CURRENT model — Φ_total (τ~N(0.40,0.12), random)") +
    heat(filter(agg, lambda == best), "phi_total",
         sprintf("Unified — Φ_total (best-fit λ = %.2f)", best)) +
    plot_layout(nrow = 1, guides = "collect")
  )
}
dev.off()
message(format(Sys.time(), "%H:%M:%S"), "  Saved ",
        file.path(PLOTS_OUT, "unified_sweep_heatmaps.pdf"), " — DONE.")
