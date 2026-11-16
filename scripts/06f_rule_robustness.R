# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SUBSTITUTES OR COMPLEMENTS? Does the premium depend on how quality and social
# proof are combined in the adoption rule? 06f_rule_robustness.R  (~3 min)
#
# The main rule ADDS them (substitutes):     Gamma + lambda * E_i >= q_i.
# Two complementary alternatives, base exposure, the networks, twins, seeds and
# frozen gate of 05:
#
# Level 1, complements in shares (no new simulation):
#     Gamma + (lambda0 + phi * Gamma) * E_i >= q_i
#   for each Gamma this is the main rule with lambda(Gamma) = lambda0 + phi * Gamma,
#   so it is a diagonal path through the (Gamma, lambda) grid that 05 already swept.
#   Adoption at lambda(Gamma) is interpolated linearly between the two nearest swept
#   lambdas; lambda(Gamma) is clamped to the swept range [0.5, 1.0].
# Level 2, Tur's gamma form (Tur et al. 2018, 2024), a new sweep:
#     Gamma * A_i^gamma >= q_i  <=>  A_i >= c*_i = (q_i / Gamma)^(1/gamma)  <=>  E_i >= c*_i / k_i
#   (A_i = adopting contacts with an open tie, k_i = degree). It is one more vector of
#   node thresholds, so run_exposure_grid() runs it unchanged with Gamma = 0, lambda = 1
#   and "q" = c*_i / k_i. Gamma = 0 is undefined in this form and left out. The simulation
#   seed of each cell is that of the main sweep at lambda = 0.7 (same gate coins).
# Checks first (stop on any difference): (1) the threshold trick reproduces the main
# rule cell by cell; (2) the gamma form of the engine equals an independent R cascade
# by counts, cell by cell; (3) six hand-solved cases on a star.
#
# Premium = mean over the plane (Gamma x h) of adoption in the network with homophily
# (NET) minus its DP twin; 95% interval by run bootstrap (B = 2000, percentile).
#
# Outputs: output/06f_rule_robustness/<SET>/gamma_<seeding>.rds   (sweep, gitignored)
#          output/06f_rule_robustness/{level1_paths,gamma_premium}.csv
#          output/06f_rule_robustness/{table_level1,table_gamma}.tex   (random seeding)
#          plots/06f_rule_robustness/{rule_robustness,gamma_heatmaps}.pdf
# Usage (repo root, after 05 for every set and seeding):
#   Rscript scripts/06f_rule_robustness.R
#   Rscript scripts/06f_rule_robustness.R --sets=GSS_x2,ATP_x2 --seedings=random
#   Rscript scripts/06f_rule_robustness.R --analysis-only     # skip the gamma sweep
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(igraph); library(parallel); library(data.table); library(ggplot2) }))
Sys.setenv(OMP_NUM_THREADS = "1")
source("scripts/05_unified_engine.R")

ARGS <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) sub(".*=", "", a[1]) else d }
SETS <- strsplit(arg("sets", "GSS_x2,ATP_x2,GSS_x1,ATP_x1"), ",")[[1]]
SEEDINGS <- strsplit(arg("seedings", "random,central,marginal,closeness,eigen"), ",")[[1]]
N_CORES <- as.integer(arg("cores", 8)); ANALYSIS_ONLY <- "--analysis-only" %in% ARGS
OUT <- "output/06f_rule_robustness"; PLOTS <- "plots/06f_rule_robustness"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
N_RUNS <- 24L; B <- 2000L
G_SWEEP <- seq(0.025, 1, by = 0.025); H_SWEEP <- seq(0, 1, by = 1/24)
GAMMAS <- c(0.2, 0.25, 1/3, 0.5, 0.75, 1)
LAMBDAS <- seq(0.5, 1, by = 0.1)
LAM_SEED <- 0.7
TOL <- 1e-9          # A_i >= c* with an integer c* must not fail by rounding

# ---- the gamma form on the engine ------------------------------------------------
c_star <- function(q, Gamma, gamma) (q / Gamma)^(1 / gamma)
tau_gamma <- function(q, deg, Gamma, gamma) { t <- pmax(0, c_star(q, Gamma, gamma) - TOL) / deg; t[!is.finite(t)] <- 2; t }
run_tau_grid <- function(top, W, ctx, tau_of_G, seeds, run, topo, H, G) {   # tau_of_G(Gamma) -> node thresholds
  S <- ctx$none; P <- ctx$pairs_none; ss <- sim_seed_matrix(run, LAM_SEED, topo, H, G)
  vapply(seq_along(G), function(gi)
    run_exposure_grid(top$src, top$dst, W, S$slot, S$slot_n, S$slot_ego, 0, 0, P$tri_ptr, P$tri_f, 0, P$eN,
                      tau_of_G(G[gi]), 0, 1, as.integer(seeds), ss[, gi, drop = FALSE], 5000L)[, "num_adopters"],
    integer(length(H)))                                                    # nH x nG
}
cascade_gamma_R <- function(top, q, n, seeds, h, Gamma, gamma, sim_seed) {   # independent reference, by counts
  cs <- c_star(q, Gamma, gamma); w <- sigmoid_gate(top$edge_d, h)
  set.seed(sim_seed); E <- length(w) %/% 2L; coin <- runif(E) <= w[seq_len(E)]; open <- c(coin, coin)
  a <- rep(FALSE, n); a[seeds] <- TRUE
  repeat {
    idx <- which(a[top$src] & !a[top$dst] & open); cnt <- tabulate(top$dst[idx], nbins = n)
    new <- !a & cnt >= 1L & cnt >= cs - TOL
    if (!any(new)) break
    a[new] <- TRUE
  }
  sum(a)
}

# ---- 0. checks ----------------------------------------------------------------------
message(format(Sys.time(), "%H:%M:%S"), "  checks")
Hc <- seq(0, 1, by = 1/8); Gc <- c(0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9); bad <- 0L
for (SET in intersect(c("GSS_x2", "ATP_x2"), SETS)) {
  nets <- readRDS(file.path("data/02_calibrated_networks", paste0(SET, ".rds")))
  for (run in c(1L, 7L)) {
    nt <- make_pair(nets[[run]], run, "random")
    for (topo in c("NET", "DP")) {
      top <- nt[[topo]]; seeds <- draw_seeds(top, run, topo, "random")
      W <- gate_matrix(top, Hc); ctx <- build_exposure_ctx(top, run, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
      lin <- run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, nt$q, seeds, run, LAM_SEED, topo, Hc, Gc)[, "num_adopters"]
      trk <- run_tau_grid(top, W, ctx, function(G) pmax(0, (nt$q - G) / LAM_SEED), seeds, run, topo, Hc, Gc)
      d1 <- sum(as.vector(t(trk)) != lin)
      ss <- sim_seed_matrix(run, LAM_SEED, topo, Hc, Gc); d2 <- 0L
      for (gam in c(0.2, 0.33, 0.5, 1)) {
        eng <- run_tau_grid(top, W, ctx, function(G) tau_gamma(nt$q, top$deg, G, gam), seeds, run, topo, Hc, Gc)
        for (hi in seq_along(Hc)) for (gi in seq_along(Gc))
          d2 <- d2 + (cascade_gamma_R(top, nt$q, nt$n, seeds, Hc[hi], Gc[gi], gam, ss[hi, gi]) != eng[hi, gi])
      }
      message(sprintf("  %s run %2d %-3s  trick vs main rule: %d/%d differ | gamma engine vs R: %d/%d differ",
                      SET, run, topo, d1, length(lin), d2, 4L * length(Hc) * length(Gc)))
      bad <- bad + d1 + d2
    }
  }
}
star <- make_star(11, mode = "undirected"); el <- as_edgelist(star, names = FALSE)   # ego 1, k = 10 leaves
st <- list(src = as.integer(c(el[, 1], el[, 2])), dst = as.integer(c(el[, 2], el[, 1])), edge_d = rep(0, 20), deg = as.integer(degree(star)))
qs <- c(0.6, rep(1e6, 10)); Ss <- make_slots(st); Ps <- no_pairs(st)
hand <- function(nseed, Gam, gam) run_exposure_grid(st$src, st$dst, matrix(1, 10, 1), Ss$slot, Ss$slot_n, Ss$slot_ego, 0, 0,
  Ps$tri_ptr, Ps$tri_f, 0, Ps$eN, pmin(2, tau_gamma(qs, st$deg, Gam, gam)), 0, 1, as.integer(1 + seq_len(nseed)),
  matrix(1L, 1, 1), 100L)[1, "num_adopters"] - nseed
h3 <- c(hand(3, 0.3, 0.5), hand(4, 0.3, 0.5), hand(1, 0.3, 1), hand(2, 0.3, 1), hand(1, 0.6, 0.5), hand(10, 0.05, 0.5))
message("  by hand (c* = 4, 4, 2, 2, q <= Gamma, 144 > k; expected 0 1 0 1 1 0): ", paste(h3, collapse = " "))
if (bad > 0 || !identical(as.integer(h3), c(0L, 1L, 0L, 1L, 1L, 0L))) stop("checks failed")
message("  all checks passed")

# ---- 1. the gamma sweep ---------------------------------------------------------------
gamma_file <- function(set, sd) file.path(OUT, set, paste0("gamma_", sd, ".rds"))
if (!ANALYSIS_ONLY) for (SET in SETS) {
  survey <- sub("_.*", "", SET); level <- sub(".*_", "", SET)
  nets <- readRDS(file.path("data/02_calibrated_networks", paste0(SET, ".rds")))
  dp_nets <- if (level == "x2") nets else readRDS(file.path("data/02_calibrated_networks", paste0(survey, "_x2.rds")))
  dir.create(file.path(OUT, SET), showWarnings = FALSE)
  for (sd in SEEDINGS) {
    t0 <- Sys.time()
    out <- mclapply(seq_len(N_RUNS), function(run) tryCatch({
      nt <- make_pair(nets[[run]], run, sd, dp_from = dp_nets[[run]])
      arr <- array(0L, dim = c(length(H_SWEEP) * length(G_SWEEP), 2L, length(GAMMAS))); nseed <- integer(2)
      for (ti in 1:2) {
        topo <- c("NET", "DP")[ti]; top <- nt[[topo]]
        seeds <- draw_seeds(top, run, topo, sd); nseed[ti] <- length(seeds)
        W <- gate_matrix(top, H_SWEEP); ctx <- build_exposure_ctx(top, run, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
        for (ki in seq_along(GAMMAS))
          arr[, ti, ki] <- as.vector(t(run_tau_grid(top, W, ctx, function(G) tau_gamma(nt$q, top$deg, G, GAMMAS[ki]),
                                                    seeds, run, topo, H_SWEEP, G_SWEEP)))     # rows h-major
      }
      imm <- outer(G_SWEEP, GAMMAS, Vectorize(function(G, gm) mean(c_star(nt$q, G, gm) - TOL > nt$NET$deg_raw)))
      list(arr = arr, nseed = nseed, n = nt$n, imm = imm)
    }, error = function(e) paste("ERROR", conditionMessage(e))), mc.cores = N_CORES, mc.preschedule = FALSE)
    err <- Filter(is.character, out); if (length(err)) stop(err[[1]])
    res <- simplify2array(lapply(out, `[[`, "arr")); dimnames(res) <- list(NULL, c("NET", "DP"), sprintf("%.3f", GAMMAS), NULL)
    saveRDS(list(set = SET, seeding = sd, H = H_SWEEP, G = G_SWEEP, gammas = GAMMAS, N = out[[1]]$n, lam_seed = LAM_SEED,
                 seeds = t(vapply(out, `[[`, integer(2), "nseed")), immune = simplify2array(lapply(out, `[[`, "imm")), res = res),
            gamma_file(SET, sd))
    message(sprintf("%s  gamma sweep %s %-9s %.1f min", format(Sys.time(), "%H:%M:%S"), SET, sd,
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

# ---- 2. analysis ------------------------------------------------------------------------
boot_ci <- function(per_run) {
  set.seed(20261110); idx <- matrix(sample.int(length(per_run), length(per_run) * B, TRUE), B)
  bs <- apply(idx, 1, function(i) mean(per_run[i])); 100 * c(quantile(bs, 0.025), quantile(bs, 0.975))
}
# level 1
PATHS <- rbind(data.frame(path = sprintf("lambda = %.1f", LAMBDAS), lambda0 = LAMBDAS, phi = 0, kind = "constant"),
               data.frame(path = c("0.5 + 0.25 G", "0.5 + 0.5 G", "0.5 + G", "0.5 + 2 G"), lambda0 = 0.5,
                          phi = c(0.25, 0.5, 1, 2), kind = "complements"),
               data.frame(path = "1 - 0.5 G", lambda0 = 1, phi = -0.5, kind = "reverse"))
lvl1 <- rbindlist(lapply(SETS, function(SET) rbindlist(lapply(SEEDINGS, function(sd) {
  d <- fread(file.path("output/05_unified_diffusion", SET, paste0("unified_premium_results_", sd, ".csv")),
             select = c("topology", "lambda", "run_id", "social_distance_h", "innovation_iul_Gamma", "num_adopters", "N_nodes"))
  w <- dcast(d, lambda + run_id + social_distance_h + innovation_iul_Gamma ~ topology, value.var = "num_adopters")
  w[, gap := (NET - DP) / d$N_nodes[1]]
  a <- xtabs(gap ~ run_id + G + lambda, w[, .(gap = mean(gap)), by = .(run_id, G = innovation_iul_Gamma, lambda)])
  Gs <- as.numeric(dimnames(a)[[2]])
  rbindlist(lapply(seq_len(nrow(PATHS)), function(p) {
    lam <- pmin(1, pmax(0.5, PATHS$lambda0[p] + PATHS$phi[p] * Gs))
    lo <- pmin(length(LAMBDAS) - 1L, floor((lam - 0.5) / 0.1 + 1e-9) + 1L); wt <- (lam - LAMBDAS[lo]) / 0.1
    per_run <- vapply(seq_len(dim(a)[1]), function(r)
      mean((1 - wt) * a[cbind(r, seq_along(Gs), lo)] + wt * a[cbind(r, seq_along(Gs), lo + 1L)]), 0)
    ci <- boot_ci(per_run)
    data.table(set = SET, seeding = sd, PATHS[p, ], gap_pp = 100 * mean(per_run), ci_lo = ci[1], ci_hi = ci[2],
               runs_positive = sum(per_run > 0))
  }))
}))))
fwrite(lvl1, file.path(OUT, "level1_paths.csv"))
# level 2
rows <- list(); cells <- list()
for (SET in SETS) for (sd in SEEDINGS) {
  x <- readRDS(gamma_file(SET, sd)); n <- x$N
  for (gi in seq_along(x$gammas)) {
    diff <- (x$res[, "NET", gi, ] - x$res[, "DP", gi, ]) / n; per_run <- colMeans(diff); ci <- boot_ci(per_run)
    rows[[length(rows) + 1]] <- data.table(set = SET, seeding = sd, gamma = x$gammas[gi],
      phi_NET = 100 * mean(x$res[, "NET", gi, ] / n), phi_DP = 100 * mean(x$res[, "DP", gi, ] / n),
      gap_pp = 100 * mean(per_run), ci_lo = ci[1], ci_hi = ci[2], runs_positive = sum(per_run > 0),
      pct_cells_DP_wins = 100 * mean(rowMeans(diff) < 0), immune_pct = 100 * mean(x$immune[, gi, ]))
    if (sd == "random") cells[[length(cells) + 1]] <- data.table(set = SET, gamma = x$gammas[gi],
      h = rep(x$H, each = length(x$G)), G = rep(x$G, times = length(x$H)), gap = 100 * rowMeans(diff))
  }
}
gm <- rbindlist(rows); fwrite(gm, file.path(OUT, "gamma_premium.csv"))

# ---- 3. tables (random seeding) and figures -------------------------------------------
fmt <- function(d) d[, sprintf("%.1f [%.1f, %.1f]", gap_pp, ci_lo, ci_hi)]
hdr <- paste0(" & ", paste(sub("_x", " $\\\\times", SETS), collapse = "$ & "), "$ \\\\")
tab <- function(first, labs, cellfun, groups) {
  body <- vapply(seq_along(labs), function(i) paste0(labs[i], " & ", paste(vapply(SETS, cellfun, "", i = i), collapse = " & "), " \\\\"), "")
  if (!is.null(groups)) body <- unlist(lapply(split(body, groups), function(b) c(b, "\\midrule")), use.names = FALSE)[-sum(lengths(split(body, groups))) - length(unique(groups))]
  c(sprintf("\\begin{tabular}{@{}l%s@{}}", strrep("c", length(SETS))), "\\toprule", paste0(first, hdr), "\\midrule", body, "\\bottomrule", "\\end{tabular}")
}
l1lab <- c(sprintf("$\\lambda = %.1f$", LAMBDAS), "$0.5 + 0.25\\,\\Gamma$", "$0.5 + 0.5\\,\\Gamma$", "$0.5 + \\Gamma$",
           "$0.5 + 2\\,\\Gamma$", "$1 - 0.5\\,\\Gamma$")
writeLines(tab("$\\lambda(\\Gamma)$", l1lab, function(s, i) fmt(lvl1[seeding == "random" & set == s & path == PATHS$path[i]]),
               factor(PATHS$kind, unique(PATHS$kind))), file.path(OUT, "table_level1.tex"))
writeLines(tab("$\\gamma$", sprintf("%.2f", GAMMAS), function(s, i) fmt(gm[seeding == "random" & set == s][i]), NULL),
           file.path(OUT, "table_gamma.tex"))

fig <- rbind(lvl1[, .(set, seeding, x = lambda0 + 0.5 * phi, gap_pp, ci_lo, ci_hi, panel = "Level 1: lambda(Gamma) paths", kind)],
             gm[, .(set, seeding, x = gamma, gap_pp, ci_lo, ci_hi, panel = "Level 2: gamma form", kind = "gamma")])
fig[, set := factor(set, SETS)]
p1 <- ggplot(lvl1[seeding == "random"][, set := factor(set, SETS)],
             aes(factor(path, PATHS$path), gap_pp, colour = kind)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.25) +
  facet_wrap(~set, nrow = 1) + coord_flip() +
  labs(x = NULL, y = "Structural premium (pp)", colour = NULL,
       title = "Level 1: complements in shares, Gamma + (lambda0 + phi Gamma) E >= q (random seeding)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9))
p2 <- ggplot(gm[, set := factor(set, SETS)], aes(gamma, gap_pp, colour = seeding)) +
  geom_hline(yintercept = 0, linewidth = 0.3) + geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.015, alpha = 0.6) +
  geom_line() + geom_point(size = 1.1) + facet_wrap(~set, nrow = 1) + scale_x_continuous(breaks = c(0.2, 0.33, 0.5, 0.75, 1)) +
  labs(x = expression(gamma), y = "Structural premium (pp)", colour = "Seeding",
       title = "Level 2: Tur's gamma form, Gamma * A^gamma >= q (95% run-bootstrap intervals)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9))
cairo_pdf(file.path(PLOTS, "rule_robustness.pdf"), width = 10, height = 4, onefile = TRUE); print(p1); print(p2); invisible(dev.off())
ggsave(file.path(PLOTS, "gamma_premium.pdf"), p2, width = 10, height = 3.6, device = cairo_pdf)
cl <- rbindlist(cells)[set %in% c("GSS_x2", "ATP_x2")]
h <- ggplot(cl, aes(G, h, fill = gap)) + geom_raster() + facet_grid(set ~ sprintf("gamma = %.2f", gamma)) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, name = "NET - DP (pp)") +
  labs(x = expression(Gamma ~ "(IUL)"), y = "h (MSP)", title = "Gamma form: premium by cell (random seeding)") +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(PLOTS, "gamma_heatmaps.pdf"), h, width = 11, height = 4.8, device = cairo_pdf)

message("\n=== premium (pp), random seeding ===")
print(dcast(lvl1[seeding == "random", .(set, path, v = fmt(.SD))], path ~ set, value.var = "v"))
print(dcast(gm[seeding == "random", .(set, gamma = sprintf("%.2f", gamma), v = fmt(.SD))], gamma ~ set, value.var = "v"))
message(sprintf("Level 1 minimum over sets, seedings and paths: %.2f pp (lowest lower bound %.2f)", min(lvl1$gap_pp), min(lvl1$ci_lo)))
message(sprintf("Level 2 minimum over sets, seedings and gammas: %.2f pp (lowest lower bound %.2f); premium positive in all 24 networks in %d of %d",
                min(gm$gap_pp), min(gm$ci_lo), sum(gm$runs_positive == N_RUNS), nrow(gm)))
