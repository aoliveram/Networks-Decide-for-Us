# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK OF THE Rcpp ENGINE (scripts/05_unified_engine.R)
#
#   1. annealed gate: run_unified_cpp vs 05's former R loop, the five outputs of
#      05, cell by cell
#   2. frozen gate:   run_unified_cpp vs an independent R cascade (adopters, steps)
#   3. the sweep engine (run_exposure_grid) with the base exposure vs
#      run_unified_cpp with the frozen gate: adopters, rational adopters, steps
#      and t50, cell by cell
#   4. the exposure formula on a star ego: the engine's exposure, recovered by
#      bisection on q, equals the values computed by hand (theta and kappa)
#   5. every condition of EXPOSURE_CONDITIONS vs an independent R
#      implementation (full recomputation every step), adopters cell by cell
#   5b. PR #5: on ELSOC network 1 (if data/02_calibrated_networks/ELSOC_x2.rds exists),
#      the conditions whose inputs come from the network's own attributes (Gower on
#      age, education, sex, religion and ideology; education x religion circles;
#      ideology circles) vs the same R implementation
#   6. time per simulation
#
# Networks GSS 001..004 of 02_GSS_network_ergm (tracked, so the check does not
# depend on the generator) and their DP twins (make_pair: set.seed(500000 + k),
# 20 x m swaps); random seeding; grid h = 0..1 by 1/8, Gamma = 0..1 by 0.1,
# lambda = 0.5 / 0.7 / 1.0 (2,376 cells per check in 1-3). One core, ~1 min.
#
# Usage (repo root): Rscript scripts/05_engine_check.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(igraph); library(intergraph) }))
source("scripts/05_unified_engine.R")

H <- seq(0, 1, by = 1/8); G <- seq(0, 1, by = 0.1); L <- c(0.5, 0.7, 1.0)
grid <- expand.grid(h = H, Gamma = G, lambda = L)
load_net <- function(k) asIgraph(readRDS(sprintf("data/02_GSS_network_ergm/GSS_net_sim_1000_%03d.rds", k)))
pairs <- lapply(1:4, function(k) make_pair(load_net(k), k))

# ---- 1-3: identities ---------------------------------------------------------------------------
res <- list(); tm <- c(R = 0, cpp = 0, frozen_R = 0, frozen_cpp = 0, grid = 0); n_cells <- 0L; n_grid <- 0L
grid_same <- c(adopters = 0L, rational = 0L, steps = 0L, t50 = 0L)
for (k in 1:4) {
  nt <- pairs[[k]]; q <- nt$q; n <- nt$n
  for (topo in c("NET", "DP")) {
    top <- nt[[topo]]
    seeds <- draw_seeds(top, k, topo)
    for (i in seq_len(nrow(grid))) {
      x <- grid[i, ]; s <- sim_seed_of(k, x$h, x$Gamma, x$lambda, topo)
      t0 <- proc.time()[[3]]; rR <- run_unified_R(top, q, n, seeds, x$h, x$Gamma, x$lambda, s)
      t1 <- proc.time()[[3]]; rC <- run_cell(top, q, seeds, x$h, x$Gamma, x$lambda, s, "annealed")
      t2 <- proc.time()[[3]]; fR <- cascade_frozen_R(top, q, n, seeds, x$h, x$Gamma, x$lambda, s)
      t3 <- proc.time()[[3]]; fC <- run_cell(top, q, seeds, x$h, x$Gamma, x$lambda, s, "frozen")
      t4 <- proc.time()[[3]]
      tm <- tm + c(t1 - t0, t2 - t1, t3 - t2, t4 - t3, 0); n_cells <- n_cells + 1L
      res[[length(res) + 1L]] <- c(
        annealed_same = all(rR == rC[names(rR)]),
        frozen_same   = all(fR == fC[c("num_adopters", "num_steps")]),
        frozen_fixed_point = fC[["stop_reason"]] %in% c(1L, 2L, 3L))
    }
    # 3: the sweep engine, whole grid per lambda, against run_cell(frozen)
    W <- gate_matrix(top, H); ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    for (lam in L) {
      t0 <- proc.time()[[3]]
      r <- run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, q, seeds, k, lam, topo, H, G)
      tm[["grid"]] <- tm[["grid"]] + proc.time()[[3]] - t0
      for (hi in seq_along(H)) for (gi in seq_along(G)) {
        fC <- run_cell(top, q, seeds, H[hi], G[gi], lam, sim_seed_of(k, H[hi], G[gi], lam, topo), "frozen")
        row <- r[(hi - 1L) * length(G) + gi, ]
        grid_same <- grid_same + c(row[["num_adopters"]] == fC[["num_adopters"]],
                                   row[["num_adopted_rational"]] == fC[["num_adopted_rational"]],
                                   row[["num_steps"]] == fC[["num_steps"]], row[["t50"]] == fC[["t50"]])
        n_grid <- n_grid + 1L
      }
    }
  }
  message(sprintf("network %d done", k))
}
ok <- colSums(do.call(rbind, res))
message(sprintf("[1] annealed: %d of %d cells identical to 05's former R loop", ok[["annealed_same"]], n_cells))
message(sprintf("[2] frozen:   %d of %d cells identical to the independent R cascade; %d of %d end at a fixed point",
                ok[["frozen_same"]], n_cells, ok[["frozen_fixed_point"]], n_cells))
message(sprintf("[3] sweep engine, base exposure: identical in %d (adopters), %d (rational), %d (steps), %d (t50) of %d cells",
                grid_same[["adopters"]], grid_same[["rational"]], grid_same[["steps"]], grid_same[["t50"]], n_grid))

# ---- 4: formula on a star ego (hand values) -------------------------------------------------------
# ego 1; contacts 2..5 = circle A, 6..9 = circle B; two adopting contacts, gate open.
star <- function(with_pairs) {
  g <- make_empty_graph(9, directed = FALSE) + edges(as.vector(rbind(1, 2:9)))
  if (with_pairs) g <- g + edges(c(combn(2:5, 2), combn(6:9, 2)))   # each circle is a clique
  el <- as_edgelist(g, names = FALSE)
  list(g = g, src = as.integer(c(el[, 1], el[, 2])), dst = as.integer(c(el[, 2], el[, 1])), deg = pmax(degree(g), 1L))
}
E_of <- function(top, lab, seeds, theta, kp) {
  slots <- make_slots(top, lab)
  P <- if (kp > 0) pairs_csr(top, pairs_real(top)) else no_pairs(top)
  W <- matrix(1, length(top$src) %/% 2L, 1L); lo <- 0; hi <- 1
  for (it in 1:40) {                      # the largest q of ego 1 that still adopts = E_1
    mid <- (lo + hi) / 2; q <- c(mid, rep(2, 8))
    r <- run_exposure_grid(as.integer(top$src), as.integer(top$dst), W, slots$slot, slots$slot_n, slots$slot_ego,
                           theta, 0, P$tri_ptr, P$tri_f, kp, P$eN, q, 0, 1, as.integer(seeds), matrix(1L, 1, 1), 50L)
    if (r[1, "num_adopters"] > length(seeds)) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}
ts <- star(FALSE); tp <- star(TRUE)
lab_star <- function(top) { l <- ifelse(top$src %in% 2:5, 1L, 2L); l[top$dst != 1] <- top$src[top$dst != 1] + 100L; l }
# hand values: two circles of 4 with A_theta(4) = 4 / (1 + 3 theta) each in the denominator;
# kappa_p: each circle has 6 pairs, so e_N = 12 and e_a = 1 if both adopters share a circle
star_cases <- list(
  list("base, same circle",            ts, lab_star(ts), c(2, 3), 0,   0,   2 / 8),
  list("theta 0.5, same circle",       ts, lab_star(ts), c(2, 3), 0.5, 0,   (2 / 1.5) / (2 * 4 / 2.5)),
  list("theta 0.5, two circles",       ts, lab_star(ts), c(2, 6), 0.5, 0,   2 / (2 * 4 / 2.5)),
  list("theta 1, same circle",         ts, lab_star(ts), c(2, 3), 1,   0,   1 / 2),
  list("kappa 0.5, connected pair",    tp, NULL,         c(2, 3), 0,   0.5, (2 + 1) / (8 + 12)),
  list("kappa 0.5, unconnected pair",  tp, NULL,         c(2, 6), 0,   0.5, 2 / (8 + 12)))
star_ok <- TRUE
for (cs in star_cases) {
  e <- E_of(cs[[2]], cs[[3]], cs[[4]], cs[[5]], cs[[6]]); star_ok <- star_ok && abs(e - cs[[7]]) < 1e-9
  message(sprintf("[4] %-30s engine %.6f  hand %.6f", cs[[1]], e, cs[[7]]))
}

# ---- 5: every condition vs the independent R implementation --------------------------------------
set.seed(12); cells <- cbind(sample(length(H), 24, TRUE), sample(length(G), 24, TRUE))
lam <- 0.7; k <- 1L; nt <- pairs[[k]]; ref_ok <- TRUE; n_ref <- 0L; n_ref_same <- 0L
for (topo in c("NET", "DP")) {
  top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo); W <- gate_matrix(top, H)
  ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP")); ssm <- sim_seed_matrix(k, lam, topo, H, G)
  E <- length(top$src) %/% 2L
  for (cn in names(EXPOSURE_CONDITIONS)) {
    cd <- EXPOSURE_CONDITIONS[[cn]]
    r <- run_condition(top, W, ctx, cd, nt$q, seeds, k, lam, topo, H, G)
    got <- ref <- integer(nrow(cells))
    for (j in seq_len(nrow(cells))) {
      hi <- cells[j, 1]; gi <- cells[j, 2]; set.seed(ssm[hi, gi]); o <- as.numeric(runif(E) <= W[, hi])
      ref[j] <- exposure_cascade_R(top, ctx[[cd$circles]], ctx[[paste0("pairs_", cd$pairs)]], cd$theta, cd$kappa_c,
                                   cd$kappa_p, nt$q, G[gi], lam, seeds, c(o, o))
      got[j] <- r[(hi - 1L) * length(G) + gi, "num_adopters"]
    }
    n_ref <- n_ref + length(got); n_ref_same <- n_ref_same + sum(got == ref)
    message(sprintf("[5] %-4s %-15s c = %+.2f: %d/%d cells identical (adopters %d-%d)",
                    topo, cn, cd$c, sum(got == ref), length(got), min(ref), max(ref)))
  }
}

# ---- 5b: ELSOC network (graph attributes blau and cell; the ideology theta arm) -----------------------
el_file <- "data/02_calibrated_networks/ELSOC_x2.rds"; n_el <- 0L; n_el_same <- 0L
if (file.exists(el_file)) {
  ge <- readRDS(el_file)[[1]]; nt <- make_pair(ge, 1L)
  stopifnot(identical(ge$blau, c("age", "educ_num", "sex", "relig", "ideol")), !is.null(nt$DP$g$cell))
  CE <- exposure_conditions(ge)[c("base", "theta_1", "theta_1_ideol", "theta_033_ideol", "theta_1_ideol_perm")]
  for (topo in c("NET", "DP")) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, 1L, topo); W <- gate_matrix(top, H)
    ctx <- build_exposure_ctx(top, 1L, as.integer(topo == "DP"), CE); ssm <- sim_seed_matrix(1L, lam, topo, H, G)
    E <- length(top$src) %/% 2L
    for (cn in names(CE)) {
      cd <- CE[[cn]]; r <- run_condition(top, W, ctx, cd, nt$q, seeds, 1L, lam, topo, H, G)
      got <- ref <- integer(nrow(cells))
      for (j in seq_len(nrow(cells))) {
        hi <- cells[j, 1]; gi <- cells[j, 2]; set.seed(ssm[hi, gi]); o <- as.numeric(runif(E) <= W[, hi])
        ref[j] <- exposure_cascade_R(top, ctx[[cd$circles]], ctx[[paste0("pairs_", cd$pairs)]], cd$theta, cd$kappa_c,
                                     cd$kappa_p, nt$q, G[gi], lam, seeds, c(o, o))
        got[j] <- r[(hi - 1L) * length(G) + gi, "num_adopters"]
      }
      n_el <- n_el + length(got); n_el_same <- n_el_same + sum(got == ref)
      message(sprintf("[5b] ELSOC %-4s %-20s c = %+.2f: %d/%d cells identical (adopters %d-%d)",
                      topo, cn, cd$c, sum(got == ref), length(got), min(ref), max(ref)))
    }
  }
  # the Gower distance of the ELSOC network uses its own Blau space (no race)
  D <- gower_of(ge); i <- 1L; j <- 2L; va <- vertex_attr(ge)
  hand <- mean(c(abs(va$age[i] - va$age[j]) / diff(range(va$age)), abs(va$educ_num[i] - va$educ_num[j]) / diff(range(va$educ_num)),
                 va$sex[i] != va$sex[j], va$relig[i] != va$relig[j], va$ideol[i] != va$ideol[j]))
  message(sprintf("[5b] ELSOC Gower d(1, 2): engine %.6f, by hand %.6f", D[i, j], hand))
  n_el <- n_el + 1L; n_el_same <- n_el_same + (abs(D[i, j] - hand) < 1e-12)
} else message("[5b] skipped: ", el_file, " not found")

message(sprintf("[6] ms per simulation (one core): R loop %.1f, Rcpp %.2f (annealed); R cascade %.1f, Rcpp %.2f (frozen); sweep engine %.2f (frozen, base)",
                1000 * tm[["R"]] / n_cells, 1000 * tm[["cpp"]] / n_cells,
                1000 * tm[["frozen_R"]] / n_cells, 1000 * tm[["frozen_cpp"]] / n_cells, 1000 * tm[["grid"]] / n_grid))
all_ok <- ok[["annealed_same"]] == n_cells && ok[["frozen_same"]] == n_cells && all(grid_same == n_grid) &&
  star_ok && n_ref_same == n_ref && n_el_same == n_el
message(if (all_ok) "ALL ENGINE CHECKS PASSED" else "SOME ENGINE CHECKS FAILED")
if (!all_ok) stop("the Rcpp engine does not match its reference")
