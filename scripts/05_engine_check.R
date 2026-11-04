# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CHECK OF THE Rcpp ENGINE (scripts/05_unified_engine.R)
#
#   1. annealed gate: Rcpp vs 05's R loop, the five outputs of 05, cell by cell
#   2. frozen gate:   Rcpp vs an independent R cascade (adopters and steps)
#   3. time per simulation of each
#
# Networks GSS 001..004 and their DP twins, built exactly as in 05 (Gower
# distances, rewire keeping_degseq with set.seed(500000 + k)); seeds = a random
# node plus 27 of its contacts; grid h = 0..1 by 1/8, Gamma = 0..1 by 0.1,
# lambda = 0.5 / 0.7 / 1.0 (2,376 cells per gate). One core, about 2 minutes.
#
# Usage (repo root): Rscript scripts/05_engine_check.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(igraph); library(intergraph); library(cluster) }))
source("scripts/05_unified_engine.R")

edge_arrays <- function(g, d_mat) {
  el <- as_edgelist(g, names = FALSE)
  src <- c(el[, 1], el[, 2]); dst <- c(el[, 2], el[, 1]); deg <- degree(g)
  list(src = src, dst = dst, edge_d = d_mat[cbind(src, dst)], deg = pmax(deg, 1L))
}
H <- seq(0, 1, by = 1/8); G <- seq(0, 1, by = 0.1); L <- c(0.5, 0.7, 1.0)
grid <- expand.grid(h = H, Gamma = G, lambda = L)

res <- list(); tm <- c(R = 0, cpp = 0, frozen_R = 0, frozen_cpp = 0); n_cells <- 0L
for (k in 1:4) {
  g <- asIgraph(readRDS(sprintf("data/02_GSS_network_ergm/GSS_net_sim_1000_%03d.rds", k)))
  attrs <- data.frame(age = V(g)$age, educ_num = V(g)$educ_num, race = as.factor(V(g)$race),
                      relig = as.factor(V(g)$relig), sex = as.factor(V(g)$sex))
  d_mat <- as.matrix(daisy(attrs, metric = "gower"))
  set.seed(500000 + k)
  g_dp <- rewire(g, with = keeping_degseq(niter = 20 * ecount(g)))
  q <- V(g)$mur_score; n <- vcount(g)
  for (topo in c("GSS", "DP")) {
    gg <- if (topo == "GSS") g else g_dp
    top <- edge_arrays(gg, d_mat)
    set.seed(k); primary <- sample.int(n, 1)
    seeds <- c(primary, head(as.integer(neighbors(gg, primary)), 27))
    for (i in seq_len(nrow(grid))) {
      x <- grid[i, ]; s <- 1e6 * k + 1e4 * (topo == "DP") + i
      t0 <- proc.time()[[3]]; rR <- run_unified_R(top, q, n, seeds, x$h, x$Gamma, x$lambda, s)
      t1 <- proc.time()[[3]]; rC <- run_cell(top, q, seeds, x$h, x$Gamma, x$lambda, s, "annealed")
      t2 <- proc.time()[[3]]; fR <- cascade_frozen_R(top, q, n, seeds, x$h, x$Gamma, x$lambda, s)
      t3 <- proc.time()[[3]]; fC <- run_cell(top, q, seeds, x$h, x$Gamma, x$lambda, s, "frozen")
      t4 <- proc.time()[[3]]
      tm <- tm + c(t1 - t0, t2 - t1, t3 - t2, t4 - t3); n_cells <- n_cells + 1L
      res[[length(res) + 1L]] <- c(
        annealed_same = all(rR == rC[names(rR)]),
        frozen_same   = all(fR == fC[c("num_adopters", "num_steps")]),
        frozen_fixed_point = fC[["stop_reason"]] %in% c(1L, 2L, 3L))
    }
  }
  message(sprintf("network %d done", k))
}
ok <- colSums(do.call(rbind, res))
message(sprintf("annealed: %d of %d cells identical to 05's R loop", ok[["annealed_same"]], n_cells))
message(sprintf("frozen:   %d of %d cells identical to the independent R cascade; %d of %d end at a fixed point",
                ok[["frozen_same"]], n_cells, ok[["frozen_fixed_point"]], n_cells))
message(sprintf("ms per simulation: R loop %.1f, Rcpp %.1f (annealed); R cascade %.1f, Rcpp %.1f (frozen)",
                1000 * tm[["R"]] / n_cells, 1000 * tm[["cpp"]] / n_cells,
                1000 * tm[["frozen_R"]] / n_cells, 1000 * tm[["frozen_cpp"]] / n_cells))
if (ok[["annealed_same"]] < n_cells || ok[["frozen_same"]] < n_cells) stop("the Rcpp engine does not match its reference")
