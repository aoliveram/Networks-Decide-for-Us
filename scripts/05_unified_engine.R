# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DIFFUSION ENGINE of 05: the Rcpp function and its R helpers.
# source("scripts/05_unified_engine.R") from the repo root.
#
#   run_unified_cpp()   Rcpp version of 05's run_unified() (05_unified_engine.cpp)
#   sigmoid_gate()      pass probability of a tie at Gower distance d
#   frozen_weights()    frozen gate: one coin per tie for the whole simulation
#   run_cell()          one simulation under either gate
#   run_unified_R()     05's R loop (reference for the annealed gate)
#   cascade_frozen_R()  independent R version of the frozen gate (reference)
#
# Two gates:
#   annealed (05): at every step each frontier tie passes with probability
#            w(d_ij) = sigmoid((h - d_ij) / 0.02), drawn anew. A tie with small w
#            still passes now and then, so with enough steps the gate only delays
#            adoption, and the result depends on the step limit (200) and the
#            patience (25 steps without adoption).
#   frozen:  one coin per tie per simulation, with the same probability w(d_ij),
#            shared by both directions of the tie. The tie is open or closed for
#            the whole run: the gate is a filter on who counts. After the draw the
#            cascade is deterministic and stops at the first step without new
#            adopters (a fixed point), so no step limit or patience is needed.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressMessages(library(Rcpp))
sourceCpp("scripts/05_unified_engine.cpp")

GATE_T <- 0.02
sigmoid_gate <- function(d, h, T = GATE_T) 1 / (1 + exp((d - h) / T))

# w: pass probabilities of the 2E directed ties, stored as in 05's edge_arrays()
# (rows 1..E are the ties, rows E+1..2E the same ties reversed). Call it right
# after set.seed(sim_seed): it draws E uniforms, one per tie.
frozen_weights <- function(w) {
  E <- length(w) %/% 2L
  open <- as.numeric(runif(E) <= w[seq_len(E)])
  c(open, open)
}

# One simulation. top: list(src, dst, deg, edge_d) as built by 05's edge_arrays().
# Returns num_adopters, num_adopted_rational, num_adopted_social, num_steps,
# initial_cluster_size (the five outputs of 05), stop_reason (1 all adopted,
# 2 frontier empty, 3 patience / fixed point, 4 step limit), t_stop, t50 and the
# adopters at each horizon (n_t25, ...).
run_cell <- function(top, q, seeds, h, Gamma, lambda, sim_seed, gate = c("annealed", "frozen"),
                     max_steps = 200L, patience = 25L, horizons = c(25L, 50L, 100L, 200L)) {
  gate <- match.arg(gate)
  w <- sigmoid_gate(top$edge_d, h)
  set.seed(sim_seed)
  if (gate == "frozen") {
    w <- frozen_weights(w)
    max_steps <- max(max_steps, length(q))   # a cascade adds >= 1 node per step
    patience <- 1L
  }
  run_unified_cpp(top$src, top$dst, top$deg, w, pmax(0, (q - Gamma) / lambda), q, Gamma,
                  as.integer(seeds), as.integer(max_steps), as.integer(patience), as.integer(horizons))
}

# 05's run_unified(), with MAX_STEPS and PATIENCE as arguments instead of globals
run_unified_R <- function(top, q, n, seeds, h, Gamma, lambda, sim_seed, MAX_STEPS = 200, PATIENCE = 25) {
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

# Frozen gate written independently: draw the open ties, then iterate the
# threshold rule over open ties until nothing changes
cascade_frozen_R <- function(top, q, n, seeds, h, Gamma, lambda, sim_seed) {
  tau_star <- pmax(0, (q - Gamma) / lambda)
  w <- sigmoid_gate(top$edge_d, h)
  set.seed(sim_seed)
  E <- length(w) %/% 2L
  coin <- runif(E) <= w[seq_len(E)]
  open <- c(coin, coin)
  a <- rep(FALSE, n); a[seeds] <- TRUE; steps <- 0L
  repeat {
    idx <- which(a[top$src] & !a[top$dst] & open)
    cnt <- tabulate(top$dst[idx], nbins = n)
    new <- !a & cnt > 0L & cnt / top$deg >= tau_star
    if (!any(new)) break
    a[new] <- TRUE; steps <- steps + 1L
  }
  c(num_adopters = sum(a), num_steps = steps)
}
