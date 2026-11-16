# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DIFFUSION ENGINE of 05: the Rcpp functions and the pieces shared by 05, its
# checks and the netdiffuseR reference. source("scripts/05_unified_engine.R")
# from the repo root.
#
# Engine (05_unified_engine.cpp)
#   run_exposure_grid()  the sweep engine: frozen gate, general exposure, one call
#                        per (network, topology, lambda, condition) for the grid
#   run_unified_cpp()    one simulation, annealed or frozen gate, base exposure
#                        (exact reference for 05's former R loop)
#
# Rule:  adopt_i  <=>  Gamma + lambda * E_i >= q_i  (and >= 1 adopting contact
# whose tie passes the gate). q_i = mur_score (aversion).
#
# Two gates. The pass probability of a tie at Gower distance d is
# w(d) = sigmoid((h - d) / 0.02), h = MSP.
#   annealed (05 until 2026-11): at every step each frontier tie passes with
#            probability w, drawn anew. A tie with small w still passes now and
#            then, so with enough steps the gate only delays adoption, and the
#            result depends on the step limit (200) and the patience (25).
#   frozen (the sweep since 2026-11): one coin per tie per simulation, with the
#            same probability w, shared by both directions. The tie is open or
#            closed for the whole run, so the gate filters who counts. After the
#            draw the cascade is deterministic and stops at a fixed point: no
#            step limit or patience.
#
# Exposure: the base model is E_i = a_i / d_i (adopting contacts behind open
# ties over all contacts). The robustness axis uses one formula with two dials
# (see run_exposure_grid() in the .cpp):
#   theta  independent reinforcement: adopters in the same circle (race x
#          religion cell of the ego's contacts) count less, Kish's discount
#   kappa  cohesive reinforcement: a bonus for every pair of adopting contacts
#          who know each other
# Common signed axis for figures, c = (value of the second adopter) - 1:
#   theta arm: c = -2 theta / (1 + theta)   (theta = -c / (2 + c),  c in [-1, 0])
#   kappa arm: c = 2 kappa                   (kappa = c / 2,         c > 0)
# EXPOSURE_CONDITIONS lists the conditions of the sweep (c = -1, -0.5, 0, 0.5,
# 1, each arm with its permuted control).
#
# Helpers
#   sigmoid_gate(), frozen_weights()       the gate
#   gower_of(), edge_arrays(), make_pair() network -> engine arrays, DP twin
#   seed_budget(), grow_cluster(),
#   draw_seeds()                           seeding (five strategies)
#   sim_seed_of(), sim_seed_matrix()       per-cell simulation seeds
#   build_exposure_ctx(), run_condition()  circles, pairs and one condition
#   run_cell(), run_unified_R(),
#   cascade_frozen_R(), exposure_cascade_R()  references used by the checks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(Rcpp); library(igraph); library(cluster) }))
sourceCpp("scripts/05_unified_engine.cpp")

GATE_T <- 0.02
TAU_REF_SEED <- 0.40          # seed budget: 0.40 x the network's maximum degree
SWAP_FACTOR <- 20             # DP twin: 20 x m degree-preserving swaps
sigmoid_gate <- function(d, h, T = GATE_T) 1 / (1 + exp((d - h) / T))

# w: pass probabilities of the 2E directed ties (rows 1..E are the ties, rows
# E+1..2E the same ties reversed). Call it right after set.seed(sim_seed): it
# draws E uniforms, one per tie.
frozen_weights <- function(w) {
  E <- length(w) %/% 2L
  open <- as.numeric(runif(E) <= w[seq_len(E)])
  c(open, open)
}

# ---- networks -> engine arrays ---------------------------------------------------------------
gower_of <- function(g) {
  attrs <- data.frame(age = V(g)$age, educ_num = V(g)$educ_num, race = as.factor(V(g)$race),
                      relig = as.factor(V(g)$relig), sex = as.factor(V(g)$sex))
  as.matrix(daisy(attrs, metric = "gower"))
}
# seeding: the primary nodes of the position-based strategies, computed once
edge_arrays <- function(g, d_mat, seeding = "random") {
  el <- as_edgelist(g, names = FALSE)
  src <- c(el[, 1], el[, 2]); dst <- c(el[, 2], el[, 1]); deg <- degree(g)
  out <- list(g = g, src = as.integer(src), dst = as.integer(dst), edge_d = d_mat[cbind(src, dst)],
              deg = as.integer(pmax(deg, 1L)), deg_raw = deg)
  if (seeding == "central")   out$primary_central   <- which.max(deg)[1]
  if (seeding == "closeness") out$primary_closeness <- which.max(closeness(g))[1]
  if (seeding == "eigen")     out$primary_eigen     <- which.max(eigen_centrality(g)$vector)[1]
  out
}
# k: network number (the DP seed is 500000 + k, as in 05 since 2026-08).
# dp_from: the network the twin is rewired from (default g). 05 passes the x2 network
# for x1 too: same nodes and degrees, so both levels share ONE twin and x2 vs x1 is a
# paired comparison. Slots: NET (the network with empirical homophily) and DP.
make_pair <- function(g, k, seeding = "random", d_mat = gower_of(g), dp_from = g) {
  stopifnot(vcount(dp_from) == vcount(g), identical(degree(dp_from), degree(g)))
  set.seed(500000 + k)
  g_dp <- rewire(dp_from, with = keeping_degseq(niter = SWAP_FACTOR * ecount(dp_from)))
  list(n = vcount(g), q = V(g)$mur_score,
       NET = edge_arrays(g, d_mat, seeding), DP = edge_arrays(g_dp, d_mat, seeding))
}

# ---- seeding -------------------------------------------------------------------------------------
seed_budget <- function(top) { k <- max(top$deg_raw); max(1, min(round(TAU_REF_SEED * k), k + 1)) }
# the primary plus its nearest nodes, shell by shell, sampled within the shell
# that overflows the budget
grow_cluster <- function(g, primary, n_seeds) {
  seeds <- primary
  if (n_seeds <= 1) return(seeds)
  d <- distances(g, v = primary)[1, ]
  for (shell_d in seq_len(max(d[is.finite(d)]))) {
    need <- n_seeds - length(seeds)
    if (need <= 0) break
    shell <- which(d == shell_d)
    seeds <- c(seeds, if (length(shell) <= need) shell else shell[sample.int(length(shell), need)])
  }
  seeds
}
# The draw is independent of lambda and of the exposure condition. The PRIMARY
# node of random, central and marginal is shared by the network (NET) and its DP
# twin: same degrees, so the same draw and the same seed budget (seed-matched
# comparison); closeness and eigen take each topology's own top node. Cluster
# members differ (the neighbourhoods differ) and use a topology stream.
draw_seeds <- function(top, run, topo, seeding = "random") {
  set.seed(run * 2000L)
  n <- length(top$deg_raw)
  primary <- switch(seeding,
    central   = top$primary_central,
    closeness = top$primary_closeness,
    eigen     = top$primary_eigen,
    random    = sample.int(n, 1),
    marginal  = as.integer(sample(order(top$deg_raw)[seq_len(ceiling(n * 0.1))], 1)))
  set.seed(run * 2000L + 1L + ifelse(topo == "DP", 777L, 0L))
  as.integer(grow_cluster(top$g, primary, seed_budget(top)))
}
# same formula as 05 since 2026-08; independent of the seeding strategy, shared by
# x1 and x2 and by every exposure condition (paired comparisons)
sim_seed_of <- function(run, h, Gamma, lam, topo)
  run * 1e6 + round(h * 24) * 1e4 + round(Gamma * 1000) + round(lam * 100) * 7L + ifelse(topo == "DP", 3e8, 0)
sim_seed_matrix <- function(run, lam, topo, H, G)
  outer(seq_along(H), seq_along(G), Vectorize(function(a, b) as.integer(sim_seed_of(run, H[a], G[b], lam, topo))))
# gate probabilities of the E ties (one column per MSP level)
gate_matrix <- function(top, H) {
  E <- length(top$src) %/% 2L
  vapply(H, function(h) sigmoid_gate(top$edge_d[seq_len(E)], h), numeric(E))
}

# ---- exposure: circles, pairs, conditions ------------------------------------------------------
c_of_theta <- function(theta) -2 * theta / (1 + theta)
c_of_kappa <- function(kappa) 2 * kappa
cond <- function(theta = 0, kappa_p = 0, circles = "none", pairs = "none", control = "")
  list(theta = theta, kappa_c = 0, kappa_p = kappa_p, circles = circles, pairs = pairs,
       c = if (theta > 0) c_of_theta(theta) else c_of_kappa(kappa_p), control = control)
# c = -1, -0.5, 0, +0.5, +1; every arm level with its permuted control (same circle
# sizes / same number of pairs per ego, random members)
EXPOSURE_CONDITIONS <- list(
  base         = cond(),
  theta_1      = cond(theta = 1,   circles = "cell"),
  theta_033    = cond(theta = 1/3, circles = "cell"),
  kappa_025    = cond(kappa_p = 0.25, pairs = "real"),
  kappa_05     = cond(kappa_p = 0.5,  pairs = "real"),
  theta_1_perm   = cond(theta = 1,   circles = "cell_perm", control = "theta_1"),
  theta_033_perm = cond(theta = 1/3, circles = "cell_perm", control = "theta_033"),
  kappa_025_perm = cond(kappa_p = 0.25, pairs = "perm", control = "kappa_025"),
  kappa_05_perm  = cond(kappa_p = 0.5,  pairs = "perm", control = "kappa_05"))

# circle label of the SOURCE of each directed tie: race x religion cell
cell_labels <- function(top) as.integer(interaction(V(top$g)$race, V(top$g)$relig, drop = TRUE))[top$src]
permute_within_ego <- function(top, lab, seed) {
  set.seed(seed); out <- lab
  for (ii in split(seq_along(top$dst), top$dst)) if (length(ii) > 1L) out[ii] <- lab[ii][sample.int(length(ii))]
  out
}
make_slots <- function(top, lab = NULL) {           # lab NULL: one circle per contact (base)
  n <- length(top$deg)
  if (is.null(lab)) {
    m <- length(top$src)
    return(list(slot = seq_len(m) - 1L, slot_n = rep.int(1L, m), slot_ego = top$dst))
  }
  key <- top$dst * (n + 1) + lab
  u <- !duplicated(key); id <- match(key, key[u])
  list(slot = as.integer(id - 1L), slot_n = tabulate(id), slot_ego = as.integer(top$dst[u]))
}
pairs_real <- function(top) {                       # records (ego, j, k): j ~ k, both ~ ego
  tr <- triangles(top$g)
  if (length(tr) == 0L) return(data.frame(ego = integer(0), j = integer(0), k = integer(0)))
  tr <- matrix(as.integer(tr), ncol = 3, byrow = TRUE)
  data.frame(ego = c(tr[, 1], tr[, 2], tr[, 3]), j = c(tr[, 2], tr[, 1], tr[, 1]), k = c(tr[, 3], tr[, 3], tr[, 2]))
}
pairs_permuted <- function(top, pr, seed) {         # same number of pairs per ego, random pairs
  set.seed(seed); adj <- adjacent_vertices(top$g, V(top$g))
  m_i <- tabulate(pr$ego, nbins = length(top$deg)); cache <- list()
  ego <- integer(0); jj <- integer(0); kk <- integer(0)
  for (i in which(m_i > 0L)) {
    v <- as.integer(adj[[i]]); d <- length(v); key <- as.character(d)
    if (is.null(cache[[key]])) cache[[key]] <- combn(d, 2)
    cm <- cache[[key]]; pick <- sample.int(ncol(cm), min(m_i[i], ncol(cm)))
    ego <- c(ego, rep.int(i, length(pick))); jj <- c(jj, v[cm[1, pick]]); kk <- c(kk, v[cm[2, pick]])
  }
  data.frame(ego = ego, j = jj, k = kk)
}
pairs_csr <- function(top, pr) {
  n <- length(top$deg); m <- length(top$src); key_e <- top$dst * (n + 1) + top$src
  ej <- match(pr$ego * (n + 1) + pr$j, key_e) - 1L; ek <- match(pr$ego * (n + 1) + pr$k, key_e) - 1L
  from <- c(ej, ek); to <- c(ek, ej); o <- order(from)
  list(tri_ptr = as.integer(c(0L, cumsum(tabulate(from + 1L, nbins = m)))), tri_f = as.integer(to[o]),
       eN = as.integer(tabulate(pr$ego, nbins = n)))
}
no_pairs <- function(top) list(tri_ptr = integer(length(top$src) + 1L), tri_f = integer(0), eN = integer(length(top$deg)))
# everything the conditions need for one topology; topo_id 0 (network) or 1 (DP twin).
# conds: the conditions that will run (only their circles and pairs are built)
build_exposure_ctx <- function(top, run, topo_id, conds = EXPOSURE_CONDITIONS) {
  need_c <- unique(vapply(conds, `[[`, "", "circles")); need_p <- unique(vapply(conds, `[[`, "", "pairs"))
  base_seed <- 910000L + run * 10L + topo_id
  ctx <- list(none = make_slots(top), pairs_none = no_pairs(top))
  if (any(c("cell", "cell_perm") %in% need_c)) {
    lab <- cell_labels(top)
    ctx$cell <- make_slots(top, lab)
    ctx$cell_perm <- make_slots(top, permute_within_ego(top, lab, base_seed + 2L))
  }
  if (any(c("real", "perm") %in% need_p)) {
    pr <- pairs_real(top)
    ctx$pairs_real <- pairs_csr(top, pr)
    if ("perm" %in% need_p) ctx$pairs_perm <- pairs_csr(top, pairs_permuted(top, pr, base_seed + 3L))
  }
  ctx
}
# one condition on one topology over the grid -> integer matrix (H-major rows)
run_condition <- function(top, W, ctx, cd, q, seeds, run, lam, topo, H, G, max_steps = 5000L) {
  S <- ctx[[cd$circles]]; P <- ctx[[paste0("pairs_", cd$pairs)]]
  run_exposure_grid(top$src, top$dst, W, S$slot, S$slot_n, S$slot_ego, cd$theta, cd$kappa_c,
                    P$tri_ptr, P$tri_f, cd$kappa_p, P$eN, q, G, lam, as.integer(seeds),
                    sim_seed_matrix(run, lam, topo, H, G), as.integer(max_steps))
}

# ---- references (scripts/05_engine_check.R) -------------------------------------------------
# One simulation with run_unified_cpp(). top: list(src, dst, deg, edge_d).
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

# 05's former R loop (annealed gate), with MAX_STEPS and PATIENCE as arguments
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

# The general exposure written independently (full recomputation every step).
# S: slots (make_slots), P: pairs (pairs_csr), open: 0/1 per directed tie.
exposure_cascade_R <- function(top, S, P, theta, kc, kp, q, Gamma, lam, seeds, open) {
  Af <- function(x, th, ka) ifelse(x <= 0, 0, x * (1 + ka * (x - 1)) / (1 + th * (x - 1)))
  n <- length(top$deg); m <- length(top$src); a <- rep(FALSE, n); a[seeds] <- TRUE
  tau <- pmax(0, (q - Gamma) / lam); eps <- if (theta == 0 && kc == 0 && kp == 0) 0 else 1e-12
  den <- as.numeric(tapply(Af(S$slot_n, theta, kc), factor(S$slot_ego, levels = seq_len(n)), sum))
  den[is.na(den)] <- 0; den <- den + 2 * kp * P$eN
  e_idx <- rep(seq_len(m), diff(P$tri_ptr)); f_idx <- P$tri_f + 1L
  repeat {
    act <- a[top$src] & open == 1 & !a[top$dst]; act_all <- a[top$src] & open == 1
    cnt <- tabulate(S$slot[act] + 1L, nbins = length(S$slot_n))
    num <- as.numeric(tapply(Af(cnt, theta, kc), factor(S$slot_ego, levels = seq_len(n)), sum)); num[is.na(num)] <- 0
    if (kp > 0 && length(e_idx)) {
      ok <- act_all[e_idx] & act_all[f_idx] & !a[top$dst[e_idx]]
      num <- num + 2 * kp * tabulate(top$dst[e_idx][ok], nbins = n) / 2
    }
    touched <- tabulate(top$dst[act], nbins = n) > 0
    new <- !a & touched & (num / den >= tau - eps)
    if (!any(new)) break
    a[new] <- TRUE
  }
  sum(a)
}
