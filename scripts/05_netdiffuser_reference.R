# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REFERENCE: the sweep engine against netdiffuseR (the reference implementation
# cited in the paper). 05_netdiffuser_reference.R
#
# Under the frozen gate a simulation is: draw the open ties once (one coin per
# tie, probability w(d_ij) = sigma((h - d_ij) / 0.02), right after
# set.seed(sim_seed)), then run a deterministic threshold cascade on the open
# ties. In netdiffuseR that is rdiffnet() on the graph of open ties with
#   - seed.nodes  = the seed cluster of 05,
#   - exposure    = number of adopting contacts behind open ties
#                   (exposure.args: normalized = FALSE),
#   - an adoption mechanism (adoption_mechanism) that divides by the FULL degree
#     (the base exposure counts open adopting contacts over all contacts) and
#     adopts when Gamma + lambda * E >= q, i.e. E >= tau = (q - Gamma) / lambda,
#     with at least one open adopting contact.
# Both are synchronous, so the final number of adopters must be IDENTICAL cell
# by cell, not only equal in distribution.
#
# Networks GSS_x2 #1-2 and their DP twins (05's recipe), random seeding,
# lambda 0.5 / 0.7 / 1.0, 30 random (h, Gamma) cells per topology and lambda.
# Usage (repo root): Rscript scripts/05_netdiffuser_reference.R        (~1-2 min)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(netdiffuseR); library(Matrix) }))
source("scripts/05_unified_engine.R")
H <- seq(0, 1, by = 1/24); G <- seq(0, 1, by = 0.025); L <- c(0.5, 0.7, 1.0)
nets <- readRDS("data/02_calibrated_networks/GSS_x2.rds")
T_MAX <- 400L      # periods given to rdiffnet; a cascade that has not stopped by then is flagged

mech <- function(expo, thresholds, not_adopted, time, pars)   # E over the full degree, >= 1 open adopter
  which(not_adopted & expo > 0 & expo / pars$deg >= thresholds)

ndr_cell <- function(top, n, q, seeds, h, Gamma, lam, sim_seed) {
  E <- length(top$src) %/% 2L
  set.seed(sim_seed)
  open <- runif(E) <= sigmoid_gate(top$edge_d[seq_len(E)], h)
  A <- sparseMatrix(i = c(top$src[seq_len(E)][open], top$dst[seq_len(E)][open]),
                    j = c(top$dst[seq_len(E)][open], top$src[seq_len(E)][open]), x = 1, dims = c(n, n))
  d <- suppressWarnings(suppressMessages(rdiffnet(
    seed.graph = A, t = T_MAX, seed.nodes = seeds, rewire = FALSE, threshold.dist = pmax(0, (q - Gamma) / lam),
    exposure.args = list(normalized = FALSE, valued = FALSE), stop.no.diff = FALSE,
    adoption_mechanism = mech, adoption_pars = list(deg = top$deg))))
  toa <- d$toa
  c(adopters = sum(!is.na(toa)), last = max(toa, na.rm = TRUE))
}

set.seed(5); rows <- list(); t_ndr <- 0; t_eng <- 0
for (run in 1:2) {
  nt <- make_pair(nets[[run]], run)
  for (topo in c("NET", "DP")) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, run, topo); W <- gate_matrix(top, H)
    ctx <- build_exposure_ctx(top, run, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    for (lam in L) {
      t0 <- proc.time()[[3]]
      eng <- run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, nt$q, seeds, run, lam, topo, H, G)
      t_eng <- t_eng + proc.time()[[3]] - t0
      cells <- cbind(sample(length(H), 30, TRUE), sample(length(G), 30, TRUE))
      for (k in seq_len(nrow(cells))) {
        hi <- cells[k, 1]; gi <- cells[k, 2]
        t0 <- proc.time()[[3]]
        r <- ndr_cell(top, nt$n, nt$q, seeds, H[hi], G[gi], lam, sim_seed_of(run, H[hi], G[gi], lam, topo))
        t_ndr <- t_ndr + proc.time()[[3]] - t0
        rows[[length(rows) + 1L]] <- data.frame(run = run, topology = topo, lambda = lam, h = H[hi], Gamma = G[gi],
                                                engine = eng[(hi - 1L) * length(G) + gi, "num_adopters"],
                                                netdiffuseR = r[["adopters"]], last_period = r[["last"]])
      }
    }
  }
  message(sprintf("network %d done", run))
}
res <- do.call(rbind, rows)
same <- res$engine == res$netdiffuseR
message(sprintf("netdiffuseR %s vs the sweep engine: %d of %d cells with the same number of adopters (range %d-%d)",
                as.character(packageVersion("netdiffuseR")), sum(same), nrow(res), min(res$engine), max(res$engine)))
message(sprintf("cascades still running at period %d in netdiffuseR: %d", T_MAX, sum(res$last_period >= T_MAX)))
message(sprintf("ms per simulation: netdiffuseR %.0f, sweep engine %.2f",
                1000 * t_ndr / nrow(res), 1000 * t_eng / (2 * 2 * length(L) * length(H) * length(G))))
if (any(!same)) { print(head(res[!same, ], 10)); stop("netdiffuseR and the engine disagree") }
message("REFERENCE CHECK PASSED")
