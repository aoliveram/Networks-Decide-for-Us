# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CONSOLIDATION: is aversion q tied to the homophily attributes, and what does
# that do to the premium? 06e_consolidation.R  (x2 networks, base exposure,
# random seeding, the seeds of 05; ~1 min)
#
# 1. Association: R2 of q on age, sex, education, race and religion over the
#    nodes, and q assortativity over the ties of the network.
# 2. Reverse corner (h >= 0.6, Gamma 0.275-0.45, lambda 0.5-0.7), 2 x 2:
#      real       q as it is (reproduces the sweep cell by cell; checked)
#      perm       q permuted across nodes: same distribution, no association
#      swap       q mapped by rank onto the OTHER survey's distribution
#      swap_perm  both
# 3. Premium over the whole plane with real vs permuted q (lambda 0.5/0.7/1.0).
# Outputs: output/06e_consolidation/{association,corner,premium_qperm}.csv
# Usage (repo root): Rscript scripts/06e_consolidation.R [--surveys=GSS,ATP]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(parallel) }))
source("scripts/05_unified_engine.R")
ARGS <- commandArgs(trailingOnly = TRUE)
SURVEYS <- { a <- grep("^--surveys=", ARGS, value = TRUE); if (length(a)) strsplit(sub(".*=", "", a), ",")[[1]] else c("GSS", "ATP") }
OUT <- "output/06e_consolidation"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NETS <- lapply(setNames(SURVEYS, SURVEYS), function(S) readRDS(sprintf("data/02_calibrated_networks/%s_x2.rds", S)))
POOL <- lapply(NETS, function(L) unlist(lapply(L, function(g) V(g)$mur_score)))

# 1. association
assoc <- rbindlist(lapply(SURVEYS, function(S) {
  va <- rbindlist(lapply(NETS[[S]], function(g) as.data.table(vertex_attr(g)[c("age", "sex", "educ_num", "race", "relig", "mur_score")])))
  data.table(survey = S, r2_q_attributes = summary(lm(mur_score ~ age + sex + educ_num + race + relig, data = va))$r.squared,
             q_assortativity = mean(sapply(NETS[[S]], function(g) assortativity(g, V(g)$mur_score, directed = FALSE))),
             q_distinct_values = uniqueN(va$mur_score), share_q1 = mean(va$mur_score == 1))
}))
fwrite(assoc, file.path(OUT, "association.csv")); print(assoc, digits = 3)

# simulations
to_marginal <- function(q, target, seed) { set.seed(seed); r <- rank(q + runif(length(q), 0, 1e-9), ties.method = "first")
  as.numeric(quantile(target, (r - 0.5) / length(q), type = 1)) }
sim <- function(S, k, conds, HS, GS, LAMS) {
  nt <- make_pair(NETS[[S]][[k]], k, "random"); set.seed(7000 + k); perm <- sample(nt$q)
  other <- if (length(POOL) > 1) POOL[[setdiff(SURVEYS, S)[1]]] else nt$q
  Q <- list(real = nt$q, perm = perm, swap = to_marginal(nt$q, other, 8000 + k), swap_perm = to_marginal(perm, other, 8000 + k))[conds]
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, HS)
    ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    rbindlist(lapply(names(Q), function(cn) rbindlist(lapply(LAMS, function(l) {
      r <- run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, Q[[cn]], seeds, k, l, topo, HS, GS)
      data.table(survey = S, run = k, topo = topo, cond = cn, lambda = l, h = rep(HS, each = length(GS)), G = rep(GS, length(HS)),
                 phi = r[, "num_adopters"] / nt$n)
    }))))
  }))
}
run_all <- function(...) rbindlist(lapply(SURVEYS, function(S) rbindlist(mclapply(seq_along(NETS[[S]]), function(k) sim(S, k, ...), mc.cores = 8))))

# 2. corner
H2 <- seq(0, 1, by = 1/24); H2 <- H2[H2 >= 0.5 - 1e-9]
res <- run_all(c("real", "perm", "swap", "swap_perm"), H2, seq(0.2, 0.5, by = 0.025), c(0.5, 0.6, 0.7))
w <- dcast(res, survey + run + cond + lambda + h + G ~ topo, value.var = "phi")[, lead := 100 * (DP - NET)]
cell <- w[, .(lead = mean(lead), z = mean(lead) / (sd(lead) / sqrt(.N))), by = .(survey, cond, lambda, h, G)]
corner <- cell[h >= 0.6 & G >= 0.275 - 1e-9 & G <= 0.45 + 1e-9,
               .(mean_gap_pp = -mean(lead), max_dp_lead_pp = max(lead), cells_dp_gt2 = sum(lead > 2), cells_dp_gt2_z3.5 = sum(lead > 2 & z > 3.5)),
               by = .(survey, cond, lambda)]
setorder(corner, survey, cond, lambda); fwrite(corner, file.path(OUT, "corner.csv")); print(corner, digits = 3)
for (S in SURVEYS) {                                  # the real condition must reproduce the sweep
  f <- sprintf("output/05_unified_diffusion/%s_x2/unified_premium_results_random.csv", S); if (!file.exists(f)) next
  sw <- dcast(fread(f)[lambda %in% c(0.5, 0.6, 0.7)], lambda + run_id + social_distance_h + innovation_iul_Gamma ~ topology, value.var = "num_adopters")
  m <- merge(w[survey == S & cond == "real"], sw, by.x = c("lambda", "run", "h", "G"), by.y = c("lambda", "run_id", "social_distance_h", "innovation_iul_Gamma"))
  message(S, ": real condition identical to the sweep in ", sum(abs(m$NET.x * 1000 - m$NET.y) < 1e-6 & abs(m$DP.x * 1000 - m$DP.y) < 1e-6), " of ", nrow(m), " cells")
}

# 3. premium with permuted q
res3 <- run_all(c("real", "perm"), seq(0, 1, by = 1/24), seq(0, 1, by = 0.025), c(0.5, 0.7, 1.0))
g <- dcast(res3[, .(phi = mean(phi)), by = .(survey, run, topo, cond, lambda)], survey + run + cond + lambda ~ topo, value.var = "phi")[, gap := 100 * (NET - DP)]
set.seed(1)
prem <- g[, .(premium_pp = mean(gap), lo = quantile(replicate(2000, mean(sample(gap, replace = TRUE))), .025),
              hi = quantile(replicate(2000, mean(sample(gap, replace = TRUE))), .975)), by = .(survey, lambda, cond)]
setorder(prem, survey, lambda, cond); fwrite(prem, file.path(OUT, "premium_qperm.csv")); print(prem, digits = 3)
