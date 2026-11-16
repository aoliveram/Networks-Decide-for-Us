# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WHICH q FOR ELSOC? The premium and the reverse corner with the networks' q (c09, willingness)
# and with the GSS format (gss4), on the SAME networks: only q changes.
# 06i_q_format.R  (Prompt 6.2, 2026-11-16; from playground/PR6.1-elsoc-mechanism/scripts/
# 04_q_sensitivity_sweep.R). Decided with the user on 2026-11-12: c09 is the main q, the GSS format
# is a declared sensitivity.
#
#   c09        1 - (sum c09 - 4) / 16, wave 2018: the q of the networks (mur_score)
#   gss4       the same respondents in the GSS format, 4 items (petition, march, strike, social
#              media): done in the past year (c08 2018) / done before (c08 in 2016-2017) / might do
#              it (c09 >= 2) / would never; q = 1 - sum(4 - x) / 12 (01_survey_data.R: elsoc_panel()).
#              Nodes without it (c08 missing) keep c09 (counted)
#   c09_perm, gss4_perm   permuted across the nodes of each network (no association with attributes)
# Setting of the sweep: ELSOC_x2, its DP twin, random seeding (0.40 x maximum degree), base exposure,
# frozen gate with the sweep's simulation seeds, plane h x Gamma = (0, 1/24, ..., 1) x (0, 0.025, ..., 1),
# lambda 0.5-1.0 and the declared sensitivity 1.25, 1.5.
# The node ids are rebuilt as 02_network_calibrated.R draws them (set.seed(280000 + 1000 + i), probability
# proportional to ponderador02) and checked against the attributes and mur_score of every node; the c09
# condition reproduces the sweep of 05 cell by cell (main sweep for 0.5-1.0; --tag=lambda_ext for 1.25, 1.5).
#
# Outputs: output/06i_q_format/q_distribution.csv (mean q, share at q = 1, R2 on the Blau attributes,
#            assortativity, per q), premium.csv (premium = mean over the plane of NET - DP, pp, by q and
#            lambda, run bootstrap; and the mean over lambda 0.5-1.0 as in the sweep), corner.csv (as
#            06e: box h >= 0.6, Gamma 0.275-0.45, and the whole upper half h >= 0.6), cells.rds
#          plots/06i_q_format/q_format_corner.pdf
# Usage (repo root): Rscript scripts/06i_q_format.R 2>&1 | tee logs/06i_q_format.log    (~1 min, 8 cores)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(igraph); library(parallel); library(ggplot2) }))
Sys.setenv(OMP_NUM_THREADS = "1")
source("scripts/05_unified_engine.R")
source("scripts/01_survey_data.R")
OUT <- "output/06i_q_format"; FIG <- "plots/06i_q_format"
for (p in c(OUT, FIG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
HS <- seq(0, 1, by = 1/24); GS <- seq(0, 1, by = 0.025); LAMS <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5)
X2 <- readRDS("data/02_calibrated_networks/ELSOC_x2.rds"); NN <- length(X2)
pop <- elsoc_nodes(3); d <- elsoc_panel()
QG <- d[year == 2018, .(id, q_gss4)]; setkey(QG, id)
node_ids <- function(i) { set.seed(280000L + 1000L + i); pop$idencuesta[sample.int(nrow(pop), 1000L, replace = TRUE, prob = pop$w)] }
t0 <- Sys.time()

# ---- the two q on the nodes; check that the ids are the nodes of the networks --------------------------------
NODES <- rbindlist(lapply(seq_len(NN), function(k) {
  g <- X2[[k]]; nid <- node_ids(k); pp <- pop[match(nid, pop$idencuesta)]
  ok <- all(pp$educ == V(g)$educ, pp$sex == V(g)$sex, pp$relig == V(g)$relig, pp$ideol == V(g)$ideol,
            abs(pp$age - V(g)$age) < 1e-9, abs(pp$mur_score - V(g)$mur_score) < 1e-12)
  data.table(run = k, node = seq_len(vcount(g)), idencuesta = nid, ok = ok, q_c09 = V(g)$mur_score, q_gss4 = QG[.(nid)]$q_gss4,
             age = V(g)$age, educ_num = V(g)$educ_num, sex = V(g)$sex, relig = V(g)$relig, ideol = V(g)$ideol)
}))
stopifnot(all(NODES$ok))
NODES[, filled := is.na(q_gss4)][is.na(q_gss4), q_gss4 := q_c09]
message("node ids rebuilt for ", NN, " networks: attributes and mur_score identical in all of them; q_gss4 missing (kept c09) for ",
        round(mean(NODES[, .(f = sum(filled)), by = run]$f), 1), " nodes per network on average")
X <- survey_config("ELSOC")$blau
QD <- rbindlist(lapply(c("q_c09", "q_gss4"), function(qv) data.table(
  q = sub("q_", "", qv), mean_q = mean(NODES[[qv]]), sd_q = sd(NODES[[qv]]), share_q1 = mean(NODES[[qv]] == 1), share_q0 = mean(NODES[[qv]] == 0),
  r2_q_attributes = summary(lm(reformulate(X, qv), data = NODES))$r.squared,
  q_assortativity = mean(sapply(seq_len(NN), function(k) assortativity(X2[[k]], NODES[run == k][[qv]], directed = FALSE))),
  corr_c09_gss4 = cor(NODES$q_c09, NODES$q_gss4))))
fwrite(QD, file.path(OUT, "q_distribution.csv"))
cat("\nq ON THE NODES OF THE 24 NETWORKS (R2 on ", paste(X, collapse = " + "), "):\n", sep = ""); print(QD, digits = 3)

# ---- simulations ------------------------------------------------------------------------------------------
sim <- function(k) {
  nt <- make_pair(X2[[k]], k, "random"); qg <- NODES[run == k]$q_gss4
  stopifnot(isTRUE(all.equal(nt$q, NODES[run == k]$q_c09)))
  set.seed(7000 + k); pc <- sample(nt$q); set.seed(7100 + k); pg <- sample(qg)
  Q <- list(c09 = nt$q, gss4 = qg, c09_perm = pc, gss4_perm = pg)
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, HS)
    ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    rbindlist(lapply(names(Q), function(cn) rbindlist(lapply(LAMS, function(l) {
      r <- run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, Q[[cn]], seeds, k, l, topo, HS, GS)
      data.table(run = k, topo = topo, q = cn, lambda = l, h = rep(HS, each = length(GS)), G = rep(GS, length(HS)), phi = r[, "num_adopters"] / nt$n)
    }))))
  }))
}
message(format(Sys.time(), "%H:%M:%S"), "  ", NN, " networks x 2 topologies x 4 q x ", length(LAMS), " lambda x ", length(HS) * length(GS), " cells")
RES <- rbindlist(mclapply(seq_len(NN), sim, mc.cores = 8L))
saveRDS(RES, file.path(OUT, "cells.rds"))
message(format(Sys.time(), "%H:%M:%S"), "  simulations done (", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
# check: the c09 condition reproduces the sweep of 05
for (f in c("output/05_unified_diffusion/ELSOC_x2/unified_premium_results_random.csv",
            "output/05_unified_diffusion/ELSOC_x2/lambda_ext/unified_premium_results_random.csv")) {
  if (!file.exists(f)) { message("not found (check skipped): ", f); next }
  sw <- fread(f); sw <- sw[lambda %in% LAMS]
  w <- dcast(RES[q == "c09" & lambda %in% unique(sw$lambda)], run + lambda + h + G ~ topo, value.var = "phi")
  m <- merge(w, dcast(sw, lambda + run_id + social_distance_h + innovation_iul_Gamma ~ topology, value.var = "num_adopters"),
             by.x = c("lambda", "run", "h", "G"), by.y = c("lambda", "run_id", "social_distance_h", "innovation_iul_Gamma"))
  message("check: c09 reproduces ", f, " in ", sum(abs(m$NET.x * 1000 - m$NET.y) < 1e-6 & abs(m$DP.x * 1000 - m$DP.y) < 1e-6), " of ", nrow(m),
          " cells (lambda ", paste(sort(unique(m$lambda)), collapse = ", "), ")")
}

# ---- premium by lambda (mean over the plane, then over networks; run bootstrap) --------------------------------
g <- dcast(RES[, .(phi = mean(phi)), by = .(run, topo, q, lambda)], run + q + lambda ~ topo, value.var = "phi")[, gap := 100 * (NET - DP)]
set.seed(1)
bs_ci <- function(x) { b <- replicate(2000, mean(sample(x, replace = TRUE))); c(lo = quantile(b, .025), hi = quantile(b, .975)) }
PR <- g[, { ci <- bs_ci(gap); .(premium_pp = mean(gap), lo = ci[[1]], hi = ci[[2]]) }, by = .(q, lambda)]
PRm <- g[lambda <= 1, .(gap = mean(gap)), by = .(run, q)][, { ci <- bs_ci(gap); .(lambda = "mean 0.5-1.0", premium_pp = mean(gap), lo = ci[[1]], hi = ci[[2]]) }, by = q]
PR <- rbind(PR[, lambda := as.character(lambda)], PRm)
fwrite(PR, file.path(OUT, "premium.csv"))
cat("\nPREMIUM (pp; mean over the plane h x Gamma; 95% run bootstrap in premium.csv):\n"); print(dcast(PR, lambda ~ q, value.var = "premium_pp"), digits = 3)

# ---- corner (as 06e) --------------------------------------------------------------------------------------------
wc <- dcast(RES, run + q + lambda + h + G ~ topo, value.var = "phi")[, lead := 100 * (DP - NET)]
cell <- wc[, .(lead = mean(lead), z = mean(lead) / (sd(lead) / sqrt(.N))), by = .(q, lambda, h, G)]
CO <- merge(cell[h >= 0.6 & G >= 0.275 - 1e-9 & G <= 0.45 + 1e-9, .(box_max_dp_lead = max(lead), box_cells_z3.5 = sum(lead > 2 & z > 3.5)), by = .(q, lambda)],
            cell[h >= 0.6, .(upper_max_dp_lead = max(lead), upper_cells_z3.5 = sum(lead > 2 & z > 3.5), G_of_max = G[which.max(lead)], h_of_max = h[which.max(lead)],
                             G_min_robust = if (any(lead > 2 & z > 3.5)) min(G[lead > 2 & z > 3.5]) else NA_real_,
                             G_max_robust = if (any(lead > 2 & z > 3.5)) max(G[lead > 2 & z > 3.5]) else NA_real_), by = .(q, lambda)], by = c("q", "lambda"))
fwrite(CO, file.path(OUT, "corner.csv"))
cat("\nREVERSE CORNER (DP ahead by > 2 pp with z > 3.5; box = h >= 0.6, Gamma 0.275-0.45; upper = h >= 0.6, every Gamma):\n")
print(CO[order(q, lambda)], digits = 3)

# ---- figure: DP lead over the plane, c09 vs gss4 ------------------------------------------------------------------
pc <- cell[q %in% c("c09", "gss4") & lambda %in% c(0.5, 0.8, 1.0, 1.5)]
pc[, q_lab := factor(c(c09 = "q willingness (c09)", gss4 = "q GSS format (4 items)")[q])]
p <- ggplot(pc, aes(G, h, fill = pmax(pmin(lead, 30), -30))) + geom_raster() +
  geom_point(data = pc[lead > 2 & z > 3.5], aes(G, h), inherit.aes = FALSE, size = 0.35, colour = "black") +
  scale_fill_gradient2(low = "#B03A2E", mid = "white", high = "#1F618D", midpoint = 0, name = "DP − network (pp)") +
  facet_grid(q_lab ~ paste0("λ = ", lambda)) + labs(x = "Γ (quality)", y = "h (similarity gate)",
                                                    title = "ELSOC: where the DP twin overtakes the network, by q",
                                                    subtitle = "blue = DP ahead (reverse corner); dots: DP > 2 pp and z > 3.5; random seeding, base exposure") +
  theme_bw(base_size = 8)
ggsave(file.path(FIG, "q_format_corner.pdf"), p, width = 7.5, height = 4.2, device = cairo_pdf)
message(format(Sys.time(), "%H:%M:%S"), "  saved ", OUT, " and ", FIG, "; total ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
