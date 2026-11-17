# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DOES THE PREMIUM SURVIVE A NOISY RULE? 06k_noisy_rule.R  (Prompt 7b, item 7 of "Qué falta";
# 2026-11-17). The paper estimates in ELSOC a rule with noise (sigma_hat ~ 0.48 in units of q, an upper
# bound: it absorbs measurement error) under which 14-15% of the less averse (q <= 0.5) take part with
# no exposure at all, and then simulates a deterministic rule that needs one adopting contact.
#
# Rule:  adopt_i  <=>  Gamma + lambda E_i - q_i + eps_i >= 0,   eps_i ~ logistic(0, sigma)
#   (the structural model of 10_ELSOC_rule.R, same scale: sd = sigma pi / sqrt(3))
# Two variants:
#   noise        (a) noise only; the condition of one adopting contact behind an open tie is kept
#   noise_spont  (b) noise + spontaneous adoption: no contact condition, so whoever has q'_i <= Gamma
#                    adopts with no exposure
# sigma in {0, 0.1, 0.25, 0.5}; sigma = 0 of (a) is the main sweep of 05 (checked cell by cell);
# sigma = 0 of (b) removes only the contact condition.
#
# Design (and why):
#   - FROZEN noise. eps_i is drawn once and kept for the whole cascade, like the coin of the gate:
#     adoption is permanent, so a fresh eps every step would let everyone adopt sooner or later.
#     A frozen eps is the same as running the engine with q'_i = q_i - eps_i: the engine reads q only
#     through tau_i = max(0, (q_i - Gamma) / lambda). 05_unified_engine.cpp is NOT touched.
#   - One eps vector per network (run) and survey: eps_i = sigma * L_i, L_i standard logistic drawn with
#     set.seed(EPS_SEED[survey] + run). It is shared by the network and its DP twin (same node, same
#     simulation: paired contrast, like the seeds), by every lambda and every cell of the plane (as 06e
#     and 06i pass one permuted q per network), and by every sigma (common random numbers: the
#     difference with sigma = 0 is paired). Each cell has the right marginal distribution of eps; the
#     run bootstrap is valid because eps is independent across runs.
#   - Spontaneous adoption WITHOUT an engine switch. With the contact condition removed, the nodes with
#     q'_i <= Gamma (tau_i = 0) adopt at step 1; for every other node tau_i > 0, so E_i >= tau_i already
#     implies an adopting contact and the condition is redundant. The cascade is monotone under the
#     frozen gate and its final set is the least fixed point containing the initial adopters, so (b) is
#     exactly the engine run with seeds = seeds of 05 + {q'_i <= Gamma}. One engine call per Gamma (the
#     extra seeds depend on Gamma), with the same simulation seed (gate coins) of each cell. Checked
#     against an independent R cascade that drops the contact condition (cascade_noisy_R below).
#   - Everything else as 05: x2 networks, base exposure, random seeding (0.40 x maximum degree), frozen
#     gate with the simulation seeds of 05, plane h x Gamma = (0, 1/24, ..., 1) x (0, 0.025, ..., 1),
#     the lambda grid of 05 (0.5-1.5 every 0.1).
#
# Outputs (output/06k_noisy_rule/; one block of rows per survey, a run replaces its surveys' rows):
#   premium.csv      premium (pp, mean over the plane of NET - DP), adoption of NET and DP, and the paired
#                    difference with sigma = 0 of the main rule ("d_vs_main", variant a) and of the same
#                    variant ("d_vs_sigma0"), by survey x variant x sigma x lambda and lambda-means
#                    (0.5-1.5, the headline; 0.5-1.0); 95% run bootstrap (B = 2000)
#   band.csv         share of the plane with gap > 10 pp, mean gap in that band, MSP_c (first h with
#                    Phi >= 0.5, per Gamma; mean over the Gamma where both topologies cross) and its shift
#                    DP - NET (as 06_premium_sensitivity.R, sections C and E)
#   corner.csv       DP ahead by > 2 pp with z > 3.5 across runs (as 06e/06i): box h >= 0.6, Gamma
#                    0.275-0.45, and the whole upper half h >= 0.6
#   corner_perm.csv  the same count with q permuted across nodes (as 06e), variant (a), lambda >= 1, h >= 0.6,
#                    Gamma <= 0.5: is the corner that noise brings back at lambda >= 1 still consolidation?
#   spontaneous.csv  share of nodes with q' <= Gamma per sigma and Gamma (all; noise-made, q' <= Gamma < q;
#                    among q <= 0.5), mean over the networks
#   elsoc_calibration.csv  which sigma matches ELSOC's 14-15% (participation with E = 0 and q <= 0.5)
#   cells_<S>.rds    every simulation (not tracked);  paragraph.txt (the paragraph for the paper, by hand)
#   plots/06k_noisy_rule/noisy_rule_premium.pdf  premium vs sigma by survey
# Usage (repo root): Rscript scripts/06k_noisy_rule.R [--surveys=GSS,ATP,GSS14,ELSOC] [--lambdas=0.7,1,1.5]
#   [--sigmas=0,0.1,0.25,0.5] [--test]   2>&1 | tee logs/06k_noisy_rule.log
#   --test: 2 networks, lambda 1.0, sandbox output/06k_noisy_rule/_test. Runtime: ~3-4 min per survey, 8 cores.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(igraph); library(parallel); library(ggplot2) }))
Sys.setenv(OMP_NUM_THREADS = "1")
source("scripts/05_unified_engine.R")
options(width = 200)
ARGS <- commandArgs(trailingOnly = TRUE)
arg_list <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) strsplit(sub(".*=", "", a), ",")[[1]] else d }
TEST <- "--test" %in% ARGS
SURVEYS <- arg_list("surveys", c("GSS", "ATP", "GSS14", "ELSOC"))
LAMS <- as.numeric(sprintf("%.1f", as.numeric(arg_list("lambdas", if (TEST) "1.0" else LAMBDA_GRID))))
SIGMAS <- as.numeric(arg_list("sigmas", c(0, 0.1, 0.25, 0.5)))
stopifnot(SIGMAS[1] == 0)                              # sigma = 0 is the reference of every paired difference
NCORE <- 8L
OUT <- if (TEST) "output/06k_noisy_rule/_test" else "output/06k_noisy_rule"
FIG <- if (TEST) file.path(OUT, "plots") else "plots/06k_noisy_rule"
for (p in c(OUT, FIG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
HS <- seq(0, 1, by = 1/24); GS <- seq(0, 1, by = 0.025)
EPS_SEED <- c(GSS = 660000L, ATP = 661000L, GSS14 = 662000L, ELSOC = 663000L)
VARIANTS <- c("noise", "noise_spont")
t0 <- Sys.time()
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
NETS <- lapply(setNames(SURVEYS, SURVEYS), function(S) readRDS(sprintf("data/02_calibrated_networks/%s_x2.rds", S)))
RUNS <- function(S) if (TEST) 1:2 else seq_along(NETS[[S]])
eps_L <- function(S, k, n) { set.seed(EPS_SEED[[S]] + k); qlogis(runif(n)) }   # standard logistic, one per node
upsert <- function(d, file) {                          # replace the rows of these surveys, keep the others
  f <- file.path(OUT, file)
  if (file.exists(f)) d <- rbind(fread(f)[!survey %in% unique(d$survey)], d, fill = TRUE)
  setorderv(d, intersect(c("survey", "variant", "sigma", "lambda", "G"), names(d))); fwrite(d, f); d
}

# ---- simulations ------------------------------------------------------------------------------------------
grid_call <- function(top, W, ctx, q, Gs, lam, seeds, SM)
  run_exposure_grid(top$src, top$dst, W, ctx$none$slot, ctx$none$slot_n, ctx$none$slot_ego, 0, 0,
                    ctx$pairs_none$tri_ptr, ctx$pairs_none$tri_f, 0, ctx$pairs_none$eN,
                    q, Gs, lam, as.integer(seeds), SM, 5000L)[, "num_adopters"]
sim <- function(S, k) {
  nt <- make_pair(NETS[[S]][[k]], k, "random"); L <- eps_L(S, k, nt$n)
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, HS)
    ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    rbindlist(lapply(LAMS, function(l) {
      SM <- sim_seed_matrix(k, l, topo, HS, GS)
      rbindlist(lapply(SIGMAS, function(s) {
        qp <- nt$q - s * L
        a <- grid_call(top, W, ctx, qp, GS, l, seeds, SM)                       # rows h-major
        b <- vapply(seq_along(GS), function(gi)                                  # columns = Gamma
          grid_call(top, W, ctx, qp, GS[gi], l, union(seeds, which(qp <= GS[gi])), SM[, gi, drop = FALSE]), numeric(length(HS)))
        data.table(survey = S, run = k, topo = topo, sigma = s, lambda = l, variant = rep(VARIANTS, each = length(HS) * length(GS)),
                   h = rep(HS, each = length(GS)), G = rep(GS, length(HS)), phi = c(a, as.vector(t(b))) / nt$n)
      }))
    }))
  }))
}
RES <- rbindlist(lapply(SURVEYS, function(S) {
  stamp(S, ": ", length(RUNS(S)), " networks x 2 topologies x ", length(LAMS), " lambda x ", length(SIGMAS), " sigma x 2 variants x ",
        length(HS) * length(GS), " cells")
  r <- rbindlist(mclapply(RUNS(S), function(k) sim(S, k), mc.cores = NCORE))
  saveRDS(r, file.path(OUT, sprintf("cells_%s.rds", S))); stamp(S, " done"); r
}))

# ---- checks -----------------------------------------------------------------------------------------------
# 1. variant (a) with sigma = 0 reproduces the sweep of 05, cell by cell
for (S in SURVEYS) {
  f <- sprintf("output/05_unified_diffusion/%s_x2/unified_premium_results_random.csv", S)
  if (!file.exists(f)) { message("not found (check skipped): ", f); next }
  sw <- fread(f)[run_id %in% RUNS(S)][, `:=`(lambda = round(lambda, 6), social_distance_h = round(social_distance_h, 6), innovation_iul_Gamma = round(innovation_iul_Gamma, 6))]
  sw <- dcast(sw[lambda %in% round(LAMS, 6)], lambda + run_id + social_distance_h + innovation_iul_Gamma ~ topology, value.var = "num_adopters")
  w <- dcast(RES[survey == S & variant == "noise" & sigma == 0][, `:=`(lambda = round(lambda, 6), h = round(h, 6), G = round(G, 6))],  # CSV keys: round
             lambda + run + h + G ~ topo, value.var = "phi")
  m <- merge(w, sw, by.x = c("lambda", "run", "h", "G"), by.y = c("lambda", "run_id", "social_distance_h", "innovation_iul_Gamma"))
  message(sprintf("check 1, %s: sigma = 0 (noise) identical to the sweep of 05 in %d of %d cells (%d expected)", S,
                  sum(abs(m$NET.x * 1000 - m$NET.y) < 1e-6 & abs(m$DP.x * 1000 - m$DP.y) < 1e-6), nrow(m), nrow(w)))
}
# 2. both variants against an independent R cascade (frozen gate, q' = q - eps; spont = TRUE drops the
#    condition of one adopting contact). Random cells of network 1 and 2, both topologies.
cascade_noisy_R <- function(top, q, n, seeds, h, Gamma, lambda, sim_seed, spont) {
  tau <- pmax(0, (q - Gamma) / lambda)
  w <- sigmoid_gate(top$edge_d, h); set.seed(sim_seed)
  E <- length(w) %/% 2L; coin <- runif(E) <= w[seq_len(E)]; open <- c(coin, coin)
  a <- rep(FALSE, n); a[seeds] <- TRUE
  repeat {
    cnt <- tabulate(top$dst[a[top$src] & open], nbins = n)
    new <- !a & (cnt > 0L | spont) & cnt / top$deg >= tau
    if (!any(new)) break
    a[new] <- TRUE
  }
  sum(a)
}
set.seed(2)
CHK <- rbindlist(lapply(SURVEYS, function(S) rbindlist(lapply(head(RUNS(S), 2), function(k) {
  nt <- make_pair(NETS[[S]][[k]], k, "random"); L <- eps_L(S, k, nt$n)
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random")
    cl <- data.table(sigma = sample(SIGMAS, 12, TRUE), lambda = sample(LAMS, 12, TRUE), h = sample(HS[HS >= 0.25], 12, TRUE),
                     G = sample(GS[GS <= 0.6], 12, TRUE), variant = rep(VARIANTS, 6))
    cl[, ref := mapply(function(s, l, h, G, v) cascade_noisy_R(top, nt$q - s * L, nt$n, seeds, h, G, l, sim_seed_of(k, h, G, l, topo), v == "noise_spont") / nt$n,
                       sigma, lambda, h, G, variant)]
    cl[, `:=`(survey = S, run = k, topo = topo)]
  }))
}))))
CHK <- merge(CHK, RES, by = c("survey", "run", "topo", "sigma", "lambda", "variant", "h", "G"))
message(sprintf("check 2: engine = independent R cascade in %d of %d random cells (noise %d, noise_spont %d)",
                sum(abs(CHK$ref - CHK$phi) < 1e-9), nrow(CHK), CHK[variant == "noise", sum(abs(ref - phi) < 1e-9)], CHK[variant == "noise_spont", sum(abs(ref - phi) < 1e-9)]))
stopifnot(all(abs(CHK$ref - CHK$phi) < 1e-9))

# ---- premium by sigma, variant and lambda (mean over the plane, then over networks; run bootstrap) ----------
set.seed(1)
bs_ci <- function(x) { b <- replicate(2000, mean(sample(x, replace = TRUE))); c(quantile(b, .025), quantile(b, .975)) }
g <- dcast(RES[, .(phi = mean(phi)), by = .(survey, run, topo, variant, sigma, lambda)], survey + run + variant + sigma + lambda ~ topo, value.var = "phi")
g[, gap := 100 * (NET - DP)]
g <- merge(g, g[variant == "noise" & sigma == 0, .(survey, run, lambda, gap_main = gap)], by = c("survey", "run", "lambda"))
g <- merge(g, g[sigma == 0, .(survey, run, lambda, variant, gap_s0 = gap)], by = c("survey", "run", "lambda", "variant"))
g[, `:=`(d_main = gap - gap_main, d_s0 = gap - gap_s0)]
summ <- function(x) x[, { a <- bs_ci(gap); b <- bs_ci(d_main); c <- bs_ci(d_s0)
  .(premium_pp = mean(gap), lo = a[[1]], hi = a[[2]], phi_NET = mean(NET), phi_DP = mean(DP),
    d_vs_main_pp = mean(d_main), d_vs_main_lo = b[[1]], d_vs_main_hi = b[[2]], d_vs_sigma0_pp = mean(d_s0), d_vs_sigma0_lo = c[[1]], d_vs_sigma0_hi = c[[2]]) },
  by = .(survey, variant, sigma, lambda)]
PR <- summ(g)[, lambda := sprintf("%.1f", as.numeric(lambda))]
for (top in c(1.5, 1.0)) if (sum(LAMS <= top + 1e-9) > 1) {
  gm <- g[lambda <= top + 1e-9, lapply(.SD, mean), by = .(survey, run, variant, sigma), .SDcols = c("gap", "NET", "DP", "d_main", "d_s0")]
  gm[, lambda := sprintf("mean %.1f-%.1f", min(LAMS), top)]
  PR <- rbind(PR, summ(gm))
}
PR <- upsert(PR, "premium.csv")
cat("\nPREMIUM (pp; mean over the plane; 95% run bootstrap in premium.csv). d = paired difference with the main rule (noise, sigma = 0):\n")
print(PR[survey %in% SURVEYS & grepl("mean|^0.7$|^1.0$|^1.5$", lambda), .(survey, variant, sigma, lambda, premium_pp, lo, hi, phi_NET, phi_DP, d_vs_main_pp, d_vs_main_lo, d_vs_main_hi)], digits = 3)

# ---- band and inflection point (as 06_premium_sensitivity.R, C and E) ---------------------------------------
cm <- RES[, .(phi = mean(phi)), by = .(survey, variant, sigma, lambda, topo, h, G)]
wc <- dcast(cm, survey + variant + sigma + lambda + h + G ~ topo, value.var = "phi")[, gap := NET - DP]
band <- wc[, .(pct_plane_gap10 = 100 * mean(gap > 0.10), gap_in_band_pp = if (any(gap > 0.10)) 100 * mean(gap[gap > 0.10]) else NA_real_,
               max_gap_pp = 100 * max(gap)), by = .(survey, variant, sigma, lambda)]
tip <- cm[, .(h_c = if (any(phi >= 0.5)) min(h[phi >= 0.5]) else NA_real_), by = .(survey, variant, sigma, lambda, topo, G)]
tip <- dcast(tip, survey + variant + sigma + lambda + G ~ topo, value.var = "h_c")
tip <- tip[, .(n_G_both_cross = sum(!is.na(NET) & !is.na(DP)), n_G_NET_only = sum(!is.na(NET) & is.na(DP)),
               MSPc_NET = mean(NET[!is.na(NET) & !is.na(DP)]), MSPc_DP = mean(DP[!is.na(NET) & !is.na(DP)]),
               MSPc_shift = mean((DP - NET)[!is.na(NET) & !is.na(DP)])), by = .(survey, variant, sigma, lambda)]
BA <- merge(band, tip, by = c("survey", "variant", "sigma", "lambda"))
BA <- upsert(BA, "band.csv")
cat("\nBAND (gap > 10 pp) AND INFLECTION POINT (MSP_c: first h with Phi >= 0.5; shift = DP - NET), lambda 0.7 / 1.0 / 1.5:\n")
print(BA[survey %in% SURVEYS & round(lambda, 1) %in% c(0.7, 1.0, 1.5)], digits = 3)

# ---- reverse corner (as 06e/06i) ----------------------------------------------------------------------------
wr <- dcast(RES, survey + run + variant + sigma + lambda + h + G ~ topo, value.var = "phi")[, lead := 100 * (DP - NET)]
cell <- wr[, .(lead = mean(lead), z = mean(lead) / (sd(lead) / sqrt(.N))), by = .(survey, variant, sigma, lambda, h, G)]
cell[is.nan(z), z := 0]
CO <- merge(cell[h >= 0.6 & G >= 0.275 - 1e-9 & G <= 0.45 + 1e-9, .(box_max_dp_lead = max(lead), box_cells_z3.5 = sum(lead > 2 & z > 3.5)), by = .(survey, variant, sigma, lambda)],
            cell[h >= 0.6, .(upper_max_dp_lead = max(lead), upper_cells_z3.5 = sum(lead > 2 & z > 3.5), G_of_max = G[which.max(lead)], h_of_max = h[which.max(lead)],
                             G_min_robust = if (any(lead > 2 & z > 3.5)) min(G[lead > 2 & z > 3.5]) else NA_real_,
                             G_max_robust = if (any(lead > 2 & z > 3.5)) max(G[lead > 2 & z > 3.5]) else NA_real_), by = .(survey, variant, sigma, lambda)],
            by = c("survey", "variant", "sigma", "lambda"))
CO <- upsert(CO, "corner.csv")
cat("\nREVERSE CORNER (DP ahead by > 2 pp with z > 3.5): robust cells in the upper half h >= 0.6, by lambda (columns) :\n")
print(dcast(CO[survey %in% SURVEYS], survey + variant + sigma ~ lambda, value.var = "upper_cells_z3.5"))
cat("largest DP lead (pp) in the upper half:\n")
print(dcast(CO[survey %in% SURVEYS], survey + variant + sigma ~ lambda, value.var = "upper_max_dp_lead"), digits = 2)

# Is the corner under noise still consolidation? Variant (a), lambda >= 1, upper half (h >= 0.6, Gamma <= 0.5),
# q as it is vs q permuted across the nodes (as 06e: set.seed(7000 + run)); the same eps is added to both.
HP <- HS[HS >= 0.6 - 1e-9]; GP <- GS[GS <= 0.5 + 1e-9]; LP <- LAMS[LAMS >= 1 - 1e-9]
if (length(LP)) {
  simp <- function(S, k) {
    nt <- make_pair(NETS[[S]][[k]], k, "random"); L <- eps_L(S, k, nt$n); set.seed(7000 + k); perm <- sample(nt$q)
    rbindlist(lapply(c("NET", "DP"), function(topo) {
      top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, HP)
      ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
      rbindlist(lapply(LP, function(l) { SM <- sim_seed_matrix(k, l, topo, HP, GP)
        rbindlist(lapply(SIGMAS, function(s) rbindlist(lapply(c("real", "perm"), function(qc)
          data.table(survey = S, run = k, topo = topo, sigma = s, lambda = l, q = qc, h = rep(HP, each = length(GP)), G = rep(GP, length(HP)),
                     phi = grid_call(top, W, ctx, (if (qc == "real") nt$q else perm) - s * L, GP, l, seeds, SM) / nt$n)))))
      }))
    }))
  }
  RP <- rbindlist(lapply(SURVEYS, function(S) rbindlist(mclapply(RUNS(S), function(k) simp(S, k), mc.cores = NCORE))))
  m <- merge(RP[q == "real"], RES[variant == "noise"], by = c("survey", "run", "topo", "sigma", "lambda", "h", "G"))
  message(sprintf("check 3: q real of the permutation test = main simulations in %d of %d cells", sum(abs(m$phi.x - m$phi.y) < 1e-12), nrow(RP[q == "real"])))
  wp <- dcast(RP, survey + run + sigma + lambda + q + h + G ~ topo, value.var = "phi")[, lead := 100 * (DP - NET)]
  cp <- wp[, .(lead = mean(lead), z = mean(lead) / (sd(lead) / sqrt(.N))), by = .(survey, sigma, lambda, q, h, G)][is.nan(z), z := 0]
  CP <- cp[, .(cells_z3.5 = sum(lead > 2 & z > 3.5), max_dp_lead = max(lead)), by = .(survey, sigma, lambda, q)]
  CP <- upsert(CP, "corner_perm.csv")
  cat("\nCORNER UNDER NOISE, q real vs permuted (variant noise; h >= 0.6, Gamma <= 0.5): robust cells (largest DP lead, pp):\n")
  print(dcast(CP[survey %in% SURVEYS][, cell := sprintf("%d (%.1f)", cells_z3.5, max_dp_lead)], survey + sigma ~ lambda + q, value.var = "cell"))
}

# ---- spontaneous adoption: share of nodes with q' <= Gamma ---------------------------------------------------
SP <- rbindlist(lapply(SURVEYS, function(S) rbindlist(lapply(RUNS(S), function(k) {
  q <- V(NETS[[S]][[k]])$mur_score; L <- eps_L(S, k, length(q)); lo <- q <= 0.5
  rbindlist(lapply(SIGMAS, function(s) { qp <- q - s * L
    data.table(survey = S, run = k, sigma = s, G = GS,
               share_all = vapply(GS, function(G) mean(qp <= G), 0), share_noise = vapply(GS, function(G) mean(qp <= G & q > G), 0),
               share_q05 = vapply(GS, function(G) mean(qp[lo] <= G), 0)) }))
}))))[, lapply(.SD, mean), by = .(survey, sigma, G), .SDcols = c("share_all", "share_noise", "share_q05")]
SP <- upsert(SP, "spontaneous.csv")
cat("\nSPONTANEOUS ADOPTION: % of nodes with q' <= Gamma (all; made by the noise, q' <= Gamma < q; among q <= 0.5), Gamma 0 / 0.1 / 0.2 / 0.3:\n")
print(dcast(SP[survey %in% SURVEYS & round(G, 3) %in% c(0, 0.1, 0.2, 0.3)][, cell := sprintf("%4.1f / %4.1f / %4.1f", 100 * share_all, 100 * share_noise, 100 * share_q05)],
            survey + sigma ~ G, value.var = "cell"))

# ---- which sigma matches ELSOC? -------------------------------------------------------------------------------
# Observed: participation with E = 0 among q <= 0.5 in the S1 sample of 10_ELSOC_rule.R (rebuilt here and
# checked against rule_E0_test.csv). Model: P(q' <= Gamma) = plogis((Gamma - q) / sigma), expected over the
# nodes of the ELSOC networks with q <= 0.5, with (i) Gamma_hat 2018 of the S1 rule (10) and (ii) the
# (Gamma, sigma) fitted on the E = 0 respondents alone (logit of y on q; noise identified by the slope in q).
if ("ELSOC" %in% SURVEYS && file.exists("output/10_ELSOC_rule/rule_ladder.csv")) {
  source("scripts/01_survey_data.R")
  d <- elsoc_panel()
  r <- d[named == TRUE & !is.na(own), .(id, year, w, y = as.integer(own >= 2), G = (eff - 1) / 4, E = ((fam - 1) / 4 + (fri - 1) / 4) / 2, q = q_c09)]
  r <- r[year %in% c(2016, 2018, 2022) & !is.na(G) & !is.na(E) & !is.na(q)]
  e0 <- r[, .(y_E0 = wmean(y[E == 0], w[E == 0])), by = .(qbin = cut(q, c(-0.01, 0.25, 0.5, 0.75, 0.999, 1), labels = c("q <= 0.25", "0.25-0.5", "0.5-0.75", "0.75-1", "q = 1")))]
  ref <- fread("output/10_ELSOC_rule/rule_E0_test.csv")
  message(sprintf("S1 sample rebuilt: %d person-waves (10: 5192); E = 0 participation by q bin equal to rule_E0_test.csv: %s", nrow(r),
                  isTRUE(all.equal(e0[match(ref$qbin, qbin)]$y_E0, ref$y_E0))))
  obs <- r[E == 0 & q <= 0.5, wmean(y, w)]; nobs <- r[E == 0 & q <= 0.5, .N]
  lad <- fread("output/10_ELSOC_rule/rule_ladder.csv")[spec == "S1"]
  G18 <- lad[stat == "Gamma_2018"]$est; sig_hat <- lad[stat == "sigma"]$est
  fit0 <- function(dd) { b <- coef(suppressWarnings(glm(y ~ q, family = binomial, data = dd, weights = w))); c(sigma = -1 / b[["q"]], Gamma = b[["(Intercept)"]] * (-1 / b[["q"]])) }
  r0 <- r[E == 0]; est0 <- fit0(r0)
  ids <- unique(r0$id); setkey(r0, id)
  bs0 <- do.call(rbind, mclapply(1:500, function(b) { set.seed(500 + b); fit0(r0[.(sample(ids, replace = TRUE)), allow.cartesian = TRUE]) }, mc.cores = NCORE))
  qe <- unlist(lapply(NETS$ELSOC, function(g) V(g)$mur_score)); qlo <- qe[qe <= 0.5]
  pred <- function(s, G) mean(plogis((G - qlo) / s))
  s_star <- uniroot(function(s) pred(s, G18) - obs, c(0.01, 5))$root
  CAL <- rbind(
    data.table(what = "observed: ELSOC S1, E = 0 and q <= 0.5 (weighted)", sigma = NA, Gamma = NA, share_pct = 100 * obs, n = nobs),
    data.table(what = "model, Gamma_hat 2018 of the S1 rule", sigma = c(SIGMAS[SIGMAS > 0], sig_hat), Gamma = G18,
               share_pct = 100 * sapply(c(SIGMAS[SIGMAS > 0], sig_hat), pred, G = G18), n = length(qlo)),
    data.table(what = "sigma that matches the observed share, Gamma_hat 2018", sigma = s_star, Gamma = G18, share_pct = 100 * obs, n = length(qlo)),
    data.table(what = sprintf("logit on the E = 0 respondents alone (sigma 95%% CI %.2f-%.2f, Gamma %.2f to %.2f; respondent bootstrap, B = 500)",
                              quantile(bs0[, "sigma"], .025), quantile(bs0[, "sigma"], .975), quantile(bs0[, "Gamma"], .025), quantile(bs0[, "Gamma"], .975)),
               sigma = est0[["sigma"]], Gamma = est0[["Gamma"]], share_pct = 100 * pred(est0[["sigma"]], est0[["Gamma"]]), n = nrow(r0)))
  fwrite(CAL, file.path(OUT, "elsoc_calibration.csv"))
  cat("\nWHICH SIGMA MATCHES ELSOC (participation with no exposure among q <= 0.5)?\n"); print(CAL, digits = 3)
}

# ---- figure: premium vs sigma by survey -------------------------------------------------------------------------
if (!TEST) {
  SL <- c(GSS = "GSS 2004", ATP = "ATP 2014", GSS14 = "GSS 2014", ELSOC = "ELSOC 2017")
  rows <- c("mean 0.5-1.5", "0.7", "1.0", "1.5")
  pf <- PR[survey %in% SURVEYS & lambda %in% rows]
  pf[, row := factor(ifelse(grepl("mean", lambda), "λ-mean 0.5–1.5", paste0("λ = ", lambda)), levels = c("λ-mean 0.5–1.5", "λ = 0.7", "λ = 1.0", "λ = 1.5"))]
  pf[, srv := factor(SL[survey], levels = SL[intersect(names(SL), SURVEYS)])]
  pf[, var_lab := factor(c(noise = "noise only (one adopting contact needed)", noise_spont = "noise + spontaneous adoption")[variant],
                         levels = c("noise only (one adopting contact needed)", "noise + spontaneous adoption"))]
  dg <- position_dodge(width = 0.025)
  p <- ggplot(pf, aes(sigma, premium_pp, colour = var_lab, shape = var_lab, linetype = var_lab)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_line(position = dg, linewidth = 0.5) + geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, position = dg, linewidth = 0.45, linetype = "solid") +
    geom_point(position = dg, size = 1.7) +
    scale_colour_manual(values = c("#2a78d6", "#eb6834"), name = NULL) + scale_shape_manual(values = c(16, 17), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    scale_x_continuous(breaks = SIGMAS, labels = format(SIGMAS)) + expand_limits(y = 0) +
    facet_grid(row ~ srv, scales = "free_y") +
    labs(x = "σ (logistic noise of the rule, units of q; ELSOC estimate 0.48 is an upper bound)", y = "premium, network − DP twin (pp)",
         title = "The structural premium under a noisy rule",
         subtitle = "x2 networks, base exposure, random seeding; mean over the (Γ, h) plane; 95% run bootstrap. σ = 0, noise only = the main result") +
    theme_bw(base_size = 8) + theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.text.y = element_text(size = 6.5))
  ggsave(file.path(FIG, "noisy_rule_premium.pdf"), p, width = 7.5, height = 6.2, device = cairo_pdf)
}
stamp("saved ", OUT, if (!TEST) paste0(" and ", FIG), "; total ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
