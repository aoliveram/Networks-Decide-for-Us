# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOW MUCH DOES EACH CHANGE OF THE NETWORK RECIPE MOVE THE PREMIUM? 06h_calibration_steps.R
#
# On 2026-11-11 the recipe of the calibrated networks changed (02_network_calibrated.R): survey
# weights in every survey, own weighted case-control coefficients for the GSS family (instead of
# Smith et al. 2014), and homophily BY GROUP for religion and sex (Smith 2012: every categorical
# dimension of the confidant data). This script takes the x2 networks of every step and computes
# the premium with base exposure, random seeding, the 41 x 25 grid, lambda 0.5-1.0 and 24
# networks, with the seeds of 05:
#   GSS    before   race by group, Smith et al. coefficients, unweighted (the networks of PR #4)
#          race_uw  race by group, own GSS 2004 case-control coefficients, unweighted
#          race_w   the same, weighted (WTSSNR nodes and targets)
#          full     + religion and sex by group (the main networks from 2026-11-11)
#   ATP    the same four steps (ATP nodes, weight_ATP; targets of GSS 2004)
#   GSS14  before   race by group, Smith et al. coefficients, weighted (WTSSNR; earlier on 2026-11-11)
#          full     the main networks
#   ELSOC  before   education and ideology by group (2026-11-11, earlier in the day)
#          eir      + religion by group
#          full     + sex (one pooled term; activity by sex)
# Steps in data/02_calibrated_networks/steps/ (not tracked; 02 --tag=steps/<name> rebuilds them,
# the 'before' files are copies of the earlier tracked networks, also in git history).
# Unit = network run; differences between consecutive steps are paired by run (same seeds; the
# node draw changes when weights change).
# Outputs: output/06h_calibration_steps/premium_by_step.csv
# Usage (repo root): Rscript scripts/06h_calibration_steps.R   (~4 min)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(parallel) }))
source("scripts/05_unified_engine.R")
OUT <- "output/06h_calibration_steps"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ST <- "data/02_calibrated_networks/steps"
STEPS <- list(
  GSS = c(before = file.path(ST, "GSS_before_x2.rds"), race_uw = file.path(ST, "race_uw/GSS_x2.rds"),
          race_w = file.path(ST, "race_w/GSS_x2.rds"), full = "data/02_calibrated_networks/GSS_x2.rds"),
  ATP = c(before = file.path(ST, "ATP_before_x2.rds"), race_uw = file.path(ST, "race_uw/ATP_x2.rds"),
          race_w = file.path(ST, "race_w/ATP_x2.rds"), full = "data/02_calibrated_networks/ATP_x2.rds"),
  GSS14 = c(before = file.path(ST, "GSS14_before_x2.rds"), full = "data/02_calibrated_networks/GSS14_x2.rds"),
  ELSOC = c(before = file.path(ST, "ELSOC_before_x2.rds"), eir = file.path(ST, "eir/ELSOC_x2.rds"),
            full = "data/02_calibrated_networks/ELSOC_x2.rds"))
H <- seq(0, 1, by = 1/24); G <- seq(0, 1, by = 0.025); LAMS <- seq(0.5, 1.0, by = 0.1)
sim <- function(nets, k) {
  nt <- make_pair(nets[[k]], k, "random")
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, H)
    ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    rbindlist(lapply(LAMS, function(l) data.table(run = k, topo = topo, lambda = l,
      phi = mean(run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, nt$q, seeds, k, l, topo, H, G)[, "num_adopters"]) / nt$n)))
  }))
}
set.seed(20261111)
BW <- t(replicate(2000, tabulate(sample.int(24, 24, TRUE), 24))) / 24
ci <- function(v) { b <- drop(BW %*% v); c(mean(v), quantile(b, c(.025, .975), names = FALSE)) }
OUT_T <- list()
for (S in names(STEPS)) {
  res <- rbindlist(lapply(names(STEPS[[S]]), function(st) { nets <- readRDS(STEPS[[S]][[st]]); stopifnot(length(nets) == 24)
    message(format(Sys.time(), "%H:%M:%S"), "  ", S, " ", st)
    rbindlist(mclapply(seq_along(nets), function(k) sim(nets, k), mc.cores = 8))[, step := st] }))
  w <- dcast(res, step + run + lambda ~ topo, value.var = "phi")[, gap := 100 * (NET - DP)]
  lm <- w[, .(gap = mean(gap), NET = mean(NET), DP = mean(DP)), by = .(step, run)]
  summ <- lm[, { p <- ci(gap[order(run)]); .(premium_pp = p[1], lo = p[2], hi = p[3], phi_net = 100 * mean(NET), phi_dp = 100 * mean(DP)) }, by = step]
  st <- names(STEPS[[S]])
  dd <- rbindlist(lapply(seq_along(st)[-1], function(j) { a <- lm[step == st[j - 1]][order(run)]$gap; b <- lm[step == st[j]][order(run)]$gap
    p <- ci(b - a); data.table(step = paste(st[j], "-", st[j - 1]), premium_pp = p[1], lo = p[2], hi = p[3]) }))
  OUT_T[[S]] <- rbind(summ, dd, fill = TRUE)[, survey := S]
}
out <- rbindlist(OUT_T); setcolorder(out, "survey")
fwrite(out, file.path(OUT, "premium_by_step.csv"))
cat("\nPremium (pp, x2, base exposure, random seeding, lambda-mean; 95% run-bootstrap CI) by step of the recipe:\n")
print(out[, .(survey, step, premium = sprintf("%.2f [%.2f, %.2f]", premium_pp, lo, hi),
              phi_net = round(phi_net, 1), phi_dp = round(phi_dp, 1))], nrows = 50)
