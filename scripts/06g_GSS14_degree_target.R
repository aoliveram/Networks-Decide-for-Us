# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSS 2014: does the degree target by race change the premium? 06g_GSS14_degree_target.R
#
# The GSS14 networks (like ATP's) take the GSS 2004 targets: relative number of confidants
# by race (White 1.09, Black 0.68, Hispanic 0.90, Asian 0.84) and same-race shares. GSS 2014
# has PEOCNTCT, contacts in a typical weekday, which is closer to the ~30 contacts of a node
# but gives flatter ratios (1.05, 0.88, 0.94, 0.76; 03_GSS14_MUR_calculation.R). This script
# compares the premium of the main GSS14 x2 networks with x2 networks built with the
# PEOCNTCT ratios (02_network_calibrated.R --survey=GSS14 --rel-deg=peocntct --tag=peocntct;
# same nodes and seeds, same same-race shares): base exposure, random seeding, the 41 x 25
# grid, lambda 0.5-1.0, 24 networks, the seeds of 05 (the main set reproduces the sweep).
# Unit = network run; the difference is paired by run (same node draw, same seeds).
# Outputs: output/06g_GSS14_degree_target/premium_by_target.csv
# Usage (repo root): Rscript scripts/06g_GSS14_degree_target.R   (~1 min)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(parallel) }))
source("scripts/05_unified_engine.R")
OUT <- "output/06g_GSS14_degree_target"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ALT <- "data/02_calibrated_networks/peocntct/GSS14_x2.rds"
if (!file.exists(ALT)) system2(file.path(R.home("bin"), "Rscript"),
                               c("scripts/02_network_calibrated.R", "--survey=GSS14", "--rel-deg=peocntct", "--tag=peocntct"))
SETS <- list(confidants_2004 = readRDS("data/02_calibrated_networks/GSS14_x2.rds"), peocntct_2014 = readRDS(ALT))
H <- seq(0, 1, by = 1/24); G <- seq(0, 1, by = 0.025); LAMS <- seq(0.5, 1.0, by = 0.1)
stopifnot(all(vapply(seq_along(SETS[[1]]), function(k) identical(V(SETS[[1]][[k]])$age, V(SETS[[2]][[k]])$age), TRUE)))   # same nodes
sim <- function(target, k) {
  nt <- make_pair(SETS[[target]][[k]], k, "random")
  rbindlist(lapply(c("NET", "DP"), function(topo) {
    top <- nt[[topo]]; seeds <- draw_seeds(top, k, topo, "random"); W <- gate_matrix(top, H)
    ctx <- build_exposure_ctx(top, k, as.integer(topo == "DP"), EXPOSURE_CONDITIONS["base"])
    rbindlist(lapply(LAMS, function(l) data.table(target = target, run = k, topo = topo, lambda = l,
      phi = mean(run_condition(top, W, ctx, EXPOSURE_CONDITIONS$base, nt$q, seeds, k, l, topo, H, G)[, "num_adopters"]) / nt$n,
      deg_black = mean(top$deg_raw[V(top$g)$race == "Black"]) / mean(top$deg_raw))))
  }))
}
res <- rbindlist(lapply(names(SETS), function(t) rbindlist(mclapply(seq_along(SETS[[t]]), function(k) sim(t, k), mc.cores = 8))))
w <- dcast(res, target + run + lambda ~ topo, value.var = "phi")[, gap := 100 * (NET - DP)]
lm_rows <- rbind(w, w[, .(NET = mean(NET), DP = mean(DP), gap = mean(gap)), by = .(target, run)][, lambda := NA_real_], fill = TRUE)
set.seed(20140601)
B <- t(replicate(2000, tabulate(sample.int(24, 24, TRUE), 24))) / 24
ci <- function(v) { b <- drop(B %*% v); c(mean(v), quantile(b, c(.025, .975), names = FALSE)) }
summ <- lm_rows[, { p <- ci(gap[order(run)]); .(premium_pp = p[1], lo = p[2], hi = p[3], phi_net = 100 * mean(NET), phi_dp = 100 * mean(DP)) },
                by = .(target, lambda)]
dd <- dcast(lm_rows, run + lambda ~ target, value.var = "gap")[, d := peocntct_2014 - confidants_2004]
diff <- dd[, { p <- ci(d[order(run)]); .(target = "peocntct - confidants", premium_pp = p[1], lo = p[2], hi = p[3]) }, by = lambda]
out <- rbind(summ, diff, fill = TRUE)
fwrite(out, file.path(OUT, "premium_by_target.csv"))
cat("\nGSS14 x2, base exposure, random seeding: premium (pp, 95% run-bootstrap CI) by degree target (lambda NA = lambda-mean):\n")
print(out[, .(target, lambda, premium = sprintf("%.2f [%.2f, %.2f]", premium_pp, lo, hi))], nrows = 50)
cat("Relative degree of Black nodes: ", paste(sprintf("%s %.3f", names(SETS), res[topo == "NET" & lambda == 0.5, mean(deg_black), by = target]$V1), collapse = ", "), "\n")
f <- "output/05_unified_diffusion/GSS14_x2/unified_premium_results_random.csv"
if (file.exists(f)) {           # the main set must reproduce the sweep
  sw <- fread(f)[, .(phi = mean(num_adopters) / 1000), by = .(run_id, topology, lambda)]
  m <- merge(res[target == "confidants_2004"], sw, by.x = c("run", "topo", "lambda"), by.y = c("run_id", "topology", "lambda"))
  message("main set vs the sweep: max |difference| in mean adoption = ", signif(max(abs(m$phi.x - m$phi.y)), 3))
}
