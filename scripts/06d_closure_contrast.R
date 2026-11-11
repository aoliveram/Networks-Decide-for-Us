# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What does closure do to the premium? x2 vs x1, paired. 06d_closure_contrast.R
#
# x1 and x2 have the same nodes, degrees and mixing; x2 has twice the triangles
# of every node. Both are compared with the SAME DP twin (rewired from x2) and
# share the seed cluster's primary and the simulation seeds, so per run
#   premium(x2) - premium(x1) = mean over the plane of Phi(x2) - Phi(x1),
# the effect of closure alone. Per survey, seeding strategy and lambda (and the
# lambda mean): the premium of each level (mean over the plane of Phi_NET -
# Phi_DP, pp), the paired difference, and the same over the selective half of
# the plane (h <= 0.5), with 95% run-bootstrap intervals (B = 2000; one set of
# run weights for every statistic).
#
# Sensitivity: if the sweep of a third level exists (x3: three times the triangles of
# every node, 02_network_calibrated.R --closure=3; 05 --level=x3, any subset of
# conditions), its premium and x3 - x1 are added (same twin, same seeds).
#
# Reads  output/05_unified_diffusion/<S>_{x2,x1[,x3]}/unified_premium_results_<seeding>.csv
# Writes output/06d_closure_contrast/<S>_closure_contrast.csv
#        plots/06d_closure_contrast/<S>_closure_contrast.pdf
# Usage: Rscript scripts/06d_closure_contrast.R --survey=GSS     (seconds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
ARGS <- commandArgs(trailingOnly = TRUE)
SURVEY <- { a <- grep("^--survey=", ARGS, value = TRUE); if (length(a)) toupper(sub(".*=", "", a[1])) else "GSS" }
B <- 2000L
STRATS <- c("central", "closeness", "eigen", "random", "marginal")
DATA_OUT <- "output/06d_closure_contrast"; PLOTS_OUT <- "plots/06d_closure_contrast"
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS_OUT, recursive = TRUE, showWarnings = FALSE)

read_level <- function(l, s) {
  f <- file.path("output/05_unified_diffusion", paste0(SURVEY, "_", l), paste0("unified_premium_results_", s, ".csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = c("topology", "lambda", "run_id", "social_distance_h", "innovation_iul_Gamma", "num_adopters", "N_nodes"))
  d[topology == "GSS", topology := "NET"][, `:=`(level = l, seeding = s)][]
}
d <- rbindlist(lapply(STRATS, function(s) rbindlist(list(read_level("x2", s), read_level("x1", s), read_level("x3", s)))))
if (!nrow(d) || !all(c("x1", "x2") %in% d$level)) stop("need the x2 and x1 sweeps of ", SURVEY)
X3 <- "x3" %in% d$level
if (X3) d <- d[seeding %in% unique(d[level == "x3", seeding])]   # strategies with all three levels
d[, phi := num_adopters / N_nodes]
w <- dcast(d, seeding + lambda + run_id + social_distance_h + innovation_iul_Gamma ~ level + topology, value.var = "phi")
same_dp <- w[, max(abs(x2_DP - x1_DP), if (X3) abs(x3_DP - x2_DP) else 0)]
message(sprintf("%s: max |Phi_DP(x2 sweep) - Phi_DP(x1 sweep)| = %.2g (the twin is shared)", SURVEY, same_dp))
stopifnot(same_dp < 1e-12)
w[, `:=`(p2 = x2_NET - x2_DP, p1 = x1_NET - x2_DP, sel = social_distance_h <= 0.5 + 1e-9)]
if (!X3) w[, p3 := NA_real_] else w[, p3 := x3_NET - x2_DP]
per_run <- w[, .(premium_x2 = mean(p2), premium_x1 = mean(p1), diff = mean(p2 - p1),
                 diff_selective = mean((p2 - p1)[sel]), premium_x3 = mean(p3), diff3 = mean(p3 - p1)),
             by = .(seeding, lambda, run_id)]
per_run <- rbind(per_run, per_run[, lapply(.SD, mean), by = .(seeding, run_id), .SDcols = -"lambda"][, lambda := NA_real_], use.names = TRUE)

runs <- sort(unique(per_run$run_id)); R <- length(runs)
set.seed(20261105)
Wt <- t(replicate(B, tabulate(sample.int(R, R, replace = TRUE), R))) / R
ci <- function(v, r) { v <- v[order(r)]; stopifnot(identical(sort(r), runs)); b <- drop(Wt %*% v); c(mean(v), quantile(b, c(.025, .975), names = FALSE)) }
summ <- per_run[, {
  a <- ci(premium_x2, run_id); b <- ci(premium_x1, run_id); c0 <- ci(diff, run_id); s <- ci(diff_selective, run_id)
  a3 <- if (X3) ci(premium_x3, run_id) else rep(NA_real_, 3); d3 <- if (X3) ci(diff3, run_id) else rep(NA_real_, 3)
  .(premium_x2 = 100 * a[1], x2_lo = 100 * a[2], x2_hi = 100 * a[3], premium_x1 = 100 * b[1], x1_lo = 100 * b[2], x1_hi = 100 * b[3],
    diff_pp = 100 * c0[1], diff_lo = 100 * c0[2], diff_hi = 100 * c0[3], rel_diff = c0[1] / b[1],
    diff_selective_pp = 100 * s[1], sel_lo = 100 * s[2], sel_hi = 100 * s[3],
    premium_x3 = 100 * a3[1], x3_lo = 100 * a3[2], x3_hi = 100 * a3[3],
    diff3_pp = 100 * d3[1], diff3_lo = 100 * d3[2], diff3_hi = 100 * d3[3], rel_diff3 = d3[1] / b[1], runs = .N)
}, by = .(seeding, lambda)]
summ[, lambda_lab := ifelse(is.na(lambda), "mean", sprintf("%.1f", lambda))]
setorder(summ, seeding, lambda, na.last = TRUE)
fwrite(summ, file.path(DATA_OUT, paste0(SURVEY, "_closure_contrast.csv")))
f3 <- function(x, lo, hi) sprintf("%.2f [%.2f, %.2f]", x, lo, hi)
cat("\n", SURVEY, ": premium (pp) of x2 and x1 against the same DP twin, and x2 - x1 (paired; 95% run-bootstrap CI)\n", sep = "")
print(summ[, .(seeding, lambda = lambda_lab, x2 = f3(premium_x2, x2_lo, x2_hi), x1 = f3(premium_x1, x1_lo, x1_hi),
               x2_minus_x1 = f3(diff_pp, diff_lo, diff_hi), relative = sprintf("%+.1f%%", 100 * rel_diff),
               selective_half = f3(diff_selective_pp, sel_lo, sel_hi))], row.names = FALSE, nrows = 100)
if (X3) {
  cat("\nSensitivity, three times the triangles of every node (x3) against the same twin:\n")
  print(summ[, .(seeding, lambda = lambda_lab, x3 = f3(premium_x3, x3_lo, x3_hi), x3_minus_x1 = f3(diff3_pp, diff3_lo, diff3_hi),
                 relative = sprintf("%+.1f%%", 100 * rel_diff3))], row.names = FALSE, nrows = 100)
}

pd <- summ[!is.na(lambda)]
p <- ggplot(pd, aes(lambda, diff_pp, color = seeding)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(ymin = diff_lo, ymax = diff_hi), size = 0.3, position = position_dodge(width = 0.05)) +
  geom_line(position = position_dodge(width = 0.05), linewidth = 0.4) +
  scale_color_manual(values = c(central = "#2a78d6", closeness = "#eb6834", eigen = "#1baf7a", random = "#eda100", marginal = "#e87ba4"), name = NULL) +
  labs(title = sprintf("%s: what closure does to the premium, x2 − x1 against the same DP twin (95%% run-bootstrap CI)", SURVEY),
       x = "λ", y = "Premium(x2) − premium(x1), pp") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(PLOTS_OUT, paste0(SURVEY, "_closure_contrast.pdf")), p, width = 7, height = 4, device = cairo_pdf)
message("Saved ", file.path(DATA_OUT, paste0(SURVEY, "_closure_contrast.csv")), " and the figure")
