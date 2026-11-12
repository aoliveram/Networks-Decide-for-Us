# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Does the premium survive the exposure axis? 06c_exposure_axis.R
#
# The main text uses the base exposure E = a / d. The supplement varies how
# adopting contacts add up with one formula and two dials (05_unified_engine.R):
#   theta  independent reinforcement (Kish discount of adopters in the same
#          race x religion circle), c = -2 theta / (1 + theta)
#   kappa  cohesive reinforcement (a bonus per pair of adopting contacts who know
#          each other), c = 2 kappa
# c = (value of the second adopter) - 1 places both arms on one signed axis:
# c = -1, -0.5 (theta = 1, 1/3), 0 (base), +0.5, +1 (kappa = 0.25, 0.5).
# Each arm level has a permuted control (same circle sizes / same number of
# pairs per ego, random members): real - permuted is the part of the effect
# that comes from WHERE the circles and pairs are, not from how many.
# Sets whose nodes carry ideology (ELSOC, GSS14; PR #5) have a second theta arm,
# arm "theta_ideol": circles = ideology groups (c = -1, -0.5, with permuted controls).
#
# Per set, seeding, lambda and condition (unit = network run, paired):
#   premium_pp     mean over the plane of Phi_network - Phi_DP (pp)
#   selective_pp   the same over h <= 0.5 (where the gate filters)
#   vs_base_pp     premium(condition) - premium(base)
#   vs_perm_pp     premium(real arm) - premium(its permuted control)
# with 95% percentile intervals from a run bootstrap (B = 2000; one set of run
# weights shared by every statistic of a panel).
#
# Reads  output/05_unified_diffusion/<set>/exposure_axis_<seeding>.rds (every
#        seeding found)
# Writes output/06c_exposure_axis/<set>/exposure_axis.csv
#        plots/06c_exposure_axis/<set>/exposure_axis.pdf
# Usage: Rscript scripts/06c_exposure_axis.R --set=GSS_x2      (seconds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
ARGS <- commandArgs(trailingOnly = TRUE)
SET <- sub("^--set=", "", grep("^--set=", ARGS, value = TRUE)); if (length(SET) == 0) SET <- "GSS_x2"
B <- 2000L
IN <- file.path("output/05_unified_diffusion", SET)
DATA_OUT <- file.path("output/06c_exposure_axis", SET); PLOTS_OUT <- file.path("plots/06c_exposure_axis", SET)
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS_OUT, recursive = TRUE, showWarnings = FALSE)
files <- list.files(IN, pattern = "^exposure_axis_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no exposure_axis_<seeding>.rds in ", IN, ": run 05 first")

per_run <- rbindlist(lapply(files, function(f) {
  x <- readRDS(f); h_cell <- rep(x$H, each = length(x$G)); sel <- h_cell <= 0.5 + 1e-9
  rbindlist(lapply(names(x$res), function(cn) {
    a <- x$res[[cn]]                                  # cells x topology x lambda x output x run
    rbindlist(lapply(seq_along(x$lambdas), function(li) {
      net <- a[, "NET", li, "num_adopters", , drop = TRUE] / x$N
      dp  <- a[, "DP",  li, "num_adopters", , drop = TRUE] / x$N
      data.table(set = x$set, seeding = x$seeding, lambda = x$lambdas[li], condition = cn, run = seq_len(ncol(net)),
                 premium = colMeans(net - dp), selective = colMeans((net - dp)[sel, , drop = FALSE]),
                 phi_net = colMeans(net), phi_dp = colMeans(dp))
    }))
  }))[, `:=`(N = x$N)][]
}))
meta <- as.data.table(readRDS(files[1])$conditions)
per_run <- merge(per_run, meta[, .(condition, theta, kappa, c, control_of, circles)], by = "condition")
per_run[, arm := fifelse(condition == "base", "base",
                         fifelse(theta > 0, fifelse(grepl("^ideol", circles), "theta_ideol", "theta"), "kappa"))]
per_run[, kind := fifelse(nzchar(control_of), "permuted", "real")]
base <- per_run[condition == "base", .(seeding, lambda, run, base = premium)]
per_run <- merge(per_run, base, by = c("seeding", "lambda", "run"))
perm <- per_run[kind == "permuted", .(seeding, lambda, run, condition = control_of, perm = premium)]
per_run <- merge(per_run, perm, by = c("seeding", "lambda", "run", "condition"), all.x = TRUE)

# run bootstrap: one weight matrix (B x runs) shared by every statistic
runs <- sort(unique(per_run$run)); R <- length(runs)
set.seed(20261105)
Wt <- t(replicate(B, tabulate(sample.int(R, R, replace = TRUE), R))) / R
ci <- function(v, r) { stopifnot(identical(sort(r), runs)); v <- v[order(r)]; b <- drop(Wt %*% v)
  c(mean(v), quantile(b, c(.025, .975), names = FALSE)) }
fmt <- function(v) sprintf("%.2f [%.2f, %.2f]", 100 * v[1], 100 * v[2], 100 * v[3])
summ <- per_run[, {
  p <- ci(premium, run); d <- ci(premium - base, run); pp <- if (all(!is.na(perm))) ci(premium - perm, run) else rep(NA_real_, 3)
  .(c = c[1], theta = theta[1], kappa = kappa[1], arm = arm[1], kind = kind[1],
    premium_pp = 100 * p[1], premium_lo = 100 * p[2], premium_hi = 100 * p[3],
    premium_ci = fmt(p), selective_pp = 100 * mean(selective),
    vs_base_ci = if (condition[1] == "base") NA_character_ else fmt(d),
    vs_perm_ci = if (any(is.na(pp))) NA_character_ else fmt(pp),
    phi_net = 100 * mean(phi_net), phi_dp = 100 * mean(phi_dp), runs = .N)
}, by = .(set, seeding, lambda, condition)]
setorder(summ, seeding, lambda, c, kind)
fwrite(summ, file.path(DATA_OUT, "exposure_axis.csv"))
cat("\nEXPOSURE AXIS,", SET, "(premium and contrasts in pp; 95% run-bootstrap intervals):\n")
print(summ[lambda %in% c(0.5, 0.7, 1.0), .(seeding, lambda, condition, c, premium_ci, vs_base_ci, vs_perm_ci)], nrows = 500)
cat("\nPremium range over every seeding, lambda and condition: ",
    sprintf("%.2f to %.2f pp; conditions with a negative premium: %d of %d", min(summ$premium_pp), max(summ$premium_pp),
            sum(summ$premium_pp < 0), nrow(summ)), "\n", sep = "")

# ---- figure: premium along c, by seeding (rows) and lambda (columns) ----
pd <- summ[lambda %in% c(0.5, 0.7, 1.0)]
pd[, lab := sprintf("λ = %.1f", lambda)]
p <- ggplot(pd, aes(c, premium_pp, color = arm, shape = kind)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_linerange(aes(ymin = premium_lo, ymax = premium_hi), linewidth = 0.4,
                 position = position_dodge(width = 0.12)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.12)) +
  facet_grid(seeding ~ lab) +
  scale_color_manual(values = c(base = "black", theta = "#2166AC", theta_ideol = "#1B7837", kappa = "#B2182B"),
                     labels = c(base = "base", theta = "θ (independent)", theta_ideol = "θ, ideology circles", kappa = "κ (cohesive)"),
                     name = NULL) +
  scale_shape_manual(values = c(real = 16, permuted = 1), name = NULL) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(title = sprintf("%s: structural premium along the exposure axis (95%% run-bootstrap CI)", SET),
       x = "c = value of the second adopter − 1   (θ arm: c < 0; κ arm: c > 0)",
       y = "Premium, Φ(network) − Φ(DP), pp") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(PLOTS_OUT, "exposure_axis.pdf"), p, width = 8, height = 1.6 + 1.7 * uniqueN(pd$seeding), device = cairo_pdf)
message("Saved ", file.path(DATA_OUT, "exposure_axis.csv"), " and ", file.path(PLOTS_OUT, "exposure_axis.pdf"))
