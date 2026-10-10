# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How certain is the structural premium? Joint BAM + RUN-LEVEL bootstrap
# 06b_premium_bootstrap.R
#
# The BAM standard errors of 05/06 (~5e-4) treat each of the 1025 (IUL, MSP)
# cells of a run as an independent observation. They are not: the cells of run
# k share the GSS network k, its DP twin and the primary seed. The unit of
# replication is the RUN (96 per topology), and run k uses the same network,
# twin and primary draw in every seeding strategy, so contrasts between
# strategies are PAIRED by run. This script resamples run_id (1..96) with
# replacement — one set of weights shared by GSS, DP and every strategy — and
# refits.
#
# Model (decided 2026-10-10, see docs/MAKING-IT-BOTTOM-UP.tex §premium), per lambda:
#
#   cbind(adopters, fail) ~ strategy + strategy:topology
#                           + te(IUL, MSP, by = strategy, k = c(12, 12))
#
#   - the baseline surface varies by STRATEGY: with one shared surface, each
#     strategy's misfit leaks into the strategy x topology terms and flattens
#     the contrasts (at lambda = 1, central - random falls from -0.20 to -0.07)
#   - it does NOT vary by topology: beta_DP is defined as the constant logit
#     shift between GSS and DP; a topology-specific surface would absorb the
#     premium itself. Its shape (a shift of the tipping frontier) is what the
#     link-free gap and the MSP_c shift of 06 describe
#   - te(12, 12) is where beta_DP converges; the s(k = 10) surface of 05/06
#     understates |beta| by up to ~13% at lambda = 1 and is kept as a
#     continuity column ("s10")
#   - pooled beta_DP = the unweighted mean of the six per-lambda betas (the
#     "pooled" fit of 05 puts one surface under all lambdas and is attenuated)
#
# Two exact shortcuts make B = 1000 cheap:
#   1. The model has no run term, so the binomial likelihood of the stacked
#      cells equals that of the run-aggregated counts: a bootstrap replicate
#      with run multiplicities w is a fit on sum_r w_r * adopters_r per cell
#      (2050 rows per strategy and lambda instead of 196,800).
#   2. The joint likelihood factorizes by strategy (own intercept, surface and
#      DP shift), so the joint fit equals separate per-strategy fits; the joint
#      model is fitted once on the full sample (and checked against the
#      separate fits), the replicates refit strategy by strategy.
#
# Variants: the five strategies with the hub budget (0.40 x the network's
# maximum degree; the default of 05 since 2026-10-10), plus random and marginal
# with the former "own" budget (0.40 x degree of the primary), read from
# _superseded/, to price the seed-size confound. Caveat: closeness and eigen
# were run on 2026-10-05 under the own rule; their primaries are near-hubs
# (28.0 vs 28.7 seeds on average), but because they are not node-matched the
# GSS and DP budgets differ by ~1 seed in ~45% of runs. Pending a rerun.
#
# Also reported:
#   - the link-free premium: mean over the plane of Phi_GSS - Phi_DP, per run,
#     with its bootstrap (and t) interval
#   - the reverse-premium corner: per cell, the paired z over the 96 runs, and
#     a one-sided simultaneous (max-T) threshold per strategy x lambda panel
#     from the same bootstrap; a cell is a ROBUST DP win only if it clears it
#
# Reads  output/05_unified_diffusion/unified_premium_results_<s>.csv
#        output/05_unified_diffusion/_superseded/unified_premium_results_{random,marginal}.csv
# Writes output/06b_premium_bootstrap/premium_bootstrap_{beta,gap,contrasts}.csv
#        output/06b_premium_bootstrap/premium_reverse_corner{,_cells}.csv
#        output/06b_premium_bootstrap/premium_bootstrap_replicates.csv
#        plots/06b_premium_bootstrap/premium_bootstrap.pdf
#
# Usage (from the repo root):
#   Rscript scripts/06b_premium_bootstrap.R               # B = 1000, 8 cores
#   Rscript scripts/06b_premium_bootstrap.R --B=10        # timing run (sandbox)
#   Rscript scripts/06b_premium_bootstrap.R --B=500 --cores=4
#   Rscript scripts/06b_premium_bootstrap.R --plots-only  # redraw from the CSVs
# With B < 200 everything goes to output/06b_premium_bootstrap/_test.
# Runtime: ~1 min of loading + joint fits, then ~8 s per replicate and core.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressWarnings(suppressMessages({
  library(data.table); library(mgcv); library(parallel)
  library(ggplot2); library(patchwork)
}))
Sys.setenv(OMP_NUM_THREADS = "1")

ARGS <- commandArgs(trailingOnly = TRUE)
arg_num <- function(name, default) {
  a <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), ARGS, value = TRUE))
  if (length(a)) as.numeric(a) else default
}
B         <- arg_num("B", 1000)
N_CORES   <- arg_num("cores", 8)
TEST_MODE <- B < 200
PLOTS_ONLY <- "--plots-only" %in% ARGS
BOOT_SEED <- 20261010

RES_DIR   <- "output/05_unified_diffusion"
DATA_OUT  <- if (TEST_MODE) "output/06b_premium_bootstrap/_test" else "output/06b_premium_bootstrap"
PLOTS_OUT <- if (TEST_MODE) "output/06b_premium_bootstrap/_test" else "plots/06b_premium_bootstrap"
dir.create(DATA_OUT,  showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_OUT, showWarnings = FALSE, recursive = TRUE)

MAIN <- c("central", "closeness", "eigen", "random", "marginal")
VARIANTS <- data.table(
  variant  = c(MAIN, "random_own", "marginal_own"),
  strategy = c(MAIN, "random", "marginal"),
  budget   = c("hub", "own", "own", "hub", "hub", "own", "own"),
  path     = c(file.path(RES_DIR, paste0("unified_premium_results_", MAIN, ".csv")),
               file.path(RES_DIR, "_superseded",
                         paste0("unified_premium_results_", c("random", "marginal"), ".csv"))))
VARIANTS <- VARIANTS[file.exists(path)]
stopifnot(all(MAIN %in% VARIANTS$variant))
V <- VARIANTS$variant

# contrasts A - B (paired by run). On beta: negative = A has the LARGER premium;
# on the gap: positive = A has the larger premium.
CONTRASTS <- data.table(
  A = c("central", "random",   "central",   "central", "closeness",
        "random", "marginal", "central", "random_own"),
  B = c("random",  "marginal", "closeness", "eigen",   "eigen",
        "random_own", "marginal_own", "random_own", "marginal_own"),
  label = c("central - random", "random - marginal", "central - closeness",
            "central - eigen", "closeness - eigen",
            "random: hub - own budget", "marginal: hub - own budget",
            "central - random (own budget)", "random - marginal (own budget)"))
CONTRASTS <- CONTRASTS[A %in% V & B %in% V]

# reverse-premium regions (the corner box is the one described in the design note)
REGIONS <- list(
  corner = function(G, h) h >= 0.6 & G >= 0.275 - 1e-9 & G <= 0.45 + 1e-9,
  lowG   = function(G, h) G <= 0.1 + 1e-9,
  other  = function(G, h) !(h >= 0.6 & G >= 0.275 - 1e-9 & G <= 0.45 + 1e-9) & G > 0.1 + 1e-9,
  plane  = function(G, h) rep(TRUE, length(G)))
BAND_PP <- 2      # the corner is counted as cells where DP wins by more than this

line <- function(x) message(paste(rep(x, 78), collapse = ""))
t_all <- Sys.time()

# --------------------------------- plots --------------------------------------
# Defined first so that --plots-only can redraw from the saved CSVs.
# Categorical slots in fixed order, one hue per strategy; budget = line type.
COLS   <- c(central = "#2a78d6", closeness = "#eb6834", eigen = "#1baf7a",
            random = "#eda100", marginal = "#e87ba4")
SHAPES <- c(central = 16, closeness = 17, eigen = 15, random = 18, marginal = 25)

draw_plots <- function(beta_tab, gap_tab, con_tab, cell_tab) {
  num_l <- function(x) as.numeric(as.character(x))
  lams  <- sort(unique(num_l(cell_tab$lambda)))
  base_theme <- theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.title = element_text(size = 10, face = "bold"))
  dodge <- position_dodge(width = 0.06)

  bm <- beta_tab[spec == "te12" & lambda != "mean" & variant %in% MAIN]
  gm <- gap_tab[lambda != "mean" & variant %in% MAIN]
  bm[, `:=`(lam = num_l(lambda), strategy = factor(strategy, MAIN))]
  gm[, `:=`(lam = num_l(lambda), strategy = factor(strategy, MAIN))]
  strat_plot <- function(df, y, title, ylab) {
    ggplot(df, aes(lam, .data[[y]], color = strategy, shape = strategy)) +
      geom_line(linewidth = 0.5, position = dodge) +
      geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.45, linewidth = 0.5,
                      fill = "white", position = dodge) +
      scale_color_manual(values = COLS, name = NULL) +
      scale_shape_manual(values = SHAPES, name = NULL) +
      scale_x_continuous(breaks = lams) +
      labs(title = title, x = "λ", y = ylab) + base_theme
  }
  p_beta <- strat_plot(bm, "beta", "Premium on the odds scale (te12), run-bootstrap 95% CI",
                       expression(beta[DP]))
  p_gap  <- strat_plot(gm, "gap_pp", "Premium on the probability scale, run-bootstrap 95% CI",
                       expression(bar(Phi)[GSS] - bar(Phi)[DP] ~ "(pp)"))

  shown <- c("central - random", "random - marginal", "central - closeness",
             "central - eigen", "closeness - eigen")
  cm <- con_tab[lambda != "mean" & scale %in% c("beta_te12", "gap_pp") & contrast %in% shown]
  cm[, `:=`(lam = num_l(lambda),
            scale = factor(scale, c("beta_te12", "gap_pp"), c("Δβ (te12)", "Δ gap (pp)")),
            contrast = factor(contrast, shown))]
  p_con <- ggplot(cm, aes(lam, est)) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "#2a78d6", alpha = 0.18) +
    geom_line(color = "#2a78d6", linewidth = 0.6) + geom_point(color = "#2a78d6", size = 2) +
    facet_grid(scale ~ contrast, scales = "free_y", switch = "y") +
    scale_x_continuous(breaks = c(0.5, 0.75, 1)) +
    labs(title = "Paired contrasts between strategies (95% CI; Δβ < 0 and Δ gap > 0: the first has the larger premium)",
         x = "λ", y = NULL) +
    base_theme + theme(strip.placement = "outside", strip.text = element_text(size = 8))

  # budget: random and marginal, hub vs own (same hue, budget = line type)
  bb <- rbind(beta_tab[spec == "te12" & lambda != "mean" & strategy %in% c("random", "marginal"),
                       .(strategy, budget, lambda, est = beta, ci_lo, ci_hi, what = "β_DP (te12)")],
              gap_tab[lambda != "mean" & strategy %in% c("random", "marginal"),
                      .(strategy, budget, lambda, est = gap_pp, ci_lo, ci_hi, what = "gap (pp)")])
  bb[, `:=`(lam = num_l(lambda), strategy = factor(strategy, c("random", "marginal")),
            what = factor(what, c("β_DP (te12)", "gap (pp)")))]
  p_bud <- ggplot(bb, aes(lam, est, color = strategy, linetype = budget, shape = budget)) +
    geom_line(linewidth = 0.5, position = dodge) +
    geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = 0.45, linewidth = 0.5,
                    fill = "white", position = dodge) +
    facet_wrap(~ what, scales = "free_y") +
    scale_color_manual(values = COLS[c("random", "marginal")], name = NULL) +
    scale_linetype_manual(values = c(hub = "solid", own = "22"), name = "budget") +
    scale_shape_manual(values = c(hub = 16, own = 21), name = "budget") +
    scale_x_continuous(breaks = lams) +
    labs(title = "Seed budget: hub (0.40 × max degree, ≈28 seeds) vs own (0.40 × degree of the primary)",
         x = "λ", y = NULL) + base_theme

  # where DP robustly beats GSS
  cc <- cell_tab[gap_pp < 0]
  cc[, `:=`(dp_pp = -gap_pp, row = factor(variant, intersect(V, unique(variant))),
            col = paste0("λ = ", num_l(lambda)))]
  p_map <- ggplot(cc, aes(G, h)) +
    geom_tile(aes(fill = dp_pp), width = 0.025, height = 1 / 24) +
    geom_point(data = cc[robust_dp == TRUE], size = 0.3, color = "grey10") +
    facet_grid(row ~ col) +
    scale_fill_gradient(low = "#C7E9C0", high = "#00441B", trans = "sqrt",
                        breaks = c(0.05, 0.5, 2, 5, 10), labels = c("0.05", "0.5", "2", "5", "10"),
                        name = "Φ(DP) − Φ(GSS), pp") +
    scale_x_continuous(limits = c(-0.0125, 1.0125), breaks = c(0, 0.5, 1),
                       labels = c("0", ".5", "1"), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-1 / 48, 1 + 1 / 48), breaks = c(0, 0.5, 1),
                       labels = c("0", ".5", "1"), expand = c(0, 0)) +
    labs(title = "Cells where DP beats GSS (fill); dot = robust DP win (one-sided max-T, familywise 5%, run bootstrap)",
         x = "IUL (Γ)", y = "MSP (h)") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "grey96", color = NA),
          panel.spacing = unit(0.6, "lines"), legend.position = "bottom",
          legend.key.width = unit(1.6, "cm"), plot.title = element_text(size = 9, face = "bold"),
          strip.text = element_text(size = 7))

  cairo_pdf(file.path(PLOTS_OUT, "premium_bootstrap.pdf"), width = 11, height = 7.5, onefile = TRUE)
  print((p_beta | p_gap) / p_con + plot_layout(heights = c(1, 1.1)))
  print(p_bud)
  print(p_map)
  invisible(dev.off())
}

if (PLOTS_ONLY) {
  rd <- function(f) fread(file.path(DATA_OUT, f), colClasses = list(character = "lambda"))
  draw_plots(rd("premium_bootstrap_beta.csv"), rd("premium_bootstrap_gap.csv"),
             rd("premium_bootstrap_contrasts.csv"), rd("premium_reverse_corner_cells.csv"))
  message("Redrew ", file.path(PLOTS_OUT, "premium_bootstrap.pdf"), " from the CSVs in ", DATA_OUT, ".")
  quit(save = "no")
}

# ----------------- load: one (2 x 1025) x 96 count matrix per panel -----------
message(format(Sys.time(), "%H:%M:%S"), "  Loading ", length(V), " variants: ",
        paste(V, collapse = ", "))
N_NODES <- NULL; LAMBDAS <- NULL; RUNS <- NULL
Y <- list()      # Y[[v]][[l]]: adopters, rows = cells (GSS then DP), cols = runs
CELLS <- NULL    # the (G, h) of the 1025 GSS rows (identical order in DP rows)
for (i in seq_len(nrow(VARIANTS))) {
  v <- VARIANTS$variant[i]
  d <- fread(VARIANTS$path[i], select = c("topology", "lambda", "run_id", "social_distance_h",
                                          "innovation_iul_Gamma", "num_adopters", "N_nodes"))
  setnames(d, c("social_distance_h", "innovation_iul_Gamma"), c("h", "G"))
  stopifnot(uniqueN(d$N_nodes) == 1,
            d[, .N, by = .(topology, lambda, run_id, G, h)][, all(N == 1)])
  if (is.null(N_NODES)) {
    N_NODES <- d$N_nodes[1]; LAMBDAS <- sort(unique(d$lambda)); RUNS <- sort(unique(d$run_id))
  }
  stopifnot(d$N_nodes[1] == N_NODES, identical(sort(unique(d$lambda)), LAMBDAS),
            identical(sort(unique(d$run_id)), RUNS))
  Y[[v]] <- list()
  for (l in LAMBDAS) {
    wd <- dcast(d[lambda == l], topology + G + h ~ run_id, value.var = "num_adopters")
    wd[, topology := factor(topology, c("GSS", "DP"))]
    setorder(wd, topology, G, h)
    m <- as.matrix(wd[, as.character(RUNS), with = FALSE])
    stopifnot(!anyNA(m), nrow(m) %% 2 == 0)
    cg <- wd[topology == "GSS", .(G, h)]
    stopifnot(isTRUE(all.equal(cg, wd[topology == "DP", .(G, h)], check.attributes = FALSE)))
    if (is.null(CELLS)) CELLS <- cg else stopifnot(isTRUE(all.equal(CELLS, cg)))
    Y[[v]][[as.character(l)]] <- m
  }
  rm(d, wd); invisible(gc())
}
N_RUNS  <- length(RUNS)
N_CELLS <- nrow(CELLS)
L <- as.character(LAMBDAS)
FRAME <- data.frame(topology = factor(rep(c("GSS", "DP"), each = N_CELLS), c("GSS", "DP")),
                    G = rep(CELLS$G, 2), h = rep(CELLS$h, 2))
message(sprintf("  %d runs x %d cells x 2 topologies x %d lambdas per variant; N = %d.",
                N_RUNS, N_CELLS, length(L), N_NODES))

# -------------------------------- the fits ------------------------------------
SPECS <- list(
  te12 = list(sep   = cbind(A, F) ~ topology + te(G, h, k = c(12, 12)),
              joint = cbind(A, F) ~ variant + variant:topology + te(G, h, by = variant, k = c(12, 12))),
  s10  = list(sep   = cbind(A, F) ~ topology + s(G, h, k = 10),
              joint = cbind(A, F) ~ variant + variant:topology + s(G, h, by = variant, k = 10)))

agg_frame <- function(m, w) {
  fr <- FRAME
  fr$A <- as.vector(m %*% w)
  fr$F <- N_NODES * sum(w) - fr$A
  fr
}
fit_beta <- function(fr, form) {
  f <- tryCatch(bam(form, family = binomial(), data = fr, discrete = TRUE),
                error = function(e) NULL)
  if (is.null(f)) NA_real_ else unname(coef(f)["topologyDP"])
}
# one replicate: beta[variant, lambda, spec] for run multiplicities w
betas_for <- function(w) {
  out <- array(NA_real_, c(length(V), length(L), length(SPECS)),
               dimnames = list(V, L, names(SPECS)))
  for (v in V) for (l in L) {
    fr <- agg_frame(Y[[v]][[l]], w)
    for (s in names(SPECS)) out[v, l, s] <- fit_beta(fr, SPECS[[s]]$sep)
  }
  out
}

# ---- point estimates: the joint BAM, checked against the separate fits -------
message(format(Sys.time(), "%H:%M:%S"), "  Joint BAMs on the full sample (", length(V),
        " variants per lambda)...")
w1 <- rep(1, N_RUNS)
joint <- list()
for (l in L) {
  jd <- do.call(rbind, lapply(V, function(v) cbind(agg_frame(Y[[v]][[l]], w1),
                                                     variant = factor(v, V))))
  for (s in names(SPECS)) {
    f  <- bam(SPECS[[s]]$joint, family = binomial(), data = jd, discrete = TRUE)
    pt <- summary(f)$p.table
    rn <- paste0("variant", V, ":topologyDP")
    joint[[length(joint) + 1]] <- data.table(variant = V, lambda = l, spec = s,
      beta = pt[rn, 1], bam_se = pt[rn, 2], dev_expl = summary(f)$dev.expl)
  }
}
joint <- rbindlist(joint)
t_rep <- Sys.time()
beta_hat <- betas_for(w1)
t_rep <- as.numeric(difftime(Sys.time(), t_rep, units = "secs"))
joint[, beta_sep := beta_hat[cbind(variant, lambda, spec)]]
dev_js <- joint[, max(abs(beta - beta_sep))]
message(sprintf("  Joint vs. separate fits: max |diff| in beta_DP = %.1e (%s).", dev_js,
                if (dev_js < 5e-3) "same model, as expected" else "CHECK"))
message(sprintf("  One replicate takes %.1f s; B = %d on %d cores -> ~%.0f min.",
                t_rep, B, N_CORES, B * t_rep / N_CORES / 60))

# ------------------------------ the bootstrap ---------------------------------
set.seed(BOOT_SEED)
W <- t(replicate(B, tabulate(sample.int(N_RUNS, N_RUNS, replace = TRUE), N_RUNS)))
message(format(Sys.time(), "%H:%M:%S"), "  Bootstrap: B = ", B, " run resamples (shared by ",
        "GSS, DP and every variant), ", N_CORES, " cores...")
t_b <- Sys.time()
reps <- mclapply(seq_len(B), function(b) betas_for(W[b, ]), mc.cores = N_CORES)
bad <- vapply(reps, function(x) inherits(x, "try-error") || is.null(x), logical(1))
if (any(bad)) stop("Bootstrap worker failure in ", sum(bad), " replicates")
BETA <- simplify2array(reps)                 # [variant, lambda, spec, b]
dimnames(BETA)[[4]] <- seq_len(B)
el_b <- as.numeric(difftime(Sys.time(), t_b, units = "mins"))
message(sprintf("  Done in %.1f min (%.2f s per replicate per core); failed fits: %d of %d.",
                el_b, 60 * el_b * N_CORES / B, sum(is.na(BETA)), length(BETA)))

ci_lo <- function(x) unname(quantile(x, 0.025, na.rm = TRUE))
ci_hi <- function(x) unname(quantile(x, 0.975, na.rm = TRUE))
p_boot <- function(x) {
  x <- x[!is.na(x)]
  max(1 / length(x), min(1, 2 * min(mean(x <= 0), mean(x >= 0))))
}
lam_label <- c(L, "mean")

# ---- beta_DP: per lambda and lambda-mean ------------------------------------
beta_tab <- rbindlist(lapply(V, function(v) rbindlist(lapply(names(SPECS), function(s) {
  rbindlist(lapply(lam_label, function(l) {
    if (l == "mean") { est <- mean(beta_hat[v, , s]); bs <- colMeans(BETA[v, , s, ]) }
    else             { est <- beta_hat[v, l, s];      bs <- BETA[v, l, s, ] }
    data.table(variant = v, lambda = l, spec = s, beta = est,
               ci_lo = ci_lo(bs), ci_hi = ci_hi(bs), boot_se = sd(bs, na.rm = TRUE),
               n_fail = sum(is.na(bs)))
  }))
}))))
beta_tab <- merge(beta_tab, joint[, .(variant, lambda, spec, bam_se, dev_expl)],
                  by = c("variant", "lambda", "spec"), all.x = TRUE, sort = FALSE)
beta_tab[, se_inflation := boot_se / bam_se]
beta_tab <- merge(VARIANTS[, .(variant, strategy, budget)], beta_tab, by = "variant", sort = FALSE)

# ---- link-free premium: mean plane gap per run ------------------------------
# D[[v]][[l]]: cells x runs paired differences Phi_GSS - Phi_DP
D <- lapply(Y, function(yv) lapply(yv, function(m)
  (m[seq_len(N_CELLS), , drop = FALSE] - m[N_CELLS + seq_len(N_CELLS), , drop = FALSE]) / N_NODES))
GAPR <- sapply(V, function(v) sapply(L, function(l) colMeans(D[[v]][[l]])), simplify = "array")
dimnames(GAPR) <- list(RUNS, L, V)           # [run, lambda, variant]
gap_boot <- function(x) as.vector(W %*% x) / N_RUNS      # x: one value per run
t_ci <- function(x) mean(x) + c(-1, 1) * qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
run_series <- function(v, l) if (l == "mean") rowMeans(GAPR[, , v]) else GAPR[, l, v]

gap_tab <- rbindlist(lapply(V, function(v) rbindlist(lapply(lam_label, function(l) {
  x <- run_series(v, l); bs <- gap_boot(x); tc <- t_ci(x)
  data.table(variant = v, lambda = l, gap_pp = 100 * mean(x),
             ci_lo = 100 * ci_lo(bs), ci_hi = 100 * ci_hi(bs),
             t_lo = 100 * tc[1], t_hi = 100 * tc[2])
}))))
gap_tab <- merge(VARIANTS[, .(variant, strategy, budget)], gap_tab, by = "variant", sort = FALSE)

# ---- contrasts (paired by run) ----------------------------------------------
con_tab <- rbindlist(lapply(seq_len(nrow(CONTRASTS)), function(i) {
  a <- CONTRASTS$A[i]; b <- CONTRASTS$B[i]
  rbindlist(lapply(lam_label, function(l) {
    rows <- lapply(names(SPECS), function(s) {
      if (l == "mean") {
        est <- mean(beta_hat[a, , s]) - mean(beta_hat[b, , s])
        bs  <- colMeans(BETA[a, , s, ]) - colMeans(BETA[b, , s, ])
      } else {
        est <- beta_hat[a, l, s] - beta_hat[b, l, s]
        bs  <- BETA[a, l, s, ] - BETA[b, l, s, ]
      }
      data.table(scale = paste0("beta_", s), est = est, ci_lo = ci_lo(bs), ci_hi = ci_hi(bs),
                 p_boot = p_boot(bs), t_lo = NA_real_, t_hi = NA_real_)
    })
    x <- run_series(a, l) - run_series(b, l); bs <- gap_boot(x); tc <- t_ci(x)
    rows[[length(rows) + 1]] <- data.table(scale = "gap_pp", est = 100 * mean(x),
      ci_lo = 100 * ci_lo(bs), ci_hi = 100 * ci_hi(bs), p_boot = p_boot(bs),
      t_lo = 100 * tc[1], t_hi = 100 * tc[2])
    cbind(data.table(contrast = CONTRASTS$label[i], A = a, B = b, lambda = l), rbindlist(rows))
  }))
}))

# ---- the reverse-premium corner: cell-level paired z + simultaneous max-T ----
message(format(Sys.time(), "%H:%M:%S"), "  Reverse-premium corner: cell z and max-T thresholds...")
cell_tab <- list(); corner_tab <- list()
for (v in V) for (l in L) {
  dm  <- D[[v]][[l]]
  gap <- rowMeans(dm)
  se  <- apply(dm, 1, sd) / sqrt(N_RUNS)
  ok  <- se > 0                               # cells that vary across runs at all
  z   <- ifelse(ok, gap / se, 0)
  gst <- dm %*% t(W) / N_RUNS                 # cells x B bootstrap gaps
  # one-sided familywise 5%: how negative can the most negative centred cell
  # statistic get by chance alone?
  min_t <- apply((gst[ok, , drop = FALSE] - gap[ok]) / se[ok], 2, min)
  q_sim <- -unname(quantile(min_t, 0.05))
  robust <- ok & z < -q_sim
  cell_tab[[length(cell_tab) + 1]] <- data.table(variant = v, lambda = l, G = CELLS$G, h = CELLS$h,
    gap_pp = 100 * gap, se_pp = 100 * se, z = z, q_sim = q_sim, robust_dp = robust)
  for (r in names(REGIONS)) {
    in_r <- REGIONS[[r]](CELLS$G, CELLS$h)
    dp   <- in_r & gap < 0
    mbs  <- colMeans(gst[in_r, , drop = FALSE])
    corner_tab[[length(corner_tab) + 1]] <- data.table(
      variant = v, lambda = l, region = r, n_cells = sum(in_r),
      n_dp = sum(dp), n_dp_band = sum(in_r & gap < -BAND_PP / 100),
      n_pointwise = sum(dp & z < -2), n_robust = sum(in_r & robust), q_sim = q_sim,
      max_dp_pp = if (any(dp)) -100 * min(gap[dp]) else 0,
      max_dp_pp_robust = if (any(in_r & robust)) -100 * min(gap[in_r & robust]) else 0,
      mean_gap_pp = 100 * mean(gap[in_r]),
      mean_ci_lo = 100 * ci_lo(mbs), mean_ci_hi = 100 * ci_hi(mbs))
  }
}
cell_tab   <- rbindlist(cell_tab)
corner_tab <- rbindlist(corner_tab)
corner_tab[, verdict := fifelse(n_robust == 0, "noise",
                        fifelse(max_dp_pp_robust >= BAND_PP, "robust (> 2 pp)",
                        fifelse(max_dp_pp_robust >= 0.2, "robust, small (0.2-2 pp)",
                                "detectable, negligible (< 0.2 pp)")))]
corner_tab <- merge(VARIANTS[, .(variant, strategy, budget)], corner_tab, by = "variant", sort = FALSE)

# --------------------------------- save ---------------------------------------
fwrite(beta_tab,   file.path(DATA_OUT, "premium_bootstrap_beta.csv"))
fwrite(gap_tab,    file.path(DATA_OUT, "premium_bootstrap_gap.csv"))
fwrite(con_tab,    file.path(DATA_OUT, "premium_bootstrap_contrasts.csv"))
fwrite(corner_tab, file.path(DATA_OUT, "premium_reverse_corner.csv"))
fwrite(cell_tab,   file.path(DATA_OUT, "premium_reverse_corner_cells.csv"))
fwrite(as.data.table(as.table(BETA))[, .(variant = V1, lambda = V2, spec = V3,
                                          b = as.integer(V4), beta = N)],
       file.path(DATA_OUT, "premium_bootstrap_replicates.csv"))

# -------------------------------- report --------------------------------------
fmt_ci <- function(e, lo, hi, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), e, lo, hi)
line("="); message("beta_DP (te12) WITH RUN-BOOTSTRAP 95% CI   (B = ", B, ")"); line("=")
print(dcast(beta_tab[spec == "te12", .(variant, budget, lambda, ci = fmt_ci(beta, ci_lo, ci_hi))],
            variant + budget ~ factor(lambda, lam_label), value.var = "ci"))
line("-"); message("Same, continuity spec s(k = 10) of 05/06"); line("-")
print(dcast(beta_tab[spec == "s10", .(variant, budget, lambda, ci = fmt_ci(beta, ci_lo, ci_hi))],
            variant + budget ~ factor(lambda, lam_label), value.var = "ci"))
line("-"); message("Bootstrap SE / naive BAM SE (te12)"); line("-")
print(dcast(beta_tab[spec == "te12" & lambda != "mean"], variant ~ lambda,
            value.var = "se_inflation"), digits = 3)
line("-"); message("LINK-FREE PREMIUM: mean plane gap (pp), bootstrap CI"); line("-")
print(dcast(gap_tab[, .(variant, budget, lambda, ci = fmt_ci(gap_pp, ci_lo, ci_hi, 2))],
            variant + budget ~ factor(lambda, lam_label), value.var = "ci"))
line("-"); message("CONTRASTS (paired by run), lambda-mean   [beta: < 0 = A larger premium]"); line("-")
print(con_tab[lambda == "mean", .(contrast, scale, est = round(est, 3),
                                  ci = sprintf("[%.3f, %.3f]", ci_lo, ci_hi), p_boot = signif(p_boot, 2))])
line("-"); message("CONTRASTS per lambda: beta (te12) and gap (pp)"); line("-")
print(dcast(con_tab[scale %in% c("beta_te12", "gap_pp"),
                    .(contrast, scale, lambda, ci = fmt_ci(est, ci_lo, ci_hi, 2))],
            contrast + scale ~ factor(lambda, lam_label), value.var = "ci"))
line("-"); message("REVERSE-PREMIUM CORNER (robust = clears the max-T threshold)"); line("-")
print(corner_tab[region != "plane", .(variant, lambda, region, n_dp, n_dp_band, n_pointwise,
                                      n_robust, max_dp_pp = round(max_dp_pp, 3),
                                      max_robust = round(max_dp_pp_robust, 3),
                                      mean_gap = fmt_ci(mean_gap_pp, mean_ci_lo, mean_ci_hi, 2),
                                      verdict)], nrows = 500)

draw_plots(beta_tab, gap_tab, con_tab, cell_tab)

message("\n", format(Sys.time(), "%H:%M:%S"), "  Saved CSVs in ", DATA_OUT, " and ",
        file.path(PLOTS_OUT, "premium_bootstrap.pdf"), " (total ",
        round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min).")
