# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIT AND FIGURES of the calibrated networks (02_network_calibrated.R)
# 02c_network_calibrated_tests.R — the counterpart of 02_*_ergm_tests.R and
# 02b_network_plots_EN.R for data/02_calibrated_networks/
#
# Per survey (GSS 2004 nodes, ATP 2014 nodes; both with GSS 2004 targets):
#   1. targets: relative degree and same-race share by race of ego, GSS 2004
#      (95% bootstrap CI over respondents) vs the 24 networks (mean, min-max).
#      x1 and the rewired levels share them by construction
#   2. race mixing: share of the ties of each ego group that go to each group,
#      GSS 2004 confidants vs the networks. Only the diagonal is a target; the
#      off-diagonal shares follow from group sizes and activity
#   3. structure of every level found (<S>_x1.rds, <S>_x2.rds and any other
#      rewired level, e.g. x2n): transitivity, local clustering by degree
#      tercile, degree assortativity, mean distance, giant component, homophily
#      on the five axes, with the spread over networks; degree distribution by
#      race and local clustering C(k)
#   4. an example: network 1, the ego network of the same node in x1 and in each
#      rewired level (same degree; the rewiring changes who the contacts are and
#      how many of them know each other), and the whole x2 network by race
#
# Outputs
#   output/02c_network_calibrated_tests/<S>_targets.csv, <S>_mixing.csv, <S>_structure.csv
#   plots/02c_network_calibrated_tests/<S>_targets_mixing.pdf, <S>_degree_Ck.pdf, <S>_example.pdf
# Usage (repo root): Rscript scripts/02c_network_calibrated_tests.R [--survey=GSS|ATP]   (~1 min each)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({
  library(igraph); library(data.table); library(haven); library(dplyr); library(ggplot2); library(patchwork)
  library(parallel); library(ggraph); library(graphlayouts)
}))
ARGS <- commandArgs(trailingOnly = TRUE)
SURVEYS <- { a <- grep("^--survey=", ARGS, value = TRUE); if (length(a)) toupper(sub(".*=", "", a)) else c("GSS", "ATP") }
NET_DIR <- "data/02_calibrated_networks"; OUT <- "output/02c_network_calibrated_tests"; PLOTS <- "plots/02c_network_calibrated_tests"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
RACES <- c("White", "Black", "Hispanic", "Asian", "Other"); CAL <- RACES[1:4]
COL_RACE <- c(White = "#7F8C8D", Black = "#C0392B", Hispanic = "#E67E22", Asian = "#2980B9", Other = "#8E44AD")
COL_LEVEL <- c(x1 = "#2980B9", x2 = "#E74C3C", x2n = "#8E44AD", x2g = "#E74C3C", x3 = "#E67E22", x3n = "#16A085")
B_BOOT <- 500L
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
theme_set(theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, face = "bold")))

# ---- GSS 2004 confidants: ego race x alter race, with a bootstrap over respondents -----------
g4 <- read_dta("data/02_GSS_network_ergm/GSS_2004_NORC.dta") |> rename_with(tolower)
race_gss <- function(d) with(d, case_when(hispanic %in% c(2:6, 10, 15, 20:24, 30, 40, 41) ~ "Hispanic",   # as in 02
                                          racecen1 == 1 ~ "White", racecen1 == 2 ~ "Black",
                                          racecen1 %in% c(3, 11, 15) ~ "Other", racecen1 %in% c(4:10) ~ "Asian",
                                          racecen1 == 16 ~ "Hispanic", TRUE ~ NA_character_))
resp <- data.table(id = seq_len(nrow(g4)), ego = race_gss(g4), numgiven = as.numeric(g4$numgiven))[!is.na(ego) & !is.na(numgiven)]
alt <- rbindlist(lapply(1:5, function(k) data.table(id = seq_len(nrow(g4)),
  alter = c("Asian", "Black", "Hispanic", "White", "Other")[as.numeric(g4[[paste0("race", k)]])])))[!is.na(alter)]
alt <- merge(alt, resp[, .(id, ego)], by = "id")
gss_mix <- function(ids) {
  a <- alt[.(ids), on = "id", nomatch = 0L, allow.cartesian = TRUE]
  a[, .N, by = .(ego, alter)][, share := N / sum(N), by = ego][CJ(ego = RACES, alter = RACES), on = .(ego, alter)][is.na(share), share := 0][]
}
mix_gss <- gss_mix(resp$id)
set.seed(2004)
mix_boot <- rbindlist(lapply(seq_len(B_BOOT), function(b) gss_mix(sample(resp$id, replace = TRUE))[, b := b]))
mix_gss <- merge(mix_gss[, .(ego, alter, gss = share)],
                 mix_boot[, .(gss_lo = quantile(share, .025), gss_hi = quantile(share, .975)), by = .(ego, alter)], by = c("ego", "alter"))

# ---- per network statistics ---------------------------------------------------------------------
net_mix <- function(g) {
  el <- as_edgelist(g, names = FALSE); r <- V(g)$race
  data.table(ego = factor(c(r[el[, 1]], r[el[, 2]]), RACES), alter = factor(c(r[el[, 2]], r[el[, 1]]), RACES))[
    , .N, by = .(ego, alter)][, share := N / sum(N), by = ego][]
}
net_struct <- function(g) {
  d <- degree(g); el <- as_edgelist(g, names = FALSE); va <- vertex_attr(g); lc <- transitivity(g, "local")
  terc <- cut(d, unique(quantile(d, c(0, 1/3, 2/3, 1))), include.lowest = TRUE, labels = FALSE)
  data.table(mean_degree = mean(d), sd_degree = sd(d), max_degree = max(d),
             transitivity = transitivity(g, "global"), local_clustering = mean(lc, na.rm = TRUE),
             local_low_deg = mean(lc[terc == 1], na.rm = TRUE), local_mid_deg = mean(lc[terc == 2], na.rm = TRUE),
             local_high_deg = mean(lc[terc == max(terc)], na.rm = TRUE),
             assort_degree = assortativity_degree(g, directed = FALSE), mean_distance = mean_distance(g),
             giant = max(components(g)$csize) / vcount(g),
             same_race = mean(va$race[el[, 1]] == va$race[el[, 2]]), same_sex = mean(va$sex[el[, 1]] == va$sex[el[, 2]]),
             same_relig = mean(va$relig[el[, 1]] == va$relig[el[, 2]]),
             age_diff = mean(abs(va$age[el[, 1]] - va$age[el[, 2]])), educ_diff = mean(abs(va$educ_num[el[, 1]] - va$educ_num[el[, 2]])))
}
ck_of <- function(g) data.table(degree = degree(g), C = transitivity(g, "local"))

for (S in SURVEYS) {
  files <- list.files(NET_DIR, pattern = paste0("^", S, "_x[0-9.]+[a-z]*\\.rds$"))
  LV <- sub(paste0("^", S, "_(.*)\\.rds$"), "\\1", files); LV <- LV[order(LV != "x1", LV != "x2", LV)]
  stamp(S, ": levels ", paste(LV, collapse = ", "))
  NETS <- lapply(setNames(LV, LV), function(l) readRDS(file.path(NET_DIR, paste0(S, "_", l, ".rds"))))
  tg <- fread(file.path("output/02_calibrated_networks", "gss2004_targets.csv"))

  # 1. targets (from x1; identical in every level: same degrees and race mixing)
  by_net <- rbindlist(lapply(seq_along(NETS$x1), function(k) {
    g <- NETS$x1[[k]]; d <- degree(g); r <- factor(V(g)$race, RACES); m <- net_mix(g)
    data.table(net = k, ego = RACES, rel_deg = as.numeric(tapply(d, r, mean)) / mean(d),
               same = m[as.character(ego) == as.character(alter)][match(RACES, as.character(ego)), share])
  }))
  for (l in setdiff(LV, "x1")) stopifnot(all(vapply(seq_along(NETS$x1), function(k)
    identical(degree(NETS$x1[[k]]), degree(NETS[[l]][[k]])), logical(1))))
  fit <- by_net[, .(rel_deg = mean(rel_deg, na.rm = TRUE), rel_min = min(rel_deg, na.rm = TRUE), rel_max = max(rel_deg, na.rm = TRUE),
                    same = mean(same, na.rm = TRUE), same_min = min(same, na.rm = TRUE), same_max = max(same, na.rm = TRUE)), by = ego]
  fit <- merge(tg[, .(ego, n_gss = n, target_rel = ifelse(ego == "Other", NA, rel_deg), rel_lo, rel_hi, target_same = same, same_lo, same_hi)],
               fit, by = "ego")[match(RACES, ego)]
  fit[, `:=`(rel_inside = rel_deg >= rel_lo & rel_deg <= rel_hi, same_inside = same >= same_lo & same <= same_hi)]
  fwrite(fit, file.path(OUT, paste0(S, "_targets.csv")))
  cat("\n", S, ": TARGETS (GSS 2004, 95% CI) vs networks (mean over ", length(NETS$x1), "; 'Other' has no degree target)\n", sep = "")
  print(fit[, .(ego, n_gss, target_rel = round(target_rel, 3), net_rel = round(rel_deg, 3), target_same = round(target_same, 3),
                net_same = round(same, 3), inside_CI = paste(rel_inside, same_inside))], row.names = FALSE)

  # 2. mixing
  mx <- rbindlist(lapply(NETS$x1, net_mix))[, .(net = mean(share)), by = .(ego, alter)]
  mx <- merge(mix_gss, mx[, .(ego = as.character(ego), alter = as.character(alter), net)], by = c("ego", "alter"), all.x = TRUE)
  mx[is.na(net), net := 0]
  fwrite(mx, file.path(OUT, paste0(S, "_mixing.csv")))

  # 3. structure by level
  st <- rbindlist(lapply(LV, function(l) rbindlist(mclapply(NETS[[l]], net_struct, mc.cores = 8))[, level := l]))
  ss <- melt(st, id.vars = "level")[, .(mean = mean(value), sd = sd(value), min = min(value), max = max(value)), by = .(level, variable)]
  fwrite(ss, file.path(OUT, paste0(S, "_structure.csv")))
  cat("\n", S, ": STRUCTURE by level (mean [min, max] over networks)\n", sep = "")
  print(dcast(ss[, .(level, variable, v = sprintf("%.3f [%.3f, %.3f]", mean, min, max))], variable ~ level, value.var = "v")[
    match(unique(ss$variable), variable)], row.names = FALSE)

  # figures: targets and mixing
  fl <- melt(fit[ego %in% CAL], id.vars = "ego", measure.vars = patterns(target = "^target_", net = "^(rel_deg|same)$",
                                                                        lo = "_lo$", hi = "_hi$", nmin = "_min$", nmax = "_max$"),
             variable.name = "stat")[, stat := c("Relative degree", "Same-race share")[stat]]
  p1 <- ggplot(fl, aes(ego)) +
    geom_pointrange(aes(y = target, ymin = lo, ymax = hi, color = "GSS 2004 (95% CI)"), position = position_nudge(x = -0.12), size = 0.3) +
    geom_pointrange(aes(y = net, ymin = nmin, ymax = nmax, color = "Networks (mean, range)"), position = position_nudge(x = 0.12), size = 0.3) +
    facet_wrap(~ stat, scales = "free_y") + scale_color_manual(values = c("black", "#E74C3C"), name = NULL) +
    labs(title = sprintf("%s networks: calibration targets by race of ego", S), x = NULL, y = NULL) + theme(legend.position = "bottom")
  p2 <- ggplot(mx[ego %in% CAL], aes(alter)) +
    geom_pointrange(aes(y = gss, ymin = gss_lo, ymax = gss_hi, color = "GSS 2004 confidants"), position = position_nudge(x = -0.12), size = 0.25) +
    geom_point(aes(y = net, color = "Networks"), position = position_nudge(x = 0.12), size = 1.4) +
    facet_wrap(~ factor(ego, CAL, paste("Ego:", CAL)), nrow = 1) + scale_color_manual(values = c("black", "#E74C3C"), name = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Race of the contacts, by race of ego (only the same-race share is a target)", x = "Race of contact", y = "Share of ties") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(PLOTS, paste0(S, "_targets_mixing.pdf")), p1 / p2, width = 8, height = 6.5, device = cairo_pdf)

  # figures: degree by race and C(k) by level
  dg <- rbindlist(lapply(NETS$x1, function(g) data.table(degree = degree(g), race = factor(V(g)$race, RACES))))
  p3 <- ggplot(dg, aes(degree, after_stat(density), fill = race)) +
    geom_histogram(binwidth = 2, position = "identity", alpha = 0.55) + facet_wrap(~ race, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = COL_RACE, guide = "none") +
    labs(title = sprintf("Degree distribution by race (%d networks; the same in every level)", length(NETS$x1)), x = "Degree", y = "Density")
  ck <- rbindlist(lapply(LV, function(l) rbindlist(lapply(NETS[[l]], ck_of))[, level := l]))
  ck[, bin := cut(degree, c(0, seq(10, 60, by = 5), Inf), right = FALSE)]
  ckb <- ck[!is.na(C), .(k = mean(degree), C = mean(C), lo = quantile(C, .25), hi = quantile(C, .75), n = .N), by = .(level, bin)][n >= 30]
  p4 <- ggplot(ckb, aes(k, C, color = level, fill = level)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) + geom_line(linewidth = 0.6) + geom_point(size = 1) +
    scale_color_manual(values = COL_LEVEL, name = NULL) + scale_fill_manual(values = COL_LEVEL, name = NULL) +
    labs(title = "Local clustering by degree, C(k) (mean and interquartile range of nodes, binned)", x = "Degree k", y = "C(k)") +
    theme(legend.position = "bottom")
  ggsave(file.path(PLOTS, paste0(S, "_degree_Ck.pdf")), p3 / p4, width = 8, height = 6, device = cairo_pdf)

  # figure: example
  set.seed(7)
  g1 <- NETS$x1[[1]]; d1 <- degree(g1)
  v0 <- which(V(g1)$race == "White" & d1 == round(median(d1)))[1]
  ego_plot <- function(g, l) {
    e <- make_ego_graph(g, 1, v0)[[1]]; V(e)$is_ego <- V(e)$name == as.character(v0)
    tri <- sum(count_triangles(e, which(V(e)$is_ego))); kk <- degree(g)[v0]
    ggraph(e, layout = "stress") + geom_edge_link(alpha = 0.35, edge_width = 0.3, color = "grey40") +
      geom_node_point(aes(color = race, size = is_ego, shape = is_ego)) +
      scale_color_manual(values = COL_RACE, name = NULL, drop = FALSE) + scale_size_manual(values = c(1.8, 3.5), guide = "none") +
      scale_shape_manual(values = c(16, 17), guide = "none") +
      labs(title = sprintf("%s: ego network of node %d (degree %d)\n%d ties among its contacts; local clustering %.3f",
                           l, v0, kk, tri, transitivity(g, "local", vids = v0))) +
      theme_void(base_size = 8) + theme(plot.title = element_text(size = 7.5, face = "bold"), legend.position = "bottom")
  }
  V(g1)$name <- as.character(seq_len(vcount(g1)))
  eps <- lapply(LV, function(l) { g <- NETS[[l]][[1]]; V(g)$name <- as.character(seq_len(vcount(g))); ego_plot(g, l) })
  gm <- NETS[[if ("x2" %in% LV) "x2" else LV[1]]][[1]]
  lay <- layout_with_stress(gm)
  pw <- ggraph(gm, layout = "manual", x = lay[, 1], y = lay[, 2]) + geom_edge_link(alpha = 0.02, edge_width = 0.15, color = "grey30") +
    geom_node_point(aes(color = race, size = degree(gm)), alpha = 0.8) +
    scale_color_manual(values = COL_RACE, name = NULL) + scale_size_continuous(range = c(0.3, 2.2), guide = "none") +
    labs(title = sprintf("%s network 1 (%s): %d nodes, %d ties, by race (size = degree)", S, if ("x2" %in% LV) "x2" else LV[1], vcount(gm), ecount(gm))) +
    theme_void(base_size = 8) + theme(plot.title = element_text(size = 8, face = "bold"), legend.position = "bottom")
  ggsave(file.path(PLOTS, paste0(S, "_example.pdf")), wrap_plots(eps, nrow = 1) / pw + plot_layout(heights = c(1, 1.6)),
         width = 3.3 * length(LV), height = 8, device = cairo_pdf)
  stamp(S, ": saved ", file.path(PLOTS, paste0(S, "_*.pdf")), " and ", file.path(OUT, paste0(S, "_*.csv")))
}
