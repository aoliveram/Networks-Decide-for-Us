# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HELD-OUT VALIDATION of the calibrated networks against GSS 2004 ego networks
# 09_network_validation.R
#
# Do the networks of 02_network_calibrated.R reproduce features of real ego
# networks that were not used to build them?
#
# 1. Split. The respondents of the GSS 2004 network module (numgiven and race
#    known) are split in two halves, stratified by race (seed 20041105).
# 2. Networks. Half A gives the generator's targets (relative degree and share of
#    same-race confidants, by race); the nodes are resampled from half B. 24
#    networks per level (02_network_calibrated.R --tag=holdout --target-ids
#    --node-ids): x1, x2 (the main level: twice the triangles of every node) and
#    three variants of the rewiring: x2g (twice the TOTAL triangles, the global
#    target used until 2026-11-05), x3 and x3g (three times).
# 3. Features of half B's confidants (up to 5 per ego) vs the networks, by race of
#    ego (White, Black, Hispanic, Asian) and for all egos:
#      rel_deg     mean numgiven of the group / overall mean   (network: degree)
#      same_race   share of confidants of ego's race: the generator's target,
#                  here out of sample (targets from A, observed in B)
#      same_sex, same_relig, age_diff, educ_diff (|difference| in the 8 GSS
#                  education categories): not targeted. Their coefficients come
#                  from Smith et al. (2014), estimated on the GSS 1985 and 2004
#                  confidants, so half B is not independent of them; what is tested
#                  is that those log-odds, with the race calibration and 30
#                  contacts, give the confidant-level shares and distances
#      blau_race, blau_relig, blau_sex  Blau heterogeneity 1 - sum_g p_g^2 of an
#                  ego's confidants (egos with >= 2), averaged over egos
#      closure_know   share of pairs of confidants who know each other
#                  ("especially close" or "know each other")
#      closure_close  share of pairs who are "especially close"
#    and closure by number of confidants (numgiven 2, 3, 4, 5, 6+).
# 4. Confidants vs contacts. A GSS ego names 2-5 confidants, the strongest of the
#    ~30 contacts of a network node. Shares and distances are compared with ALL
#    the contacts: the generator assumes that confidants have the homophily of
#    the whole personal network (that is how the targets enter). Blau diversity
#    and closure depend on which k contacts are taken. Every half-B ego is
#    matched to a network node of the same race and the same degree rank
#    (numgiven as a proxy of network size, as in rel_deg), and for egos with
#    k = min(numgiven, 5) >= 2, k of the node's contacts are taken by three rules:
#      random    k at random: the network's average contact
#      nearest   the k closest in Gower distance (the model's d_ij, as the gate)
#      embedded  the k with the most common neighbours with the ego (the overlap
#                of a tie grows with its strength)
#    Shares and distances are also reported for the three rules.
# 5. Intervals. Stratified bootstrap of the half-B respondents (within race),
#    jointly with a bootstrap of the 24 networks (B = 1000; one set of weights for
#    every statistic). Relative error = predicted / observed - 1, with its 95%
#    percentile interval.
#
# Outputs
#   output/09_network_validation/split.csv, ids_A.csv, ids_B.csv
#   output/09_network_validation/validation_table.csv          every feature, group, level, rule
#   output/09_network_validation/validation_closure_by_size.csv
#   output/09_network_validation/validation_structure.csv     clustering, assortativity per level
#   output/09_network_validation/validation_table.tex         compact table for the paper
#   plots/09_network_validation/network_validation.pdf
# Usage (repo root): Rscript scripts/09_network_validation.R      (~3 min the first
#   time, when it builds the networks; --regen rebuilds them)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({
  library(haven); library(dplyr); library(data.table); library(igraph); library(cluster)
  library(ggplot2); library(patchwork); library(parallel)
}))
ARGS <- commandArgs(trailingOnly = TRUE)
B_BOOT <- 1000L; SPLIT_SEED <- 20041105L; BOOT_SEED <- 20041106L; TAG <- "holdout"; N_CORES <- 8L
RACES <- c("White", "Black", "Hispanic", "Asian", "Other"); GROUPS <- c(RACES[1:4], "All")
RELIG <- c("Protestant", "Catholic", "Jewish", "None", "OtherRelig")
OUT <- "output/09_network_validation"; PLOTS <- "plots/09_network_validation"
NET_DIR <- file.path("data/02_calibrated_networks", TAG)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
LEVELS <- data.table(level = c("x1", "x2", "x2g", "x3", "x3g"), closure = c(NA, 2, 2, 3, 3),
                     mode = c(NA, "node", "global", "node", "global"))
RULES <- c("random", "nearest", "embedded")
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
num <- function(x) as.numeric(x)

# ---- 1. GSS 2004: egos, confidants and pairs of confidants ----------------------------------
g4 <- read_dta("data/02_GSS_network_ergm/GSS_2004_NORC.dta") |> rename_with(tolower)
race_gss <- function(d) with(d, case_when(hispanic %in% c(2:6, 10, 15, 20:24, 30, 40, 41) ~ "Hispanic",   # as in 02
                                          racecen1 == 1 ~ "White", racecen1 == 2 ~ "Black",
                                          racecen1 %in% c(3, 11, 15) ~ "Other", racecen1 %in% c(4:10) ~ "Asian",
                                          racecen1 == 16 ~ "Hispanic", TRUE ~ NA_character_))
relig5 <- function(x) RELIG[pmin(x, 5)]                    # ego codes 5-13 -> OtherRelig (as in 02)
# ego education in the 8 categories used for the confidants (educ1-5): years + highest degree
# 0 = 1-6 years, 1 = 7-9, 2 = 10-12 without diploma, 3 = high school, 4 = some college,
# 5 = associate, 6 = bachelor, 7 = graduate (degree: 0 < high school ... 4 graduate)
educ_cat_of <- function(y, d) {
  by_years <- ifelse(y <= 6, 0, ifelse(y <= 9, 1, ifelse(y <= 11, 2, ifelse(y == 12, 3, ifelse(y <= 15, 4, ifelse(y == 16, 6, 7))))))
  out <- by_years                                          # degree unknown: years only
  out[d %in% 4] <- 7; out[d %in% 3] <- 6; out[d %in% 2] <- 5
  i1 <- d %in% 1; out[i1] <- ifelse(!is.na(y[i1]) & y[i1] >= 13, 4, 3)
  i0 <- d %in% 0; out[i0] <- ifelse(y[i0] >= 10, 2, ifelse(y[i0] >= 7, 1, 0))
  out
}
yrs <- num(g4$educ); yrs[yrs %in% c(98, 99)] <- NA
EG <- data.table(gss_id = num(g4$id), race = race_gss(g4), numgiven = num(g4$numgiven),
                 sex = num(g4$sex), relig = relig5(num(g4$relig)),
                 age = ifelse(num(g4$age) %in% c(98, 99), NA_real_, num(g4$age)),
                 educ_cat = educ_cat_of(yrs, num(g4$degree)))
EDUC_CAT <- setNames(EG$educ_cat, EG$gss_id)                 # every respondent, for the network nodes
EG <- EG[!is.na(numgiven) & !is.na(race)]
AL <- rbindlist(lapply(1:5, function(k) data.table(
  gss_id = num(g4$id), k = k, race = c("Asian", "Black", "Hispanic", "White", "Other")[num(g4[[paste0("race", k)]])],
  sex = num(g4[[paste0("sex", k)]]), relig = RELIG[num(g4[[paste0("relig", k)]])], age = num(g4[[paste0("age", k)]]),
  educ_cat = num(g4[[paste0("educ", k)]]))))
AL <- AL[gss_id %in% EG$gss_id & !(is.na(race) & is.na(sex) & is.na(relig) & is.na(age) & is.na(educ_cat))]
PR <- rbindlist(lapply(combn(5, 2, simplify = FALSE), function(p) data.table(
  gss_id = num(g4$id), close = num(g4[[sprintf("close%d%d", p[1], p[2])]]))))[gss_id %in% EG$gss_id & close %in% 1:3]

# split, stratified by race
set.seed(SPLIT_SEED)
EG[, half := { n <- .N; h <- rep(c("A", "B"), each = n %/% 2); if (n %% 2) h <- c(h, sample(c("A", "B"), 1)); sample(h) },
   by = race]
fwrite(EG[, .(gss_id, race, numgiven, half)], file.path(OUT, "split.csv"))
IDS_A <- file.path(OUT, "ids_A.csv"); IDS_B <- file.path(OUT, "ids_B.csv")
fwrite(EG[half == "A", .(id = gss_id)], IDS_A); fwrite(EG[half == "B", .(id = gss_id)], IDS_B)
stamp("GSS 2004 network module: ", nrow(EG), " respondents; half A ", sum(EG$half == "A"), ", half B ", sum(EG$half == "B"),
      " (", paste(EG[, .(n = .N), by = .(race, half)][order(race, half), sprintf("%s %s %d", race, half, n)], collapse = ", "), ")")

# ---- 2. networks: targets from A, nodes from B -------------------------------------------------
rscript <- file.path(R.home("bin"), "Rscript")
for (i in which(!is.na(LEVELS$closure))) {
  f <- file.path(NET_DIR, paste0("GSS_", LEVELS$level[i], ".rds"))
  if (file.exists(f) && !"--regen" %in% ARGS) next
  stamp("building the ", LEVELS$level[i], " networks (02_network_calibrated.R, tag ", TAG, ")")
  st <- system2(rscript, c("scripts/02_network_calibrated.R", "--survey=GSS", paste0("--tag=", TAG),
                           paste0("--target-ids=", IDS_A), paste0("--node-ids=", IDS_B),
                           paste0("--closure=", LEVELS$closure[i]), paste0("--closure-mode=", LEVELS$mode[i])))
  if (st != 0) stop("02_network_calibrated.R failed for ", LEVELS$level[i])
}
NETS <- lapply(setNames(LEVELS$level, LEVELS$level), function(l) readRDS(file.path(NET_DIR, paste0("GSS_", l, ".rds"))))
N_NETS <- length(NETS$x1)
# the networks must come from this split: nodes from half B, targets equal to half A's
stopifnot(all(unlist(lapply(NETS, function(L) lapply(L, function(g) V(g)$gss_id))) %in% EG[half == "B", gss_id]))
tgA <- fread(file.path("output/02_calibrated_networks", TAG, "gss2004_targets.csv"))
chkA <- merge(EG[half == "A", .(n = .N, deg = mean(numgiven)), by = .(ego = race)],
              merge(AL[!is.na(race)], EG[half == "A", .(gss_id, ego = race)], by = "gss_id")[, .(same = mean(race == ego)), by = ego],
              by = "ego")
chkA <- merge(chkA, tgA[, .(ego, n_t = n, same_t = same)], by = "ego")
stopifnot(all(chkA$n == chkA$n_t), isTRUE(all.equal(chkA$same, chkA$same_t)))
stamp("networks: ", N_NETS, " per level (", paste(LEVELS$level, collapse = ", "), "); targets of half A verified")

# ---- 3. observed: one row of sufficient statistics per half-B ego ---------------------------------
blau <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else 1 - sum(prop.table(table(x))^2) }
EB <- EG[half == "B"]
EB[, k := pmin(numgiven, 5)]
a <- merge(AL, EB[, .(gss_id, e_race = race, e_sex = sex, e_relig = relig, e_age = age, e_educ = educ_cat)], by = "gss_id")
ob <- a[, .(a_race = sum(!is.na(race)), s_race = sum(race == e_race, na.rm = TRUE),
            a_sex = sum(!is.na(sex) & !is.na(e_sex)), s_sex = sum(sex == e_sex, na.rm = TRUE),
            a_relig = sum(!is.na(relig) & !is.na(e_relig)), s_relig = sum(relig == e_relig, na.rm = TRUE),
            a_age = sum(!is.na(age) & !is.na(e_age)), d_age = sum(abs(age - e_age), na.rm = TRUE),
            a_educ = sum(!is.na(educ_cat) & !is.na(e_educ)), d_educ = sum(abs(educ_cat - e_educ), na.rm = TRUE),
            b_race = blau(race), b_relig = blau(relig), b_sex = blau(sex)), by = gss_id]
pr <- PR[gss_id %in% EB$gss_id, .(p_all = .N, p_know = sum(close %in% 1:2), p_close = sum(close == 1)), by = gss_id]
EB <- merge(merge(EB, ob, by = "gss_id", all.x = TRUE), pr, by = "gss_id", all.x = TRUE)
for (v in c("a_race", "s_race", "a_sex", "s_sex", "a_relig", "s_relig", "a_age", "d_age", "a_educ", "d_educ", "p_all", "p_know", "p_close"))
  set(EB, which(is.na(EB[[v]])), v, 0)
for (v in c("race", "relig", "sex")) {
  EB[, paste0("i_", v) := as.numeric(!is.na(get(paste0("b_", v))))]
  set(EB, which(is.na(EB[[paste0("b_", v)]])), paste0("b_", v), 0)
}
EB[, size := factor(pmin(numgiven, 6), levels = 0:6, labels = c(0:5, "6+"))]
setorder(EB, gss_id)
N_E <- nrow(EB)

# features as ratios of weighted sums (numerator column / denominator column)
FEATS <- list(same_race = c("s_race", "a_race"), same_sex = c("s_sex", "a_sex"), same_relig = c("s_relig", "a_relig"),
              age_diff = c("d_age", "a_age"), educ_diff = c("d_educ", "a_educ"),
              blau_race = c("b_race", "i_race"), blau_relig = c("b_relig", "i_relig"), blau_sex = c("b_sex", "i_sex"),
              closure_know = c("p_know", "p_all"), closure_close = c("p_close", "p_all"))
MATCHED <- c("blau_race", "blau_relig", "blau_sex", "closure_know", "closure_close")
gmask <- function(gr) if (gr == "All") rep(1, N_E) else as.numeric(EB$race == gr)
SIZES <- c("2", "3", "4", "5", "6+")

# bootstrap weights: respondents within race (W, B x egos) and networks (V, B x nets); row 1 = point estimate
set.seed(BOOT_SEED)
W <- matrix(0, B_BOOT, N_E)
for (r in unique(EB$race)) { ii <- which(EB$race == r)
  W[, ii] <- t(vapply(seq_len(B_BOOT), function(b) tabulate(sample.int(length(ii), length(ii), TRUE), length(ii)), numeric(length(ii)))) }
V <- t(vapply(seq_len(B_BOOT), function(b) tabulate(sample.int(N_NETS, N_NETS, TRUE), N_NETS), numeric(N_NETS)))
W <- rbind(1, W); V <- rbind(1, V)

obs_feature <- function(f, gr, mask = gmask(gr)) {
  cols <- FEATS[[f]]
  drop(W %*% (EB[[cols[1]]] * mask)) / drop(W %*% (EB[[cols[2]]] * mask))
}
obs_reldeg <- function(gr) (drop(W %*% (EB$numgiven * gmask(gr))) / drop(W %*% gmask(gr))) / (drop(W %*% EB$numgiven) / rowSums(W))
OBS <- list()
for (gr in GROUPS) {
  if (gr != "All") OBS[[paste("rel_deg", gr)]] <- obs_reldeg(gr)
  for (f in names(FEATS)) OBS[[paste(f, gr)]] <- obs_feature(f, gr)
}
for (s in SIZES) for (f in c("closure_know", "closure_close"))
  OBS[[paste(f, s)]] <- obs_feature(f, "All", as.numeric(EB$size == s))

# ---- 4. predicted -----------------------------------------------------------------------------------
gower_of <- function(g)                                     # the model's d_ij (05_unified_engine.R)
  as.matrix(daisy(data.frame(age = V(g)$age, educ_num = V(g)$educ_num, race = as.factor(V(g)$race),
                             relig = as.factor(V(g)$relig), sex = as.factor(V(g)$sex)), metric = "gower"))
# every half-B ego -> a node of the same race and degree rank (ties at random)
match_nodes <- function(race, deg) {
  node <- integer(N_E)
  for (r in unique(EB$race)) {
    ei <- which(EB$race == r); ni <- which(race == r)
    if (!length(ni)) stop("no ", r, " nodes in a network")
    ni <- ni[order(deg[ni] + runif(length(ni)))]; eo <- ei[order(EB$numgiven[ei] + runif(length(ei)))]
    node[eo] <- ni[pmax(1L, ceiling((seq_along(eo) - 0.5) / length(eo) * length(ni)))]
  }
  node
}
# statistics of one network: all-contacts features by group, structure, and per-ego matched statistics
net_stats <- function(g, seed) {
  set.seed(seed)
  va <- vertex_attr(g); n <- vcount(g); deg <- degree(g); va$educ_cat <- unname(EDUC_CAT[as.character(va$gss_id)])
  el <- as_edgelist(g, names = FALSE); i <- c(el[, 1], el[, 2]); j <- c(el[, 2], el[, 1])
  te <- data.table(r = va$race[i], same_race = va$race[i] == va$race[j], same_sex = va$sex[i] == va$sex[j],
                   same_relig = va$relig[i] == va$relig[j], age_diff = abs(va$age[i] - va$age[j]),
                   educ_diff = abs(va$educ_cat[i] - va$educ_cat[j]))
  all_c <- rbind(te[, lapply(.SD, mean), by = r, .SDcols = -"r"], te[, lapply(.SD, mean), .SDcols = -"r"][, r := "All"])
  all_c <- merge(all_c, data.table(r = c(RACES, "All"), rel_deg = c(tapply(deg, factor(va$race, RACES), mean) / mean(deg), 1)), by = "r")
  lc <- transitivity(g, "local"); q5 <- cut(deg, unique(quantile(deg, 0:5 / 5)), include.lowest = TRUE, labels = FALSE)
  str <- as.data.table(c(list(transitivity = transitivity(g, "global"), local_clustering = mean(lc, na.rm = TRUE),
                              assort_degree = assortativity_degree(g, directed = FALSE)),
                         as.list(setNames(tapply(lc, q5, mean, na.rm = TRUE), paste0("local_q", sort(unique(q5)))))))
  # matched egos
  A <- as_adjacency_matrix(g, sparse = FALSE) > 0; D <- gower_of(g); nb <- adjacent_vertices(g, V(g))
  node <- match_nodes(va$race, deg)
  res <- lapply(RULES, function(rule) {
    M <- matrix(0, N_E, 16, dimnames = list(NULL, c("a_race", "s_race", "a_sex", "s_sex", "a_relig", "s_relig", "a_age", "d_age",
                                                    "a_educ", "d_educ", "b_race", "b_relig", "b_sex", "i_any", "p_all", "p_know")))
    for (e in which(EB$k >= 2)) {
      v <- node[e]; nv <- as.integer(nb[[v]]); k <- min(EB$k[e], length(nv)); if (k < 2) next
      S <- switch(rule,
        random   = nv[sample.int(length(nv), k)],
        nearest  = nv[order(D[v, nv], runif(length(nv)))[seq_len(k)]],
        embedded = nv[order(-rowSums(A[nv, nv, drop = FALSE]), runif(length(nv)))[seq_len(k)]])
      M[e, ] <- c(k, sum(va$race[S] == va$race[v]), k, sum(va$sex[S] == va$sex[v]), k, sum(va$relig[S] == va$relig[v]),
                  k, sum(abs(va$age[S] - va$age[v])), k, sum(abs(va$educ_cat[S] - va$educ_cat[v])),
                  blau(va$race[S]), blau(va$relig[S]), blau(va$sex[S]), 1, k * (k - 1) / 2, sum(A[S, S]) / 2)
    }
    M
  })
  names(res) <- RULES
  list(all = all_c, structure = str, matched = res)
}
ST <- list()
for (l in LEVELS$level) {
  stamp("statistics of ", l)
  ST[[l]] <- mclapply(seq_len(N_NETS), function(k) net_stats(NETS[[l]][[k]], 9000L + k), mc.cores = N_CORES)
}

# predicted, with the network weights V (and, for matched statistics, the respondent weights W)
pred_all <- function(l, f, gr) {
  x <- vapply(ST[[l]], function(s) s$all[r == gr][[f]], numeric(1))
  drop(V %*% x) / rowSums(V)
}
# matched: sum over networks and egos, num / den; observed availability (i_*, p_all > 0) is imposed
pred_matched <- function(l, rule, f, mask) {
  num_col <- c(same_race = "s_race", same_sex = "s_sex", same_relig = "s_relig", age_diff = "d_age", educ_diff = "d_educ",
               blau_race = "b_race", blau_relig = "b_relig", blau_sex = "b_sex", closure_know = "p_know")[[f]]
  den_col <- c(same_race = "a_race", same_sex = "a_sex", same_relig = "a_relig", age_diff = "a_age", educ_diff = "a_educ",
               blau_race = "i_any", blau_relig = "i_any", blau_sex = "i_any", closure_know = "p_all")[[f]]
  avail <- switch(f, blau_race = EB$i_race, blau_relig = EB$i_relig, blau_sex = EB$i_sex,
                  closure_know = as.numeric(EB$p_all > 0), rep(1, N_E))
  m <- mask * avail
  Xn <- vapply(ST[[l]], function(s) s$matched[[rule]][, num_col] * m, numeric(N_E))
  Xd <- vapply(ST[[l]], function(s) s$matched[[rule]][, den_col] * m, numeric(N_E))
  rowSums((W %*% Xn) * V) / rowSums((W %*% Xd) * V)
}

# ---- 5. table -----------------------------------------------------------------------------------------
q3 <- function(x) c(x[1], quantile(x[-1], c(.025, .975), na.rm = TRUE, names = FALSE))
row_of <- function(feature, group, level, source, pred, obs) {
  o <- q3(obs); p <- q3(pred); re <- q3(pred / obs - 1)
  data.table(feature = feature, group = group, level = level, source = source,
             observed = o[1], obs_lo = o[2], obs_hi = o[3], predicted = p[1], pred_lo = p[2], pred_hi = p[3],
             rel_err = re[1], rel_lo = re[2], rel_hi = re[3], fits = re[2] <= 0 & re[3] >= 0)
}
TAB <- list(); CBS <- list()
for (l in LEVELS$level) for (gr in GROUPS) {
  if (gr != "All") TAB[[length(TAB) + 1]] <- row_of("rel_deg", gr, l, "all contacts", pred_all(l, "rel_deg", gr), OBS[[paste("rel_deg", gr)]])
  for (f in c("same_race", "same_sex", "same_relig", "age_diff", "educ_diff"))
    TAB[[length(TAB) + 1]] <- row_of(f, gr, l, "all contacts", pred_all(l, f, gr), OBS[[paste(f, gr)]])
  # shares and distances under a rule: observed restricted to the same egos (k >= 2)
  for (rule in RULES) for (f in c("same_race", "same_sex", "same_relig", "age_diff", "educ_diff"))
    TAB[[length(TAB) + 1]] <- row_of(f, gr, l, rule, pred_matched(l, rule, f, gmask(gr)), obs_feature(f, gr, gmask(gr) * (EB$k >= 2)))
  for (rule in RULES) for (f in MATCHED[MATCHED != "closure_close"])
    TAB[[length(TAB) + 1]] <- row_of(f, gr, l, rule, pred_matched(l, rule, f, gmask(gr)), OBS[[paste(f, gr)]])
  # "especially close" is compared with the same connected pairs
  for (rule in RULES)
    TAB[[length(TAB) + 1]] <- row_of("closure_close", gr, l, rule, pred_matched(l, rule, "closure_know", gmask(gr)), OBS[[paste("closure_close", gr)]])
}
for (l in LEVELS$level) for (rule in RULES) for (s in SIZES) for (f in c("closure_know", "closure_close"))
  CBS[[length(CBS) + 1]] <- row_of(f, s, l, rule, pred_matched(l, rule, "closure_know", as.numeric(EB$size == s)), OBS[[paste(f, s)]])
TAB <- rbindlist(TAB); CBS <- rbindlist(CBS)
# profile of closure with ego-network size: closure among egos with 5 or more confidants minus
# closure among egos with 2 or 3 (GSS: flat to falling); the network side uses the rank-matched nodes
big <- as.numeric(EB$size %in% c("5", "6+")); small <- as.numeric(EB$size %in% c("2", "3"))
TREND <- rbindlist(lapply(LEVELS$level, function(l) rbindlist(lapply(RULES, function(rule) rbindlist(lapply(c("closure_know", "closure_close"), function(f) {
  o <- obs_feature(f, "All", big) - obs_feature(f, "All", small)
  pp <- pred_matched(l, rule, "closure_know", big) - pred_matched(l, rule, "closure_know", small)
  oq <- q3(o); pq <- q3(pp); dq <- q3(pp - o)
  data.table(feature = f, level = l, source = rule, observed = oq[1], obs_lo = oq[2], obs_hi = oq[3],
             predicted = pq[1], pred_lo = pq[2], pred_hi = pq[3], diff = dq[1], diff_lo = dq[2], diff_hi = dq[3])
}))))))
TAB[, n_egos := vapply(group, function(gr) if (gr == "All") N_E else sum(EB$race == gr), numeric(1))]
fwrite(TAB, file.path(OUT, "validation_table.csv")); fwrite(CBS, file.path(OUT, "validation_closure_by_size.csv"))
fwrite(TREND, file.path(OUT, "validation_closure_trend.csv"))
STR <- rbindlist(lapply(LEVELS$level, function(l) rbindlist(lapply(ST[[l]], `[[`, "structure"), fill = TRUE)[, lapply(.SD, mean)][, level := l]))
fwrite(STR, file.path(OUT, "validation_structure.csv"))

# ---- 6. report ------------------------------------------------------------------------------------------
fmt <- function(x, lo, hi, d = 3) sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), x, lo, hi)
cat("\n==== HELD-OUT VALIDATION (targets: half A; nodes and observed: half B; ", N_E, " egos, ", N_NETS,
    " networks per level) ====\n", sep = "")
show <- TAB[level %in% c("x1", "x2", "x2g") & (source == "all contacts" | feature %in% MATCHED)]
show[, obs := fmt(observed, obs_lo, obs_hi)][, pred := sprintf("%.3f", predicted)][, err := sprintf("%+.0f%% [%+.0f, %+.0f]%s", 100 * rel_err, 100 * rel_lo, 100 * rel_hi, ifelse(fits, " ok", ""))]
for (f in unique(show$feature)) {
  cat("\n--", f, "\n")
  w <- dcast(show[feature == f], group + obs ~ level + source, value.var = "err")
  print(w[match(intersect(GROUPS, w$group), group)], row.names = FALSE)
}
cat("\nClosure by number of confidants (all egos): observed (know / especially close) and predicted share of connected pairs:\n")
cb <- dcast(CBS[feature == "closure_know"], group ~ level + source, value.var = "predicted")
cb <- merge(dcast(unique(CBS[level == "x1" & source == "random", .(group, feature, observed)]), group ~ feature, value.var = "observed"), cb, by = "group")
print(cb[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
cat("\nClosure trend: (numgiven >= 5) minus (numgiven 2-3), pp; observed vs predicted; diff = predicted - observed [95% CI]:\n")
print(TREND[, .(feature, level, source, observed = sprintf("%+.1f [%+.1f, %+.1f]", 100 * observed, 100 * obs_lo, 100 * obs_hi),
                predicted = sprintf("%+.1f", 100 * predicted), diff = sprintf("%+.1f [%+.1f, %+.1f]", 100 * diff, 100 * diff_lo, 100 * diff_hi))],
      row.names = FALSE, nrows = 200)
cat("\nStructure of the validation networks (means over networks):\n"); print(STR, digits = 3)

# compact table for the paper: all egos (and same-race by group), x1 / x2 / x2g
pick <- rbind(TAB[group == "All" & source == "all contacts" & feature %in% c("same_race", "same_sex", "same_relig", "age_diff", "educ_diff")],
              TAB[group %in% RACES[1:4] & source == "all contacts" & feature %in% c("rel_deg", "same_race")],
              TAB[group == "All" & feature %in% MATCHED & source %in% RULES])[level %in% c("x1", "x2", "x2g")]
lab_f <- c(rel_deg = "Relative degree", same_race = "Same race", same_sex = "Same sex", same_relig = "Same religion",
           age_diff = "Age difference (years)", educ_diff = "Education difference (categories)",
           blau_race = "Blau diversity, race", blau_relig = "Blau diversity, religion", blau_sex = "Blau diversity, sex",
           closure_know = "Pairs who know each other", closure_close = "Pairs especially close")
pick[, row := paste0(lab_f[feature], ifelse(group == "All", "", paste0(", ", group)), ifelse(source == "all contacts", "", paste0(" (", source, ")")))]
tex <- dcast(pick, row + feature + group + source + observed + obs_lo + obs_hi ~ level, value.var = c("predicted", "rel_err", "fits"))
tex <- tex[order(match(feature, names(lab_f)), match(group, GROUPS), match(source, c("all contacts", RULES)))]
cell <- function(p, e, ok) sprintf("%.3f & %s%+.0f\\%%%s", p, if (ok) "" else "\\textbf{", 100 * e, if (ok) "" else "}")
lines <- c("% generated by scripts/09_network_validation.R; bold = 95% interval of the relative error excludes 0",
           "\\begin{tabular}{lcrrrrrr}", "\\hline",
           " & GSS 2004, half B & \\multicolumn{2}{c}{x1} & \\multicolumn{2}{c}{x2 (main)} & \\multicolumn{2}{c}{x2g (global)} \\\\",
           "Feature & observed [95\\% CI] & pred. & error & pred. & error & pred. & error \\\\", "\\hline",
           vapply(seq_len(nrow(tex)), function(i) with(tex[i], paste0(row, " & ", sprintf("%.3f [%.3f, %.3f]", observed, obs_lo, obs_hi), " & ",
                                                                     cell(predicted_x1, rel_err_x1, fits_x1), " & ", cell(predicted_x2, rel_err_x2, fits_x2), " & ",
                                                                     cell(predicted_x2g, rel_err_x2g, fits_x2g), " \\\\")), ""),
           "\\hline",
           " & \\multicolumn{7}{l}{Closure, egos with 5+ confidants minus egos with 2--3 (pp): pred. \\& pred.$-$obs.} \\\\", "\\hline",
           vapply(RULES, function(rule) {
             tr <- TREND[feature == "closure_know" & source == rule]
             cl <- function(l) with(tr[level == l], sprintf("%+.1f & %s%+.1f%s", 100 * predicted, if (diff_lo <= 0 & diff_hi >= 0) "" else "\\textbf{",
                                                            100 * diff, if (diff_lo <= 0 & diff_hi >= 0) "" else "}"))
             o <- tr[level == "x1"]
             paste0("Pairs who know each other (", rule, ") & ", sprintf("%+.1f [%+.1f, %+.1f]", 100 * o$observed, 100 * o$obs_lo, 100 * o$obs_hi),
                    " & ", cl("x1"), " & ", cl("x2"), " & ", cl("x2g"), " \\\\")
           }, ""),
           "\\hline", "\\end{tabular}")
writeLines(lines, file.path(OUT, "validation_table.tex"))

# ---- 7. figure ------------------------------------------------------------------------------------------
COLS <- c(x1 = "#2980B9", x2 = "#E74C3C", x2g = "#7F8C8D", x3 = "#E67E22", x3g = "#BDC3C7")
pd <- TAB[group == "All" & level %in% c("x1", "x2", "x2g") & (source == "all contacts" | feature %in% MATCHED)]
pd[, lab := paste0(lab_f[feature], ifelse(source == "all contacts", "", paste0(" — ", source)))]
pd[, lab := factor(lab, rev(unique(lab[order(match(feature, names(lab_f)), match(source, c("all contacts", RULES)))])))]
p1 <- ggplot(pd, aes(100 * rel_err, lab, color = level)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_linerange(aes(xmin = 100 * rel_lo, xmax = 100 * rel_hi), position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.6) +
  scale_color_manual(values = COLS, name = NULL) +
  labs(title = "Relative error, networks vs GSS 2004 half B (all egos; 95% bootstrap CI)",
       x = "Predicted / observed − 1 (%)", y = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 8, face = "bold"))
gd <- TAB[group %in% RACES[1:4] & source == "all contacts" & feature %in% c("rel_deg", "same_race", "same_sex", "same_relig") & level == "x2"]
p2 <- ggplot(gd, aes(group)) +
  geom_pointrange(aes(y = observed, ymin = obs_lo, ymax = obs_hi), color = "black", size = 0.25, position = position_nudge(x = -0.12)) +
  geom_pointrange(aes(y = predicted, ymin = pred_lo, ymax = pred_hi), color = COLS[["x2"]], size = 0.25, position = position_nudge(x = 0.12)) +
  facet_wrap(~ lab_f[feature], scales = "free_y", nrow = 1) +
  labs(title = "By race of ego: GSS 2004 half B (black) vs x2 networks (red)", x = NULL, y = NULL) +
  theme_bw(base_size = 8) + theme(plot.title = element_text(size = 8, face = "bold"))
cd <- CBS[feature == "closure_know" & level %in% c("x1", "x2", "x2g", "x3", "x3g")]
co <- unique(CBS[feature %in% c("closure_know", "closure_close") & level == "x1" & source == "random", .(group, feature, observed, obs_lo, obs_hi)])
p3 <- ggplot(cd, aes(group, predicted, color = level, group = level)) +
  geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
  geom_pointrange(data = co, aes(group, observed, ymin = obs_lo, ymax = obs_hi, shape = feature), inherit.aes = FALSE, size = 0.25) +
  scale_shape_manual(values = c(closure_know = 16, closure_close = 1), labels = c(closure_know = "GSS: know each other", closure_close = "GSS: especially close"), name = NULL) +
  facet_wrap(~ source, nrow = 1) + scale_color_manual(values = COLS, name = NULL) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Closure among k = min(numgiven, 5) contacts, by number of confidants (egos matched by race and degree rank)",
       x = "Number of confidants named (numgiven)", y = "Share of pairs connected") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 8, face = "bold"))
ggsave(file.path(PLOTS, "network_validation.pdf"), p1 / p2 / p3 + plot_layout(heights = c(1.5, 0.8, 1)),
       width = 8.5, height = 11, device = cairo_pdf)
stamp("saved ", file.path(OUT, "validation_table.csv"), ", validation_table.tex and ", file.path(PLOTS, "network_validation.pdf"))
