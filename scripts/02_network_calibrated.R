# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NETWORKS FOR THE SWEEP: calibrated to GSS 2004 (x1), then rewired to twice the
# triangles of every node, with the same degrees and the same mixing (x2, the
# main networks).
#
# Replaces the ERGM networks of 02_GSS_network_ergm.R / 02_ATP_network_ergm.R
# for 05-07 (those stay in data/ as the earlier reference). Three changes:
#   1. Nodes are RESAMPLED for every network (with replacement, own seed) from
#      the survey's complete cases. The 100 ERGM networks shared the same 1000.
#   2. Race enters by group. The tie model is dyad independent:
#        logit p_ij = delta_r(i) + delta_r(j) + gamma_r [race_i = race_j = r] + rest
#        rest = sex, religion, age and education terms of 02 (Smith et al. 2014,
#               2004 effects; 'Other' keeps 02's common race coefficient)
#      delta (activity, 5 groups) and gamma (homophily, White/Black/Hispanic/
#      Asian) are solved by damped Newton so that the EXPECTED statistics of each
#      network equal GSS 2004 targets for its own nodes (moment matching = the
#      MLE of this dyad-independent ERGM; what ergm.ego estimates for
#      nodefactor + nodematch(race, diff = TRUE)):
#        sum of degrees of group r = n_r x 30 x rel_deg_r
#        ties inside group r       = n_r x d_r x same_r / 2
#      rel_deg_r: mean number of confidants named by group-r respondents over
#      the overall mean; same_r: share of their confidants of race r. Mean
#      degree 30 (Boase et al. 2006). Ties are then drawn as independent coins
#      (exact, no MCMC). This is level x1.
#   3. Clustering (level x2). Degree-preserving swaps (a,b),(c,d) -> (a,d),(c,b)
#      that keep the race mixing matrix and the same-sex and same-religion counts
#      EXACT and the total age and education differences within 1%, accepted
#      toward twice the triangles of EVERY node of x1 (02_closure_rewire.cpp,
#      rewire_closure_node(): Metropolis on sum_i |t_i - 2 t_i(x1)|, stopped at a
#      5% relative error), then 5e7 more proposals that keep the triangles near
#      the target while the network mixes. x1 and x2 have the same nodes,
#      degrees and mixing; they differ only in clustering (transitivity ~0.075 ->
#      ~0.15, local clustering doubled at every degree: the clustering-degree
#      profile of x1 is kept). The age and education constraints bind: both
#      totals end ~1% higher.
#      x2 is not a draw from a probability model; it is the main set because
#      real personal networks are clustered (81-84% of GSS 2004 confidant pairs
#      know each other) and a Bernoulli graph has almost no closure.
#      The global variant (--closure-mode=global, level x2g: twice the TOTAL
#      triangles) puts the new triangles among high-degree nodes, so closure rises
#      with degree and degree assortativity jumps; the held-out GSS 2004 ego
#      networks reject that profile (09_network_validation.R), so x2g is kept for
#      the validation only. (x1 and x2g reproduce bit for bit the exploratory
#      networks of playground/PR4-rcpp-engine.)
#
# ATP 2014 has no ego networks: ATP nodes get the SAME GSS 2004 targets, as the
# ATP ERGM of 02 used the GSS-based coefficients of Smith et al.
#
# Outputs (<S> = GSS or ATP):
#   data/02_calibrated_networks/<S>_x1.rds, <S>_x2.rds   lists of 24 igraph objects
#       vertex attributes: age, sex, educ_num, race, relig, mur_score, propensity_score
#   data/02_calibrated_networks/<S>_calibration.rds      parameters, Newton and rewiring info
#   output/02_calibrated_networks/gss2004_targets.csv    targets with bootstrap CI
#   output/02_calibrated_networks/<S>_by_group.csv       degree and same-race share by group vs targets
#   output/02_calibrated_networks/<S>_structure.csv      structure of x1 and x2 (means over networks)
#   output/02_calibrated_networks/<S>_checks.csv         per network: degrees, mixing, simple graph
#
# Seeds: nodes set.seed(80000 + 1000 + i) for GSS and 180000 + 1000 + i for ATP;
# rewiring 70000 + i and 170000 + i (both modes).
#
# Usage (repo root):
#   Rscript scripts/02_network_calibrated.R --survey=GSS     # ~15 s on 8 cores
#   Rscript scripts/02_network_calibrated.R --survey=ATP
#   Rscript scripts/02_network_calibrated.R --survey=GSS --test   # 2 networks, sandbox
# Options: --nets=24 --closure=2 --cores=8
#   --closure-mode=node    every node moves to closure x its own triangles (default;
#                          keeps the clustering-degree profile of x1)
#   --closure-mode=global  the target is the total number of triangles (level x<closure>g)
# Validation runs (scripts/09_network_validation.R; none of these changes the default output):
#   --target-ids=<csv>  GSS 2004 respondents (column id) from which the targets are computed
#   --node-ids=<csv>    GSS node pool restricted to these respondents (column id)
#   --tag=<name>        write to data/02_calibrated_networks/<name>/ and
#                       output/02_calibrated_networks/<name>/, and keep the respondent id
#                       of every GSS node as the vertex attribute gss_id
#   --name=<level>      file name of the rewired level (default x<closure>, plus "g" in global mode)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({
  library(igraph); library(data.table); library(haven); library(dplyr); library(parallel); library(Rcpp)
}))
Sys.setenv(OMP_NUM_THREADS = "1")
sourceCpp("scripts/02_closure_rewire.cpp")

ARGS <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) sub(".*=", "", a[1]) else d }
TEST    <- "--test" %in% ARGS
SURVEY  <- toupper(arg("survey", "GSS")); stopifnot(SURVEY %in% c("GSS", "ATP"))
N_NETS  <- as.integer(arg("nets", if (TEST) 2 else 24))
CLOSURE <- as.numeric(arg("closure", 2))
N_CORES <- as.integer(arg("cores", 8))
MODE    <- arg("closure-mode", "node"); stopifnot(MODE %in% c("node", "global"))
E_TOL   <- 0.05          # node mode: phase 1 ends when sum_i |t_i - t_i*| <= 5% of sum_i t_i*
LEVEL2  <- arg("name", paste0("x", format(CLOSURE), if (MODE == "global") "g" else ""))
TAG     <- arg("tag", "")
TARGET_IDS <- arg("target-ids", ""); NODE_IDS <- arg("node-ids", "")
read_ids <- function(f) { x <- read.csv(f); stopifnot("id" %in% names(x)); as.numeric(x$id) }
N <- 1000L; K_MEAN <- 30
sub_dir <- function(d) if (nzchar(TAG)) file.path(d, TAG) else d
NET_OUT <- sub_dir(if (TEST) "data/02_calibrated_networks/_test" else "data/02_calibrated_networks")
OUT     <- sub_dir(if (TEST) "output/02_calibrated_networks/_test" else "output/02_calibrated_networks")
dir.create(NET_OUT, recursive = TRUE, showWarnings = FALSE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED_NODES <- if (SURVEY == "GSS") 80000L else 180000L
SEED_REWIRE <- if (SURVEY == "GSS") 70000L else 170000L
RACES <- c("White", "Black", "Hispanic", "Asian", "Other"); CAL <- 1:4      # groups with their own targets
# 02's homophily coefficients (Smith et al. 2014, Table 3 model 2, effects for 2004), nodematch scale
COEF <- c(race = 1.959 - 0.264, sex = 0.373 - 0.144, relig = 1.270 + 0.215,
          age = -0.047 - 0.005, educ = -0.157 - 0.044)
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
t_start <- Sys.time()

# ---- 1. node pool: complete cases, recoded as in 02 and 03 ----------------------------------
GSS_DTA <- "data/02_GSS_network_ergm/GSS_2004_NORC.dta"
race_gss <- function(d) with(d, case_when(hispanic %in% c(2:6, 10, 15, 20:24, 30, 40, 41) ~ "Hispanic",
                                          racecen1 == 1 ~ "White", racecen1 == 2 ~ "Black",
                                          racecen1 %in% c(3, 11, 15) ~ "Other", racecen1 %in% c(4:10) ~ "Asian",
                                          racecen1 == 16 ~ "Hispanic", TRUE ~ NA_character_))
g4 <- read_dta(GSS_DTA) |> rename_with(tolower)
if (SURVEY == "GSS") {
  items <- c("signdpet", "avoidbuy", "joindem", "attrally", "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")
  pop <- g4 |>
    mutate(gss_id = as.numeric(id),
           sex = as.character(factor(na_if(as.numeric(sex), 9), levels = 1:2, labels = c("Male", "Female"))),
           race = race_gss(g4),
           educ_num = ifelse(educ %in% c(98, 99), NA_real_, as.numeric(educ)),
           age = ifelse(age %in% c(98, 99), NA_real_, as.numeric(age)),
           relig = case_when(relig == 1 ~ "Protestant", relig == 2 ~ "Catholic", relig == 3 ~ "Jewish",
                             relig == 4 ~ "None", relig %in% 5:13 ~ "OtherRelig", TRUE ~ NA_character_))
  for (v in items) { x <- as.numeric(pop[[v]]); x[x %in% c(8, 9)] <- NA; pop[[v]] <- x }
  pop <- pop[complete.cases(pop[, c("age", "sex", "educ_num", "race", "relig", items)]), ]
  pop$propensity_score <- rowSums(4 - as.matrix(pop[, items])) / 27        # 03_GSS_MUR_calculation.R
} else {
  items <- paste0("METECH_", LETTERS[1:6])
  pop <- readRDS("data/01_ATP_GSS_imputation/ATP_W3_W4_input.rds")
  for (v in items) { x <- as.numeric(zap_labels(pop[[v]])); x[!x %in% c(0, 1)] <- NA; pop[[v]] <- x }
  pop$sex <- as.character(pop$sex); pop$race <- as.character(pop$race); pop$relig <- as.character(pop$relig)
  pop <- pop[complete.cases(pop[, c("age", "sex", "educ_num", "race", "relig", items)]), ]
  pop$propensity_score <- with(pop, METECH_A + (1 - METECH_B) + METECH_C + METECH_D +   # 03_ATP_MUR_calculation.R
                                    (1 - METECH_E) + (1 - METECH_F)) / 6
}
pop$mur_score <- 1 - pop$propensity_score                                    # q = aversion
if (nzchar(NODE_IDS)) {
  stopifnot(SURVEY == "GSS"); pop <- pop[pop$gss_id %in% read_ids(NODE_IDS), ]
}
VCOLS <- c("age", "sex", "educ_num", "race", "relig", "mur_score", "propensity_score",
           if (nzchar(TAG) && SURVEY == "GSS") "gss_id")
pop <- as.data.frame(pop)[, VCOLS]
stopifnot(all(pop$race %in% RACES))
stamp(SURVEY, " node pool: ", nrow(pop), " complete cases", if (nzchar(NODE_IDS)) paste0(" (restricted to ", NODE_IDS, ")"),
      "; mean mur_score = ", round(mean(pop$mur_score), 3))

# ---- 2. GSS 2004 targets by race of ego (with bootstrap CI) ---------------------------------
alt_code <- c("1" = "Asian", "2" = "Black", "3" = "Hispanic", "4" = "White", "5" = "Other")
resp <- data.table(id = seq_len(nrow(g4)), gss_id = as.numeric(g4$id), ego = race_gss(g4), numgiven = as.numeric(g4$numgiven))
resp[numgiven > 6, numgiven := NA]
if (nzchar(TARGET_IDS)) resp <- resp[gss_id %in% read_ids(TARGET_IDS)]
alt <- rbindlist(lapply(1:5, function(k) data.table(id = seq_len(nrow(g4)),
                                                    alter = alt_code[as.character(as.numeric(g4[[paste0("race", k)]]))])))[!is.na(alter)]
resp <- resp[!is.na(ego) & !is.na(numgiven)]
alt <- merge(alt, resp[, .(id, ego)], by = "id")
gss_targets <- function(ids) {
  r <- resp[.(ids), on = "id"]; a <- alt[.(ids), on = "id", nomatch = 0L, allow.cartesian = TRUE]
  d <- r[, .(n = .N, deg = mean(numgiven)), by = ego]; d[, rel_deg := deg / weighted.mean(deg, n)]
  s <- a[, .(alters = .N, same = mean(alter == ego)), by = ego]
  merge(d, s, by = "ego", all = TRUE)[match(RACES, ego)]
}
tg <- gss_targets(resp$id)
set.seed(2004)
boot <- rbindlist(lapply(1:1000, function(b) gss_targets(sample(resp$id, replace = TRUE))[, b := b]))
tg <- merge(tg, boot[, .(rel_lo = quantile(rel_deg, .025, na.rm = TRUE), rel_hi = quantile(rel_deg, .975, na.rm = TRUE),
                         same_lo = quantile(same, .025, na.rm = TRUE), same_hi = quantile(same, .975, na.rm = TRUE)), by = ego],
            by = "ego")[match(RACES, ego)]
REL  <- setNames(ifelse(RACES == "Other", 1, tg$rel_deg), RACES)      # 'Other' (18 respondents) not calibrated
SAME <- setNames(tg$same, RACES)
fwrite(tg, file.path(OUT, "gss2004_targets.csv"))
cat("\nGSS 2004 TARGETS (", nrow(resp), " respondents with numgiven", if (nzchar(TARGET_IDS)) paste0(", from ", TARGET_IDS),
    "; 95% bootstrap CI):\n", sep = ""); print(tg, digits = 3)

# ---- 3. level x1: calibrated dyad-independent draw -------------------------------------------
gen_cal <- function(i) {
  set.seed(SEED_NODES + N + i)
  VA <- pop[sample.int(nrow(pop), N, replace = TRUE), ]; rownames(VA) <- NULL
  ij <- which(upper.tri(matrix(FALSE, N, N)), arr.ind = TRUE)
  gi <- match(VA$race, RACES); g1 <- gi[ij[, 1]]; g2 <- gi[ij[, 2]]; ng <- tabulate(gi, 5)
  same <- g1 == g2
  rest <- with(VA, COEF["sex"] * (sex[ij[, 1]] == sex[ij[, 2]]) + COEF["relig"] * (relig[ij[, 1]] == relig[ij[, 2]]) +
                   COEF["age"] * abs(age[ij[, 1]] - age[ij[, 2]]) + COEF["educ"] * abs(educ_num[ij[, 1]] - educ_num[ij[, 2]])) +
          COEF["race"] * (same & g1 == 5L)
  d_r <- K_MEAN * REL * N / sum(ng * REL); deg_t <- ng * d_r
  W_t <- (deg_t * SAME / 2)[CAL]
  cap <- 0.9 * ng[CAL] * (ng[CAL] - 1) / 2
  capped <- RACES[CAL][W_t > cap]; W_t <- pmin(W_t, cap)     # within-group density above 90%: capped
  Tt <- c(deg_t, W_t)
  pk <- (pmin(g1, g2) - 1L) * 5L + pmax(g1, g2)              # group pair of each dyad
  X <- matrix(0, 25, 9)
  for (a in 1:5) for (b in a:5) { k <- (a - 1) * 5 + b; X[k, a] <- X[k, a] + 1; X[k, b] <- X[k, b] + 1
                                   if (a == b && a %in% CAL) X[k, 5 + a] <- 1 }
  live <- c(ng > 0, ng[CAL] > 1)                             # parameters of absent groups stay at 0
  th <- c(rep(uniroot(function(c) sum(plogis(c + rest + COEF["race"] * same)) - K_MEAN * N / 2, c(-20, 5))$root / 2, 5),
          rep(COEF["race"], 4))
  for (it in 1:300) {
    p <- plogis(rest + drop(X %*% th)[pk])
    sp <- tabulate(pk, 25) * 0; r1 <- rowsum(p, pk); sp[as.integer(rownames(r1))] <- r1[, 1]
    sq <- sp * 0; r2 <- rowsum(p * (1 - p), pk); sq[as.integer(rownames(r2))] <- r2[, 1]
    E <- drop(crossprod(X, sp)); gap <- (Tt - E)[live]
    if (max(abs(gap / Tt[live])) < 1e-6) break
    H <- crossprod(X, X * sq)[live, live]
    step <- solve(H, gap); if (max(abs(step)) > 0.5) step <- 0.5 * step / max(abs(step))   # damped Newton
    th[live] <- th[live] + step
  }
  keep <- runif(length(p)) < p
  g <- make_empty_graph(N, directed = FALSE) + edges(as.vector(t(ij[keep, ])))
  vertex_attr(g) <- as.list(VA)
  list(g = g, theta = setNames(th, c(paste0("delta_", RACES), paste0("gamma_", RACES[CAL]))), iters = it,
       max_rel_gap = max(abs(gap / Tt[live])), capped = capped)
}

# ---- 4. level x2: constrained rewiring toward CLOSURE x the triangles of every node ----------
code <- function(x) as.integer(factor(x)) - 1L
rewire_one <- function(g, i) {
  set.seed(SEED_REWIRE + i)
  args <- list(as_edgelist(g, names = FALSE), vcount(g), code(V(g)$race), code(V(g)$sex), code(V(g)$relig),
               as.numeric(V(g)$age), as.numeric(V(g)$educ_num), target_factor = CLOSURE, tol = 0.01,
               max_props = 2e8, temp = 0.5, mix_props = 5e7)
  r <- if (MODE == "global") do.call(rewire_closure, args) else do.call(rewire_closure_node, c(args, e_tol = E_TOL))
  g2 <- make_empty_graph(vcount(g), directed = FALSE) + edges(as.vector(t(r$edges)))
  vertex_attr(g2) <- vertex_attr(g)
  list(g = g2, info = data.table(net = i, T0 = r$T0, T = r$T, accepted = r$accepted, proposals = r$proposals,
                                 age_change = r$age_rel_change, educ_change = r$educ_rel_change,
                                 node_error = if (is.null(r$E_rel)) NA_real_ else r$E_rel))
}

one_net <- function(i) { x1 <- gen_cal(i); x2 <- rewire_one(x1$g, i); list(x1 = x1, x2 = x2) }
res <- mclapply(seq_len(N_NETS), one_net, mc.cores = N_CORES)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("generator failure: ", res[[which(bad)[1]]])
g1s <- lapply(res, function(x) x$x1$g); g2s <- lapply(res, function(x) x$x2$g)
saveRDS(g1s, file.path(NET_OUT, paste0(SURVEY, "_x1.rds")))
saveRDS(g2s, file.path(NET_OUT, paste0(SURVEY, "_", LEVEL2, ".rds")))
calib <- list(survey = SURVEY, targets = tg, coef_rest = COEF, closure = CLOSURE,
              closure_mode = MODE, level2 = LEVEL2, tag = TAG, target_ids = TARGET_IDS, node_ids = NODE_IDS,
              theta = rbindlist(lapply(res, function(x) as.data.table(as.list(x$x1$theta)))),
              newton = data.table(net = seq_len(N_NETS), iters = vapply(res, function(x) x$x1$iters, integer(1)),
                                  max_rel_gap = vapply(res, function(x) x$x1$max_rel_gap, numeric(1)),
                                  capped = vapply(res, function(x) paste(x$x1$capped, collapse = ","), character(1))),
              rewiring = rbindlist(lapply(res, function(x) x$x2$info)))
saveRDS(calib, file.path(NET_OUT, paste0(SURVEY, if (LEVEL2 == "x2") "" else paste0("_", LEVEL2), "_calibration.rds")))
stamp(N_NETS, " networks x 2 levels in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
cat("\nNEWTON: iterations ", paste(range(calib$newton$iters), collapse = "-"), "; max relative gap ",
    signif(max(calib$newton$max_rel_gap), 2), "; capped within-group targets: ",
    if (any(nzchar(calib$newton$capped))) paste(unique(calib$newton$capped[nzchar(calib$newton$capped)]), collapse = " ") else "none", "\n", sep = "")
cat("REWIRING (", MODE, ", ", LEVEL2, "): triangles x", round(mean(calib$rewiring$T / calib$rewiring$T0), 3), "; age diff change max ",
    sprintf("%+.4f", max(abs(calib$rewiring$age_change))), "; educ ", sprintf("%+.4f", max(abs(calib$rewiring$educ_change))),
    if (MODE == "node") sprintf("; per-node L1 error %.3f (max %.3f)", mean(calib$rewiring$node_error), max(calib$rewiring$node_error)),
    "\n", sep = "")

# ---- 5. checks ---------------------------------------------------------------------------------
mix_of <- function(g) { el <- as_edgelist(g, names = FALSE); r <- V(g)$race
  table(factor(pmin(r[el[, 1]], r[el[, 2]]), RACES), factor(pmax(r[el[, 1]], r[el[, 2]]), RACES)) }
same_count <- function(g, a) { el <- as_edgelist(g, names = FALSE); x <- vertex_attr(g, a); sum(x[el[, 1]] == x[el[, 2]]) }
checks <- rbindlist(lapply(seq_len(N_NETS), function(i) data.table(
  net = i, same_degrees = identical(degree(g1s[[i]]), degree(g2s[[i]])),
  same_race_mixing = identical(as.vector(mix_of(g1s[[i]])), as.vector(mix_of(g2s[[i]]))),
  same_sex_count = same_count(g1s[[i]], "sex") == same_count(g2s[[i]], "sex"),
  same_relig_count = same_count(g1s[[i]], "relig") == same_count(g2s[[i]], "relig"),
  simple_x1 = is_simple(g1s[[i]]), simple_x2 = is_simple(g2s[[i]]))))
SFX <- if (LEVEL2 == "x2") "" else paste0("_", LEVEL2)      # output tables of other levels get a suffix
fwrite(checks, file.path(OUT, paste0(SURVEY, SFX, "_checks.csv")))
cat("\nCHECKS x2 vs x1 (networks passing / ", N_NETS, "): degrees ", sum(checks$same_degrees), ", race mixing ",
    sum(checks$same_race_mixing), ", same-sex ties ", sum(checks$same_sex_count), ", same-religion ties ",
    sum(checks$same_relig_count), ", simple graphs ", sum(checks$simple_x1 & checks$simple_x2), "\n", sep = "")

by_group <- function(g) {
  d <- degree(g); r <- V(g)$race; el <- as_edgelist(g, names = FALSE)
  e2 <- data.table(ego = c(r[el[, 1]], r[el[, 2]]), alter = c(r[el[, 2]], r[el[, 1]]))
  ng <- as.numeric(table(factor(r, RACES)))
  w <- e2[, .(same = mean(ego == alter), within = sum(ego == alter) / 2), by = ego][match(RACES, ego)]
  data.table(ego = RACES, n = ng, rel_deg = as.numeric(tapply(d, factor(r, RACES), mean)) / mean(d),
             same = w$same, within_density = w$within / (ng * (ng - 1) / 2))
}
grp <- rbindlist(lapply(g1s, by_group))[, .(n = mean(n), rel_deg = mean(rel_deg, na.rm = TRUE), same = mean(same, na.rm = TRUE),
                                             within_density = mean(within_density, na.rm = TRUE)), by = ego]
grp <- merge(tg[, .(ego, target_rel_deg = REL[ego], target_same = same, same_lo, same_hi)], grp, by = "ego")[match(RACES, ego)]
fwrite(grp, file.path(OUT, paste0(SURVEY, "_by_group.csv")))      # x1 only: identical for every level
cat("\nBY RACE OF EGO, x1 (= x2 by construction) vs GSS 2004 targets (means over networks):\n"); print(grp, digits = 3)

structure_of <- function(g) {
  d <- degree(g); el <- as_edgelist(g, names = FALSE); va <- vertex_attr(g); lc <- transitivity(g, "local")
  terc <- cut(d, unique(quantile(d, c(0, 1/3, 2/3, 1))), include.lowest = TRUE, labels = FALSE)
  data.table(mean_degree = mean(d), sd_degree = sd(d), max_degree = max(d),
             transitivity = transitivity(g, "global"), local_clustering = mean(lc, na.rm = TRUE),
             local_clust_low_deg = mean(lc[terc == 1], na.rm = TRUE), local_clust_high_deg = mean(lc[terc == max(terc)], na.rm = TRUE),
             same_race = mean(va$race[el[, 1]] == va$race[el[, 2]]), same_sex = mean(va$sex[el[, 1]] == va$sex[el[, 2]]),
             same_relig = mean(va$relig[el[, 1]] == va$relig[el[, 2]]),
             age_diff = mean(abs(va$age[el[, 1]] - va$age[el[, 2]])), educ_diff = mean(abs(va$educ_num[el[, 1]] - va$educ_num[el[, 2]])),
             giant = max(components(g)$csize) / vcount(g), mean_distance = mean_distance(g),
             mean_mur = mean(va$mur_score))
}
st <- rbind(rbindlist(mclapply(g1s, structure_of, mc.cores = N_CORES))[, lapply(.SD, mean)][, level := "x1"],
            rbindlist(mclapply(g2s, structure_of, mc.cores = N_CORES))[, lapply(.SD, mean)][, level := LEVEL2])
fwrite(st, file.path(OUT, paste0(SURVEY, SFX, "_structure.csv")))
cat("\nSTRUCTURE (means over networks):\n"); print(st, digits = 3)
ok <- all(checks$same_degrees, checks$same_race_mixing, checks$same_sex_count, checks$same_relig_count,
          checks$simple_x1, checks$simple_x2)
stamp(if (ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED", "; total ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
if (!ok) quit(status = 1)
