# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NETWORKS FOR THE SWEEP: calibrated to GSS 2004 (x1), then rewired to twice the
# triangles of every node, with the same degrees and the same mixing (x2, the
# main networks).
#
# Replaces the ERGM networks of 02_GSS_network_ergm.R / 02_ATP_network_ergm.R
# for 05-07 (those stay in data/ as the earlier reference). Three changes:
#   1. Nodes are RESAMPLED for every network (with replacement, own seed) from
#      the survey's complete cases. The 100 ERGM networks shared the same 1000.
#   2. Race enters by group (PR #4; extended to every categorical dimension below). Dyad independent:
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
# THE RECIPE SINCE 2026-11-11 (decided with the user): every categorical dimension of the
# confidant data is calibrated BY GROUP, jointly (Smith 2012: "include homophily terms for each
# demographic dimension available in the sampled data"), and survey weights are used from the start.
#   logit p_ij = sum over calibrated dimensions k of
#                  [activity a_k(i) + a_k(j)] + [homophily h_kl  if i and j share level l of k]
#                + common terms: age and education distance (and the common coefficient inside
#                  levels left uncalibrated)
#   Dimensions (01_survey_data.R, CFG$dims):
#     GSS, ATP, GSS14  race (5; 'Other', 18 respondents, keeps the common race term), religion (5:
#                      Protestant, Catholic, Jewish, None, Other), sex (2)
#     ELSOC            education (5), ideology (3; Neither has no homophily term: residual, and
#                      with it the system is nearly flat), religion (5), sex (2)
#   Sex has activity by group and ONE pooled homophily term (target: all within-sex ties).
#   With two groups, degrees and both within-group shares cannot all hold in an undirected
#   network, because confidant reports are directed and do not balance (ELSOC: men's confidants
#   41% women, women's 32% men); the pooled term keeps the degrees and the overall share.
#   Newton runs on the cells (combinations of the calibrated levels; 50 in the GSS, 150 in ELSOC):
#   one activity set carries the constant (the first dimension); the others drop their first level.
#   Targets: relative number of confidants and same-group share by group of ego, weighted (GSS 2004:
#   WTSSNR; ELSOC 2017: ponderador02). Within-group densities above 90% are capped (Jewish,
#   ELSOC 'Other' religion, GSS14 Asian). Common coefficients: weighted case-control
#   (02_GSS_homophily.R: GSS 2004, close to Smith et al.'s; 02_ELSOC_homophily.R: ELSOC 2017).
#   Nodes are drawn with probability proportional to the survey weight (GSS 2004 WTSSNR, ATP
#   weight_ATP, GSS 2014 WTSSNR, ELSOC ponderador02).
#   How much each change moves the premium: 06h_calibration_steps.R. Until 2026-11-11 only race
#   (GSS family) or education and ideology (ELSOC) were calibrated, with Smith et al.'s
#   coefficients for the rest of the GSS family and no weights for GSS 2004 and ATP.
#   Options: --dims=race,relig (a subset of the dimensions; the rest take the common coefficient),
#   --weights=0 (uniform node draws, unweighted targets); both with --tag for the steps of 06h.
#
# ATP 2014 has no ego networks: ATP nodes get the SAME GSS 2004 targets.
# GSS14  GSS 2014 nodes (collective action, 8 of the 9 items), built as ATP; POLVIEWS kept as
#        attribute 'ideol' (3 groups) for the circles of the theta arm (not in the tie model).
# ELSOC  ELSOC 2018 nodes (willingness c09); targets and case-control from the 2017 confidants.
#        No alter-alter ties: x2 is a choice there too. ELSOC and GSS14 networks carry the graph
#        attributes survey, blau and cell (the Gower inputs and theta circles of 05's engine).
# The x2 rewiring keeps the mixing matrix of every categorical attribute exact.
#
# Outputs (<S> = GSS, ATP, GSS14 or ELSOC):
#   data/02_calibrated_networks/<S>_x1.rds, <S>_x2.rds   lists of 24 igraph objects
#       vertex attributes: age, sex, educ_num, race, relig, mur_score, propensity_score
#       (GSS14: + ideol; ELSOC: age, sex, educ, educ_num, relig, ideol, mur_score, propensity_score)
#   data/02_calibrated_networks/<S>_calibration.rds      parameters, Newton and rewiring info
#   output/02_calibrated_networks/gss2004_targets.csv    targets with bootstrap CI (elsoc_targets.csv)
#   output/02_calibrated_networks/<S>_by_group.csv       degree and same-group share by group vs targets
#                                                        (ELSOC also <S>_by_ideol.csv, elsoc_targets_ideol.csv)
#   output/02_calibrated_networks/<S>_structure.csv      structure of x1 and x2 (means over networks)
#   output/02_calibrated_networks/<S>_checks.csv         per network: degrees, mixing, simple graph
#
# Seeds: nodes set.seed(SEED + 1000 + i) and rewiring REWIRE + i, with (SEED, REWIRE) =
# (80000, 70000) GSS, (180000, 170000) ATP, (280000, 270000) ELSOC, (380000, 370000) GSS14.
#
# Usage (repo root):
#   Rscript scripts/02_network_calibrated.R --survey=GSS     # ~15 s on 8 cores
#   Rscript scripts/02_network_calibrated.R --survey=ATP      # also GSS14, ELSOC
#   Rscript scripts/02_network_calibrated.R --survey=GSS --test   # 2 networks, sandbox
# Options: --nets=24 --closure=2 --cores=8
#   --closure-mode=node    every node moves to closure x its own triangles (default;
#                          keeps the clustering-degree profile of x1)
#   --closure-mode=global  the target is the total number of triangles (level x<closure>g)
# Validation runs (scripts/09_network_validation.R; none of these changes the default output):
#   --target-ids=<csv>  GSS 2004 (ELSOC) respondents (column id) from which the targets are computed
#   --node-ids=<csv>    GSS (ELSOC) node pool restricted to these respondents (column id)
#   --tag=<name>        write to data/02_calibrated_networks/<name>/ and
#                       output/02_calibrated_networks/<name>/, and keep the respondent id
#                       of every GSS (ELSOC) node as the vertex attribute gss_id (elsoc_id)
#   --name=<level>      file name of the rewired level (default x<closure>, plus "g" in global mode)
# Alternative degree target (GSS14 only; use with --tag, e.g. --tag=peocntct):
#   --rel-deg=peocntct  relative degree by race from PEOCNTCT (contacts in a typical weekday, GSS 2014
#                       cross-section weighted by WTSSNR, category midpoints 2, 7, 14.5, 34.5, 60) instead of the GSS 2004
#                       relative number of confidants; same-race shares unchanged (03_GSS14_MUR_calculation.R)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({
  library(igraph); library(data.table); library(haven); library(dplyr); library(parallel); library(Rcpp)
}))
Sys.setenv(OMP_NUM_THREADS = "1")
sourceCpp("scripts/02_closure_rewire.cpp")
source("scripts/01_survey_data.R")

ARGS <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) sub(".*=", "", a[1]) else d }
TEST    <- "--test" %in% ARGS
SURVEY  <- toupper(arg("survey", "GSS")); CFG <- survey_config(SURVEY)
N_NETS  <- as.integer(arg("nets", if (TEST) 2 else 24))
CLOSURE <- as.numeric(arg("closure", 2))
N_CORES <- as.integer(arg("cores", 8))
MODE    <- arg("closure-mode", "node"); stopifnot(MODE %in% c("node", "global"))
E_TOL   <- 0.05          # node mode: phase 1 ends when sum_i |t_i - t_i*| <= 5% of sum_i t_i*
LEVEL2  <- arg("name", paste0("x", format(CLOSURE), if (MODE == "global") "g" else ""))
TAG     <- arg("tag", "")
TARGET_IDS <- arg("target-ids", ""); NODE_IDS <- arg("node-ids", "")
REL_DEG <- arg("rel-deg", "confidants"); stopifnot(REL_DEG == "confidants" || (REL_DEG == "peocntct" && SURVEY == "GSS14" && nzchar(arg("tag", ""))))
read_ids <- function(f) { x <- read.csv(f); stopifnot("id" %in% names(x)); as.numeric(x$id) }
N <- 1000L; K_MEAN <- 30
sub_dir <- function(d) if (nzchar(TAG)) file.path(d, TAG) else d
NET_OUT <- sub_dir(if (TEST) "data/02_calibrated_networks/_test" else "data/02_calibrated_networks")
OUT     <- sub_dir(if (TEST) "output/02_calibrated_networks/_test" else "output/02_calibrated_networks")
dir.create(NET_OUT, recursive = TRUE, showWarnings = FALSE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
SEED_NODES  <- c(GSS = 80000L, ATP = 180000L, ELSOC = 280000L, GSS14 = 380000L)[[SURVEY]]
SEED_REWIRE <- c(GSS = 70000L, ATP = 170000L, ELSOC = 270000L, GSS14 = 370000L)[[SURVEY]]
GROUP <- CFG$group; GROUPS <- CFG$groups; G <- length(GROUPS)
CATS <- CFG$cats; NUMS <- CFG$nums
# common coefficients (nodematch scale; age and education per year), weighted case-control of 2004 / 2017:
#   GSS family: 02_GSS_homophily.R (GSS 2004, all ties, WTSSNR);  ELSOC: 02_ELSOC_homophily.R (2017, ponderador02)
# A categorical attribute calibrated by group uses COEF only inside its uncalibrated levels (uncal_common);
# one left out with --dims uses it for every pair (the common term of the recipe before 2026-11-11).
COEF <- if (SURVEY == "ELSOC") c(sex = 0.467, relig = 0.990, ideol = 0.528, age = -0.030, educ = -0.195) else
                                c(race = 1.887, relig = 1.498, sex = 0.267, age = -0.0521, educ = -0.1866)
HOM_ATTRS <- if (SURVEY == "ELSOC") c("sex", "relig", "ideol") else c("race", "relig", "sex")
NUM_TERMS <- c(age = "age", educ = "educ_num")
DIM_SEL <- strsplit(arg("dims", paste(vapply(CFG$dims, `[[`, "", "attr"), collapse = ",")), ",")[[1]]
DIMS <- Filter(function(d) d$attr %in% DIM_SEL, CFG$dims); D <- length(DIMS)
stopifnot(D >= 1, DIMS[[1]]$attr == GROUP)
COMMON <- setdiff(HOM_ATTRS, vapply(DIMS, `[[`, "", "attr"))          # categorical terms with a common coefficient
USE_W <- isTRUE(CFG$weights) && arg("weights", "1") != "0"
CAL <- DIMS[[1]]$hom; UNCAL <- setdiff(seq_len(G), CAL)
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
                             relig == 4 ~ "None", relig %in% 5:13 ~ "OtherRelig", TRUE ~ NA_character_),
           w = as.numeric(zap_labels(wtssnr)))
  for (v in items) { x <- as.numeric(pop[[v]]); x[x %in% c(8, 9)] <- NA; pop[[v]] <- x }
  pop <- pop[complete.cases(pop[, c("age", "sex", "educ_num", "race", "relig", "w", items)]), ]
  pop$propensity_score <- rowSums(4 - as.matrix(pop[, items])) / 27        # 03_GSS_MUR_calculation.R
} else if (SURVEY == "ATP") {
  items <- paste0("METECH_", LETTERS[1:6])
  pop <- readRDS("data/01_ATP_GSS_imputation/ATP_W3_W4_input.rds")
  for (v in items) { x <- as.numeric(zap_labels(pop[[v]])); x[!x %in% c(0, 1)] <- NA; pop[[v]] <- x }
  pop$sex <- as.character(pop$sex); pop$race <- as.character(pop$race); pop$relig <- as.character(pop$relig)
  pop$w <- as.numeric(pop$weight_ATP)
  pop <- pop[complete.cases(pop[, c("age", "sex", "educ_num", "race", "relig", "w", items)]), ]
  pop$propensity_score <- with(pop, METECH_A + (1 - METECH_B) + METECH_C + METECH_D +   # 03_ATP_MUR_calculation.R
                                    (1 - METECH_E) + (1 - METECH_F)) / 6
} else if (SURVEY == "GSS14") {
  pop <- gss14_nodes("ca")                    # 2014 cross-section, 8 items: sum(4 - x) / 24 (03_GSS14_MUR_calculation.R)
} else {
  pop <- elsoc_nodes(3)[, elsoc_id := idencuesta]   # 2018, willingness c09 (03_ELSOC_MUR_calculation.R)
}
pop$mur_score <- 1 - pop$propensity_score                                    # q = aversion
if (nzchar(NODE_IDS)) {
  stopifnot(SURVEY %in% c("GSS", "ELSOC"))
  pop <- if (SURVEY == "GSS") pop[pop$gss_id %in% read_ids(NODE_IDS), ] else pop[pop$elsoc_id %in% read_ids(NODE_IDS), ]
}
VCOLS <- switch(SURVEY,
  ELSOC = c("age", "sex", "educ", "educ_num", "relig", "ideol", "mur_score", "propensity_score", if (nzchar(TAG)) "elsoc_id"),
  GSS14 = c("age", "sex", "educ_num", "race", "relig", "ideol", "mur_score", "propensity_score"),
  c("age", "sex", "educ_num", "race", "relig", "mur_score", "propensity_score", if (nzchar(TAG) && SURVEY == "GSS") "gss_id"))
PROB <- if (USE_W) pop$w else NULL                                       # survey weights: node draw ~ w
pop <- as.data.frame(pop)[, VCOLS]
stopifnot(all(pop[[GROUP]] %in% GROUPS))
stamp(SURVEY, " node pool: ", nrow(pop), " complete cases", if (nzchar(NODE_IDS)) paste0(" (restricted to ", NODE_IDS, ")"),
      "; mean mur_score = ", round(mean(pop$mur_score), 3),
      if (!is.null(PROB)) sprintf(" (weighted %.3f; Kish effective n %d)", weighted.mean(pop$mur_score, PROB), round(sum(PROB)^2 / sum(PROB^2))))

# ---- 2. targets by group of ego, for every calibrated dimension (weighted, bootstrap CI) ------------
# EGO: one row per respondent (id, weight w, number of confidants deg, attributes); ALT: confidants
if (GROUP == "race") {               # GSS 2004 confidants (GSS, ATP, GSS14), weighted by WTSSNR
  numz <- function(x) as.numeric(zap_labels(x))
  EGO <- data.table(id = numz(g4$id), w = numz(g4$wtssnr), deg = numz(g4$numgiven), race = race_gss(g4),
                    relig = RELIG_US[pmin(numz(g4$relig), 5)], sex = SEX2[numz(g4$sex)])
  EGO[deg > 6, deg := NA]
  ALT <- rbindlist(lapply(1:5, function(k) data.table(id = numz(g4$id),
    a_race = c("Asian", "Black", "Hispanic", "White", "Other")[numz(g4[[paste0("race", k)]])],
    a_relig = RELIG_US[numz(g4[[paste0("relig", k)]])], a_sex = SEX2[numz(g4[[paste0("sex", k)]])])))
  TG_PREFIX <- "gss2004_targets"; TG_SRC <- if (USE_W) "GSS 2004 (weighted)" else "GSS 2004"
} else {                             # ELSOC confidants (r13, 2017), weighted by ponderador02
  EGO <- elsoc_egos(2)[, .(id = idencuesta, w, deg = nredes, educ, ideol, relig, sex)]
  ALT <- elsoc_alters(2)[, .(id = idencuesta, a_educ = educ, a_ideol = ideol, a_relig = relig, a_sex = sex)]
  TG_PREFIX <- "elsoc_targets"; TG_SRC <- if (USE_W) "ELSOC 2017 (weighted)" else "ELSOC 2017"
}
EGO <- EGO[!is.na(deg) & !is.na(w)]
if (nzchar(TARGET_IDS)) EGO <- EGO[id %in% read_ids(TARGET_IDS)]
if (!USE_W) EGO[, w := 1]
EGO[, w := w / mean(w)]
DATTR <- vapply(DIMS, `[[`, "", "attr")
ALT <- merge(ALT, setnames(EGO[, c("id", "w", DATTR), with = FALSE], DATTR, paste0("e_", DATTR)), by = "id")
setkey(EGO, id); setkey(ALT, id)
targets_of <- function(ids, v, lv) {     # ids: respondents (a bootstrap draw repeats them)
  r <- EGO[.(ids), nomatch = 0L][!is.na(get(v))]
  a <- ALT[.(ids), nomatch = 0L, allow.cartesian = TRUE][!is.na(get(paste0("a_", v))) & !is.na(get(paste0("e_", v)))]
  d <- r[, .(n = .N, wn = sum(w), deg = weighted.mean(deg, w)), by = .(ego = get(v))][, rel_deg := deg / weighted.mean(deg, wn)]
  s <- a[, .(alters = .N, same = weighted.mean(get(paste0("a_", v)) == get(paste0("e_", v)), w)), by = .(ego = get(paste0("e_", v)))]
  merge(d[, .(ego, n, deg, rel_deg)], s, by = "ego", all = TRUE)[match(lv, ego)]
}
add_ci <- function(tg, boot, lv) merge(tg, boot[, .(rel_lo = quantile(rel_deg, .025, na.rm = TRUE), rel_hi = quantile(rel_deg, .975, na.rm = TRUE),
                                                    same_lo = quantile(same, .025, na.rm = TRUE), same_hi = quantile(same, .975, na.rm = TRUE)), by = ego],
                                      by = "ego")[match(lv, ego)]
uid <- unique(EGO$id); N_RESP <- length(uid)
TG <- list()
for (k in seq_len(D)) {
  dd <- DIMS[[k]]; set.seed(c(race = 2004, educ = 2017, ideol = 2018, relig = 2019, sex = 2020)[[dd$attr]])
  boot <- rbindlist(lapply(1:1000, function(b) targets_of(sample(uid, replace = TRUE), dd$attr, dd$levels)[, b := b]))
  TG[[dd$attr]] <- add_ci(targets_of(uid, dd$attr, dd$levels), boot, dd$levels)
}
if (REL_DEG == "peocntct") {         # GSS14 sensitivity: relative daily contacts by race (GSS 2014 cross-section)
  g14 <- gss_load(2014); pc <- data.table(race = race_gss(g14), x = as.numeric(zap_labels(g14$peocntct)), w = as.numeric(zap_labels(g14$wtssnr)))[
    as.numeric(g14$samptype) == 2014 & x %in% 1:5 & !is.na(race)][, x := c(2, 7, 14.5, 34.5, 60)[x]]
  pcr <- pc[, .(wn = sum(w), m = weighted.mean(x, w)), by = race][, rel := m / weighted.mean(m, wn)]     # weighted by WTSSNR
  TG$race[, rel_deg := pcr$rel[match(ego, pcr$race)]][, `:=`(rel_lo = NA_real_, rel_hi = NA_real_)]
  cat("\nRELATIVE DEGREE TARGETS REPLACED BY PEOCNTCT (GSS 2014):\n"); print(TG$race[, .(ego, rel_deg, same)], digits = 3)
}
for (k in seq_len(D)) {               # activity targets (uncalibrated levels: 1) and within-group shares
  dd <- DIMS[[k]]; tgk <- TG[[dd$attr]]; H <- length(dd$levels)
  DIMS[[k]]$REL <- setNames(ifelse(seq_len(H) %in% dd$act_lv, tgk$rel_deg, 1), dd$levels)
  DIMS[[k]]$SAME <- setNames(tgk$same, dd$levels)
  DIMS[[k]]$g0 <- if (dd$attr %in% HOM_ATTRS) COEF[[dd$attr]] else 0.5          # Newton start of its homophily terms
  fwrite(tgk, file.path(OUT, paste0(TG_PREFIX, if (k == 1) "" else paste0("_", dd$attr), ".csv")))
  cat("\n", TG_SRC, " TARGETS BY ", toupper(dd$attr), " (", N_RESP, " respondents", if (nzchar(TARGET_IDS)) paste0(", from ", TARGET_IDS),
      "; 95% bootstrap CI; activity calibrated for ", paste(dd$levels[dd$act_lv], collapse = ", "), "; homophily for ",
      paste(dd$levels[dd$hom], collapse = ", "), "):\n", sep = ""); print(tgk, digits = 3)
}
tg <- TG[[GROUP]]; REL <- DIMS[[1]]$REL; SAME <- DIMS[[1]]$SAME; TG_FILE <- paste0(TG_PREFIX, ".csv")
cat("\nCommon terms: ", paste(c(COMMON, names(NUM_TERMS)), collapse = ", "), "; weights ", if (USE_W) "on" else "off", "\n", sep = "")

# ---- 3. level x1: calibrated dyad-independent draw -------------------------------------------
gen_cal <- function(i) {
  set.seed(SEED_NODES + N + i)
  VA <- pop[sample.int(nrow(pop), N, replace = TRUE, prob = PROB), ]; rownames(VA) <- NULL
  ij <- which(upper.tri(matrix(FALSE, N, N)), arr.ind = TRUE)
  Hs <- vapply(DIMS, function(d) length(d$levels), 1L)
  idx <- lapply(DIMS, function(d) match(VA[[d$attr]], d$levels)); nn <- lapply(seq_len(D), function(k) tabulate(idx[[k]], Hs[k]))
  i1 <- lapply(idx, function(x) x[ij[, 1]]); same <- lapply(seq_len(D), function(k) i1[[k]] == idx[[k]][ij[, 2]])
  # rest of the tie model: common categorical terms, continuous distances, common term inside uncalibrated levels
  terms <- c(lapply(COMMON, function(t) COEF[[t]] * (VA[[t]][ij[, 1]] == VA[[t]][ij[, 2]])),
             lapply(names(NUM_TERMS), function(t) COEF[[t]] * abs(VA[[NUM_TERMS[[t]]]][ij[, 1]] - VA[[NUM_TERMS[[t]]]][ij[, 2]])))
  rest <- Reduce(`+`, terms)
  for (k in seq_len(D)) { un <- setdiff(seq_len(Hs[k]), DIMS[[k]]$hom)
    if (isTRUE(DIMS[[k]]$uncal_common) && length(un)) rest <- rest + COEF[[DIMS[[k]]$attr]] * (same[[k]] & i1[[k]] %in% un) }
  # cells = combinations of the calibrated attributes; X maps cell pairs to parameters
  stride <- as.integer(rev(cumprod(rev(c(Hs[-1], 1L))))); NC <- as.integer(prod(Hs))
  ci <- 1L + Reduce(`+`, lapply(seq_len(D), function(k) (idx[[k]] - 1L) * stride[k]))
  LV <- matrix(vapply(seq_len(D), function(k) ((seq_len(NC) - 1L) %/% stride[k]) %% Hs[k] + 1L, integer(NC)), ncol = D)
  pr <- which(upper.tri(matrix(TRUE, NC, NC), diag = TRUE), arr.ind = TRUE); kk <- (pr[, 1] - 1L) * NC + pr[, 2]
  la <- LV[pr[, 1], , drop = FALSE]; lb <- LV[pr[, 2], , drop = FALSE]
  # one design column per parameter; Newton matches the same statistics (degree sums, within-group ties).
  # 'pooled' dimensions (sex) have one homophily term whose target is the TOTAL within-group ties
  cols <- list(); Tt <- numeric(0); live <- logical(0); th0 <- numeric(0); nm <- character(0); capped <- character(0)
  for (k in seq_len(D)) {
    dd <- DIMS[[k]]; H <- Hs[k]; nk <- nn[[k]]
    deg_k <- nk * K_MEAN * dd$REL * N / sum(nk * dd$REL)                       # degree sum by group
    if (k == 1 || dd$act) {
      L <- if (k == 1) seq_len(H) else seq_len(H)[-1]                          # one activity set carries the constant
      for (l in L) cols[[length(cols) + 1]] <- (la[, k] == l) + (lb[, k] == l)
      Tt <- c(Tt, deg_k[L]); live <- c(live, nk[L] > 0); th0 <- c(th0, rep(0, length(L)))
      nm <- c(nm, if (k == 1) paste0("delta_", dd$levels) else paste0("eta_", dd$attr, "_", dd$levels[L]))
    }
    W <- (deg_k * dd$SAME / 2)[dd$hom]; cap <- 0.9 * nk[dd$hom] * (nk[dd$hom] - 1) / 2
    capped <- c(capped, dd$levels[dd$hom][W > cap]); W <- pmin(W, cap)        # within-group density above 90%: capped
    if (isTRUE(dd$pooled)) {
      cols[[length(cols) + 1]] <- as.numeric(la[, k] == lb[, k] & la[, k] %in% dd$hom)
      Tt <- c(Tt, sum(W)); live <- c(live, TRUE); th0 <- c(th0, dd$g0); nm <- c(nm, paste0("kappa_", dd$attr))
    } else {
      for (l in dd$hom) cols[[length(cols) + 1]] <- as.numeric(la[, k] == l & lb[, k] == l)
      Tt <- c(Tt, W); live <- c(live, nk[dd$hom] > 1); th0 <- c(th0, rep(dd$g0, length(dd$hom)))
      nm <- c(nm, paste0(if (k == 1) "gamma_" else paste0("kappa_", dd$attr, "_"), dd$levels[dd$hom]))
    }
  }
  X <- matrix(0, NC * NC, length(cols)); X[kk, ] <- do.call(cbind, cols)
  pk <- (pmin(ci[ij[, 1]], ci[ij[, 2]]) - 1L) * NC + pmax(ci[ij[, 1]], ci[ij[, 2]])     # cell pair of each dyad
  hom0 <- Reduce(`+`, lapply(seq_len(D), function(k) DIMS[[k]]$g0 * (same[[k]] & i1[[k]] %in% DIMS[[k]]$hom)))
  th <- th0; th[seq_len(G)] <- uniroot(function(c) sum(plogis(c + rest + hom0)) - K_MEAN * N / 2, c(-20, 5))$root / 2
  for (it in 1:300) {
    p <- plogis(rest + drop(X %*% th)[pk])
    sp <- numeric(NC * NC); r1 <- rowsum(p, pk); sp[as.integer(rownames(r1))] <- r1[, 1]
    sq <- numeric(NC * NC); r2 <- rowsum(p * (1 - p), pk); sq[as.integer(rownames(r2))] <- r2[, 1]
    E <- drop(crossprod(X, sp)); gap <- (Tt - E)[live]
    if (max(abs(gap / Tt[live])) < 1e-6) break
    Hm <- crossprod(X, X * sq)[live, live]
    step <- solve(Hm, gap); if (max(abs(step)) > 0.5) step <- 0.5 * step / max(abs(step))   # damped Newton
    th[live] <- th[live] + step
  }
  keep <- runif(length(p)) < p
  g <- make_empty_graph(N, directed = FALSE) + edges(as.vector(t(ij[keep, ])))
  vertex_attr(g) <- as.list(VA)
  if (!SURVEY %in% c("GSS", "ATP")) { g$survey <- SURVEY; g$blau <- CFG$blau; g$cell <- CFG$cell }   # read by 05's engine
  list(g = g, theta = setNames(th, nm), iters = it, max_rel_gap = max(abs(gap / Tt[live])), capped = unique(capped))
}

# ---- 4. level x2: constrained rewiring toward CLOSURE x the triangles of every node ----------
code <- function(x) as.integer(factor(x)) - 1L
rewire_one <- function(g, i) {
  set.seed(SEED_REWIRE + i)
  cats <- sapply(CATS, function(a) code(vertex_attr(g, a))); storage.mode(cats) <- "integer"
  nums <- sapply(NUMS, function(a) as.numeric(vertex_attr(g, a)))
  args <- list(as_edgelist(g, names = FALSE), vcount(g), cats, nums, target_factor = CLOSURE, tol = 0.01,
               max_props = 2e8, temp = 0.5, mix_props = 5e7)
  r <- if (MODE == "global") do.call(rewire_closure, args) else do.call(rewire_closure_node, c(args, e_tol = E_TOL))
  g2 <- make_empty_graph(vcount(g), directed = FALSE) + edges(as.vector(t(r$edges)))
  vertex_attr(g2) <- vertex_attr(g)
  if (!SURVEY %in% c("GSS", "ATP")) graph_attr(g2) <- graph_attr(g)
  list(g = g2, info = data.table(net = i, T0 = r$T0, T = r$T, accepted = r$accepted, proposals = r$proposals,
                                 age_change = r$num_rel_change[1], educ_change = r$num_rel_change[2],
                                 node_error = if (is.null(r$E_rel)) NA_real_ else r$E_rel))
}

one_net <- function(i) { x1 <- gen_cal(i); x2 <- rewire_one(x1$g, i); list(x1 = x1, x2 = x2) }
res <- mclapply(seq_len(N_NETS), one_net, mc.cores = N_CORES)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("generator failure: ", res[[which(bad)[1]]])
g1s <- lapply(res, function(x) x$x1$g); g2s <- lapply(res, function(x) x$x2$g)
saveRDS(g1s, file.path(NET_OUT, paste0(SURVEY, "_x1.rds")))
saveRDS(g2s, file.path(NET_OUT, paste0(SURVEY, "_", LEVEL2, ".rds")))
calib <- list(survey = SURVEY, targets = tg, targets_dims = TG, dims = DIM_SEL, common = COMMON, weighted = USE_W, coef_rest = COEF, closure = CLOSURE,
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
mix_of <- function(g) { el <- as_edgelist(g, names = FALSE); r <- vertex_attr(g, GROUP)
  table(factor(pmin(r[el[, 1]], r[el[, 2]]), GROUPS), factor(pmax(r[el[, 1]], r[el[, 2]]), GROUPS)) }
same_count <- function(g, a) { el <- as_edgelist(g, names = FALSE); x <- vertex_attr(g, a); sum(x[el[, 1]] == x[el[, 2]]) }
OTHER_CATS <- setdiff(CATS, GROUP)
checks <- rbindlist(lapply(seq_len(N_NETS), function(i) {
  x <- data.table(net = i, same_degrees = identical(degree(g1s[[i]]), degree(g2s[[i]])),
                  same_group_mixing = identical(as.vector(mix_of(g1s[[i]])), as.vector(mix_of(g2s[[i]]))))
  for (a in OTHER_CATS) set(x, j = paste0("same_", a, "_count"), value = same_count(g1s[[i]], a) == same_count(g2s[[i]], a))
  x[, `:=`(simple_x1 = is_simple(g1s[[i]]), simple_x2 = is_simple(g2s[[i]]))][]
}))
setnames(checks, "same_group_mixing", paste0("same_", GROUP, "_mixing"))
SFX <- if (LEVEL2 == "x2") "" else paste0("_", LEVEL2)      # output tables of other levels get a suffix
fwrite(checks, file.path(OUT, paste0(SURVEY, SFX, "_checks.csv")))
cat("\nCHECKS x2 vs x1 (networks passing / ", N_NETS, "): degrees ", sum(checks$same_degrees), ", ", GROUP, " mixing ",
    sum(checks[[paste0("same_", GROUP, "_mixing")]]),
    paste0(", same-", OTHER_CATS, " ties ", vapply(OTHER_CATS, function(a) sum(checks[[paste0("same_", a, "_count")]]), 0), collapse = ""),
    ", simple graphs ", sum(checks$simple_x1 & checks$simple_x2), "\n", sep = "")

by_group <- function(g, GROUP = get("GROUP", globalenv()), GROUPS = get("GROUPS", globalenv())) {
  d <- degree(g); r <- vertex_attr(g, GROUP); el <- as_edgelist(g, names = FALSE)
  e2 <- data.table(ego = c(r[el[, 1]], r[el[, 2]]), alter = c(r[el[, 2]], r[el[, 1]]))
  ng <- as.numeric(table(factor(r, GROUPS)))
  w <- e2[, .(same = mean(ego == alter), within = sum(ego == alter) / 2), by = ego][match(GROUPS, ego)]
  data.table(ego = GROUPS, n = ng, rel_deg = as.numeric(tapply(d, factor(r, GROUPS), mean)) / mean(d),
             same = w$same, within_density = w$within / (ng * (ng - 1) / 2))
}
grp <- rbindlist(lapply(g1s, by_group))[, .(n = mean(n), rel_deg = mean(rel_deg, na.rm = TRUE), same = mean(same, na.rm = TRUE),
                                             within_density = mean(within_density, na.rm = TRUE)), by = ego]
grp <- merge(tg[, .(ego, target_rel_deg = REL[ego], target_same = same, same_lo, same_hi)], grp, by = "ego")[match(GROUPS, ego)]
fwrite(grp, file.path(OUT, paste0(SURVEY, "_by_group.csv")))      # x1 only: identical for every level
cat("\nBY ", toupper(GROUP), " OF EGO, x1 (= x2 by construction) vs ", TG_SRC, " targets (means over networks):\n", sep = "")
print(grp, digits = 3)
for (k in seq_len(D)[-1]) {
  dd <- DIMS[[k]]
  grp2 <- rbindlist(lapply(g1s, by_group, GROUP = dd$attr, GROUPS = dd$levels))[
    , .(n = mean(n), rel_deg = mean(rel_deg, na.rm = TRUE), same = mean(same, na.rm = TRUE), within_density = mean(within_density, na.rm = TRUE)), by = ego]
  grp2 <- merge(TG[[dd$attr]][, .(ego, target_rel_deg = dd$REL[ego], target_same = same, same_lo, same_hi)], grp2, by = "ego")[match(dd$levels, ego)][
    , `:=`(rel_targeted = dd$act & ego %in% dd$levels[dd$act_lv], same_targeted = ego %in% dd$levels[dd$hom] & !isTRUE(dd$pooled))]
  fwrite(grp2, file.path(OUT, paste0(SURVEY, "_by_", dd$attr, ".csv")))
  cat("\nBY ", toupper(dd$attr), " OF EGO, x1 vs targets:\n", sep = ""); print(grp2, digits = 3)
}

structure_of <- function(g) {
  d <- degree(g); el <- as_edgelist(g, names = FALSE); va <- vertex_attr(g); lc <- transitivity(g, "local")
  terc <- cut(d, unique(quantile(d, c(0, 1/3, 2/3, 1))), include.lowest = TRUE, labels = FALSE)
  as.data.table(c(list(mean_degree = mean(d), sd_degree = sd(d), max_degree = max(d),
                       transitivity = transitivity(g, "global"), local_clustering = mean(lc, na.rm = TRUE),
                       local_clust_low_deg = mean(lc[terc == 1], na.rm = TRUE), local_clust_high_deg = mean(lc[terc == max(terc)], na.rm = TRUE)),
                  setNames(lapply(CATS, function(a) mean(va[[a]][el[, 1]] == va[[a]][el[, 2]])), paste0("same_", CATS)),
                  list(age_diff = mean(abs(va$age[el[, 1]] - va$age[el[, 2]])), educ_diff = mean(abs(va$educ_num[el[, 1]] - va$educ_num[el[, 2]])),
                       giant = max(components(g)$csize) / vcount(g), mean_distance = mean_distance(g),
                       mean_mur = mean(va$mur_score))))
}
st <- rbind(rbindlist(mclapply(g1s, structure_of, mc.cores = N_CORES))[, lapply(.SD, mean)][, level := "x1"],
            rbindlist(mclapply(g2s, structure_of, mc.cores = N_CORES))[, lapply(.SD, mean)][, level := LEVEL2])
fwrite(st, file.path(OUT, paste0(SURVEY, SFX, "_structure.csv")))
cat("\nSTRUCTURE (means over networks):\n"); print(st, digits = 3)
ok <- all(unlist(checks[, -"net"]))
stamp(if (ok) "ALL CHECKS PASSED" else "SOME CHECKS FAILED", "; total ",
      round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
if (!ok) quit(status = 1)
