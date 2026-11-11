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
# PR #5 adds two surveys (configuration and recodes in 01_survey_data.R):
#   GSS14  GSS 2014 nodes (collective action, 8 of the 9 items), built as ATP: GSS 2004
#          race targets and Smith et al. coefficients; POLVIEWS kept as attribute 'ideol'
#          (3 groups) for the circles of the theta arm (it is not in the tie model).
#   ELSOC  ELSOC 2018 nodes (willingness c09). TWO dimensions are calibrated by group,
#          jointly: EDUCATION (5 groups) and IDEOLOGY (3 groups: Left, Neither, Right):
#            logit p_ij = delta_e(i) + delta_e(j) + gamma_e [educ_i = educ_j = e]
#                       + eta_l(i) + eta_l(j) + kappa_l [ideol_i = ideol_j = l, l = Left, Right] + rest
#          (eta_Left = 0: the two activity sets share one constant). Newton matches, for the
#          nodes of each network, the degree sum and the within-group ties of every
#          education group, the degree sum of every ideology group and the within-group ties
#          of Left and Right to the ELSOC confidants of 2017 (relative number of confidants
#          and same-group share by group of ego, weighted by ponderador02;
#          02_ELSOC_homophily.R). Neither has no homophily term (as 'Other' in GSS has no
#          target): it is a residual category (no left-right identification, same / chance
#          1.14), its same-group share follows from its size, and with its own term the
#          Newton system is nearly flat (the left-right ties implied by the targets are close
#          to zero; on the held-out halves Newton diverges). The rest of the tie model comes from the ELSOC
#          case-control estimates (all ties, 2017, weighted): same sex, same religion,
#          |age| and |education in years|. Nodes are drawn with probability proportional to
#          ponderador02. The x2 rewiring keeps the education, sex, religion and ideology
#          mixing matrices exact. ELSOC has no alter-alter ties: x2 is a choice there too.
#          Networks carry the graph attributes survey, blau and cell (the Gower inputs and
#          the theta circles of 05_unified_engine.R); GSS and ATP keep the defaults.
#   GSS14 nodes are drawn with probability proportional to WTSSNR (GSS 2004 and ATP: uniform).
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
GROUP <- CFG$group; GROUPS <- CFG$groups; G <- length(GROUPS); CAL <- CFG$cal      # groups with their own targets
UNCAL <- setdiff(seq_len(G), CAL)                                                  # GSS: 'Other'
CATS <- CFG$cats; NUMS <- CFG$nums
if (SURVEY == "ELSOC") {
  # ELSOC case-control (02_ELSOC_homophily.R: all ties, 2017, weighted, B = 1000), nodematch scale;
  # educ per year of the 5 education categories (6, 10, 12, 14, 16). Same ideology (0.528 there) is
  # not a common term: ideology is calibrated by group (start value GAMMA2_0)
  COEF <- c(sex = 0.467, relig = 0.990, age = -0.030, educ = -0.195)
  CAT_TERMS <- c(sex = "sex", relig = "relig"); GAMMA0 <- 0.5; GAMMA2_0 <- 0.528
} else {
  # 02's homophily coefficients (Smith et al. 2014, Table 3 model 2, effects for 2004), nodematch scale
  COEF <- c(race = 1.959 - 0.264, sex = 0.373 - 0.144, relig = 1.270 + 0.215,
            age = -0.047 - 0.005, educ = -0.157 - 0.044)
  CAT_TERMS <- c(sex = "sex", relig = "relig"); GAMMA0 <- COEF[["race"]]
}
NUM_TERMS <- c(age = "age", educ = "educ_num")
TWO <- !is.null(CFG$group2)                                   # second calibrated dimension (ELSOC: ideology)
if (TWO) { GROUP2 <- CFG$group2; GROUPS2 <- CFG$groups2; H2 <- length(GROUPS2); CAL2 <- CFG$cal2 }
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
} else if (SURVEY == "ATP") {
  items <- paste0("METECH_", LETTERS[1:6])
  pop <- readRDS("data/01_ATP_GSS_imputation/ATP_W3_W4_input.rds")
  for (v in items) { x <- as.numeric(zap_labels(pop[[v]])); x[!x %in% c(0, 1)] <- NA; pop[[v]] <- x }
  pop$sex <- as.character(pop$sex); pop$race <- as.character(pop$race); pop$relig <- as.character(pop$relig)
  pop <- pop[complete.cases(pop[, c("age", "sex", "educ_num", "race", "relig", items)]), ]
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
PROB <- if (isTRUE(CFG$weights)) pop$w else NULL                         # survey weights: node draw ~ w
pop <- as.data.frame(pop)[, VCOLS]
stopifnot(all(pop[[GROUP]] %in% GROUPS))
stamp(SURVEY, " node pool: ", nrow(pop), " complete cases", if (nzchar(NODE_IDS)) paste0(" (restricted to ", NODE_IDS, ")"),
      "; mean mur_score = ", round(mean(pop$mur_score), 3),
      if (!is.null(PROB)) sprintf(" (weighted %.3f; Kish effective n %d)", weighted.mean(pop$mur_score, PROB), round(sum(PROB)^2 / sum(PROB^2))))

# ---- 2. targets by group of ego (with bootstrap CI) -----------------------------------------
if (GROUP == "race") {               # GSS 2004 confidants, by race of ego (GSS, ATP, GSS14)
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
TG_FILE <- "gss2004_targets.csv"; TG_SRC <- "GSS 2004"; N_RESP <- nrow(resp)
} else {                             # ELSOC confidants (r13, 2017), by education and by ideology of ego, weighted
eg_r <- elsoc_egos(2)[!is.na(nredes)]
if (nzchar(TARGET_IDS)) eg_r <- eg_r[idencuesta %in% read_ids(TARGET_IDS)]
al_r <- merge(elsoc_alters(2)[, .(idencuesta, a_educ = educ, a_ideol = ideol)],
              eg_r[, .(idencuesta, w, e_educ = educ, e_ideol = ideol)], by = "idencuesta")
setkey(eg_r, idencuesta); setkey(al_r, idencuesta)
elsoc_targets <- function(ids, v, lv) {     # ids: respondents (a bootstrap draw repeats them)
  r <- eg_r[.(ids), nomatch = 0L][!is.na(get(v))]
  a <- al_r[.(ids), nomatch = 0L, allow.cartesian = TRUE][!is.na(get(paste0("a_", v))) & !is.na(get(paste0("e_", v)))]
  d <- r[, .(n = .N, wn = sum(w), deg = weighted.mean(nredes, w)), by = .(ego = get(v))][, rel_deg := deg / weighted.mean(deg, wn)]
  s <- a[, .(alters = .N, same = weighted.mean(get(paste0("a_", v)) == get(paste0("e_", v)), w)), by = .(ego = get(paste0("e_", v)))]
  merge(d[, .(ego, n, deg, rel_deg)], s, by = "ego", all = TRUE)[match(lv, ego)]
}
uid <- unique(eg_r$idencuesta)
tg <- elsoc_targets(uid, GROUP, GROUPS)
set.seed(2017)
boot <- rbindlist(lapply(1:1000, function(b) elsoc_targets(sample(uid, replace = TRUE), GROUP, GROUPS)[, b := b]))
tg2 <- elsoc_targets(uid, GROUP2, GROUPS2)
set.seed(2018)
boot2 <- rbindlist(lapply(1:1000, function(b) elsoc_targets(sample(uid, replace = TRUE), GROUP2, GROUPS2)[, b := b]))
TG_FILE <- "elsoc_targets.csv"; TG_SRC <- "ELSOC 2017 (weighted)"; N_RESP <- length(uid)
}
add_ci <- function(tg, boot, lv) merge(tg, boot[, .(rel_lo = quantile(rel_deg, .025, na.rm = TRUE), rel_hi = quantile(rel_deg, .975, na.rm = TRUE),
                                                    same_lo = quantile(same, .025, na.rm = TRUE), same_hi = quantile(same, .975, na.rm = TRUE)), by = ego],
                                      by = "ego")[match(lv, ego)]
tg <- add_ci(tg, boot, GROUPS)
if (REL_DEG == "peocntct") {         # GSS14 sensitivity: relative daily contacts by race (GSS 2014 cross-section)
  g14 <- gss_load(2014); pc <- data.table(race = race_gss(g14), x = as.numeric(zap_labels(g14$peocntct)), w = as.numeric(zap_labels(g14$wtssnr)))[
    as.numeric(g14$samptype) == 2014 & x %in% 1:5 & !is.na(race)][, x := c(2, 7, 14.5, 34.5, 60)[x]]
  pcr <- pc[, .(wn = sum(w), m = weighted.mean(x, w)), by = race][, rel := m / weighted.mean(m, wn)]     # weighted by WTSSNR
  tg[, rel_deg := pcr$rel[match(ego, pcr$race)]][, `:=`(rel_lo = NA_real_, rel_hi = NA_real_)]
  cat("\nRELATIVE DEGREE TARGETS REPLACED BY PEOCNTCT (GSS 2014):\n"); print(tg[, .(ego, rel_deg, same)], digits = 3)
}
REL  <- setNames(ifelse(seq_len(G) %in% CAL, tg$rel_deg, 1), GROUPS)  # uncalibrated ('Other', 18 GSS respondents): 1
SAME <- setNames(tg$same, GROUPS)
fwrite(tg, file.path(OUT, TG_FILE))
cat("\n", TG_SRC, " TARGETS (", N_RESP, " respondents", if (nzchar(TARGET_IDS)) paste0(", from ", TARGET_IDS),
    "; 95% bootstrap CI):\n", sep = ""); print(tg, digits = 3)
if (TWO) {
  tg2 <- add_ci(tg2, boot2, GROUPS2); REL2 <- setNames(tg2$rel_deg, GROUPS2); SAME2 <- setNames(tg2$same, GROUPS2)
  fwrite(tg2, file.path(OUT, sub("\\.csv$", paste0("_", GROUP2, ".csv"), TG_FILE)))
  cat("\n", TG_SRC, " TARGETS BY ", toupper(GROUP2), ":\n", sep = ""); print(tg2, digits = 3)
}

# ---- 3. level x1: calibrated dyad-independent draw -------------------------------------------
gen_cal <- function(i) {
  set.seed(SEED_NODES + N + i)
  VA <- pop[sample.int(nrow(pop), N, replace = TRUE, prob = PROB), ]; rownames(VA) <- NULL
  ij <- which(upper.tri(matrix(FALSE, N, N)), arr.ind = TRUE)
  gi <- match(VA[[GROUP]], GROUPS); g1 <- gi[ij[, 1]]; g2 <- gi[ij[, 2]]; ng <- tabulate(gi, G)
  same <- g1 == g2
  # cells: the groups (one dimension) or the education x ideology combinations (two)
  if (TWO) { hi <- match(VA[[GROUP2]], GROUPS2); nh <- tabulate(hi, H2); same2 <- hi[ij[, 1]] == hi[ij[, 2]]
             ci <- (gi - 1L) * H2 + hi; NC <- G * H2 } else { ci <- gi; NC <- G }
  cg <- function(a) if (TWO) (a - 1L) %/% H2 + 1L else a          # first-dimension group of cell a
  ch <- function(a) (a - 1L) %% H2 + 1L                              # second-dimension group of cell a
  # rest of the tie model, in the order sex, relig (, ideol), age, educ (GSS: then the common race
  # coefficient inside 'Other'), so GSS and ATP sums are bit-identical to the earlier version
  terms <- c(lapply(names(CAT_TERMS), function(t) COEF[t] * (VA[[CAT_TERMS[[t]]]][ij[, 1]] == VA[[CAT_TERMS[[t]]]][ij[, 2]])),
             lapply(names(NUM_TERMS), function(t) COEF[t] * abs(VA[[NUM_TERMS[[t]]]][ij[, 1]] - VA[[NUM_TERMS[[t]]]][ij[, 2]])))
  rest <- Reduce(`+`, terms)
  if (length(UNCAL)) rest <- rest + COEF[GROUP] * (same & g1 %in% UNCAL)
  d_r <- K_MEAN * REL * N / sum(ng * REL); deg_t <- ng * d_r
  W_t <- (deg_t * SAME / 2)[CAL]
  cap <- 0.9 * ng[CAL] * (ng[CAL] - 1) / 2
  capped <- GROUPS[CAL][W_t > cap]; W_t <- pmin(W_t, cap)    # within-group density above 90%: capped
  Tt <- c(deg_t, W_t)
  c1 <- ci[ij[, 1]]; c2 <- ci[ij[, 2]]
  pk <- (pmin(c1, c2) - 1L) * NC + pmax(c1, c2)              # cell pair of each dyad
  P1 <- G + length(CAL)                                      # columns: delta (G), gamma (CAL) [, eta (H2 - 1), kappa (CAL2)]
  X <- matrix(0, NC * NC, P1 + if (TWO) H2 - 1L + length(CAL2) else 0L)
  for (a in 1:NC) for (b in a:NC) { k <- (a - 1) * NC + b; ga <- cg(a); gb <- cg(b)
    X[k, ga] <- X[k, ga] + 1; X[k, gb] <- X[k, gb] + 1
    if (ga == gb && ga %in% CAL) X[k, G + match(ga, CAL)] <- 1
    if (TWO) { ha <- ch(a); hb <- ch(b)
      if (ha > 1) X[k, P1 + ha - 1] <- X[k, P1 + ha - 1] + 1; if (hb > 1) X[k, P1 + hb - 1] <- X[k, P1 + hb - 1] + 1
      if (ha == hb && ha %in% CAL2) X[k, P1 + H2 - 1 + match(ha, CAL2)] <- 1 } }
  live <- c(ng > 0, ng[CAL] > 1)                             # parameters of absent groups stay at 0
  hom0 <- GAMMA0 * same
  if (TWO) {
    d_h <- K_MEAN * REL2 * N / sum(nh * REL2); deg_h <- nh * d_h
    W_h <- (deg_h * SAME2 / 2)[CAL2]; cap2 <- 0.9 * nh[CAL2] * (nh[CAL2] - 1) / 2
    capped <- c(capped, GROUPS2[CAL2][W_h > cap2]); W_h <- pmin(W_h, cap2)
    Tt <- c(Tt, deg_h[-1], W_h); live <- c(live, nh[-1] > 0, nh[CAL2] > 1)
    hom0 <- hom0 + GAMMA2_0 * (same2 & hi[ij[, 1]] %in% CAL2)
  }
  th <- c(rep(uniroot(function(c) sum(plogis(c + rest + hom0)) - K_MEAN * N / 2, c(-20, 5))$root / 2, G),
          rep(GAMMA0, length(CAL)), if (TWO) c(rep(0, H2 - 1), rep(GAMMA2_0, length(CAL2))))
  for (it in 1:300) {
    p <- plogis(rest + drop(X %*% th)[pk])
    sp <- tabulate(pk, NC * NC) * 0; r1 <- rowsum(p, pk); sp[as.integer(rownames(r1))] <- r1[, 1]
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
  if (!SURVEY %in% c("GSS", "ATP")) { g$survey <- SURVEY; g$blau <- CFG$blau; g$cell <- CFG$cell }   # read by 05's engine
  list(g = g, theta = setNames(th, c(paste0("delta_", GROUPS), paste0("gamma_", GROUPS[CAL]),
                                     if (TWO) c(paste0("eta_", GROUPS2[-1]), paste0("kappa_", GROUPS2[CAL2])))), iters = it,
       max_rel_gap = max(abs(gap / Tt[live])), capped = capped)
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
calib <- list(survey = SURVEY, targets = tg, targets2 = if (TWO) tg2, weighted = !is.null(PROB), coef_rest = COEF, closure = CLOSURE,
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
if (TWO) {
  grp2 <- rbindlist(lapply(g1s, by_group, GROUP = GROUP2, GROUPS = GROUPS2))[
    , .(n = mean(n), rel_deg = mean(rel_deg, na.rm = TRUE), same = mean(same, na.rm = TRUE), within_density = mean(within_density, na.rm = TRUE)), by = ego]
  grp2 <- merge(tg2[, .(ego, target_rel_deg = rel_deg, target_same = same, same_lo, same_hi)], grp2, by = "ego")[match(GROUPS2, ego)][
    , same_targeted := ego %in% GROUPS2[CAL2]]
  fwrite(grp2, file.path(OUT, paste0(SURVEY, "_by_", GROUP2, ".csv")))
  cat("\nBY ", toupper(GROUP2), " OF EGO, x1 vs targets:\n", sep = ""); print(grp2, digits = 3)
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
