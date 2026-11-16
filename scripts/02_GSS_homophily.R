# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMOPHILY AMONG GSS 2004 CONFIDANTS, WEIGHTED: the common coefficients of the GSS, ATP
# and GSS 2014 networks of 02_network_calibrated.R. 02_GSS_homophily.R
#
# Until 2026-11-11 the networks took Smith, McPherson & Smith-Lovin (2014, Table 3,
# 2004 effects) for the terms that are not calibrated by group. Race, religion and sex are
# now calibrated by group (02), so only age and education keep a common coefficient, plus
# the common same-race and same-religion terms of the small groups left uncalibrated
# ('Other' race, Jewish). They are estimated here with the design of Smith et al., on GSS
# 2004 only and weighted (decided with the user: survey weights from the start; WTSSNR, the
# weight of GSS 2014 too):
#   cases    ego-confidant pairs (up to 5 per ego), every pair weighted by the ego's weight
#   controls 20 random respondents per ego, drawn with probability proportional to their
#            weight, no self-pairs, weighted by the ego's weight
#   tie ~ same_race + same_relig + same_sex + |age difference| + |education difference|
#   race: White, Black, Hispanic, Asian, Other (as in 02); religion: Protestant, Catholic,
#   Jewish, None, Other; education in years (ego: years; confidant: categories 0-7 mapped
#   to 3.5, 8, 11, 12, 13.5, 14, 16, 18, as in the replication of Smith et al.)
# All ties (main) and non-kin ties (spouse, parent, sibling, child and other family dropped
# from the cases). Bootstrap of respondents that rebuilds cases and controls, B = 1000.
# Unweighted, the all-ties model should land near Smith et al.'s 2004 effects.
#
# Outputs: output/02_GSS_homophily/case_control.csv
# Usage (repo root): Rscript scripts/02_GSS_homophily.R   (~1 min on 8 cores)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(haven); library(parallel) }))
source("scripts/01_survey_data.R")
ARGS <- commandArgs(trailingOnly = TRUE)
B <- as.integer(sub(".*=", "", grep("^--B=", ARGS, value = TRUE)[1])); if (is.na(B)) B <- 1000L
M_CONTROLS <- 20L; N_CORES <- 8L
OUT <- "output/02_GSS_homophily"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
num <- function(x) as.numeric(zap_labels(x))

g4 <- gss_load(2004)
RELIG <- c("Protestant", "Catholic", "Jewish", "None", "OtherRelig")
ALT_EDUC_YEARS <- c(3.5, 8, 11, 12, 13.5, 14, 16, 18)          # categories 0-7
E <- data.table(gss_id = num(g4$id), numgiven = num(g4$numgiven), race = race_gss(g4),
                relig = RELIG[pmin(num(g4$relig), 5)], sex = c("Male", "Female")[num(g4$sex)],
                age = ifelse(num(g4$age) %in% c(98, 99), NA, num(g4$age)),
                educ = ifelse(num(g4$educ) %in% c(98, 99), NA, num(g4$educ)), w = num(g4$wtssnr))
E <- E[!is.na(numgiven) & numgiven <= 6 & complete.cases(E)][, w := w / mean(w)][, eid := .I]
A <- rbindlist(lapply(1:5, function(k) {
  kin <- Reduce(`|`, lapply(c("spouse", "parent", "sibling", "child", "othfam"), function(v) num(g4[[paste0(v, k)]]) %in% 1))
  data.table(gss_id = num(g4$id), a_race = c("Asian", "Black", "Hispanic", "White", "Other")[num(g4[[paste0("race", k)]])],
             a_relig = RELIG[num(g4[[paste0("relig", k)]])], a_sex = c("Male", "Female")[num(g4[[paste0("sex", k)]])],
             a_age = num(g4[[paste0("age", k)]]), a_educ = ALT_EDUC_YEARS[num(g4[[paste0("educ", k)]]) + 1], kin = kin) }))
A[a_age > 98, a_age := NA]
cases <- merge(A, E[, .(gss_id, eid, w, e_race = race, e_relig = relig, e_sex = sex, e_age = age, e_educ = educ)], by = "gss_id")
cases <- cases[complete.cases(cases)]; setkey(cases, eid)
stamp("GSS 2004 egos with complete attributes: ", nrow(E), " (Kish effective n ", round(sum(E$w)^2 / sum(E$w^2)), "); cases ",
      nrow(cases), " (", sum(!cases$kin), " non-kin)")

dyads <- function(d) d[, .(tie, w, same_race = as.numeric(a_race == e_race), same_relig = as.numeric(a_relig == e_relig),
                           same_sex = as.numeric(a_sex == e_sex), age_diff = abs(a_age - e_age), educ_diff = abs(a_educ - e_educ))]
controls_of <- function(e, weighted) {
  n <- nrow(e); i <- rep(seq_len(n), each = M_CONTROLS)
  j <- sample.int(n, length(i), replace = TRUE, prob = if (weighted) e$w else NULL)
  keep <- e$eid[i] != e$eid[j]; i <- i[keep]; j <- j[keep]
  data.table(tie = 0, w = e$w[i], e_race = e$race[i], e_relig = e$relig[i], e_sex = e$sex[i], e_age = e$age[i], e_educ = e$educ[i],
             a_race = e$race[j], a_relig = e$relig[j], a_sex = e$sex[j], a_age = e$age[j], a_educ = e$educ[j]) }
FORM <- tie ~ same_race + same_relig + same_sex + age_diff + educ_diff
fit <- function(ids, nonkin, weighted) {
  cs <- cases[.(ids), nomatch = 0L, allow.cartesian = TRUE]; if (nonkin) cs <- cs[kin == FALSE]
  d <- rbind(dyads(cs[, tie := 1]), dyads(controls_of(E[ids], weighted)))
  X <- model.matrix(FORM, d)
  setNames(glm.fit(X, d$tie, weights = if (weighted) d$w else NULL, family = quasibinomial())$coefficients, colnames(X))
}
SPECS <- data.table(ties = c("all", "nonkin", "all"), weighted = c(TRUE, TRUE, FALSE))
cc <- rbindlist(lapply(seq_len(nrow(SPECS)), function(s) {
  sp <- SPECS[s]; ids <- E$eid
  stamp("case-control: ", sp$ties, " ties, ", if (sp$weighted) "weighted" else "unweighted", " (B = ", B, ")")
  set.seed(7); pt <- fit(ids, sp$ties == "nonkin", sp$weighted)
  bb <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(20000 + b)
    fit(sample(ids, replace = TRUE), sp$ties == "nonkin", sp$weighted) }, mc.cores = N_CORES))
  data.table(ties = sp$ties, weighted = sp$weighted, term = names(pt), coef = pt, se = apply(bb, 2, sd),
             lo = apply(bb, 2, quantile, .025), hi = apply(bb, 2, quantile, .975))
}))[term != "(Intercept)"]
smith <- c(same_race = 1.959 - 0.264, same_relig = 1.270 + 0.215, same_sex = 0.373 - 0.144, age_diff = -0.047 - 0.005, educ_diff = -0.157 - 0.044)
cc[ties == "all", smith_2004 := smith[term]]
fwrite(cc, file.path(OUT, "case_control.csv"))
cat("\nGSS 2004 CASE-CONTROL (log-odds; bootstrap 95% CI, B = ", B, "); Smith et al. 2014, Table 3 model 2, 2004 effects:\n", sep = "")
print(cc[, .(ties, weighted, term, coef = sprintf("%.3f [%.3f, %.3f]", coef, lo, hi), smith_2004)])
m <- cc[ties == "all" & weighted == TRUE]
cat(sprintf("\nCOEF for 02_network_calibrated.R (GSS family, all ties, weighted):\n  c(race = %.3f, relig = %.3f, sex = %.3f, age = %.4f, educ = %.4f)\n",
            m[term == "same_race", coef], m[term == "same_relig", coef], m[term == "same_sex", coef], m[term == "age_diff", coef], m[term == "educ_diff", coef]))
