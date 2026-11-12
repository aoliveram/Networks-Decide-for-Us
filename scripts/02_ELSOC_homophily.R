# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMOPHILY AMONG ELSOC CONFIDANTS (r13, wave 2017): the inputs of the ELSOC networks
# of 02_network_calibrated.R. 02_ELSOC_homophily.R
#
# Data: the ELSOC respondents of wave 2 (2017) named 1-5 confidants (none named 0; mean
# 3.3), with the confidant's sex, age, education, religion, perceived ideology and relation.
# No alter-alter ties. Recodes in 01_survey_data.R (ideology: Left 0-3, Neither 4-6/11/12,
# Right 7-10). Only 2017, the wave before the outburst of October 2019 (decided with the
# user, 2026-11-11: the 2019 confidants carry the contingency; the shift 2017 -> 2019 is a
# separate question). Every estimate is weighted by the ego's ponderador02 (mean 1); the
# unweighted main model is kept as a check.
#
# 1. Marginals of egos and confidants (education, religion, ideology, sex).
# 2. By group of ego (education 5, ideology 3, religion 5, sex 2):
#      rel_deg   weighted mean number of confidants of the group / overall weighted mean
#      same      share of the group's confidants in the same group
#      pop       weighted share of the group among egos;  same / pop = homophily relative to chance
#    95% intervals from a bootstrap of respondents.
# 3. Case-control logit (Smith, McPherson & Smith-Lovin 2014): cases = ego-confidant pairs,
#    controls = random pairs of respondents (20 per ego, the second drawn with probability
#    proportional to its weight, no self-pairs); every pair weighted by the ego's weight:
#      tie ~ same_sex + same_relig + same_ideol + |age difference| + |education difference|
#    education difference in approximate years of the 5 categories (6, 10, 12, 14, 16), so
#    the slope is per year as in Smith et al. All ties (main: the GSS 2004 targets also use
#    all confidants) and non-kin ties (spouse, child and other relatives dropped from the
#    cases). Standard errors and 95% intervals from a bootstrap of respondents that rebuilds
#    cases and controls in every iteration (B = 1000), as in Smith et al.
# 4. Homophily by group, in log-odds: the same model with one same-group dummy per
#    education group (and, separately, per ideology group).
# 5. Comparison with Smith et al. (2014), Table 3, 2004 effects (model 2 all ties; model 4
#    non-kin), on the dimensions both have (sex, religion, age, education).
#
# The coefficients of the main model (all ties, weighted) for sex, religion, age and education
# are copied into 02_network_calibrated.R (COEF of ELSOC); education and ideology are
# calibrated there by group (activity and homophily) from the targets of section 2, so the
# common same-ideology coefficient is reported here but not used.
#
# Outputs: output/02_ELSOC_homophily/{marginals,group_targets,case_control,case_control_by_group,
#          smith_comparison}.csv;  plots/02_ELSOC_homophily/elsoc_homophily.pdf
# Usage (repo root): Rscript scripts/02_ELSOC_homophily.R   (~1 min on 8 cores)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(parallel); library(ggplot2); library(patchwork) }))
source("scripts/01_survey_data.R")
ARGS <- commandArgs(trailingOnly = TRUE)
B <- as.integer(sub(".*=", "", grep("^--B=", ARGS, value = TRUE)[1])); if (is.na(B)) B <- 1000L
M_CONTROLS <- 20L; N_CORES <- 8L
OUT <- "output/02_ELSOC_homophily"; PLOTS <- "plots/02_ELSOC_homophily"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)

eg <- elsoc_egos(2); al <- elsoc_alters(2)
a <- merge(al, eg[, .(idencuesta, ola, w, e_sex = sex, e_age = age, e_educ = educ, e_educ_num = educ_num, e_relig = relig, e_ideol = ideol)],
           by = c("idencuesta", "ola"))
stamp("egos 2017: ", nrow(eg), "; confidants: ", nrow(a), "; mean number named ", round(mean(eg$nredes), 2),
      " (weighted ", round(weighted.mean(eg$nredes, eg$w), 2), "); Kish effective n of the weights ", round(sum(eg$w)^2 / sum(eg$w^2)))

# ---- 1. marginals ------------------------------------------------------------------------------------
LV <- list(educ = EDU5, relig = REL5_CL, ideol = IDE3, sex = c("Male", "Female"))
wshare <- function(x, w, lv) { k <- !is.na(x); v <- tapply(w[k], factor(x[k], lv), sum); v[is.na(v)] <- 0; as.numeric(v / sum(v)) }
marg <- rbindlist(lapply(names(LV), function(v)
  data.table(attribute = v, level = LV[[v]], egos = wshare(eg[[v]], eg$w, LV[[v]]), confidants = wshare(a[[v]], a$w, LV[[v]]),
             egos_unweighted = wshare(eg[[v]], rep(1, nrow(eg)), LV[[v]]), confidants_missing = mean(is.na(a[[v]])))))
fwrite(marg, file.path(OUT, "marginals.csv"))
cat("\nMARGINALS 2017 (weighted shares among egos and among their confidants):\n"); print(marg, digits = 3)

# ---- 2. targets by group of ego ----------------------------------------------------------------------
setkey(eg, idencuesta); setkey(a, idencuesta)
grp_stats <- function(ids, v) {
  Eb <- eg[.(ids), nomatch = 0L][!is.na(get(v))]
  Ab <- a[.(ids), nomatch = 0L, allow.cartesian = TRUE][!is.na(get(v)) & !is.na(get(paste0("e_", v)))]
  d <- Eb[, .(n = .N, wn = sum(w), deg = weighted.mean(nredes, w)), by = .(group = get(v))][
    , rel_deg := deg / weighted.mean(deg, wn)][, pop := wn / sum(wn)]
  s <- Ab[, .(alters = .N, same = weighted.mean(get(v) == get(paste0("e_", v)), w)), by = .(group = get(paste0("e_", v)))]
  merge(d, s, by = "group")[, ratio := same / pop][]
}
ids0 <- unique(eg$idencuesta)
tg <- rbindlist(lapply(names(LV), function(v) {
  pt <- grp_stats(ids0, v)
  set.seed(2017)
  bb <- rbindlist(lapply(seq_len(min(B, 500L)), function(b) grp_stats(sample(ids0, replace = TRUE), v)))
  ci <- bb[, .(rel_lo = quantile(rel_deg, .025), rel_hi = quantile(rel_deg, .975), same_lo = quantile(same, .025),
               same_hi = quantile(same, .975), ratio_lo = quantile(ratio, .025), ratio_hi = quantile(ratio, .975)), by = group]
  merge(pt, ci, by = "group")[match(LV[[v]], group)][, attribute := v]
}))
setcolorder(tg, c("attribute", "group"))
fwrite(tg, file.path(OUT, "group_targets.csv"))
cat("\nBY GROUP OF EGO (2017, weighted): relative number of confidants, same-group share, share among egos, same / pop (95% CI):\n")
print(tg[, .(attribute, group, n, rel_deg = sprintf("%.3f [%.3f, %.3f]", rel_deg, rel_lo, rel_hi),
             same = sprintf("%.3f [%.3f, %.3f]", same, same_lo, same_hi), pop = round(pop, 3),
             ratio = sprintf("%.2f [%.2f, %.2f]", ratio, ratio_lo, ratio_hi))], nrows = 200)

# ---- 3-4. case-control -------------------------------------------------------------------------------
E <- eg[complete.cases(eg[, .(sex, age, educ, relig, ideol)])][, eid := .I]
cases <- merge(a[, .(idencuesta, ola, kin, a_sex = sex, a_age = age, a_educ = educ, a_relig = relig, a_ideol = ideol)],
               E[, .(idencuesta, ola, eid, w, e_sex = sex, e_age = age, e_educ = educ, e_relig = relig, e_ideol = ideol)],
               by = c("idencuesta", "ola"))
cases <- cases[complete.cases(cases)]; setkey(cases, eid)
stamp("case-control: ", nrow(E), " egos with complete attributes; ", nrow(cases), " cases with complete confidant attributes (",
      sum(!cases$kin), " non-kin)")
dyads <- function(d) {
  out <- d[, .(tie, w, same_sex = as.numeric(a_sex == e_sex), same_relig = as.numeric(a_relig == e_relig),
               same_ideol = as.numeric(a_ideol == e_ideol), age_diff = abs(a_age - e_age),
               educ_diff = abs(EDU_YEARS[a_educ] - EDU_YEARS[e_educ]))]
  for (g in EDU5) set(out, j = paste0("same_educ_", g), value = as.numeric(d$a_educ == g & d$e_educ == g))
  for (g in IDE3) set(out, j = paste0("same_ideol_", g), value = as.numeric(d$a_ideol == g & d$e_ideol == g))
  out
}
controls_of <- function(e, weighted) {
  n <- nrow(e); i <- rep(seq_len(n), each = M_CONTROLS)
  j <- sample.int(n, length(i), replace = TRUE, prob = if (weighted) e$w else NULL)
  keep <- e$eid[i] != e$eid[j]; i <- i[keep]; j <- j[keep]
  data.table(tie = 0, w = e$w[i], e_sex = e$sex[i], e_age = e$age[i], e_educ = e$educ[i], e_relig = e$relig[i], e_ideol = e$ideol[i],
             a_sex = e$sex[j], a_age = e$age[j], a_educ = e$educ[j], a_relig = e$relig[j], a_ideol = e$ideol[j]) }
MODELS <- list(
  main = "same_sex + same_relig + same_ideol + age_diff + educ_diff",
  educ_by_group = paste("same_sex + same_relig + same_ideol + age_diff + educ_diff +", paste0("same_educ_", EDU5, collapse = " + ")),
  ideol_by_group = paste("same_sex + same_relig + age_diff + educ_diff +", paste0("same_ideol_", IDE3, collapse = " + ")))
fit <- function(ids, nonkin, model, weighted) {
  Eb <- E[ids]; cs <- cases[.(ids), nomatch = 0L, allow.cartesian = TRUE]
  if (nonkin) cs <- cs[kin == FALSE]
  d <- rbind(dyads(cs[, tie := 1]), dyads(controls_of(Eb, weighted)), fill = TRUE)
  X <- model.matrix(as.formula(paste("tie ~", MODELS[[model]])), d)
  setNames(glm.fit(X, d$tie, weights = if (weighted) d$w else NULL, family = quasibinomial())$coefficients, colnames(X))
}
SPECS <- rbind(data.table(ties = c("all", "nonkin", "all"), weighted = c(TRUE, TRUE, FALSE), model = "main"),
               data.table(ties = "all", weighted = TRUE, model = c("educ_by_group", "ideol_by_group")))
cc <- rbindlist(lapply(seq_len(nrow(SPECS)), function(s) {
  sp <- SPECS[s]; ids <- E$eid
  stamp("case-control: ", sp$ties, " ties, ", if (sp$weighted) "weighted" else "unweighted", ", model ", sp$model, " (B = ", B, ")")
  set.seed(7); pt <- fit(ids, sp$ties == "nonkin", sp$model, sp$weighted)
  bb <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(10000 + b)
    fit(sample(ids, replace = TRUE), sp$ties == "nonkin", sp$model, sp$weighted) }, mc.cores = N_CORES))
  data.table(ties = sp$ties, weighted = sp$weighted, model = sp$model, term = names(pt), coef = pt, se = apply(bb, 2, sd),
             lo = apply(bb, 2, quantile, .025), hi = apply(bb, 2, quantile, .975))
}))
cc <- cc[!grepl("Intercept", term)]
fwrite(cc[model == "main"], file.path(OUT, "case_control.csv"))
fwrite(cc[model != "main"], file.path(OUT, "case_control_by_group.csv"))
cat("\nCASE-CONTROL 2017 (log-odds; same_* = same category, *_diff per year; bootstrap SE and 95% CI, B = ", B, "):\n", sep = "")
print(cc[model == "main", .(ties, weighted, term, coef = sprintf("%.3f [%.3f, %.3f]", coef, lo, hi), se = round(se, 3))])
cat("\nHOMOPHILY BY GROUP (all ties, weighted; same-group dummies on top of the other terms):\n")
print(cc[model != "main", .(model, term, coef = sprintf("%.3f [%.3f, %.3f]", coef, lo, hi))])

# ---- 5. Smith et al. (2014), Table 3: 2004 effects = 1985 coefficient + year interaction ------------------
smith <- data.table(term = c("same_sex", "same_relig", "age_diff", "educ_diff"),
                    smith_all = c(0.373 - 0.144, 1.270 + 0.215, -0.047 - 0.005, -0.157 - 0.044),       # model 2 (nodematch scale)
                    smith_nonkin = c(1.101 - 0.035, 0.818 + 0.239, -0.107 - 0.019, -0.199 - 0.026))    # model 4
cmp <- merge(smith, dcast(cc[model == "main" & weighted == TRUE], term ~ ties, value.var = c("coef", "lo", "hi")), by = "term", all.y = TRUE)
fwrite(cmp, file.path(OUT, "smith_comparison.csv"))
cat("\nELSOC 2017 (weighted) vs Smith et al. (2014), 2004 effects (nodematch scale: same_* > 0 = homophily; *_diff < 0):\n")
print(cmp[, .(term, ELSOC_all = sprintf("%.3f [%.3f, %.3f]", coef_all, lo_all, hi_all), Smith_all = smith_all,
              ELSOC_nonkin = sprintf("%.3f [%.3f, %.3f]", coef_nonkin, lo_nonkin, hi_nonkin), Smith_nonkin = smith_nonkin)])
cat("\nCOEF for 02_network_calibrated.R (ELSOC, all ties, 2017, weighted; same_ideol not used, ideology is calibrated):\n")
m <- cc[model == "main" & ties == "all" & weighted == TRUE]
cat(sprintf("  c(sex = %.3f, relig = %.3f, age = %.3f, educ = %.3f)   [same_ideol %.3f]\n",
            m[term == "same_sex", coef], m[term == "same_relig", coef], m[term == "age_diff", coef], m[term == "educ_diff", coef],
            m[term == "same_ideol", coef]))

# ---- figure --------------------------------------------------------------------------------------------
LAB <- c(same_sex = "Same sex", same_relig = "Same religion", same_ideol = "Same ideology", age_diff = "Age difference (year)",
         educ_diff = "Education difference (year)")
pd <- rbind(cc[model == "main" & weighted == TRUE, .(term, ties, source = "ELSOC 2017 (weighted)", coef, lo, hi)],
            melt(smith, id.vars = "term", variable.name = "ties", value.name = "coef")[
              , .(term, ties = sub("smith_", "", ties), source = "Smith et al. 2014 (GSS 2004)", coef, lo = NA_real_, hi = NA_real_)])
pd[, lab := factor(LAB[term], rev(LAB))][, ties := factor(ties, c("all", "nonkin"), c("All confidants", "Non-kin confidants"))]
p1 <- ggplot(pd, aes(coef, lab, color = source)) + geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), na.rm = TRUE) +
  geom_point(position = position_dodge(width = 0.5), size = 1.8) + facet_wrap(~ ties, scales = "free_x") +
  scale_color_manual(values = c("#C0392B", "grey30"), name = NULL) +
  labs(title = "Case-control log-odds of a confidant tie (95% bootstrap CI)", x = "Log-odds", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(size = 9, face = "bold"))
gd <- tg[attribute %in% c("educ", "ideol")][, group := factor(group, c(EDU5, IDE3))]
fl <- facet_wrap(~ attribute, scales = "free_x", labeller = as_labeller(c(educ = "Education", ideol = "Ideology")))
p2 <- ggplot(gd, aes(group, rel_deg)) + geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(ymin = rel_lo, ymax = rel_hi), size = 0.25, color = "#2980B9") + fl +
  labs(title = "Relative number of confidants by group of ego (2017, weighted)", x = NULL, y = "Relative degree") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
p3 <- ggplot(gd, aes(group, ratio)) + geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(ymin = ratio_lo, ymax = ratio_hi), size = 0.25, color = "#2980B9") + fl +
  labs(title = "Same-group share of confidants over the group's share among egos", x = NULL, y = "Same / population") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(PLOTS, "elsoc_homophily.pdf"), p1 / p2 / p3, width = 8, height = 9.5, device = cairo_pdf)
stamp("saved ", OUT, "/*.csv and ", file.path(PLOTS, "elsoc_homophily.pdf"))
