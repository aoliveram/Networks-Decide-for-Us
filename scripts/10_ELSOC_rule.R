# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE RULE IN THE DATA: Gamma + lambda E >= q estimated in the same people (ELSOC panel).
# 10_ELSOC_rule.R  (Prompt 6.2, 2026-11-16; from playground/PR6.1-elsoc-mechanism/scripts/
# 02_rule_estimation.R and 02b_structural_exposure.R, specifications S1-S7 unchanged)
#
# ELSOC is the only one of the four surveys that measures q, E and y in the same people: the
# willingness items c09 (the q of the networks), how often family and friends took part in the
# movement the respondent values most (c23, c24: E) and whether the respondent took part (c22: y).
#
# Structural model: y*_i = Gamma_i + lambda E_i - q_i + e_i, e logistic with scale sigma; y = 1[y* >= 0]
#   => logit P(y = 1) = b0 + bG Gamma_i + bE E_i + bq q_i;  the rule: bG = -bq = 1/sigma, bE = lambda/sigma
#   lambda_hat = bE / (-bq)    exposure in units of q (the engine's mur_score)
#   wG_hat     = bG / (-bq)    weight of perceived efficacy in units of q (the rule: 1)
#   sigma_hat  = 1 / (-bq)     noise of the rule in units of q (the engine has none)
#   Gamma_t    = (b0 + b_year_t + bG mean(Gamma)) / (-bq): the q at which a person of average efficacy
#                and E = 0 takes part with probability 1/2, by year
# Measurement: y = c22 >= 2 (took part in the movement valued most, last 12 months); Gamma = (c21_05 - 1)/4
#   ("the actions and protests of this movement can bring social change"); E = mean of (c23 - 1)/4 and
#   (c24 - 1)/4 (family, friends; same movement); q = q_c09 (01_survey_data.R: elsoc_panel()). Only
#   respondents who name a movement answer c21-c24 (selection; IPW sensitivity below).
# Identification ladder (none causal):
#   S1  pooled 2016, 2018, 2022, q of the same wave (association; perceived E: UPPER BOUND)
#   S2  S1 + attributes, movement, identification
#   S3  lagged q (2016/2018/2022 -> y of 2017/2019/2023); S3b 2018 and 2019 with q of 2018 (the jump of
#       Gamma across the October 2019 outburst, E held fixed)
#   S4  lagged E and lagged y, consecutive waves, same movement (removes contemporaneous projection)
#   S5  person fixed effects (linear probability)
#   S6  family and friends separately
#   S7  STRUCTURAL exposure: for every 2017 confidant (r13), the probability that a person with its
#       attributes takes part (weighted logit on the 2017 cross-section), averaged over the ego's
#       confidants; free of projection, not of homophily (controls: the ego's own attributes and own
#       predicted probability). It is what the imputed networks feed the engine. q of 2016.
# Also: robustness of S1 (measurement, outcome, selection), functional forms (additive, complements,
# pure product; AIC, BIC, 10-fold CV by respondent), the E = 0 test, the slope of exposure near and far
# from the threshold and the shape of exposure (GAMs), and which q fits the rule (c09, GSS format 4 and 3
# items).
# Correction to the Prompt 6.1 numbers (playground 02_rule_estimation.R): there, an in-place data.table
# assignment (d1 <- SPECS$S1$data; d1[, y := own >= 3]) changed the S1 sample, so the S1 forms, the E = 0
# test, the threshold bins, the GAM and the IPW row used y = c22 >= 3. Here every sample is copied. The
# ladder, S7, the other robustness rows, the S3 forms and the q comparison are identical to 6.1.
# Weights: ponderador02 of each wave (w, mean 1 per wave). Intervals: bootstrap over respondents (B).
#
# Outputs:
#   output/10_ELSOC_rule/rule_ladder.csv, rule_struct_exposure.csv, rule_robustness.csv, rule_forms.csv,
#     rule_E0_test.csv, rule_threshold_bins.csv, rule_threshold_gam.csv, rule_q_compare.csv, rule_fits.rds
#   output/10_ELSOC_rule/tab_rule_{ladder,forms,threshold,exposure_bins,E0,robustness}.tex  (English, booktabs; \input by
#     paper-NDFU/sections/rule_in_the_data.tex and docs/MAKING-IT-BOTTOM-UP.tex; no number typed by hand)
#   plots/10_ELSOC_rule/rule_lambda_ladder.pdf, rule_threshold.pdf
# Usage (repo root): Rscript scripts/10_ELSOC_rule.R [--B=500] [--cores=8] 2>&1 | tee logs/10_ELSOC_rule.log
# Runtime ~4 min with B = 500 on 8 cores.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(fixest); library(mgcv); library(MASS); library(ggplot2); library(parallel) }))
source("scripts/01_survey_data.R")
OUT <- "output/10_ELSOC_rule"; FIG <- "plots/10_ELSOC_rule"
for (p in c(OUT, FIG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
options(width = 180)
ARGS <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE); if (length(a)) sub(".*=", "", a[1]) else d }
B <- as.integer(arg("B", 500)); NC <- as.integer(arg("cores", 8))
t0 <- Sys.time()
d <- elsoc_panel()

# ---- analysis frame: respondents who named a movement -----------------------------------------------
r <- d[named == TRUE & !is.na(own), .(id, year, w, muestra, c20, educ, ideol, agegrp, sex, relig,
                                     y = as.integer(own >= 2), own, G = (eff - 1) / 4,
                                     E = ((fam - 1) / 4 + (fri - 1) / 4) / 2, Efam = (fam - 1) / 4, Efri = (fri - 1) / 4,
                                     G3 = ((c21_03 + c21_04 + c21_05) / 3 - 1) / 4, ident = ((c21_01 + c21_02) / 2 - 1) / 4,
                                     q = q_c09, q_gss4, q_gss3, march = as.integer(p2 >= 2))]
# q of the closest EARLIER measurement (c09 waves 2016, 2018, 2022) for the lagged designs
qm <- d[!is.na(q_c09), .(id, qyear = year, q_lag = q_c09, qg4_lag = q_gss4, qg3_lag = q_gss3)]
lag_of <- c(`2017` = 2016, `2019` = 2018, `2023` = 2022, `2018` = 2016, `2022` = 2018)
r[, qyear := lag_of[as.character(year)]]
r <- merge(r, qm, by = c("id", "qyear"), all.x = TRUE)
r <- merge(r, d[year == 2018 & !is.na(q_c09), .(id, q18 = q_c09, qg4_18 = q_gss4, qg3_18 = q_gss3)], by = "id", all.x = TRUE)
r[, Ebin := as.integer(E > 0)]
cat("Analysis frame: person-waves of respondents naming a movement with c22:", nrow(r), "\n")
print(r[, .(n = .N, y = round(wmean(y, w), 3), E = round(wmean(E, w), 3), G = round(wmean(G, w), 3), q = round(wmean(q, w), 3)), by = year][order(year)])

# ---- fitting helpers ----------------------------------------------------------------------------------
wlogit <- function(f, dat) suppressWarnings(glm(f, family = binomial, data = dat, weights = w))   # non-integer weights
ratios <- function(m, Gbar = NULL, qname = "q", Ename = "E", Gname = "G") {
  b <- coef(m); bq <- b[[qname]]
  out <- c(bG = if (Gname %in% names(b)) b[[Gname]] else NA, bE = if (Ename %in% names(b)) b[[Ename]] else NA, bq = bq,
           lambda = if (Ename %in% names(b)) b[[Ename]] / (-bq) else NA, wG = if (Gname %in% names(b)) b[[Gname]] / (-bq) else NA,
           sigma = 1 / (-bq))
  if (!is.null(Gbar)) {                                   # Gamma by year: b0 + b_year + bG * mean(G)
    yrs <- grep("^factor\\(year\\)", names(b), value = TRUE)
    base <- b[["(Intercept)"]] + (if (Gname %in% names(b)) b[[Gname]] * Gbar else 0)
    gt <- c(base, base + b[yrs]) / (-bq); names(gt) <- paste0("Gamma_", c(attr(Gbar, "base_year"), sub("factor\\(year\\)", "", yrs)))
    out <- c(out, gt)
  }
  out
}
boot_ids <- function(dat, fun, B, seed) {               # respondent bootstrap of a vector-valued statistic
  ids <- unique(dat$id); setkey(dat, id)
  res <- mclapply(seq_len(B), function(b) { set.seed(seed + b); s <- sample(ids, replace = TRUE)
    bd <- dat[.(s), allow.cartesian = TRUE]; tryCatch(fun(bd), error = function(e) NULL) }, mc.cores = NC)
  do.call(rbind, Filter(Negate(is.null), res))
}
summ <- function(est, bs, spec, n, nid) {
  data.table(spec = spec, stat = names(est), est = unname(est),
             lo = apply(bs, 2, quantile, .025, na.rm = TRUE)[names(est)], hi = apply(bs, 2, quantile, .975, na.rm = TRUE)[names(est)],
             n = n, n_id = nid)
}
Gbar_of <- function(dat) { g <- wmean(dat$G, dat$w); attr(g, "base_year") <- min(dat$year); g }

# ---- the identification ladder, S1-S6 -------------------------------------------------------------------
SPECS <- list(
  S1 = list(data = r[year %in% c(2016, 2018, 2022) & !is.na(G) & !is.na(E) & !is.na(q)],
            f = y ~ G + E + q + factor(year), what = "pooled 2016, 2018, 2022; q of the same wave"),
  S2 = list(data = r[year %in% c(2016, 2018, 2022) & !is.na(G) & !is.na(E) & !is.na(q) & !is.na(educ) & !is.na(ideol) & !is.na(relig) & !is.na(ident) & !is.na(agegrp)],
            f = y ~ G + E + q + factor(year) + educ + ideol + agegrp + sex + relig + factor(c20) + ident, what = "S1 + attributes, movement, identification"),
  S3 = list(data = r[year %in% c(2017, 2019, 2023) & !is.na(G) & !is.na(E) & !is.na(q_lag)][, q := q_lag],
            f = y ~ G + E + q + factor(year), what = "lagged q: q of 2016/2018/2022 -> y of 2017/2019/2023"),
  S3b = list(data = r[year %in% c(2018, 2019) & !is.na(G) & !is.na(E) & !is.na(q18)][, q := q18],
             f = y ~ G + E + q + factor(year), what = "2018 and 2019 with q of 2018 (the 18-O jump)"),
  S6 = list(data = r[year %in% c(2016, 2018, 2022) & !is.na(G) & !is.na(Efam) & !is.na(Efri) & !is.na(q)],
            f = y ~ G + Efam + Efri + q + factor(year), what = "family and friends separately"))
LAD <- list(); FITS <- list()
cat("\n==== IDENTIFICATION LADDER (lambda, wG, sigma in units of q; 95% respondent bootstrap, B = ", B, ") ====\n", sep = "")
for (s in names(SPECS)) {
  sp <- SPECS[[s]]; dat <- copy(sp$data)
  fun <- function(bd) { m <- wlogit(sp$f, bd); v <- ratios(m, Gbar_of(bd))
    if (s == "S3b") v <- c(v, Gamma_jump_2019 = v[["Gamma_2019"]] - v[["Gamma_2018"]])          # the 18-O jump, E held fixed
    if (s == "S6") { b <- coef(m); v <- c(v, bFam = b[["Efam"]], bFri = b[["Efri"]], lambda_fam = b[["Efam"]] / (-b[["q"]]),
                                          lambda_fri = b[["Efri"]] / (-b[["q"]]), lambda_both = (b[["Efam"]] + b[["Efri"]]) / 2 / (-b[["q"]])) }
    v }
  m <- wlogit(sp$f, dat); FITS[[s]] <- m
  est <- fun(dat); bs <- boot_ids(dat, fun, B, seed = 1000 * match(s, names(SPECS)))
  LAD[[s]] <- summ(est, bs, s, nrow(dat), uniqueN(dat$id))[, what := sp$what]
  t <- bs[, "bG"] + bs[, "bq"]                         # the rule's restriction bG = -bq (bootstrap of bG + bq)
  LAD[[s]] <- rbind(LAD[[s]], data.table(spec = s, stat = "bG_plus_bq", est = est[["bG"]] + est[["bq"]], lo = quantile(t, .025, na.rm = TRUE),
                                         hi = quantile(t, .975, na.rm = TRUE), n = nrow(dat), n_id = uniqueN(dat$id), what = sp$what))
  cat(sprintf("%-4s n = %5d (%4d ids)  lambda = %.2f [%.2f, %.2f]  wG = %.2f  sigma = %.2f   %s\n", s, nrow(dat), uniqueN(dat$id),
              est[["lambda"]], quantile(bs[, "lambda"], .025, na.rm = TRUE), quantile(bs[, "lambda"], .975, na.rm = TRUE), est[["wG"]], est[["sigma"]], sp$what))
  if (s == "S3b") {
    j <- LAD[[s]][stat == "Gamma_jump_2019"]
    cat(sprintf("     Gamma 2018 -> 2019 with E held fixed (the 18-O jump, units of q): %+.3f [%+.3f, %+.3f]  (Gamma 2018 %.3f, 2019 %.3f)\n",
                j$est, j$lo, j$hi, est[["Gamma_2018"]], est[["Gamma_2019"]]))
  }
}
# S4: lagged E and lagged y, consecutive waves, same movement in both
pairs <- list(c(2016, 2017), c(2017, 2018), c(2018, 2019), c(2022, 2023))
L <- rbindlist(lapply(pairs, function(p) {
  a <- r[year == p[1], .(id, c20_l = c20, E_l = E, y_l = y, G_l = G)]; b <- r[year == p[2]]
  merge(b, a, by = "id")[c20 == c20_l]
}))
L[, q_any := fifelse(!is.na(q), q, fifelse(!is.na(q_lag), q_lag, q18))]
S4d <- L[!is.na(E_l) & !is.na(G) & !is.na(q_any) & !is.na(y_l)][, q := q_any]
f4 <- y ~ E_l + y_l + G + q + factor(year)
fun4 <- function(bd) { m <- wlogit(f4, bd); b <- coef(m); c(bE_lag = b[["E_l"]], by_lag = b[["y_l"]], bG = b[["G"]], bq = b[["q"]],
                                                            lambda_lag = b[["E_l"]] / (-b[["q"]]), sigma = 1 / (-b[["q"]])) }
est4 <- fun4(S4d); bs4 <- boot_ids(S4d, fun4, B, 9000)
LAD$S4 <- summ(est4, bs4, "S4", nrow(S4d), uniqueN(S4d$id))[, what := "lagged E and lagged y, same movement, consecutive waves"]
cat(sprintf("S4   n = %5d (%4d ids)  lambda(lagged E | lagged y) = %.2f [%.2f, %.2f]\n", nrow(S4d), uniqueN(S4d$id), est4[["lambda_lag"]],
            quantile(bs4[, "lambda_lag"], .025), quantile(bs4[, "lambda_lag"], .975)))
m4b <- wlogit(y ~ E_l + G + q + factor(year), S4d); cat(sprintf("     without lagged y: lambda = %.2f\n", coef(m4b)[["E_l"]] / -coef(m4b)[["q"]]))
# S5: person fixed effects, linear probability (fixest), all waves with c22-c24
F5 <- r[!is.na(E) & !is.na(G)]
fe_pool <- feols(y ~ E + G | year, data = F5, weights = ~w, cluster = ~id)
fe_id   <- feols(y ~ E + G | id + year, data = F5, weights = ~w, cluster = ~id, notes = FALSE)
F5q <- r[year %in% c(2016, 2018, 2022) & !is.na(E) & !is.na(G) & !is.na(q)]
lpm_pool_q <- feols(y ~ E + G + q | year, data = F5q, weights = ~w, cluster = ~id)
lpm_id_q   <- feols(y ~ E + G + q | id + year, data = F5q, weights = ~w, cluster = ~id, notes = FALSE)
ratio_fe <- coef(fe_id)[["E"]] / coef(fe_pool)[["E"]]
S5 <- data.table(spec = "S5", stat = c("bE_lpm_pooled", "bE_lpm_fe", "ratio_fe_pooled", "bE_lpm_pooled_q", "bq_lpm_pooled_q", "bE_lpm_fe_q", "bq_lpm_fe_q",
                                       "lambda_lpm_pooled_q", "lambda_lpm_fe_q", "lambda_S1_x_ratio"),
                 est = c(coef(fe_pool)[["E"]], coef(fe_id)[["E"]], ratio_fe, coef(lpm_pool_q)[["E"]], coef(lpm_pool_q)[["q"]], coef(lpm_id_q)[["E"]], coef(lpm_id_q)[["q"]],
                         coef(lpm_pool_q)[["E"]] / -coef(lpm_pool_q)[["q"]], coef(lpm_id_q)[["E"]] / -coef(lpm_id_q)[["q"]],
                         LAD$S1[stat == "lambda"]$est * ratio_fe),
                 lo = NA_real_, hi = NA_real_, n = c(rep(nobs(fe_pool), 3), rep(nobs(lpm_pool_q), 6), NA), n_id = uniqueN(F5$id),
                 what = "person fixed effects, linear probability")
se <- function(m, v) sqrt(vcov(m)[v, v])
S5[stat == "bE_lpm_pooled", `:=`(lo = est - 1.96 * se(fe_pool, "E"), hi = est + 1.96 * se(fe_pool, "E"))]
S5[stat == "bE_lpm_fe", `:=`(lo = est - 1.96 * se(fe_id, "E"), hi = est + 1.96 * se(fe_id, "E"))]
S5[stat == "bq_lpm_fe_q", `:=`(lo = est - 1.96 * se(lpm_id_q, "q"), hi = est + 1.96 * se(lpm_id_q, "q"))]
LAD$S5 <- S5
cat(sprintf("S5   LPM bE pooled %.3f, with person FE %.3f (ratio %.2f); with q: pooled lambda %.2f, FE lambda %.2f (bq FE %.3f)\n",
            coef(fe_pool)[["E"]], coef(fe_id)[["E"]], ratio_fe, S5[stat == "lambda_lpm_pooled_q"]$est, S5[stat == "lambda_lpm_fe_q"]$est, coef(lpm_id_q)[["q"]]))
LADDER <- rbindlist(LAD, fill = TRUE)
fwrite(LADDER, file.path(OUT, "rule_ladder.csv"))

# ---- S7: structural exposure from the 2017 confidants (r13) -----------------------------------------------
ATTR <- c("sex", "agegrp", "educ", "relig", "ideol")
p17 <- d[year == 2017][complete.cases(d[year == 2017, ..ATTR])][, ypop := as.integer(named & !is.na(own) & own >= 2)]
mp <- suppressWarnings(glm(ypop ~ sex + agegrp + educ + relig + ideol, family = binomial, data = p17, weights = w))
cat("\nS7 population model of participation (2017, weighted): n =", nrow(p17), "; mean =", round(wmean(p17$ypop, p17$w), 3), "\n")
al <- elsoc_alters(2)[ola == 2][, agegrp := age_group(age)]
al <- al[complete.cases(al[, ..ATTR])]
al[, p := predict(mp, al, type = "response")]
ES <- al[, .(E_struct = mean(p), n_alt = .N, kin_share = mean(kin)), by = .(id = idencuesta)]
r7 <- d[year == 2017 & named == TRUE & !is.na(own), .(id, w, y = as.integer(own >= 2), E = ((fam - 1) / 4 + (fri - 1) / 4) / 2,
                                                        Efri = (fri - 1) / 4, G = (eff - 1) / 4, sex, agegrp, educ, relig, ideol)]
r7 <- merge(r7, ES, by = "id")
r7 <- merge(r7, d[year == 2016 & !is.na(q_c09), .(id, q = q_c09)], by = "id")
r7 <- r7[complete.cases(r7[, c("y", "E", "G", "q", ATTR), with = FALSE])]
r7[, p_own := predict(mp, r7, type = "response")]
cat("Egos of 2017 naming a movement with confidants, c22-c24, efficacy and q of 2016:", nrow(r7), "\n")
cat(sprintf("E_struct: mean %.3f, sd %.3f, range %.3f-%.3f; corr with perceived E %.3f, with own predicted p %.3f\n",
            wmean(r7$E_struct, r7$w), sd(r7$E_struct), min(r7$E_struct), max(r7$E_struct), cor(r7$E_struct, r7$E), cor(r7$E_struct, r7$p_own)))
FORM7 <- list(P2 = y ~ E_struct + G + q, P3 = y ~ E_struct + G + q + sex + agegrp + educ + relig + ideol,
              P4 = y ~ E_struct + p_own + G + q, P5 = y ~ E + E_struct + G + q + sex + agegrp + educ + relig + ideol)
stat7 <- function(dat) {
  l1 <- lm(E ~ E_struct, data = dat, weights = w); l1b <- lm(E ~ E_struct + sex + agegrp + educ + relig + ideol, data = dat, weights = w)
  out <- c(P1_slope = coef(l1)[["E_struct"]], P1_slope_ctrl = coef(l1b)[["E_struct"]])
  for (nm in names(FORM7)) { b <- coef(wlogit(FORM7[[nm]], dat))
    out[paste0(nm, "_bEs")] <- b[["E_struct"]]; out[paste0(nm, "_bq")] <- b[["q"]]; out[paste0(nm, "_lambda_struct")] <- b[["E_struct"]] / -b[["q"]]
    if ("E" %in% names(b)) out[paste0(nm, "_lambda_perceived")] <- b[["E"]] / -b[["q"]] }
  out }
est7 <- stat7(r7)
ids7 <- unique(r7$id); setkey(r7, id)
bs7 <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(5000 + b); s <- sample(ids7, replace = TRUE)
  tryCatch(stat7(r7[.(s), allow.cartesian = TRUE]), error = function(e) NULL) }, mc.cores = NC))
S7 <- data.table(stat = names(est7), est = est7, lo = apply(bs7, 2, quantile, .025)[names(est7)], hi = apply(bs7, 2, quantile, .975)[names(est7)], n = nrow(r7))
fwrite(S7, file.path(OUT, "rule_struct_exposure.csv"))
cat("S7 STRUCTURAL EXPOSURE (P1: perceived E on E_struct; P2 no controls; P3 + own attributes; P4 + own predicted p; P5 + perceived E):\n")
print(S7, digits = 3)

# ---- robustness of S1: measurement, outcome, selection ---------------------------------------------------
d1 <- copy(SPECS$S1$data)
rob_fit <- function(dat, f, lab, qn = "q", En = "E", Gn = "G") { m <- wlogit(f, dat); b <- coef(m)
  data.table(variant = lab, n = nrow(dat), bG = b[[Gn]], bE = b[[En]], bq = b[[qn]], lambda = b[[En]] / -b[[qn]], wG = b[[Gn]] / -b[[qn]], sigma = 1 / -b[[qn]]) }
ROB <- list(rob_fit(d1, y ~ G + E + q + factor(year), "main: y = c22 >= 2, E family+friends, Gamma c21_05"),
            rob_fit(d1[!is.na(G3)], y ~ G3 + E + q + factor(year), "Gamma = 3-item efficacy (c21_03-05)", Gn = "G3"),
            rob_fit(d1, y ~ G + Efri + q + factor(year), "E = friends only", En = "Efri"),
            rob_fit(d1, y ~ G + Efam + q + factor(year), "E = family only", En = "Efam"),
            rob_fit(d1, y ~ G + Ebin + q + factor(year), "E binary (any participation of family or friends)", En = "Ebin"),
            rob_fit(copy(d1)[, y := as.integer(own >= 3)], y ~ G + E + q + factor(year), "y = c22 >= 3 (at least 'a veces')"),
            rob_fit(d1[!is.na(march)][, y := march], y ~ G + E + q + factor(year), "y = march (c08_02 >= 2)"))
po <- suppressWarnings(polr(factor(own) ~ G + E + q + factor(year), data = d1, weights = w, Hess = TRUE))
ROB[[length(ROB) + 1]] <- data.table(variant = "ordinal c22 (1-5), ordered logit", n = nrow(d1), bG = coef(po)[["G"]], bE = coef(po)[["E"]], bq = coef(po)[["q"]],
                                     lambda = coef(po)[["E"]] / -coef(po)[["q"]], wG = coef(po)[["G"]] / -coef(po)[["q"]], sigma = 1 / -coef(po)[["q"]])
sel <- d[year %in% c(2016, 2018, 2022) & !is.na(q_c09) & !is.na(educ) & !is.na(ideol) & !is.na(agegrp) & !is.na(sex) & !is.na(relig)]
ms <- suppressWarnings(glm(named ~ q_c09 + educ + ideol + agegrp + sex + relig + factor(year), family = binomial, data = sel, weights = w))
sel[, p_named := fitted(ms)]
d1i <- merge(d1, sel[, .(id, year, p_named)], by = c("id", "year"))[, w := w / p_named][, w := w / mean(w), by = year]
ROB[[length(ROB) + 1]] <- rob_fit(d1i, y ~ G + E + q + factor(year), "IPW for naming a movement (attributes + q)")
ROB <- rbindlist(ROB); fwrite(ROB, file.path(OUT, "rule_robustness.csv"))
cat("\nROBUSTNESS OF S1 (point estimates):\n"); print(ROB, digits = 3)

# ---- additive vs multiplicative (S1 and S3 samples) ------------------------------------------------------
cv_logloss <- function(dat, f, K = 10, seed = 7) {     # folds by respondent, weighted log-loss
  set.seed(seed); ids <- unique(dat$id); fold <- setNames(sample(rep_len(1:K, length(ids))), ids); dat$fold <- fold[as.character(dat$id)]
  ll <- 0; ww <- 0
  for (k in 1:K) { tr <- dat[fold != k]; te <- dat[fold == k]; m <- wlogit(f, tr)
    p <- pmin(pmax(predict(m, te, type = "response"), 1e-6), 1 - 1e-6)
    ll <- ll - sum(te$w * (te$y * log(p) + (1 - te$y) * log(1 - p))); ww <- ww + sum(te$w) }
  ll / ww }
FORMS <- list(`A additive` = y ~ G + E + q + factor(year), `C complements (A + G:E)` = y ~ G + E + G:E + q + factor(year),
              `P pure product (G*E)` = y ~ I(G * E) + q + factor(year), `no exposure` = y ~ G + q + factor(year),
              `no efficacy` = y ~ E + q + factor(year), `A + E:q` = y ~ G + E + E:q + q + factor(year))
FM <- rbindlist(lapply(c("S1", "S3"), function(s) { dat <- copy(SPECS[[s]]$data)
  rbindlist(lapply(names(FORMS), function(fn) { m <- wlogit(FORMS[[fn]], dat); b <- coef(m)
    ll <- as.numeric(logLik(m)); k <- length(b)
    data.table(sample = s, form = fn, n = nrow(dat), loglik = ll, aic = -2 * ll + 2 * k, bic = -2 * ll + log(nrow(dat)) * k,
               cv_logloss = cv_logloss(dat, FORMS[[fn]]), b_GxE = if ("G:E" %in% names(b)) b[["G:E"]] else NA_real_,
               p_GxE = if ("G:E" %in% names(b)) coef(summary(m))["G:E", 4] else NA_real_,
               b_Exq = if ("E:q" %in% names(b)) b[["E:q"]] else NA_real_) })) }))
FM[, d_aic := aic - min(aic), by = sample]
fwrite(FM, file.path(OUT, "rule_forms.csv"))
cat("\nFUNCTIONAL FORMS (weighted pseudo-likelihood; cv = 10-fold log-loss by respondent):\n")
print(FM[, .(sample, form, n, aic = round(aic, 1), d_aic = round(d_aic, 1), bic = round(bic, 1), cv_logloss = round(cv_logloss, 4), b_GxE = round(b_GxE, 2), p_GxE = signif(p_GxE, 2), b_Exq = round(b_Exq, 2))])
cat(sprintf("Pure product minus additive: dAIC %.1f (S1), %.1f (S3)\n",
            FM[sample == "S1" & form == "P pure product (G*E)"]$aic - FM[sample == "S1" & form == "A additive"]$aic,
            FM[sample == "S3" & form == "P pure product (G*E)"]$aic - FM[sample == "S3" & form == "A additive"]$aic))
# the E = 0 test of the pure product: who takes part without any exposure?
E0 <- d1[, .(n = .N, share_E0 = NA_real_, y_E0 = wmean(y[E == 0], w[E == 0]), n_E0 = sum(E == 0), y_Epos = wmean(y[E > 0], w[E > 0])),
         by = .(qbin = cut(q, c(-0.01, 0.25, 0.5, 0.75, 0.999, 1), labels = c("q <= 0.25", "0.25-0.5", "0.5-0.75", "0.75-1", "q = 1")))][order(qbin)]
E0[, share_E0 := n_E0 / n]
fwrite(E0, file.path(OUT, "rule_E0_test.csv"))
cat("\nE = 0 TEST (S1 sample): participation with no exposure (family and friends 'never'), by q:\n"); print(E0, digits = 3)

# ---- non-linearity near the threshold ---------------------------------------------------------------------
mA <- FITS$S1; bA <- coef(mA)
d1[, m_hat := predict(mA, d1, type = "link") - bA[["E"]] * E]           # index without exposure (the margin)
d1[, mbin := cut(m_hat, quantile(m_hat, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))]
d1[, p_hat := predict(mA, d1, type = "response")]
TH <- d1[, { lp <- lm(y ~ E, weights = w); pe <- lm(p_hat ~ E, weights = w)   # empirical vs model-implied slope of P on E
  .(n = .N, m_mean = wmean(m_hat, w), m_in_q = wmean(m_hat, w) / -bA[["q"]], y_mean = wmean(y, w),
    slope_emp = coef(lp)[["E"]], slope_emp_se = sqrt(vcov(lp)["E", "E"]), slope_model = coef(pe)[["E"]],
    y_E0 = wmean(y[E == 0], w[E == 0]), y_Ehi = wmean(y[E >= 0.5], w[E >= 0.5])) }, by = mbin][order(mbin)]
fwrite(TH, file.path(OUT, "rule_threshold_bins.csv"))
cat("\nNEAR THE THRESHOLD: slope of P(y) on E by quintile of the margin (index without E), empirical vs additive logit:\n"); print(TH, digits = 3)
# where is the non-linearity? te(margin, E) = curvature and interaction; s(margin) + s(E) = curvature in E only
gm <- gam(y ~ te(m_hat, E, k = c(6, 4)), family = binomial, data = d1, weights = w, method = "REML")
ga <- gam(y ~ s(m_hat, k = 6) + s(E, k = 4), family = binomial, data = d1, weights = w, method = "REML")
gl <- gam(y ~ m_hat + E, family = binomial, data = d1, weights = w, method = "REML")
cat(sprintf("GAM te(margin, E) %.1f, additive s(margin) + s(E) %.1f, linear index %.1f (AIC); edf %.1f and %.1f\n",
            AIC(gm), AIC(ga), AIC(gl), sum(gm$edf), sum(ga$edf)))
fwrite(data.table(aic_gam = AIC(gm), aic_gam_additive = AIC(ga), aic_lin = AIC(gl), edf = sum(gm$edf), edf_additive = sum(ga$edf)),
       file.path(OUT, "rule_threshold_gam.csv"))
# observed vs rule by exposure bin (pooled over the margin)
d1[, Ebin5 := cut(E, c(-0.01, 0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.125-0.25", "0.375-0.5", "0.625-0.75", "0.875-1"))]
EB <- d1[, .(n = .N, observed = wmean(y, w), rule = wmean(p_hat, w), gam = wmean(predict(ga, .SD, type = "response"), w)), by = Ebin5][order(Ebin5)]
fwrite(EB, file.path(OUT, "rule_exposure_bins.csv"))
cat("\nOBSERVED vs RULE (additive logit) vs additive GAM, by exposure bin (S1 sample):\n"); print(EB, digits = 3)
grid <- CJ(m_hat = seq(quantile(d1$m_hat, .02), quantile(d1$m_hat, .98), length.out = 60), E = c(0, 0.25, 0.5, 0.75, 1))
grid[, `:=`(p_gam = predict(ga, grid, type = "response"), p_lin = predict(gl, grid, type = "response"))]   # te() adds nothing to ga
p <- ggplot(grid, aes(m_hat / -bA[["q"]], colour = factor(E))) + geom_line(aes(y = p_gam)) + geom_line(aes(y = p_lin), linetype = 2) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "margin without exposure, in units of q:  (index without E) / (−βq)", y = "P(takes part)", colour = "E",
       title = "Exposure is concave, at every distance from the threshold",
       subtitle = sprintf("solid: GAM s(margin) + s(E) (te(margin, E): AIC %+.1f); dashed: linear index (the rule). ELSOC, S1, weighted",
                          AIC(gm) - AIC(ga))) +
  theme_bw(base_size = 9)
ggsave(file.path(FIG, "rule_threshold.pdf"), p, width = 6, height = 3.6, device = cairo_pdf)

# ---- which q fits the rule? -----------------------------------------------------------------------------
QC <- list()
qfit <- function(dat, qv, lab, f_rhs = "G + E + qq") {
  dat <- copy(dat)[, qq := get(qv)]
  f <- as.formula(paste("y ~", f_rhs, if (uniqueN(dat$year) > 1) "+ factor(year)" else ""))
  m <- wlogit(f, dat); b <- coef(m); ll <- as.numeric(logLik(m))
  data.table(design = lab, q = qv, n = nrow(dat), aic = -2 * ll + 2 * length(b), cv_logloss = cv_logloss(dat, f),
             bq = b[["qq"]], lambda = if ("E" %in% names(b)) b[["E"]] / -b[["qq"]] else NA_real_, sigma = 1 / -b[["qq"]]) }
dC <- r[year %in% c(2018, 2022) & !is.na(G) & !is.na(E) & !is.na(q) & !is.na(q_gss4) & !is.na(q_gss3)]
dL <- r[year == 2019 & !is.na(G) & !is.na(E) & !is.na(q18) & !is.na(qg4_18) & !is.na(qg3_18)]
for (qv in c("q", "q_gss4", "q_gss3")) QC[[length(QC) + 1]] <- qfit(dC, qv, "contemporaneous 2018, 2022 (q_gss contains own c08)")
for (qv in c("q18", "qg4_18", "qg3_18")) QC[[length(QC) + 1]] <- qfit(dL, qv, "lagged: q of 2018 -> y of 2019")
for (qv in c("q18", "qg4_18", "qg3_18")) QC[[length(QC) + 1]] <- qfit(dL, qv, "lagged, q only (no E, no Gamma)", "qq")
dm <- merge(d[year == 2019 & !is.na(p2), .(id, y = as.integer(p2 >= 2), w)], d[year == 2018, .(id, q18 = q_c09, qg4_18 = q_gss4, qg3_18 = q_gss3)], by = "id")[
  !is.na(q18) & !is.na(qg4_18) & !is.na(qg3_18)][, year := 2019]
for (qv in c("q18", "qg4_18", "qg3_18")) QC[[length(QC) + 1]] <- qfit(dm, qv, "march 2019 (whole panel) from q of 2018 alone", "qq")
QC <- rbindlist(QC)
auc_w <- function(y, p, w) { o <- order(p); y <- y[o]; w <- w[o]; cw1 <- cumsum(w * y); tot1 <- sum(w * y); tot0 <- sum(w * (1 - y))
  sum(w * (1 - y) * (tot1 - cw1)) / (tot1 * tot0) }    # P(p1 > p0), weighted (ties ignored)
QC[, auc := NA_real_]
for (i in seq_len(nrow(QC))) {
  dat <- switch(QC$design[i], `contemporaneous 2018, 2022 (q_gss contains own c08)` = dC, `march 2019 (whole panel) from q of 2018 alone` = dm, dL)
  QC$auc[i] <- auc_w(dat$y, -dat[[QC$q[i]]], dat$w)                # AUC of q alone (lower q -> takes part)
}
fwrite(QC, file.path(OUT, "rule_q_compare.csv"))
cat("\nWHICH q? (same sample within each design; auc = q alone):\n"); print(QC, digits = 3)
saveRDS(list(S1 = FITS$S1, S3 = FITS$S3, S3b = FITS$S3b, ladder = LADDER, struct = S7), file.path(OUT, "rule_fits.rds"))

# ---- figure: lambda across the ladder, against the swept range ----------------------------------------------
st <- function(s, k) LADDER[spec == s & stat == k]
s7 <- function(k) S7[stat == k]
LF <- rbind(
  data.table(lab = "S1 perceived E, same wave", st("S1", "lambda")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S2 S1 + own attributes", st("S2", "lambda")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S3 lagged q", st("S3", "lambda")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S3b 2018-2019, q of 2018", st("S3b", "lambda")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S6 family only", st("S6", "lambda_fam")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S6 friends only", st("S6", "lambda_fri")[, .(est, lo, hi)], kind = "perceived (upper bound)"),
  data.table(lab = "S4 lagged E and lagged y", st("S4", "lambda_lag")[, .(est, lo, hi)], kind = "projection removed"),
  data.table(lab = "S7 structural E, no controls", s7("P2_lambda_struct")[, .(est, lo, hi)], kind = "projection removed"),
  data.table(lab = "S7 structural E + own attributes", s7("P3_lambda_struct")[, .(est, lo, hi)], kind = "projection removed"),
  data.table(lab = "S7 structural E + own predicted p", s7("P4_lambda_struct")[, .(est, lo, hi)], kind = "projection removed"))
LF[, lab := factor(lab, rev(lab))]
p2 <- ggplot(LF, aes(est, lab, colour = kind)) +
  annotate("rect", xmin = 0.5, xmax = 1, ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey40") +
  annotate("rect", xmin = 1, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "grey40") +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.25, orientation = "y", na.rm = TRUE) + geom_point() +
  scale_colour_manual(values = c(`perceived (upper bound)` = "#B03A2E", `projection removed` = "#1F618D"), name = NULL) +
  coord_cartesian(xlim = c(-0.5, 9)) +
  labs(x = expression(hat(lambda) == beta[E] / (-beta[q]) ~ "(units of q)"), y = NULL,
       title = "λ estimated in ELSOC along the identification ladder",
       subtitle = "bands: swept range 0.5–1.0 (dark), declared sensitivity to 1.5 (light); 95% respondent bootstrap") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "rule_lambda_ladder.pdf"), p2, width = 6.4, height = 3.8, device = cairo_pdf)

# ---- LaTeX tables (English, booktabs), written from the CSVs above ------------------------------------------
neg <- function(s) gsub("(^|[ \\[(])-", "\\1$-$", s)       # typographic minus
f2 <- function(x, dd = 2) ifelse(is.na(x), "---", neg(formatC(x, format = "f", digits = dd)))
ci <- function(e, lo, hi, dd = 2) paste0(f2(e, dd), " [", f2(lo, dd), ", ", f2(hi, dd), "]")
pc <- function(x, dd = 1) ifelse(is.na(x), "---", neg(formatC(100 * x, format = "f", digits = dd)))
thou <- function(n) gsub(",", "\\\\,", format(n, big.mark = ","))
row <- function(...) paste0(paste(..., sep = " & "), " \\\\")
wr <- function(lines, f) { writeLines(c("% generated by scripts/10_ELSOC_rule.R from output/10_ELSOC_rule/*.csv; do not edit by hand", lines), file.path(OUT, f)) }
# ladder
SP <- list(c("S1", "lambda", "Pooled 2016, 2018, 2022; q of the same wave"), c("S2", "lambda", "S1 + own attributes, movement, identification"),
           c("S3", "lambda", "Lagged q (2016/18/22 $\\to$ y of 2017/19/23)"), c("S3b", "lambda", "2018 and 2019, q of 2018"),
           c("S6", "lambda_fam", "Family only (friends in the model)"), c("S6", "lambda_fri", "Friends only (family in the model)"),
           c("S4", "lambda_lag", "Lagged E and lagged y, same movement"))
L <- c("\\begin{tabular}{@{}llrrlll@{}}", "\\toprule",
       " & Specification & Obs. & Resp. & $\\hat\\lambda$ [95\\% CI] & $\\hat w_\\Gamma$ [95\\% CI] & $\\hat\\sigma$ \\\\", "\\midrule",
       "\\multicolumn{7}{@{}l}{\\emph{Perceived exposure (family and friends in the same movement)}} \\\\")
for (s in SP) { x <- st(s[1], s[2]); wg <- st(s[1], "wG"); sg <- st(s[1], "sigma")
  if (s[1] == "S4") L <- c(L, "\\multicolumn{7}{@{}l}{\\emph{Projection removed}} \\\\")
  L <- c(L, row(s[1], s[3], thou(x$n), thou(x$n_id), ci(x$est, x$lo, x$hi), if (nrow(wg)) ci(wg$est, wg$lo, wg$hi) else "---", if (nrow(sg)) f2(sg$est) else "---")) }
S7ROWS <- list(c("P2", "Structural E of the 2017 confidants"), c("P3", "\\quad + own attributes"), c("P4", "\\quad + own predicted probability"))
for (s in S7ROWS) { x <- s7(paste0(s[1], "_lambda_struct")); L <- c(L, row("S7", s[2], thou(x$n), thou(x$n), ci(x$est, x$lo, x$hi), "---", "---")) }
s5 <- LADDER[spec == "S5"]
L <- c(L, "\\midrule",
       row("S5", "Person fixed effects (LPM): $\\beta_E^{FE}/\\beta_E^{pooled}$", thou(s5[stat == "bE_lpm_fe"]$n), thou(s5$n_id[1]),
           paste0(f2(s5[stat == "ratio_fe_pooled"]$est), " ($\\to\\hat\\lambda\\approx$ ", f2(s5[stat == "lambda_S1_x_ratio"]$est), ")"), "---", "---"),
       "\\bottomrule", "\\end{tabular}")
wr(L, "tab_rule_ladder.tex")
# functional forms (S1 and S3 samples) with the E = 0 test as a note line
FN <- c("A additive" = "Additive (the rule)", "C complements (A + G:E)" = "Additive + $\\Gamma E$ (complements)", "P pure product (G*E)" = "Pure product $\\Gamma E$",
        "no exposure" = "No exposure", "no efficacy" = "No efficacy", "A + E:q" = "Additive + $E q$")
L <- c("\\begin{tabular}{@{}lrrrl@{}}", "\\toprule", "Form & $\\Delta$AIC & BIC & CV log-loss & Interaction \\\\", "\\midrule")
for (s in c("S1", "S3")) { L <- c(L, paste0("\\multicolumn{5}{@{}l}{\\emph{Sample ", s, if (s == "S1") " (q of the same wave), n = " else " (lagged q), n = ",
                                            thou(FM[sample == s]$n[1]), "}} \\\\"))
  for (fn in names(FN)) { x <- FM[sample == s & form == fn]
    L <- c(L, row(paste0("\\quad ", FN[[fn]]), f2(x$d_aic, 1), f2(x$bic, 0), f2(x$cv_logloss, 4),
                  if (!is.na(x$b_GxE)) paste0(f2(x$b_GxE), if (x$p_GxE < 0.001) " ($p < 0.001$)" else paste0(" ($p$ = ", formatC(x$p_GxE, format = "f", digits = 3), ")")) else if (!is.na(x$b_Exq)) f2(x$b_Exq) else "---")) } }
wr(c(L, "\\bottomrule", "\\end{tabular}"), "tab_rule_forms.tex")
# E = 0 test
L <- c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule", "q & Obs. & With $E = 0$ (\\%) & Takes part if $E = 0$ (\\%) & Takes part if $E > 0$ (\\%) \\\\", "\\midrule",
       vapply(seq_len(nrow(E0)), function(i) row(c("$\\le 0.25$", "0.25--0.5", "0.5--0.75", "0.75--1", "$= 1$")[i], thou(E0$n[i]), pc(E0$share_E0[i]), pc(E0$y_E0[i]), pc(E0$y_Epos[i])), ""),
       "\\bottomrule", "\\end{tabular}")
wr(L, "tab_rule_E0.tex")
# threshold
L <- c("\\begin{tabular}{@{}lrrlrrr@{}}", "\\toprule",
       "Quintile of the margin & Margin (q) & Takes part (\\%) & Slope (SE) & Slope, rule & $E = 0$ (\\%) & $E \\ge 0.5$ (\\%) \\\\", "\\midrule",
       vapply(seq_len(nrow(TH)), function(i) row(as.character(TH$mbin[i]), f2(TH$m_in_q[i]), pc(TH$y_mean[i]), paste0(f2(TH$slope_emp[i]), " (", f2(TH$slope_emp_se[i]), ")"),
                                                f2(TH$slope_model[i]), pc(TH$y_E0[i]), pc(TH$y_Ehi[i])), ""),
       "\\bottomrule", "\\end{tabular}")
wr(L, "tab_rule_threshold.tex")
GG <- fread(file.path(OUT, "rule_threshold_gam.csv"))
L <- c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule", "E & Obs. & Observed (\\%) & Rule (\\%) & GAM $s(m) + s(E)$ (\\%) \\\\", "\\midrule",
       vapply(seq_len(nrow(EB)), function(i) row(sub("-", "--", as.character(EB$Ebin5[i])), thou(EB$n[i]), pc(EB$observed[i]), pc(EB$rule[i]), pc(EB$gam[i])), ""),
       "\\midrule", paste0("\\multicolumn{5}{@{}l}{AIC: linear index ", f2(GG$aic_lin, 0), "; $s(m) + s(E)$ ", f2(GG$aic_gam_additive, 0),
                          "; $te(m, E)$ ", f2(GG$aic_gam, 0), "} \\\\"), "\\bottomrule", "\\end{tabular}")
wr(L, "tab_rule_exposure_bins.tex")
# robustness of S1
RB <- c("main: y = c22 >= 2, E family+friends, Gamma c21_05" = "Main: y = c22 $\\ge 2$; E family and friends; $\\Gamma$ = c21\\_05",
        "Gamma = 3-item efficacy (c21_03-05)" = "$\\Gamma$ = 3-item efficacy (c21\\_03--05)", "E = friends only" = "E = friends only", "E = family only" = "E = family only",
        "E binary (any participation of family or friends)" = "E binary (any family member or friend took part)", "y = c22 >= 3 (at least 'a veces')" = "y = c22 $\\ge 3$ (at least sometimes)",
        "y = march (c08_02 >= 2)" = "y = march (c08\\_02 $\\ge 2$)", "ordinal c22 (1-5), ordered logit" = "c22 ordinal (1--5), ordered logit",
        "IPW for naming a movement (attributes + q)" = "IPW for naming a movement (attributes + q)")
L <- c("\\begin{tabular}{@{}lrrrr@{}}", "\\toprule", "Variant & Obs. & $\\hat\\lambda$ & $\\hat w_\\Gamma$ & $\\hat\\sigma$ \\\\", "\\midrule",
       vapply(names(RB), function(v) { x <- ROB[variant == v]; row(RB[[v]], thou(x$n), f2(x$lambda), f2(x$wG), f2(x$sigma)) }, ""),
       "\\bottomrule", "\\end{tabular}")
wr(L, "tab_rule_robustness.tex")

# ---- summary --------------------------------------------------------------------------------------------------
j <- st("S3b", "Gamma_jump_2019")
cat("\n==== SUMMARY ====\n")
cat(sprintf("lambda (units of q): S1 %s (perceived, upper bound); S4 %s (lagged); S7 %s (structural, net of own attributes)\n",
            ci(st("S1", "lambda")$est, st("S1", "lambda")$lo, st("S1", "lambda")$hi), ci(st("S4", "lambda_lag")$est, st("S4", "lambda_lag")$lo, st("S4", "lambda_lag")$hi),
            ci(s7("P3_lambda_struct")$est, s7("P3_lambda_struct")$lo, s7("P3_lambda_struct")$hi)))
cat(sprintf("weight of perceived efficacy wG (S1): %s; sigma (S1) %.2f\n", ci(st("S1", "wG")$est, st("S1", "wG")$lo, st("S1", "wG")$hi), st("S1", "sigma")$est))
cat(sprintf("Gamma jump 2018 -> 2019 with E held fixed (S3b): %+.3f [%+.3f, %+.3f] (units of q)\n", j$est, j$lo, j$hi))
cat(sprintf("pure product vs additive: dAIC +%.0f (S1), +%.0f (S3); participation with E = 0 and q <= 0.5: %.1f%%\n",
            FM[sample == "S1" & form == "P pure product (G*E)"]$d_aic - FM[sample == "S1" & form == "A additive"]$d_aic,
            FM[sample == "S3" & form == "P pure product (G*E)"]$d_aic - FM[sample == "S3" & form == "A additive"]$d_aic,
            100 * d1[q <= 0.5 & E == 0, wmean(y, w)]))
message(format(Sys.time(), "%H:%M:%S"), "  saved ", OUT, "/rule_*.csv, tab_rule_*.tex and ", FIG, "/rule_*.pdf (",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min)")
