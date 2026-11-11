# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GSS 2014: aversion indices, POLVIEWS and PEOCNTCT. 03_GSS14_MUR_calculation.R
#
# 1. Collective action, 8 of the 9 ISSP items of GSS 2004 (INTERPOL was not asked in 2014):
#    seven acts (1 done in the past year, 2 in the more distant past, 3 might do it, 4 would
#    never do it) and ACTLAW (1 very likely ... 4 not at all likely):
#      propensity = sum(4 - x) / 24,  q = mur_score = 1 - propensity
#    Compared with the 9 items of GSS 2004 (sum / 27, 03_GSS_MUR_calculation.R) and with the
#    same 8 items in 2004: alpha and distribution. These are the q of the GSS14 networks
#    (01_survey_data.R, gss14_nodes("ca")): the 2014 cross-section (the ISSP module was asked
#    only there), Blau attributes and POLVIEWS complete.
# 2. Attitude toward science (innovation), 4 items: NEXTGEN and ADVFRONT (agree = pro),
#    TOOFAST (agree = con), SCIBNFTS (benefits greater / about equal / harmful greater), each
#    rescaled to [0, 1] toward pro-science and averaged: propensity, q = 1 - propensity; alpha.
#    Attitudes, not reported adoption (unlike ATP's METECH).
# 3. POLVIEWS (1 extremely liberal ... 7 extremely conservative) -> ideol: Liberal (1-3),
#    Moderate (4), Conservative (5-7), the circles of the ideology theta arm. Its association
#    with q, and its sorting in the networks: GSS 2014 has no ego networks, so the networks
#    sort on ideology only through its correlation with the five Blau attributes.
# 4. PEOCNTCT (contacts in a typical weekday: 0-4, 5-9, 10-19, 20-49, 50+) by race, as an
#    alternative to the GSS 2004 relative number of confidants (the targets of 02).
# 5. Association of each q with the Blau attributes (R2), as in 06e.
# Weights: the GSS 2014 nodes are drawn with probability proportional to WTSSNR (01_survey_data.R);
# their mean q, the POLVIEWS shares and PEOCNTCT by race are also reported weighted (WTSSNR of
# the 2014 cross-section). GSS 2004 stays unweighted, as in its networks.
# Outputs: output/03_GSS14_MUR/{alpha,distribution,peocntct_by_race,polviews,association}.csv
#          plots/03_GSS14_MUR/gss14_q.pdf
# Usage (repo root): Rscript scripts/03_GSS14_MUR_calculation.R   (seconds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(psych); library(ggplot2); library(patchwork) }))
source("scripts/01_survey_data.R")
OUT <- "output/03_GSS14_MUR"; PLOTS <- "plots/03_GSS14_MUR"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
alpha_of <- function(X) { X <- X[complete.cases(X), , drop = FALSE]
  a <- suppressWarnings(psych::alpha(as.data.frame(X), warnings = FALSE, check.keys = FALSE))
  data.table(n = nrow(X), alpha_raw = a$total$raw_alpha, alpha_std = a$total$std.alpha,
             r_drop_min = min(a$item.stats$r.drop), r_drop_max = max(a$item.stats$r.drop)) }
dist_of <- function(p) { p <- p[!is.na(p)]
  data.table(n = length(p), mean_propensity = mean(p), sd = sd(p), mean_q = 1 - mean(p), share_q1 = mean(p == 0),
             share_q0 = mean(p == 1), distinct = uniqueN(p)) }

g14 <- gss_load(2014); g04 <- gss_load(2004); cs14 <- as.numeric(g14$samptype) == 2014
A14 <- gss_attributes(g14); A04 <- gss_attributes(g04)
blau_ok <- function(A) complete.cases(A[, .(age, sex, educ_num, race, relig)])
X14 <- gss_items(g14, GSS_CA8)[cs14, ]; X04_9 <- gss_items(g04, GSS_CA9); X04_8 <- X04_9[, GSS_CA8]
S14 <- sci_items(g14)[cs14, ]

# ---- 1-2. alpha and distributions -------------------------------------------------------------------
AL <- rbind(cbind(data.table(index = "collective action, 9 items", survey = "GSS 2004"), alpha_of(X04_9)),
            cbind(data.table(index = "collective action, same 8 items", survey = "GSS 2004"), alpha_of(X04_8)),
            cbind(data.table(index = "collective action, 8 items", survey = "GSS 2014"), alpha_of(X14)),
            cbind(data.table(index = "science attitudes, 4 items", survey = "GSS 2014"), alpha_of(S14)))
fwrite(AL, file.path(OUT, "alpha.csv"))
cat("\nCRONBACH ALPHA:\n"); print(AL, digits = 3)
p04_9 <- ca_propensity(X04_9)[blau_ok(A04)]; p04_8 <- ca_propensity(X04_8)[blau_ok(A04)]
nodes_ca <- gss14_nodes("ca"); nodes_sci <- gss14_nodes("sci")
cat(sprintf("GSS 2014 nodes, weighted (WTSSNR): mean q %.3f (unweighted %.3f); Kish effective n %d of %d\n",
            weighted.mean(nodes_ca$mur_score, nodes_ca$w), mean(nodes_ca$mur_score), round(sum(nodes_ca$w)^2 / sum(nodes_ca$w^2)), nrow(nodes_ca)))
DI <- rbind(cbind(data.table(index = "collective action, 9 items", survey = "GSS 2004 (Blau complete)"), dist_of(p04_9)),
            cbind(data.table(index = "collective action, same 8 items", survey = "GSS 2004 (Blau complete)"), dist_of(p04_8)),
            cbind(data.table(index = "collective action, 8 items", survey = "GSS 2014 nodes"), dist_of(nodes_ca$propensity_score)),
            cbind(data.table(index = "science attitudes, 4 items", survey = "GSS 2014 (Blau + polviews complete)"), dist_of(nodes_sci$propensity_score)))
fwrite(DI, file.path(OUT, "distribution.csv"))
cat("\nDISTRIBUTION of the propensity (q = 1 - propensity):\n"); print(DI, digits = 3)
cat("\nItem means (1 done in the past year ... 4 would never; actlaw 1 very likely ... 4 not at all):\n")
print(rbind(data.table(survey = "GSS 2004", t(round(colMeans(X04_9, na.rm = TRUE), 2))),
            data.table(survey = "GSS 2014", t(round(colMeans(X14, na.rm = TRUE), 2))), fill = TRUE))
both <- complete.cases(X14) & complete.cases(S14)
cat(sprintf("\nGSS 2014 cross-section with both blocks: %d; correlation of the two propensities: %.3f\n", sum(both),
            cor(ca_propensity(X14)[both], rowMeans(S14)[both])))
cat("Kolmogorov-Smirnov, GSS 2014 (8) vs GSS 2004 (same 8): D =",
    round(suppressWarnings(ks.test(nodes_ca$propensity_score, p04_8))$statistic, 3), "\n")

# ---- 3. POLVIEWS -------------------------------------------------------------------------------------
pv <- nodes_ca[, .(n = .N, share = .N / nrow(nodes_ca), share_w = sum(w) / sum(nodes_ca$w), mean_q = mean(mur_score),
                   mean_q_w = weighted.mean(mur_score, w)), by = ideol][order(match(ideol, c("Liberal", "Moderate", "Conservative")))]
fwrite(pv, file.path(OUT, "polviews.csv"))
cat("\nPOLVIEWS in the GSS 2014 nodes (3 groups) and mean q:\n"); print(pv, digits = 3)
cat("R2 of polviews (3 groups) on the Blau attributes (multinomial pseudo-R2 via lm on each dummy):\n")
for (k in c("Liberal", "Moderate", "Conservative"))
  cat("  ", k, ":", round(summary(lm(I(ideol == k) ~ age + sex + educ_num + race + relig, data = nodes_ca))$r.squared, 3), "\n")

# ---- 4. PEOCNTCT by race ------------------------------------------------------------------------------
mid <- c(2, 7, 14.5, 34.5, 60)                     # midpoints of 0-4, 5-9, 10-19, 20-49, 50+ (50+ at 60)
pc <- data.table(race = A14$race, peocntct = as.numeric(zap_labels(g14$peocntct)), w = as.numeric(zap_labels(g14$wtssnr)))[
  cs14 & peocntct %in% 1:5 & !is.na(race)]
pc[, contacts := mid[peocntct]]
pcr <- pc[, .(n = .N, wn = sum(w), mean_contacts = weighted.mean(contacts, w), mean_category = weighted.mean(peocntct, w)), by = race]
pcr[, rel_contacts := mean_contacts / weighted.mean(mean_contacts, wn)][, rel_category := mean_category / weighted.mean(mean_category, wn)]
tg04 <- fread("output/02_calibrated_networks/gss2004_targets.csv")
pcr <- merge(pcr, tg04[, .(race = ego, gss2004_rel_confidants = rel_deg, rel_lo, rel_hi)], by = "race")[match(RACES, race)]
set.seed(2014)
bt <- rbindlist(lapply(1:1000, function(b) { x <- pc[sample.int(.N, replace = TRUE)][, .(wn = sum(w), m = weighted.mean(contacts, w)), by = race]
  x[, .(race, rel = m / weighted.mean(m, wn))] }))[, .(peo_lo = quantile(rel, .025), peo_hi = quantile(rel, .975)), by = race]
pcr <- merge(pcr, bt, by = "race")[match(RACES, race)]
fwrite(pcr, file.path(OUT, "peocntct_by_race.csv"))
cat("\nPEOCNTCT (GSS 2014 cross-section, weighted by WTSSNR) vs GSS 2004 confidants: relative degree by race (95% CI):\n")
print(pcr[, .(race, n, mean_contacts = round(mean_contacts, 1), peocntct = sprintf("%.3f [%.3f, %.3f]", rel_contacts, peo_lo, peo_hi),
              by_category = round(rel_category, 3), gss2004_confidants = sprintf("%.3f [%.3f, %.3f]", gss2004_rel_confidants, rel_lo, rel_hi))])

# ---- 5. association with the Blau attributes ------------------------------------------------------------
r2 <- function(d) summary(lm(mur_score ~ age + sex + educ_num + race + relig, data = d))$r.squared
n04 <- cbind(A04, mur_score = 1 - ca_propensity(X04_9))[complete.cases(cbind(A04[, .(age, sex, educ_num, race, relig)], ca_propensity(X04_9)))]
assoc <- data.table(index = c("collective action, 9 items (GSS 2004)", "collective action, 8 items (GSS 2014 nodes)",
                              "science attitudes (GSS 2014)", "collective action, 8 items (GSS 2014 nodes) + polviews"),
                    r2_q_blau = c(r2(n04), r2(nodes_ca), r2(nodes_sci),
                                  summary(lm(mur_score ~ age + sex + educ_num + race + relig + ideol, data = nodes_ca))$r.squared))
fwrite(assoc, file.path(OUT, "association.csv"))
cat("\nR2 of q on age, sex, education, race and religion (respondents):\n"); print(assoc, digits = 3)

# ---- figure -------------------------------------------------------------------------------------------
pd <- rbind(data.table(set = "GSS 2004, 9 items", q = 1 - p04_9), data.table(set = "GSS 2004, same 8 items", q = 1 - p04_8),
            data.table(set = "GSS 2014, 8 items (nodes)", q = nodes_ca$mur_score),
            data.table(set = "GSS 2014, science attitudes", q = nodes_sci$mur_score))
pd[, set := factor(set, unique(set))]
p1 <- ggplot(pd, aes(q)) + geom_histogram(aes(y = after_stat(density)), bins = 25, boundary = 0, fill = "#2980B9", color = "white") +
  facet_wrap(~ set, nrow = 1) + labs(title = "Aversion q = 1 - propensity", x = "q", y = "Density") +
  theme_bw(base_size = 8) + theme(plot.title = element_text(size = 8, face = "bold"))
pp <- melt(pcr[race != "Other" | TRUE, .(race, `PEOCNTCT 2014` = rel_contacts, `Confidants 2004` = gss2004_rel_confidants)], id.vars = "race")
p2 <- ggplot(pp, aes(factor(race, RACES), value, color = variable)) + geom_hline(yintercept = 1, color = "grey60", linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.4), size = 2) + scale_color_manual(values = c("#C0392B", "black"), name = NULL) +
  labs(title = "Relative degree by race: daily contacts (GSS 2014) vs confidants (GSS 2004, the targets)", x = NULL, y = "Relative to the mean") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 8, face = "bold"))
ggsave(file.path(PLOTS, "gss14_q.pdf"), p1 / p2, width = 9, height = 5.5, device = cairo_pdf)
message("Saved ", OUT, "/*.csv and ", file.path(PLOTS, "gss14_q.pdf"))
