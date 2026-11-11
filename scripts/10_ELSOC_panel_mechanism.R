# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# A DESCRIPTIVE TEST OF THE MECHANISM WITH THE ELSOC PANEL. 10_ELSOC_panel_mechanism.R
#
# Did participation rise more in 2019 (the wave was fielded after the social outburst of
# 18 October 2019: 670 respondents name it as the movement they value most; attending a
# march at least "almost never" went from 12% in 2017 to 26% in 2019) among respondents
# whose 2017 confidants were SIMILAR to them or among those whose confidants were DIVERSE?
# It does not identify influence (confidants are chosen, and the 2017 network is measured
# once); it gives the sign of the association, net of the obvious confounders.
#
# Sample: respondents with confidants in 2017 (r13, wave 2) who named at least 2 (so that a
# diversity index exists) and answered c08 in 2017 and 2019.
# Outcomes (2019):
#   c08 index     mean of the four acts (petition, march, strike, social media), 1-5 -> [0, 1]
#   march         attended a march or peaceful demonstration at least "almost never" (c08_02 >= 2)
#   cacerolazo    took part in cacerolazos at least "almost never" (c08_05 >= 2; asked from 2019 on,
#                 so its model has no lagged outcome, only the 2017 c08 index)
# 2017 network measures (confidants with the attribute known; egos with >= 2 such confidants):
#   blau_ideol, blau_relig   Blau heterogeneity 1 - sum_g p_g^2 of the confidants (ideology in 3
#                            groups: Left, CentreNone, Right; religion in 5)
#   same_ideol, same_relig   share of confidants in the ego's own group (similarity)
#   left_share               share of confidants on the left (political composition; model 3)
# Controls (2017): the outcome in 2017 (lagged dependent variable), sex, age, education,
# religion and ideology of ego, number of confidants named. OLS (index) and logit (binary),
# HC1 robust standard errors; standardized network measures (1 sd).
# Outputs: output/10_ELSOC_panel_mechanism/{models,descriptives}.csv;
#          plots/10_ELSOC_panel_mechanism/panel_mechanism.pdf
# Usage (repo root): Rscript scripts/10_ELSOC_panel_mechanism.R   (seconds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(sandwich); library(lmtest); library(ggplot2) }))
source("scripts/01_survey_data.R")
OUT <- "output/10_ELSOC_panel_mechanism"; PLOTS <- "plots/10_ELSOC_panel_mechanism"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
df <- elsoc_load()

# ---- data ----------------------------------------------------------------------------------------
part <- function(w) { d <- df[ola == w]; X <- sapply(c(ELSOC_C08, "c08_05"), function(v) num_pos(d[[v]]))
  data.table(idencuesta = as.numeric(d$idencuesta), c08 = (rowMeans(X[, 1:4]) - 1) / 4, march = as.numeric(X[, 2] >= 2),
             caceroleo = as.numeric(X[, 5] >= 2)) }
p17 <- part(2); p19 <- part(4)
eg <- elsoc_egos(2); al <- elsoc_alters(2)
blau <- function(x) { x <- x[!is.na(x)]; if (length(x) < 2) NA_real_ else 1 - sum(prop.table(table(x))^2) }
a <- merge(al, eg[, .(idencuesta, e_ideol = ideol, e_relig = relig)], by = "idencuesta")
nw <- a[, .(n_ideol = sum(!is.na(ideol)), n_relig = sum(!is.na(relig)), blau_ideol = blau(ideol), blau_relig = blau(relig),
            same_ideol = mean(ideol == e_ideol, na.rm = TRUE), same_relig = mean(relig == e_relig, na.rm = TRUE),
            left_share = mean(ideol == "Left", na.rm = TRUE)), by = idencuesta]
d <- Reduce(function(x, y) merge(x, y, by = "idencuesta"),
            list(eg[, .(idencuesta, sex, age, educ, relig, ideol, nredes)], nw,
                 setnames(copy(p17), c("idencuesta", "c08_17", "march_17", "cac_17")),
                 setnames(copy(p19), c("idencuesta", "c08_19", "march_19", "cac_19"))))
d <- d[n_ideol >= 2 & n_relig >= 2 & complete.cases(d[, .(sex, age, educ, relig, ideol, nredes, blau_ideol, blau_relig, same_ideol,
                                                           same_relig, c08_17, c08_19, march_17, march_19)])]
d[, educ := factor(educ, EDU5)][, relig := factor(relig, REL5_CL)][, ideol := factor(ideol, IDE3)]
Z <- c("blau_ideol", "blau_relig", "same_ideol", "same_relig", "left_share")
for (v in Z) d[, paste0("z_", v) := (get(v) - mean(get(v))) / sd(get(v))]
message(format(Sys.time(), "%H:%M:%S"), "  panel sample: ", nrow(d), " respondents (confidants in 2017, >= 2 with ideology and religion known, c08 in 2017 and 2019)")

desc <- d[, .(n = .N, c08_17 = mean(c08_17), c08_19 = mean(c08_19), march_17 = mean(march_17), march_19 = mean(march_19),
              cacerolazo_19 = mean(cac_19, na.rm = TRUE), blau_ideol = mean(blau_ideol), blau_relig = mean(blau_relig),
              same_ideol = mean(same_ideol), same_relig = mean(same_relig), left_share = mean(left_share), nredes = mean(nredes))]
fwrite(desc, file.path(OUT, "descriptives.csv")); cat("\nDESCRIPTIVES:\n"); print(desc, digits = 3)
cat("\nCorrelations among the 2017 network measures:\n"); print(round(cor(d[, ..Z]), 2))

# ---- models --------------------------------------------------------------------------------------
CTRL <- "sex + age + educ + relig + ideol + nredes"
SPECS <- list(
  diversity  = "z_blau_ideol + z_blau_relig",
  similarity = "z_same_ideol + z_same_relig",
  diversity_composition = "z_blau_ideol + z_blau_relig + z_left_share")
OUTCOMES <- list(c08 = list(y = "c08_19", lag = "c08_17", fam = "ols"), march = list(y = "march_19", lag = "march_17", fam = "logit"),
                 cacerolazo = list(y = "cac_19", lag = "c08_17", fam = "logit"))
res <- rbindlist(lapply(names(OUTCOMES), function(o) rbindlist(lapply(names(SPECS), function(s) {
  oc <- OUTCOMES[[o]]; f <- as.formula(paste(oc$y, "~", oc$lag, "+", SPECS[[s]], "+", CTRL))
  dd <- d[!is.na(get(oc$y))]
  m <- if (oc$fam == "ols") lm(f, data = dd) else glm(f, data = dd, family = binomial())
  ct <- coeftest(m, vcov = vcovHC(m, type = "HC1"))
  keep <- grep("^z_", rownames(ct), value = TRUE)
  data.table(outcome = o, spec = s, term = keep, estimate = ct[keep, 1], se = ct[keep, 2], p = ct[keep, 4], n = nobs(m),
             scale = if (oc$fam == "ols") "c08 index points (0-1) per sd" else "log-odds per sd")
}))))
res[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
fwrite(res, file.path(OUT, "models.csv"))
cat("\nMODELS (2019 outcome on 2017 network measures, standardized; lagged outcome and controls; HC1 SE):\n")
print(res[, .(outcome, spec, term, est = sprintf("%+.4f [%+.4f, %+.4f]", estimate, lo, hi), p = signif(p, 2), n)], nrows = 100)

# ---- figure --------------------------------------------------------------------------------------
pd <- copy(res)[, term := factor(sub("^z_", "", term), rev(Z))]
p <- ggplot(pd, aes(estimate, term, color = spec)) + geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.6) + facet_wrap(~ outcome, scales = "free_x") +
  scale_color_manual(values = c(diversity = "#2166AC", similarity = "#B2182B", diversity_composition = "#1B7837"), name = NULL) +
  labs(title = "ELSOC: 2019 participation (post 18-O) on the 2017 confidants (1 sd; 95% CI; lagged outcome and controls)",
       x = "Effect per sd (c08: index points; march, cacerolazo: log-odds)", y = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(size = 8, face = "bold"))
ggsave(file.path(PLOTS, "panel_mechanism.pdf"), p, width = 8.5, height = 3.6, device = cairo_pdf)
message("Saved ", OUT, "/*.csv and ", file.path(PLOTS, "panel_mechanism.pdf"))
