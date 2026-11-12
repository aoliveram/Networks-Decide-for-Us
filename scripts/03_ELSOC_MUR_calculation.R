# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ELSOC: the aversion q of the nodes. 03_ELSOC_MUR_calculation.R
#
# Index (decided with the user, 2026-11-11): willingness c09, four items (sign a petition,
# attend a peaceful march, take part in a strike, use social media to give opinions;
# 1 = not at all willing ... 5 = totally willing), wave 2018:
#   propensity_score = (sum of the 4 items - 4) / 16  in [0, 1] (high = willing)
#   mur_score = q = 1 - propensity_score            (aversion; the engine reads only q)
# Why c09 and not the frequency c08: c09 measures aversion directly, as the "would never do
# it" of the ISSP items of GSS 2004; c08 is participation history (reliability 0.56-0.70,
# half of the respondents at "never" on every act). Why 2018: the largest wave (original +
# refresh sample), the year after the network wave of 2017. The nodes of the networks are
# the 2018 respondents with c09 and the five Blau attributes complete, drawn with probability
# proportional to ponderador02; homophily comes from the confidants of 2017
# (02_ELSOC_homophily.R). The index is built in 01_survey_data.R (elsoc_nodes) and read by
# 02_network_calibrated.R. Shares and means are reported unweighted and weighted
# (ponderador02 of the wave; the networks see the weighted distribution).
#
# Reports: Cronbach's alpha of c09 (2016, 2018, 2022) and c08 (4 items every wave; 5 with the
# cacerolazo in 2019 and 2022); distribution (mean, sd, floor, ceiling); c08-c09 correlation
# in the same wave; panel stability; coverage of the network respondents; and the share of
# the variance of q explained by the Blau attributes in the node pool (consolidation, see 06e).
# Outputs: output/03_ELSOC_MUR/{alpha,distribution,node_pool}.csv; plots/03_ELSOC_MUR/elsoc_q.pdf
# Usage (repo root): Rscript scripts/03_ELSOC_MUR_calculation.R   (seconds)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(psych); library(ggplot2) }))
source("scripts/01_survey_data.R")
OUT <- "output/03_ELSOC_MUR"; PLOTS <- "plots/03_ELSOC_MUR"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)
df <- elsoc_load()
for (v in c(ELSOC_C08, "c08_05", ELSOC_C09)) df[[v]] <- num_pos(df[[v]])

alpha_of <- function(X) { X <- X[complete.cases(X), , drop = FALSE]
  a <- suppressWarnings(psych::alpha(as.data.frame(X), warnings = FALSE, check.keys = FALSE))
  data.table(n = nrow(X), alpha_raw = a$total$raw_alpha, alpha_std = a$total$std.alpha,
             r_drop_min = min(a$item.stats$r.drop), r_drop_max = max(a$item.stats$r.drop)) }
prop_of <- function(d, items) (rowMeans(as.matrix(d[, ..items])) - 1) / 4      # = (sum - k) / (4 k)
AL <- list(); DI <- list()
for (w in 1:7) {
  d <- df[ola == w]; yr <- ELSOC_YEAR[[as.character(w)]]
  sets <- list(`c09 (4)` = ELSOC_C09, `c08 (4)` = ELSOC_C08, `c08 (5)` = c(ELSOC_C08, "c08_05"))
  for (nm in names(sets)) {
    it <- sets[[nm]]; if (any(colSums(!is.na(as.matrix(d[, ..it]))) < 100)) next
    AL[[length(AL) + 1]] <- cbind(data.table(index = nm, year = yr), alpha_of(as.matrix(d[, ..it])))
    p <- prop_of(d, it); w <- as.numeric(d$ponderador02)[!is.na(p)]; p <- p[!is.na(p)]
    DI[[length(DI) + 1]] <- data.table(index = nm, year = yr, n = length(p), mean_propensity = mean(p), sd = sd(p),
                                       mean_q = 1 - mean(p), floor_q1 = mean(p == 0), ceiling_q0 = mean(p == 1), distinct = uniqueN(p),
                                       mean_q_w = 1 - weighted.mean(p, w), floor_q1_w = weighted.mean(p == 0, w))
  }
}
AL <- rbindlist(AL); DI <- rbindlist(DI)
fwrite(AL, file.path(OUT, "alpha.csv")); fwrite(DI, file.path(OUT, "distribution.csv"))
cat("\nCRONBACH ALPHA (ELSOC, by index and wave):\n"); print(AL, digits = 3)
cat("\nDISTRIBUTION of the propensity (items rescaled to [0, 1]); floor_q1 = share with q = 1:\n"); print(DI, digits = 3)

cat("\nCorrelation of c08 (4) and c09 in the same wave:\n")
for (w in c(1, 3, 6)) { d <- df[ola == w]; cat(" ", ELSOC_YEAR[[as.character(w)]], ": r =",
                                               round(cor(prop_of(d, ELSOC_C08), prop_of(d, ELSOC_C09), use = "complete.obs"), 3), "\n") }
wide <- function(w, it) { d <- df[ola == w]; data.table(idencuesta = d$idencuesta, p = prop_of(d, it)) }
s <- merge(wide(1, ELSOC_C09), wide(3, ELSOC_C09), by = "idencuesta")
cat("Panel stability, c09 2016 -> 2018: r =", round(cor(s$p.x, s$p.y, use = "complete.obs"), 3), " (n =", sum(complete.cases(s)), ")\n")
s <- merge(wide(3, ELSOC_C09), wide(6, ELSOC_C09), by = "idencuesta")
cat("Panel stability, c09 2018 -> 2022: r =", round(cor(s$p.x, s$p.y, use = "complete.obs"), 3), " (n =", sum(complete.cases(s)), ")\n")

# coverage: respondents of the network waves with c09 in 2018
nodes <- elsoc_nodes(3); eg <- elsoc_egos(2)
cat(sprintf("\nNetwork respondents (2017) with a 2018 node record: %d of %d\n",
            sum(unique(eg$idencuesta) %in% nodes$idencuesta), uniqueN(eg$idencuesta)))

# node pool and association of q with the Blau attributes
r2 <- summary(lm(mur_score ~ age + sex + educ_num + relig + ideol, data = nodes))$r.squared
r2e <- summary(lm(mur_score ~ age + sex + educ + relig + ideol, data = nodes))$r.squared
r2w <- summary(lm(mur_score ~ age + sex + educ_num + relig + ideol, data = nodes, weights = w))$r.squared
np <- data.table(survey = "ELSOC", wave = 2018, n = nrow(nodes), mean_q = mean(nodes$mur_score), sd_q = sd(nodes$mur_score),
                 share_q1 = mean(nodes$mur_score == 1), distinct = uniqueN(nodes$mur_score),
                 r2_q_blau = r2, r2_q_blau_educ_groups = r2e,
                 mean_q_w = weighted.mean(nodes$mur_score, nodes$w), share_q1_w = weighted.mean(nodes$mur_score == 1, nodes$w),
                 r2_q_blau_w = r2w, kish_n = sum(nodes$w)^2 / sum(nodes$w^2))
fwrite(np, file.path(OUT, "node_pool.csv"))
cat("\nNODE POOL (ELSOC 2018, c09 and Blau attributes complete):\n"); print(np, digits = 3)
cat("\nMean q by group of the node pool:\n")
for (v in c("educ", "ideol", "relig", "sex")) print(nodes[, .(n = .N, share_w = round(sum(w) / sum(nodes$w), 3), mean_q = round(mean(mur_score), 3),
                                                              mean_q_w = round(weighted.mean(mur_score, w), 3)), by = c(v)][order(get(v))])
print(round(coef(summary(lm(mur_score ~ age + sex + educ_num + relig + ideol, data = nodes))), 4))

# figure: distribution of q, 2018 (nodes) and the c08 alternative
pd <- rbind(data.table(index = "c09 willingness, 2018 (used)", q = nodes$mur_score, w = nodes$w),
            data.table(index = "c08 frequency, 2018 (not used)", q = 1 - prop_of(df[ola == 3], ELSOC_C08), w = as.numeric(df[ola == 3]$ponderador02))[!is.na(q)])
p <- ggplot(pd, aes(q, weight = w)) + geom_histogram(binwidth = 1/16, boundary = 0, fill = "#2980B9", color = "white") +
  facet_wrap(~ index, scales = "free_y") +
  labs(title = "ELSOC: aversion q = 1 - propensity (weighted by ponderador02)", x = "q", y = "Respondents (weighted)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(size = 9, face = "bold"))
ggsave(file.path(PLOTS, "elsoc_q.pdf"), p, width = 7, height = 3, device = cairo_pdf)
message("Saved ", OUT, "/*.csv and ", file.path(PLOTS, "elsoc_q.pdf"))
