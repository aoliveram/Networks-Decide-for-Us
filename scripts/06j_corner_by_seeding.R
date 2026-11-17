# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THE REVERSE CORNER BY SEEDING AND LAMBDA, up to lambda = 1.5
# 06j_corner_by_seeding.R  (2026-11-16, after Prompt 6.2; reads the main sweep since Prompt 6.3)
#
# Does the region where the DP twin beats the network survive lambda >= 1 under every
# seeding strategy, not only random? Reads the base-exposure results of the main sweep of 05
# (x2, lambda 0.5-1.5 every 0.1) for every seeding found, and counts, per
# survey x seeding x lambda, the cells where the DP twin is ahead (criterion of 06e: mean lead
# > 2 pp over the 24 networks and z > 3.5), over the whole plane and the upper half (h >= 0.6),
# with the largest DP lead and where it sits; also the premium (mean NET - DP over the plane).
# 06b's one-sided max-T test is stricter; this is the screening criterion of 06e and 06i.
#
# Inputs:  output/05_unified_diffusion/<S>_x2/unified_premium_results_<seeding>.csv (05 --level=x2)
# Outputs: output/06j_corner_by_seeding/corner_by_seeding.csv
#          plots/06j_corner_by_seeding/corner_by_seeding.pdf
# Usage (repo root): Rscript scripts/06j_corner_by_seeding.R    (~1 min)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
OUT <- "output/06j_corner_by_seeding"; FIG <- "plots/06j_corner_by_seeding"
for (p in c(OUT, FIG)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
SURVEYS <- c("GSS", "ATP", "GSS14", "ELSOC"); SEEDINGS <- c("random", "central", "marginal", "closeness", "eigen")
LAB <- c(GSS = "GSS 2004", ATP = "ATP 2014", GSS14 = "GSS 2014", ELSOC = "ELSOC 2018")
RES <- rbindlist(lapply(SURVEYS, function(S) rbindlist(lapply(SEEDINGS, function(sd) {
  f <- sprintf("output/05_unified_diffusion/%s_x2/unified_premium_results_%s.csv", S, sd)
  if (!file.exists(f)) { message("missing (skipped): ", f); return(NULL) }
  x <- fread(f)
  w <- dcast(x, lambda + run_id + social_distance_h + innovation_iul_Gamma ~ topology, value.var = "num_adopters")
  w[, lead := 100 * (DP - NET) / x$N_nodes[1]]
  cell <- w[, .(lead = mean(lead), z = mean(lead) / (sd(lead) / sqrt(.N))), by = .(lambda, h = social_distance_h, G = innovation_iul_Gamma)]
  cell[is.na(z), z := 0]
  cell[, .(survey = S, seeding = sd, premium_pp = -mean(lead),
           robust_cells = sum(lead > 2 & z > 3.5), robust_cells_upper = sum(lead > 2 & z > 3.5 & h >= 0.6 - 1e-9),
           cells_dp_ahead = sum(lead > 0), max_dp_lead_pp = max(lead),
           G_of_max = G[which.max(lead)], h_of_max = h[which.max(lead)]), by = lambda]
}))))
setcolorder(RES, c("survey", "seeding", "lambda")); setorder(RES, survey, seeding, lambda)
fwrite(RES, file.path(OUT, "corner_by_seeding.csv"))
cat("ROBUST DP-AHEAD CELLS over the whole plane (lead > 2 pp, z > 3.5), by survey x seeding x lambda:\n")
print(dcast(RES, survey + seeding ~ lambda, value.var = "robust_cells"))
cat("\nLARGEST DP LEAD anywhere on the plane (pp):\n")
print(dcast(RES[, .(survey, seeding, lambda, m = round(max_dp_lead_pp, 2))], survey + seeding ~ lambda, value.var = "m"))
cat("\nPREMIUM (pp, mean over the plane):\n")
print(dcast(RES[, .(survey, seeding, lambda, p = round(premium_pp, 2))], survey + seeding ~ lambda, value.var = "p"))
cat(sprintf("\nlambda >= 1: %d robust cells and a largest DP lead of %.2f pp over %d survey x seeding combinations\n",
            sum(RES[lambda >= 1 - 1e-9]$robust_cells), max(RES[lambda >= 1 - 1e-9]$max_dp_lead_pp), uniqueN(RES[, .(survey, seeding)])))
cat(sprintf("lambda > 1 (1.1-1.5): %d robust cells and a largest DP lead of %.2f pp\n",
            sum(RES[lambda > 1 + 1e-9]$robust_cells), max(RES[lambda > 1 + 1e-9]$max_dp_lead_pp)))
p <- ggplot(RES, aes(lambda, robust_cells, colour = seeding)) + geom_line() + geom_point(size = 1.2) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  facet_wrap(~ factor(LAB[survey], LAB), nrow = 1) +
  labs(x = "λ", y = "robust cells with the DP twin ahead", colour = "seeding",
       title = "The reverse corner by seeding and λ",
       subtitle = "x2 networks, base exposure; cells of the Γ × h plane with DP ahead by > 2 pp and z > 3.5 over 24 networks") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "corner_by_seeding.pdf"), p, width = 9, height = 3.4, device = cairo_pdf)
message("saved ", OUT, "/corner_by_seeding.csv and ", FIG, "/corner_by_seeding.pdf")
