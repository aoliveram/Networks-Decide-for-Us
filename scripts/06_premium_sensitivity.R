# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# How sensitive is the structural premium to lambda? (decile-by-decile)
# 07_premium_sensitivity.R
#
# Reads output/05_unified_diffusion/unified_premium_results.csv (produced by
# 05_unified_diffusion_sweep.R) and asks, for every lambda on the grid:
#
#   A. the premium on the ODDS scale (BAM beta_DP) and its step-to-step change
#   B. the premium on the PROBABILITY scale (mean Phi_GSS - Phi_DP) — link-free
#   C. WHERE the premium lives on the (IUL, MSP) plane, and how wide that band is
#   D. the CHANNEL decomposition (rational vs social) of the premium
#   E. the TIPPING SHIFT: MSP at which each topology crosses Phi = 0.5
#   F. saturation diagnostics (how much of the plane is already at Phi ~ 1)
#
# The point of B-F is that beta_DP is a single constant on the logit scale, so
# it can move for two very different reasons: because the two surfaces really
# separate more, or merely because more of the plane sits in the saturated
# region where the logit link inflates differences. B-F tell these apart.
#
# Outputs:
#   output/06_premium_sensitivity/premium_lambda_sensitivity.csv
#   plots/06_premium_sensitivity/premium_lambda_sensitivity.pdf
#
# Runtime: ~1-2 min.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(mgcv)
})

ARGS    <- commandArgs(trailingOnly = TRUE)
SEEDING <- sub("^--seeding=", "", grep("^--seeding=", ARGS, value = TRUE))
if (length(SEEDING) == 0) SEEDING <- "random"

DATA_OUT  <- "output/06_premium_sensitivity"
PLOTS_OUT <- "plots/06_premium_sensitivity"
RES_CSV   <- file.path("output/05_unified_diffusion",
                       paste0("unified_premium_results_", SEEDING, ".csv"))

stopifnot(file.exists(RES_CSV))
raw <- read.csv(RES_CSV)
LAMBDAS <- sort(unique(raw$lambda))
message("Loaded ", nrow(raw), " rows; lambdas: ", paste(LAMBDAS, collapse = ", "))

raw <- raw |>
  mutate(topology = factor(topology, levels = c("GSS", "DP")),
         phi   = num_adopters / N_nodes,
         phi_r = num_adopted_rational / N_nodes,
         phi_s = num_adopted_social  / N_nodes,
         fail  = N_nodes - num_adopters)

# ---- cell means: one value per (topology, lambda, IUL, MSP) ------------------
cells <- raw |>
  group_by(topology, lambda, G = innovation_iul_Gamma, h = social_distance_h) |>
  summarise(phi = mean(phi), phi_r = mean(phi_r), phi_s = mean(phi_s),
            .groups = "drop")

wide <- cells |>
  pivot_wider(names_from = topology, values_from = c(phi, phi_r, phi_s)) |>
  mutate(gap = phi_GSS - phi_DP, gap_r = phi_r_GSS - phi_r_DP,
         gap_s = phi_s_GSS - phi_s_DP)

# ---- A. odds-scale premium, per lambda and per channel -----------------------
fit_beta <- function(df, yes, no) {
  m <- bam(cbind(df[[yes]], df[[no]]) ~ s(innovation_iul_Gamma, social_distance_h, k = 10)
           + topology, family = binomial(), data = df, discrete = TRUE)
  summary(m)$p.table["topologyDP", 1]
}

message("Fitting BAMs (total / rational / social) per lambda...")
betas <- lapply(LAMBDAS, function(l) {
  d <- filter(raw, lambda == l)
  d$no_r <- d$N_nodes - d$num_adopted_rational
  d$no_s <- d$N_nodes - d$num_adopted_social
  data.frame(lambda = l,
             beta_total    = fit_beta(d, "num_adopters", "fail"),
             beta_rational = fit_beta(d, "num_adopted_rational", "no_r"),
             beta_social   = fit_beta(d, "num_adopted_social",   "no_s"))
}) |> bind_rows()

# ---- B-F. link-free and structural diagnostics -------------------------------
diag <- wide |>
  group_by(lambda) |>
  summarise(
    # B. probability-scale premium
    mean_phi_GSS  = mean(phi_GSS),
    mean_phi_DP   = mean(phi_DP),
    mean_gap      = mean(gap),
    max_gap       = max(gap),
    # C. where the premium lives
    pct_plane_gap10 = 100 * mean(gap > 0.10),
    gap_in_band     = ifelse(any(gap > 0.10), mean(gap[gap > 0.10]), NA_real_),
    # D. channel split (probability scale)
    mean_gap_rational = mean(gap_r),
    mean_gap_social   = mean(gap_s),
    # F. saturation
    pct_GSS_sat = 100 * mean(phi_GSS > 0.95),
    pct_DP_sat  = 100 * mean(phi_DP  > 0.95),
    .groups = "drop")

# ---- E. tipping shift: first MSP where Phi crosses 0.5, per IUL --------------
tip <- cells |>
  group_by(topology, lambda, G) |>
  summarise(h_c = if (any(phi >= 0.5)) min(h[phi >= 0.5]) else NA_real_, .groups = "drop") |>
  pivot_wider(names_from = topology, values_from = h_c) |>
  filter(!is.na(GSS), !is.na(DP)) |>
  group_by(lambda) |>
  summarise(mean_MSPc_GSS = mean(GSS), mean_MSPc_DP = mean(DP),
            mean_MSPc_shift = mean(DP - GSS), .groups = "drop")

sens <- betas |> left_join(diag, by = "lambda") |> left_join(tip, by = "lambda") |>
  mutate(OR = exp(beta_total), odds_drop_pct = 100 * (exp(beta_total) - 1),
         d_beta_per_0.1 = c(NA, diff(beta_total)),
         d_gap_per_0.1  = c(NA, diff(mean_gap)))

write.csv(sens, file.path(DATA_OUT, paste0("premium_lambda_sensitivity_", SEEDING, ".csv")), row.names = FALSE)

# ------------------------------- report ---------------------------------------
line <- function(x) message(paste(rep(x, 78), collapse = ""))
line("="); message("A. ODDS-SCALE PREMIUM (BAM beta_DP) AND ITS DECILE STEPS"); line("=")
print(sens |> select(lambda, beta_total, OR, odds_drop_pct, d_beta_per_0.1) |>
        as.data.frame(), digits = 3, row.names = FALSE)

line("-"); message("B. PROBABILITY-SCALE PREMIUM (link-free)"); line("-")
print(sens |> select(lambda, mean_phi_GSS, mean_phi_DP, mean_gap, d_gap_per_0.1, max_gap) |>
        as.data.frame(), digits = 3, row.names = FALSE)

line("-"); message("C. WHERE THE PREMIUM LIVES + F. SATURATION"); line("-")
print(sens |> select(lambda, pct_plane_gap10, gap_in_band, pct_GSS_sat, pct_DP_sat) |>
        as.data.frame(), digits = 3, row.names = FALSE)

line("-"); message("D. CHANNEL DECOMPOSITION"); line("-")
print(sens |> select(lambda, beta_rational, beta_social, mean_gap_rational, mean_gap_social) |>
        as.data.frame(), digits = 3, row.names = FALSE)

line("-"); message("E. TIPPING SHIFT (mean MSP at which Phi crosses 0.5)"); line("-")
print(sens |> select(lambda, mean_MSPc_GSS, mean_MSPc_DP, mean_MSPc_shift) |>
        as.data.frame(), digits = 3, row.names = FALSE)

line("="); message("Legacy engine reference: beta_DP = -0.850 (OR 0.43, -57% odds)"); line("=")

# --------------------------------- plots --------------------------------------
p1 <- ggplot(sens, aes(lambda, beta_total)) +
  geom_hline(yintercept = -0.850, linetype = "22", color = "#B0161C") +
  annotate("text", x = min(LAMBDAS), y = -0.815, hjust = 0, size = 3.2,
           color = "#B0161C", label = "legacy engine (-0.850)") +
  geom_line(linewidth = 0.9, color = "#1E5AC8") +
  geom_point(size = 2.4, color = "#1E5AC8") +
  labs(title = "Premium on the odds scale", x = "λ", y = expression(beta[DP])) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

p2 <- ggplot(sens, aes(lambda, mean_gap)) +
  geom_line(linewidth = 0.9, color = "#007636") +
  geom_point(size = 2.4, color = "#007636") +
  labs(title = "Premium on the probability scale",
       x = "λ", y = expression(bar(Phi)[GSS] - bar(Phi)[DP])) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

p3 <- sens |> select(lambda, GSS = mean_phi_GSS, DP = mean_phi_DP) |>
  pivot_longer(-lambda) |>
  ggplot(aes(lambda, value, color = name, linetype = name)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  scale_color_manual(values = c(GSS = "#B0161C", DP = "#007636"), name = NULL) +
  scale_linetype_manual(values = c(GSS = "solid", DP = "22"), name = NULL) +
  labs(title = "Mean adoption over the whole plane", x = "λ", y = expression(bar(Phi))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

cairo_pdf(file.path(PLOTS_OUT, paste0("premium_lambda_sensitivity_", SEEDING, ".pdf")),
          width = 11, height = 3.6, onefile = TRUE)
print(p1 + p2 + p3 + plot_layout(nrow = 1))

# gap surfaces, one panel per lambda
gp <- ggplot(wide, aes(G, h, fill = gap)) +
  geom_tile() +
  scale_fill_gradient2(low = "#007636", mid = "white", high = "#B0161C",
                       midpoint = 0, limits = c(-1, 1), name = "Φ(GSS) − Φ(DP)") +
  facet_wrap(~ paste0("λ = ", lambda), nrow = 2) +
  labs(title = "Where the structural premium lives", x = "IUL (Γ)", y = "MSP (h)") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
print(gp)
dev.off()
message("\nSaved ", file.path(DATA_OUT, paste0("premium_lambda_sensitivity_", SEEDING, ".csv")), " and ",
        file.path(PLOTS_OUT, paste0("premium_lambda_sensitivity_", SEEDING, ".pdf")))
