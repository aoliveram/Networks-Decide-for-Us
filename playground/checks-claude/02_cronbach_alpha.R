# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cronbach's alpha for MUR constructs — from ORIGINAL survey data
# 02_cronbach_alpha.R
#
# IMPORTANT: alpha is computed on the ORIGINAL respondents in data/, NOT on the
# imputed N=1000 networks (which resample respondents to fill nodes and would
# distort the reliability estimate). Sources:
#   - GSS collective action: data/02_GSS_network_ergm/GSS_2004_NORC.dta  (2812 resp.)
#   - ATP innovation:        data/01_ATP_GSS_imputation/ATP_W3_W4.rds     (3620 resp.)
#
# RESPONSE DIRECTIONS are handled explicitly (see below). Reverse-coding only
# affects alpha when items run in OPPOSITE directions (the ATP pro/anti case);
# alpha is invariant to a uniform linear recoding (the GSS case).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(haven)
  library(psych)
})

OUT_DIR <- "playground/checks-claude/data/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

report_alpha <- function(label, item_df, note = "") {
  complete <- item_df[complete.cases(item_df), , drop = FALSE]
  a <- psych::alpha(complete, warnings = FALSE)
  cat("\n========== CRONBACH'S ALPHA:", label, "==========\n")
  if (nzchar(note)) cat(note, "\n")
  cat("Items (", ncol(item_df), "):", paste(names(item_df), collapse = ", "), "\n")
  cat("Complete-case respondents:", nrow(complete), "\n")
  cat(sprintf("Raw alpha        = %.4f\n", a$total$raw_alpha))
  cat(sprintf("Standardized     = %.4f\n", a$total$std.alpha))
  cat("Item-rest correlations (r.drop):\n")
  ir <- round(a$item.stats[, "r.drop"], 3); names(ir) <- names(item_df); print(ir)
  verdict <- if (a$total$raw_alpha >= 0.70) "GOOD (>= 0.70)" else
             if (a$total$raw_alpha >= 0.60) "ACCEPTABLE (0.60-0.70)" else "LOW (< 0.60)"
  cat("Verdict:", verdict, "\n==================================================\n")
  data.frame(construct = label, source = "original survey data (data/)",
             n_items = ncol(item_df), n_respondents = nrow(complete),
             raw_alpha = a$total$raw_alpha, std_alpha = a$total$std.alpha)
}

results <- list()

# ---------------------------------------------------------------------------
# (1) GSS collective-action propensity  -- GSS_2004_NORC.dta
#     9 items, all SAME direction: 1 = "did in the past year" ... 4 = "would
#     never". Construct = propensity to act, so recode 4 - x (1->3 ... 4->0).
#     All items point the same way => no reverse items; alpha is identical
#     with or without this uniform recoding (kept for construct fidelity).
# ---------------------------------------------------------------------------
GSS_DTA <- "data/02_GSS_network_ergm/GSS_2004_NORC.dta"
gss_items <- c("signdpet", "avoidbuy", "joindem", "attrally",
               "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")
if (file.exists(GSS_DTA)) {
  g <- read_dta(GSS_DTA)
  nm_map <- setNames(names(g), tolower(names(g)))
  gss_df <- as.data.frame(lapply(gss_items, function(it) {
    x <- as.numeric(g[[nm_map[[it]]]])
    ifelse(x %in% 1:4, 4 - x, NA)          # uniform recode, same direction
  }))
  names(gss_df) <- gss_items
  results$gss <- report_alpha(
    "GSS Collective-Action Propensity", gss_df,
    note = "Source: GSS_2004_NORC.dta | direction: all items 1='past year'->high; recode 4-x.")
} else cat("WARNING: missing", GSS_DTA, "\n")

# ---------------------------------------------------------------------------
# (2) ATP innovation propensity  -- ATP_W3_W4.rds
#     6 binary items. Directions DIFFER:
#       pro-innovation : METECH_A, _C, _D  (1 = pro)
#       anti-innovation: METECH_B, _E, _F  (1 = anti)  -> REVERSE so 1 = pro.
#     Raw file carries 99 sentinels / NA; keep only respondents with all 0/1.
# ---------------------------------------------------------------------------
ATP_RDS <- "data/01_ATP_GSS_imputation/ATP_W3_W4.rds"
atp_pro  <- c("METECH_A", "METECH_C", "METECH_D")
atp_anti <- c("METECH_B", "METECH_E", "METECH_F")
if (file.exists(ATP_RDS)) {
  d <- readRDS(ATP_RDS)
  atp_df <- d[, c(atp_pro, atp_anti)]
  keep <- apply(atp_df, 1, function(r) all(r %in% c(0, 1)))   # drop 99/NA
  atp_df <- atp_df[keep, , drop = FALSE]
  for (v in atp_anti) atp_df[[v]] <- 1 - atp_df[[v]]          # reverse-score
  results$atp <- report_alpha(
    "ATP Innovation Propensity", atp_df,
    note = "Source: ATP_W3_W4.rds | direction: A/C/D pro (1=pro), B/E/F anti reversed to 1=pro.")
} else cat("WARNING: missing", ATP_RDS, "\n")

if (length(results)) {
  summ <- do.call(rbind, results)
  write.csv(summ, file.path(OUT_DIR, "cronbach_alpha_summary.csv"), row.names = FALSE)
  cat("\nSaved", file.path(OUT_DIR, "cronbach_alpha_summary.csv"), "\n")
}
