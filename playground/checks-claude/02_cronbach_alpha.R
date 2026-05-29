# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cronbach's alpha for MUR constructs (non-destructive standalone)
# 02_cronbach_alpha.R
#
# Reads node attributes from one representative imputed network for GSS and ATP
# and computes internal-consistency (Cronbach's alpha) of the MUR item batteries.
# Does NOT overwrite any network file (unlike scripts/03_*_MUR_calculation.R).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressMessages({
  library(network)
  library(psych)
})

GSS_NET <- "data/02_GSS_network_ergm/GSS_net_sim_1000_001.rds"
ATP_NET <- "data/02_ATP_network_ergm/ATP_net_sim_1000_001.rds"
OUT_DIR <- "playground/checks-claude/data/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Item batteries -------------------------------------------------------------
# GSS collective-action items. Original GSS coding: 1=did in past year ...
# 4=would never. Recoded for propensity as 4 - x  (higher = more prone).
gss_items <- c("signdpet", "avoidbuy", "joindem", "attrally",
               "cntctgov", "polfunds", "usemedia", "interpol", "actlaw")

# ATP innovation items (binary). Anti-innovation items are reverse-scored so
# that 1 always means "pro-innovation".
atp_items_pro  <- c("metech_a", "metech_c", "metech_d")   # pro-innovation
atp_items_anti <- c("metech_b", "metech_e", "metech_f")   # anti-innovation (reverse)

extract_items <- function(net, items) {
  df <- as.data.frame(lapply(items, function(v) get.vertex.attribute(net, v)))
  names(df) <- items
  df
}

report_alpha <- function(label, item_df, n_items) {
  complete <- item_df[complete.cases(item_df), , drop = FALSE]
  a <- psych::alpha(complete, warnings = FALSE)
  raw_alpha <- a$total$raw_alpha
  std_alpha <- a$total$std.alpha
  cat("\n========== CRONBACH'S ALPHA:", label, "==========\n")
  cat("Items (", n_items, "):", paste(names(item_df), collapse = ", "), "\n")
  cat("Complete cases (nodes):", nrow(complete), "\n")
  cat(sprintf("Raw alpha       = %.4f\n", raw_alpha))
  cat(sprintf("Standardized    = %.4f\n", std_alpha))
  cat("Item-rest correlations:\n")
  print(round(a$item.stats[, c("r.drop")], 3))
  verdict <- if (raw_alpha >= 0.70) "PASS (>= 0.70)" else
             if (raw_alpha >= 0.60) "ACCEPTABLE (0.60-0.70)" else
             "LOW (< 0.60)"
  cat("Verdict:", verdict, "\n")
  cat("==================================================\n")
  data.frame(construct = label, n_items = n_items,
             n_cases = nrow(complete),
             raw_alpha = raw_alpha, std_alpha = std_alpha)
}

results <- list()

# --- GSS (collective action) ---
if (file.exists(GSS_NET)) {
  g <- readRDS(GSS_NET)
  gss_raw <- extract_items(g, gss_items)
  # Recode 1..4 -> 3..0 (higher = more prone to act). Out-of-range -> NA.
  gss_recoded <- as.data.frame(lapply(gss_raw, function(x) {
    x <- ifelse(x %in% 1:4, 4 - x, NA)
    x
  }))
  names(gss_recoded) <- gss_items
  results$gss <- report_alpha("GSS Collective-Action Propensity", gss_recoded, length(gss_items))
} else {
  cat("WARNING: GSS network not found at", GSS_NET, "\n")
}

# --- ATP (innovation) ---
if (file.exists(ATP_NET)) {
  a_net <- readRDS(ATP_NET)
  atp_raw <- extract_items(a_net, c(atp_items_pro, atp_items_anti))
  # Reverse-score anti-innovation items so 1 = pro-innovation throughout.
  atp_scored <- atp_raw
  for (v in atp_items_anti) atp_scored[[v]] <- 1 - atp_scored[[v]]
  results$atp <- report_alpha("ATP Innovation Propensity", atp_scored, ncol(atp_scored))
} else {
  cat("WARNING: ATP network not found at", ATP_NET, "\n")
}

if (length(results) > 0) {
  summary_df <- do.call(rbind, results)
  write.csv(summary_df, file.path(OUT_DIR, "cronbach_alpha_summary.csv"), row.names = FALSE)
  cat("\nSaved summary to", file.path(OUT_DIR, "cronbach_alpha_summary.csv"), "\n")
}
