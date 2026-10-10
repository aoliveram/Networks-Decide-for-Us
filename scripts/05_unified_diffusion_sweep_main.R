# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Runner: the unified diffusion sweep across ALL seeding strategies
# 05_unified_diffusion_sweep_main.R
#
# Calls scripts/05_unified_diffusion_sweep.R once per seeding strategy, each in
# its own R session, and then the sensitivity analysis for that strategy. Every
# strategy writes to its own checkpoint directory and its own result files, so
# nothing overwrites anything.
#
# Fully resumable: the sweep skips (topology, lambda) combinations whose
# checkpoint already exists, so re-running this script after an interruption
# only computes what is missing. The `random` strategy is already complete, so
# it is skipped in seconds (it still re-fits its BAM and re-draws its plots).
#
# Cost: ~50 min per strategy on 8 cores; ~3.5 h for the four that remain.
#
# Usage:
#   Rscript scripts/05_unified_diffusion_sweep_main.R                  # all five
#   Rscript scripts/05_unified_diffusion_sweep_main.R central marginal # a subset
#   NDFU_SKIP_SENSITIVITY=1 Rscript scripts/05_unified_diffusion_sweep_main.R
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ALL_SEEDINGS <- c("random", "central", "marginal", "closeness", "eigen")

args     <- commandArgs(trailingOnly = TRUE)
seedings <- if (length(args) > 0) args else ALL_SEEDINGS
unknown  <- setdiff(seedings, ALL_SEEDINGS)
if (length(unknown) > 0)
  stop("Unknown seeding strategy: ", paste(unknown, collapse = ", "),
       "\nValid: ", paste(ALL_SEEDINGS, collapse = ", "))

RUN_SENSITIVITY <- !nzchar(Sys.getenv("NDFU_SKIP_SENSITIVITY"))

rscript <- file.path(R.home("bin"), "Rscript")
banner  <- function(...) {
  message("\n", paste(rep("=", 78), collapse = ""))
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
  message(paste(rep("=", 78), collapse = ""))
}

t_all <- Sys.time()
for (s in seedings) {
  banner("SEEDING: ", toupper(s), "   (", match(s, seedings), "/", length(seedings), ")")

  st <- system2(rscript, c("scripts/05_unified_diffusion_sweep.R",
                           paste0("--seeding=", s)))
  if (st != 0) stop("Sweep failed for seeding '", s, "' (exit ", st, ")")

  if (RUN_SENSITIVITY) {
    message("\n--- lambda sensitivity for seeding '", s, "' ---")
    st <- system2(rscript, c("scripts/06_premium_sensitivity.R",
                             paste0("--seeding=", s)))
    if (st != 0) warning("Sensitivity analysis failed for seeding '", s,
                         "' (exit ", st, ") — the sweep results are still saved.")
  }
}

banner("ALL DONE in ",
       round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min")
message("Results:  output/05_unified_diffusion/unified_premium_{results,bam}_<seeding>.csv")
message("Premium:  output/06_premium_sensitivity/premium_lambda_sensitivity_<seeding>.csv")
message("Figures:  plots/05_unified_diffusion/ and plots/06_premium_sensitivity/")
