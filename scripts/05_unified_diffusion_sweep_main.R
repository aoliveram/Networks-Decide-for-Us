# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Runner: the whole relaunch, from networks to diagnostics
# 05_unified_diffusion_sweep_main.R
#
# Each step runs in its own R session; every step is resumable or cheap:
#   1. 02_network_calibrated.R  per survey (skipped if the networks exist)
#   2. 05_engine_check.R        once (exactness of the engine; ~1 min) and
#      05_netdiffuser_reference.R once (same adopters as netdiffuseR; ~1 min)
#   3. 05 sweep                 per survey x level x seeding, all exposure conditions
#                               (checkpoints per condition and network: an interrupted
#                               run resumes where it stopped)
#   4. 06                       per set and seeding (lambda sensitivity)
#   5. 06b                      per set (joint BAM + run bootstrap, B = 1000)
#   6. 06c                      per set (exposure axis)
#   7. 06d                      per survey, with both levels (x2 - x1 against the same DP twin)
#   8. 07                       per survey (structural diagnosis: x2, x1, DP)
# Sets: <SURVEY>_<LEVEL>, SURVEY in GSS, ATP; LEVEL in x2 (main), x1.
#
# Cost on 8 performance cores (2026-11-05: 99 min in all): ~1.2 min per survey x
# level x seeding with all 9 conditions (~25 min for 2 x 2 x 5); 06b ~17 min per set
# (B = 1000); the rest a few minutes.
#
# Usage (repo root; log with 2>&1 | tee logs/05_main.log):
#   Rscript scripts/05_unified_diffusion_sweep_main.R
#   Rscript scripts/05_unified_diffusion_sweep_main.R --surveys=GSS --levels=x2 --seedings=random,central
#   Rscript scripts/05_unified_diffusion_sweep_main.R --skip-analysis      # only steps 1-3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ARGS <- commandArgs(trailingOnly = TRUE)
arg_list <- function(k, d) { a <- grep(paste0("^--", k, "="), ARGS, value = TRUE)
  if (length(a)) strsplit(sub(".*=", "", a[1]), ",")[[1]] else d }
SURVEYS  <- toupper(arg_list("surveys", c("GSS", "ATP")))
LEVELS   <- tolower(arg_list("levels", c("x2", "x1")))
SEEDINGS <- arg_list("seedings", c("random", "central", "marginal", "closeness", "eigen"))
stopifnot(all(SURVEYS %in% c("GSS", "ATP")), all(LEVELS %in% c("x1", "x2")),
          all(SEEDINGS %in% c("random", "central", "marginal", "closeness", "eigen")))
SKIP_ANALYSIS <- "--skip-analysis" %in% ARGS
ALL_FIVE <- setequal(SEEDINGS, c("random", "central", "marginal", "closeness", "eigen"))

rscript <- file.path(R.home("bin"), "Rscript")
banner  <- function(...) {
  message("\n", paste(rep("=", 78), collapse = ""))
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
  message(paste(rep("=", 78), collapse = ""))
}
run <- function(script, args = character(0), fatal = TRUE) {
  st <- system2(rscript, c(file.path("scripts", script), args))
  if (st != 0) { msg <- paste0(script, " ", paste(args, collapse = " "), " failed (exit ", st, ")")
                 if (fatal) stop(msg) else warning(msg) }
  invisible(st)
}
t_all <- Sys.time()

for (s in SURVEYS) {
  if (!all(file.exists(file.path("data/02_calibrated_networks", paste0(s, c("_x1.rds", "_x2.rds")))))) {
    banner("NETWORKS: ", s); run("02_network_calibrated.R", paste0("--survey=", s))
  }
}
banner("ENGINE CHECK"); run("05_engine_check.R")
banner("netdiffuseR REFERENCE"); run("05_netdiffuser_reference.R")

for (s in SURVEYS) for (l in LEVELS) for (sd in SEEDINGS) {
  banner("SWEEP ", s, "_", l, " | seeding ", toupper(sd))
  run("05_unified_diffusion_sweep.R", c(paste0("--survey=", s), paste0("--level=", l), paste0("--seeding=", sd)))
}

if (!SKIP_ANALYSIS) {
  for (s in SURVEYS) for (l in LEVELS) {
    set <- paste0(s, "_", l)
    for (sd in SEEDINGS) { banner("06 ", set, " | ", sd); run("06_premium_sensitivity.R", c(paste0("--set=", set), paste0("--seeding=", sd)), fatal = FALSE) }
    if (ALL_FIVE) { banner("06b ", set); run("06b_premium_bootstrap.R", paste0("--set=", set), fatal = FALSE) }
    banner("06c ", set); run("06c_exposure_axis.R", paste0("--set=", set), fatal = FALSE)
  }
  if (all(c("x1", "x2") %in% LEVELS))
    for (s in SURVEYS) { banner("06d ", s); run("06d_closure_contrast.R", paste0("--survey=", s), fatal = FALSE) }
  for (s in SURVEYS) { banner("07 ", s); run("07_structural_diagnosis.R", paste0("--survey=", s), fatal = FALSE) }
}

banner("ALL DONE in ", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), " min")
message("Sweep:    output/05_unified_diffusion/<set>/  (unified_premium_results_<seeding>.csv, exposure_axis_<seeding>.rds)")
message("Premium:  output/06_premium_sensitivity/<set>/, output/06b_premium_bootstrap/<set>/, output/06c_exposure_axis/<set>/,")
message("          output/06d_closure_contrast/")
message("Figures:  plots/05_unified_diffusion/<set>/, plots/06*/<set>/, plots/07_structural_diagnosis/<survey>/")
