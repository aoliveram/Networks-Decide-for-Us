# Universality Checks for NDFU Critical Exponent γ

This folder contains comprehensive checks to validate the universality claim (γ ≈ 1) in the main NDFU manuscript.

## Motivation

The paper claims that empirical (plausible) networks exhibit **Mean-Field universality** (γ ≈ 0.977 ≈ 1), while degree-preserving randomized (DP) baselines do not (γ ≈ 3.7–4.7). This folder documents robustness and identifies potential issues.

## Scripts and Checks

### 00_exponent_stability_checks.R
**Primary check:** Compares γ_GSS vs. γ_DP directly across all tested IUL values.

**Key outputs:**
- `gamma_comparison_gss_vs_dp.csv`: Side-by-side γ values for GSS and DP networks.
- `gamma_comparison.pdf`: Visualization of separation.
- Summary statistics (mean γ, SD, ratio).

**Interpretation:**
- If γ_GSS ≈ 1 ± 0.2 and γ_DP > 2.5, universality claim is supported.
- If they overlap, claim is weakened.

---

### 01_gamma_stability_across_msp.R
**Question:** For a fixed IUL, is γ stable as we slide around the critical MSP?

**Key outputs:**
- `gamma_stability_gss_iul*.csv`: γ computed in sliding windows of MSP.
- `gamma_stability_gss.pdf`, `gamma_stability_dp.pdf`: Line plots showing γ vs. MSP window center.

**Interpretation:**
- True universality should show γ ≈ constant in the critical regime.
- High variance in γ suggests the exponent is unstable or the critical regime is narrow.

---

### 02_averaging_schemes.R (To be written)
**Question:** If we average γ over multiple (IUL, MSP) pairs, does γ_avg converge to 1?

**Proposed outputs:**
- γ_avg_dense: Average γ over 8–10 IUL values × all MSP values.
- γ_avg_critical: Average γ only over IUL/MSP pairs near criticality.
- Comparison table: raw γ values vs. averages.

**Interpretation:**
- Averaging smooths out outliers like the IUL=0.275 case (γ=1.539).
- If γ_avg_dense ≈ 1.0 ± 0.15, universality is robust.

---

### 03_finite_size_scaling.R (To be written)
**Question:** Does γ converge to 1 as network size N increases?

**Required data:** Rerun simulations at N=500, 1000, 2000.

**Key outputs:**
- Table: γ(N) for each N.
- Plot: γ vs. log(N), with fit to see if γ → 1 as N → ∞.

**Interpretation:**
- Finite-size effects at N=1000 could account for γ ≈ 0.977 < 1.
- If γ(N→∞) → 1, Mean-Field universality is confirmed.

---

### 04_fit_quality_assessment.R (To be written)
**Question:** Are the log-log power-law fits actually good?

**Proposed outputs:**
- Table: (IUL, MSP_c, γ, r²) for each exponent estimate.
- Histogram of r² values.
- Scatter: r² vs. γ (to detect if bad fits give suspicious γ values).

**Interpretation:**
- If r² < 0.85 for most fits, γ estimates are unreliable.
- If r² > 0.90, power-law is justified.

---

### 05_bootstrap_confidence_intervals.R (To be written)
**Question:** What is the uncertainty in γ?

**Proposed approach:**
- Resample simulation runs (with replacement) for each (IUL, MSP, topology) combo.
- Recompute γ for each resample.
- Report 95% CI.

**Outputs:**
- Table: γ ± 95% CI for GSS and DP.
- Plot: Forest plot comparing intervals.

**Interpretation:**
- If CI(γ_GSS) = [0.85, 1.15] and CI(γ_DP) = [3.2, 5.0], separation is robust.
- If CI's overlap, claim is weakened.

---

## Running the Checks

```r
# Run all available checks:
source("00_exponent_stability_checks.R")
source("01_gamma_stability_across_msp.R")
# source("02_averaging_schemes.R")  # Once written
# source("03_finite_size_scaling.R")  # Requires new data
# source("04_fit_quality_assessment.R")  # Once written
# source("05_bootstrap_confidence_intervals.R")  # Once written
```

Or use the master script (when available):
```r
source("_run_all_checks.R")
```

---

## Expected Outputs

**If universality is robust:**
- γ_GSS ≈ 1.0 ± 0.2 across all checks.
- γ_DP ≈ 4.0 ± 1.0, clearly separated.
- r² > 0.90 for most GSS fits.
- γ_GSS is stable across MSP windows.

**If universality is fragile:**
- γ_GSS varies wildly (e.g., 0.8–1.5) depending on IUL or window.
- γ_DP sometimes dips below 2.0, overlapping with GSS.
- r² is poor, suggesting power-law is not appropriate.
- Confidence intervals are wide and overlapping.

---

## Notes

- All scripts assume you have run the main diffusion simulations (04_GSS_diffusion_sims.R, etc.) and the phase transition analysis (07_phase_transition_deeper.R).
- Output files go to `playground/checks-claude/plots/` and `playground/checks-claude/data/`.
- Adapt paths and filenames as needed if your directory structure differs.

---

**Created by:** Claude (AI assistant for Aníbal Olivera)  
**Date:** May 2026  
**Purpose:** Robustness validation of NDFU universality claims before manuscript resubmission.
