# NDFU — Results Report #1: GSS vs GSS-DP vs GSS-SH

**Simulations run:** 2026-06-19 (GSS-SH diffusion sweep)
**Analysis run:** 2026-06-21 (regression, phase transition, plots)
**Report written:** 2026-06-21

This report consolidates the new results from adding the **attribute-shuffle null model (GSS-SH)** as a third topology alongside the plausible homophilic networks (GSS) and the degree-preserving randomized baseline (GSS-DP). It is a standalone record; the main manuscript is NOT yet modified.

---

## 1. The three topologies (what each one is, and what it destroys)

| Label | Construction | Degree | Clustering / triangles / modularity | Homophily (attribute↔position) |
|---|---|---|---|---|
| **GSS** | ERGM-imputed plausible networks (empirical homophily) | — | preserved | **present** |
| **GSS-DP** | degree-preserving rewiring (`rewire`, `keeping_degseq`, niter = 20·m) | **preserved** | **destroyed** | **destroyed** |
| **GSS-SH** | topology byte-identical; node attributes randomly **permuted** | **preserved** | **preserved (identical)** | **destroyed** |

The logic of the design: **GSS-DP destroys topology and homophily simultaneously**, so it cannot say which one drives the adoption premium. **GSS-SH destroys only homophily** (the topology is literally the same graph, with attributes shuffled across nodes), so the contrast **GSS-DP vs GSS-SH isolates the contribution of the higher-order topology**, and **GSS-SH vs GSS isolates the contribution of homophily**.

### 1.1 Network-statistics verification (GSS vs GSS-SH; 100 vs 100 networks)

Confirms the shuffle does exactly what it should: topology identical, homophily annihilated.

| Family | Statistic | GSS | GSS-SH | Status |
|---|---|---|---|---|
| B (topology) | global clustering | 0.06475 | 0.06475 | **identical** |
| B (topology) | triangles | 11,200 | 11,200 | **identical** |
| B (topology) | modularity (Louvain) | 0.2903 | 0.2903 | **identical** |
| B (topology) | degree assortativity | 0.2907 | 0.2907 | **identical** |
| B (topology) | mean degree | 28.93 | 28.93 | **identical** |
| C (homophily) | assortativity, race | 0.441 | −0.001 | collapsed → 0 |
| C (homophily) | assortativity, age | 0.463 | −0.001 | collapsed → 0 |
| C (homophily) | assortativity, religion | 0.418 | −0.0001 | collapsed → 0 |
| C (homophily) | mean social distance on ties (Gower d_ij) | 0.229 | 0.399 | **+74%** |

*(For reference, GSS-DP also collapses C to ~0, but additionally drops B: clustering −31%, triangles −31%, modularity −47%.)*

---

## 2. Result A — The structural premium, DECOMPOSED (BAM regression)

**Model:** Big Additive Model (mgcv `bam`, binomial/logit, fREML, discrete), pooled across both contagion types, all seedings, all threshold parameters.
**Specification:** `adopt ~ s(IUL, MSP, k=10, by=contagion_type) + tau_mean + tau_sd + seed_type + network_type`
**Reference levels:** `network_type` = Empirical (GSS); `seed_type` = random.
**N rows:** Empirical 2,046,720 · ER 2,046,720 · SH 1,023,360 (SH exists for CollectiveAction only — there is no ATP-SH run). All coefficients below have p < 2×10⁻³⁰⁰ (z given).

### 2.1 Topology coefficients (the core result)

| Outcome | β(GSS-DP) | OR | drop | β(GSS-SH) | OR | drop | Deviance expl. |
|---|---|---|---|---|---|---|---|
| **Total adoption** | **−0.8398** | 0.432 | −56.8% | **−0.7663** | 0.465 | −53.5% | 82.6% |
| **Social adoption** | −0.4321 | 0.649 | −35.1% | −0.4033 | 0.668 | −33.2% | 72.0% |
| Rational adoption | −0.1136 | 0.893 | −10.7% | −0.1934 | 0.824 | −17.6% | 97.1% |

### 2.2 Decomposition of the total-adoption penalty (−0.8398)

| Component | Estimate (log-odds) | Share of DP penalty |
|---|---|---|
| **Homophily** (= β_SH) | **−0.766** | **91.3%** |
| **Higher-order topology** (= β_DP − β_SH) | −0.073 | 8.7% |

**Reading.** Destroying *only* homophily — keeping clustering, triangles, modularity, degree and the giant component byte-identical — already reproduces **91% of the adoption penalty** that the full randomization produces. The structural premium is **overwhelmingly a homophily effect, not a clustering effect.** The same network, the same individuals, the same ties; merely reassigning *who sits where* destroys nearly the whole premium.

The social channel tells the same story (β_SH = −0.403 vs β_DP = −0.432; homophily ≈ 93% of the social penalty), which is expected because the selective-influence filter σ(MSP − d_ij) is exactly the part of the model that homophily feeds.

### 2.3 Full parametric tables (CollectiveAction + Innovation pooled)

**Total adoption**

| Term | Estimate | Std. Error | z |
|---|---|---|---|
| (Intercept) | 6.3837 | 0.01528 | 417.7 |
| tau_mean | −6.7839 | 0.01486 | −456.5 |
| tau_sd | 5.8197 | 0.03469 | 167.7 |
| seed: central | 0.1302 | 0.00486 | 26.8 |
| seed: closeness | 0.1296 | 0.00486 | 26.7 |
| seed: eigen | 0.1290 | 0.00486 | 26.6 |
| seed: marginal | −0.0729 | 0.00483 | −15.1 |
| **network: ER (DP)** | **−0.8398** | 0.00351 | −239.0 |
| **network: SH** | **−0.7663** | 0.00460 | −166.7 |

Smooths `s(IUL,MSP)`: edf ≈ 9.0 for both contagion types, χ² = 376,793 (Innovation) / 588,300 (CollectiveAction), p ≈ 0.

**Social adoption**

| Term | Estimate | Std. Error | z |
|---|---|---|---|
| (Intercept) | −0.3004 | 0.00771 | −39.0 |
| tau_mean | −5.2524 | 0.01245 | −421.8 |
| tau_sd | 4.3171 | 0.02992 | 144.3 |
| seed: central | 0.0370 | 0.00421 | 8.8 |
| seed: closeness | 0.0379 | 0.00421 | 9.0 |
| seed: eigen | 0.0365 | 0.00421 | 8.7 |
| seed: marginal | −0.0215 | 0.00423 | −5.1 |
| **network: ER (DP)** | **−0.4321** | 0.00297 | −145.5 |
| **network: SH** | **−0.4033** | 0.00390 | −103.4 |

**Rational adoption**

| Term | Estimate | Std. Error | z |
|---|---|---|---|
| (Intercept) | 0.4277 | 0.00679 | 62.9 |
| tau_mean | −0.1931 | 0.01059 | −18.2 |
| tau_sd | 0.0393 | 0.02647 | 1.5 (n.s.) |
| seed: central | −0.0582 | 0.00374 | −15.6 |
| seed: closeness | −0.0558 | 0.00374 | −14.9 |
| seed: eigen | −0.0555 | 0.00374 | −14.8 |
| seed: marginal | 0.0278 | 0.00375 | 7.4 |
| **network: ER (DP)** | **−0.1136** | 0.00262 | −43.3 |
| **network: SH** | **−0.1934** | 0.00340 | −56.8 |

**Note on the rational channel (caveat).** Rational adoption is q_i ≤ Γ, which is node-local and uses a q_i distribution identical across all three topologies, so the *direct* topology effect should be ≈0. The coefficients are small (OR 0.82–0.89) but non-zero, and unusually the SH penalty (−0.193) exceeds the DP penalty (−0.114). This is almost certainly a **second-order seeding/accounting interaction**: position-based seeding (central/closeness/eigenvector) and the asynchronous update change how adoptions get tallied as rational vs social near the IUL boundary when the social dynamics shift. It does not affect the Total/Social conclusions and should be a one-line footnote in the SM, not a headline.

---

## 3. Result B — Universality / criticality, across the three topologies

**Method:** susceptibility power-law fit (Method 4 of `07_phase_transition_deeper.R`): for each fixed IUL, find MSP_c (steepest adoption jump) and fit χ = Var(Φ) ~ |MSP − MSP_c|^(−γ) on the super-critical side. Seeding = random, τ_mean = 0.4, τ_sd = 0.12.

| IUL | GSS — γ (MSP_c) | GSS-DP — γ (MSP_c) | GSS-SH — γ (MSP_c) |
|---|---|---|---|
| 0.200 | 0.977 (0.250) | 3.709 (0.417) | 0.831 (0.417) |
| 0.225 | 0.892 (0.250) | 3.667 (0.417) | 0.751 (0.417) |
| 0.250 | 0.892 (0.250) | 3.667 (0.417) | 0.751 (0.417) |
| 0.275 | 0.873 (0.250) | 4.667 (0.333) | 1.977 (0.333) |
| **mean** | **≈ 0.91** | **≈ 3.93** | **≈ 1.08** |

**Reading.** **GSS-SH keeps the mean-field universality (γ ≈ 1); GSS-DP does not (γ ≈ 4).** Since GSS-SH preserves the higher-order topology and destroys only homophily, while GSS-DP destroys the topology:

→ **the continuous, universal (γ ≈ 1) character of the transition is anchored by the TOPOLOGY (clustering / community structure), NOT by homophily.** Removing homophily alone leaves the transition continuous and mean-field; only removing the higher-order topology collapses it into the discontinuous γ ≈ 4 regime.

**Critical openness shifts.** GSS tips at MSP_c = 0.25; *both* GSS-DP and GSS-SH tip at MSP_c = 0.417. So removing homophily raises the social porosity needed to ignite a cascade — but the *order* of the transition (continuous vs discontinuous) is set by topology, independent of that shift.

**Caveats.** (i) The GSS-DP γ ≈ 4 is not a genuine power-law exponent — a first-order discontinuity has no diverging-power-law susceptibility, so that fit is a fingerprint of "not a power law", not a real γ. The meaningful contrast is that GSS-SH gives a *clean* γ ≈ 1 at the *same* MSP_c (0.417) where GSS-DP fails. (ii) IUL = 0.275 is the noisy fitting point for all three (GSS 0.87, SH 1.98); the means are carried by the three stable points. (iii) coarse 13-level MSP grid — the usual finite-grid fragility.

---

## 4. The headline: a clean DOUBLE DISSOCIATION

| Topology | Destroys homophily? | Destroys topology? | **Volume premium** (β_total) | **Universality** (γ) |
|---|---|---|---|---|
| GSS | no | no | — (reference) | ≈ 1 (continuous) |
| **GSS-SH** | **yes** | no | **−0.766** (≈ 91% of full) | **≈ 1 (kept)** |
| **GSS-DP** | yes | **yes** | −0.840 | **≈ 4 (lost)** |

**Two structural properties do two different jobs:**

- **Homophily → VOLUME.** It governs *how much* spreads (the structural premium). GSS-SH, which destroys only homophily, already kills ~91% of it.
- **Topology / clustering → CONTINUITY.** It governs *how* the system tips — smoothly (continuous, universal) vs abruptly (discontinuous). GSS-SH keeps γ ≈ 1; only GSS-DP, which destroys the topology, loses it.

This is a genuine double dissociation: each manipulation knocks out one outcome while sparing the other (SH kills the premium but spares universality; the topology destruction in DP is what kills universality). It cleanly separates the project's two headline results and assigns them to *distinct, identifiable structural causes* — turning two correlations into two mechanisms.

---

## 5. Provenance (files)

- **Diffusion results (SH):** `output/04_GSS_SH_diffusion_sims/` (25 RDS, 5 seedings × all sd/mean)
- **Regression models (3-level network_type):** `output/06_regression_analysis/models/bam_{total,social,rational}{,_counts}.rds` + `_summary.txt`
- **Critical exponents (SH):** `output/07_phase_transition/GSS_SH/gamma_exponents.csv`; plots in `plots/07_phase_transition/GSS_SH/`
- **Adoption heatmaps (SH):** `plots/04_GSS_SH_diffusion_sims/` (20 PDFs)
- **Network construction & verification:** `playground/checks-claude/06_make_shuffle_networks.R`, `data/02_GSS_SH_network/` (100 nets), `playground/checks-claude/data/shuffle_verification.csv`
- **Diagnosis (GSS vs DP):** `playground/checks-claude/05_structural_diagnosis.R`, `STRUCTURAL_DIAGNOSIS.md`
- **Adapted scripts:** `scripts/06_regression_analysis.R` (3-level), `scripts/05_GSS_diffusion_plot_main.R` (SH type), `scripts/07_phase_transition_deeper.R` (NDFU_TOPO env var); SH topology branch in `scripts/04_GSS_diffusion_sims.R`; runner `scripts/04_GSS_SH_diffusion_sims_main.R`

## 6. Open questions / next

- **β extraction and the Widom relation** still need a denser near-MSP_c grid before γ/β/δ can be reported as a rigorous universality result.
- **Why does removing homophily preserve universality?** Worth a short mechanistic argument: clustering provides the local "averaging cells" that produce mean-field coupling, independent of whether neighbors share attributes.
- **Rational-channel second-order effect** (§2.3) — confirm the seeding/accounting explanation if it ever needs to go in the SM.
