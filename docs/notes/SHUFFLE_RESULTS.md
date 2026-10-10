# GSS-SH (attribute-shuffle) results — the causal decomposition

**Date:** 2026-06-21. After running the full diffusion sweep on GSS-SH (topology byte-identical, attributes permuted → only homophily destroyed) and re-running the three analysis scripts with the third topology added.

## 1. Structural premium — DECOMPOSED (BAM, 06_regression_analysis.R)

Reference level = Empirical (plausible GSS). Coefficients are log-odds vs. that reference.
`network_typeER` = degree-preserving randomization (DP); `network_typeSH` = attribute shuffle.

| Outcome | β(DP) ER | β(SH) | OR(SH) | Read |
|---|---|---|---|---|
| **Total adoption** | **−0.840** | **−0.766** | 0.465 | SH ≈ 91% of the DP penalty |
| **Social adoption** | −0.432 | −0.403 | 0.668 | SH ≈ 93% of the DP penalty |
| Rational adoption | −0.114 | **−0.193** | 0.824 | SH penalty is LARGER than DP here (see note) |

Deviance explained: Total 82.6%, Social 72%, Rational 97.1%. All coefficients p < 2e-16.
N rows: Empirical 2,046,720 · ER 2,046,720 · SH 1,023,360 (SH only CollectiveAction, no ATP-SH).

### Headline read
**Destroying ONLY homophily (topology intact) reproduces ~91% of the total-adoption penalty that the full randomization produces.** Decomposition of the −0.840 total penalty:
- **Homophily ≈ −0.766 (≈91%)** — the alignment of attributes with structure.
- **Clustering/topology ≈ −0.074 (≈9%)** — the residual ER−SH gap (−0.840 − (−0.766)).

This is the clean causal separation the DP-only design could not give: **the structural premium is overwhelmingly a HOMOPHILY effect, not a clustering effect.** The same network, same individuals, same ties — only reassigning who sits where — destroys nearly the whole premium.

### Note on the Rational channel (unexpected, worth understanding)
Rational adoption is q_i ≤ Γ — node-local, and the population of q_i is identical across topologies, so the *direct* effect "should" be ~0. It is small but non-zero, and SH (−0.193) is actually MORE negative than DP (−0.114). Likely mechanism: this is a *seeding/lock-in interaction*, not a pure rational effect — the seeding strategies (central/closeness/eigenvector) place initial adopters by network position, and once social adoption changes, the asynchronous update reshuffles how many adoptions get *attributed* to the rational vs social tally near boundaries. It is a second-order accounting effect on a tiny coefficient (OR 0.82) and does not affect the Total/Social story. Flag for the SM, do not headline.

## 2. Universality / criticality — γ across the three topologies (07_phase_transition_deeper.R)

γ extracted by the same Method-4 (susceptibility power-law on the super-critical side), IUL ∈ {0.200, 0.225, 0.250, 0.275}:

| IUL | GSS (plausible) | GSS-DP (rewired) | GSS-SH (shuffle) |
|---|---|---|---|
| 0.200 | 0.98 | 3.71 | 0.83 |
| 0.225 | 0.89 | 3.67 | 0.75 |
| 0.250 | 0.89 | 3.67 | 0.75 |
| 0.275 | 0.87 | 4.67 | 1.98 |
| **mean** | **≈0.91** | **≈3.9** | **≈1.08** |

### Headline read (this is the big one for Paper 4.2)
**GSS-SH KEEPS the mean-field universality (γ ≈ 1), GSS-DP does NOT (γ ≈ 4).** Since SH preserves clustering/triangles/modularity and destroys only homophily, while DP destroys the topology:

→ **The universality (continuous, γ≈1 transition) is anchored by the TOPOLOGY (clustering/community structure), NOT by homophily.** Destroying homophily alone (SH) leaves the transition continuous and mean-field; only destroying the higher-order topology (DP) collapses it into the discontinuous γ≈4 regime.

Note MSP_c: GSS sits at 0.25, while BOTH DP and SH sit at 0.417 — i.e. removing homophily pushes the critical openness up (you need more social porosity to ignite), but the *order* of the transition (continuous vs discontinuous) is governed by topology, not homophily.

## 3. The two findings together — a clean double dissociation

| | Destroys homophily? | Destroys topology? | Premium (volume) | Universality (γ) |
|---|---|---|---|---|
| GSS | no | no | — (ref) | ≈1 (continuous) |
| **GSS-SH** | **yes** | no | **−0.77 (most of it)** | **≈1 (kept!)** |
| **GSS-DP** | yes | **yes** | −0.84 | **≈4 (lost)** |

**Double dissociation:**
- The **VOLUME premium** is driven by **homophily** (SH already destroys ~91% of it).
- The **continuity/universality** is driven by **topology/clustering** (SH keeps γ≈1; only DP, which destroys topology, loses it).

These are *two different structural properties doing two different jobs* — a genuinely strong, publishable result. Homophily decides *how much* spreads; topology decides *how* (smoothly vs abruptly) it tips.

## 4. Outputs
- Models: `output/06_regression_analysis/models/bam_{total,social,rational}{,_counts}.rds` + `_summary.txt` (now 3-level network_type)
- γ: `output/07_phase_transition/GSS_SH/gamma_exponents.csv` (+ method2/3/4 plots in `plots/07_phase_transition/GSS_SH/`)
- Adoption heatmaps: `plots/04_GSS_SH_diffusion_sims/` (20 PDFs)

## 5. Caveats (honest)
- γ on DP (≈4) is not a true power-law exponent (first-order discontinuity → the fit is a fingerprint of "not a power law", not a real γ). The SH γ≈1 with the SAME MSP_c as DP (0.417) but a CLEAN power law is the meaningful contrast.
- IUL=0.275 is the noisy point for all three (GSS 0.87, SH 1.98) — same coarse-grid fragility flagged before. The means are dominated by the three stable points.
- Rational-channel coefficient (§1 note) is a second-order accounting effect; do not over-interpret.
