# Structural Diagnosis: what does GSS → GSS-DP destroy? (Paper 4.1.1)

**Script:** `05_structural_diagnosis.R` · **Data:** 100 GSS plausible nets + 100 DP counterparts (same `rewire(keeping_degseq)` as the diffusion sims). 100 vs 100 → Mann-Whitney U + Cohen's d.

## Results (all p = 2.6e-34 except the controls)

| Family | Statistic | GSS | DP | % change | Cohen's d |
|---|---|---|---|---|---|
| **A control** | mean degree | 28.9 | 28.9 | 0% | 0 (n.s.) |
| **A control** | density | 0.029 | 0.029 | 0% | 0 (n.s.) |
| **A control** | giant component frac. | 0.999 | 0.999 | 0% | 0 (n.s.) |
| **B topology** | global clustering | 0.065 | 0.045 | −31% | 22.9 |
| **B topology** | avg local clustering | 0.056 | 0.045 | −20% | 10.3 |
| **B topology** | triangles | 11,200 | 7,690 | −31% | 12.9 |
| **B topology** | degree assortativity | 0.291 | −0.007 | −102% | 36.4 |
| **B topology** | mean path length | 2.48 | 2.41 | −3% | 12.4 |
| **B topology** | modularity (Louvain) | 0.290 | 0.153 | −47% | 29.6 |
| **C homophily** | assort. race | 0.441 | 0.000 | −101% | 57.3 |
| **C homophily** | assort. religion | 0.418 | 0.000 | −100% | 69.1 |
| **C homophily** | assort. age | 0.463 | 0.000 | −100% | 58.8 |
| **C homophily** | assort. education | 0.283 | 0.000 | −100% | 34.9 |
| **C homophily** | assort. sex | 0.105 | 0.000 | −102% | 12.4 |
| **C homophily** | **mean Gower dist. on edges** | 0.229 | 0.325 | **+42%** | **−62.8** |

## Reading

1. **Controls behave (A).** Degree, density, giant component identical (d=0) — the DP randomization is correct; whatever changes is NOT a size/density artifact.
2. **Homophily is annihilated (C).** Every Blau-axis assortativity collapses from strongly positive (0.10–0.46) to ~0 (random mixing). Effect sizes are enormous (d up to 69). The model's own tie distance — **mean Gower distance on edges — rises 42%** (0.23 → 0.33): in plausible nets you are tied to *similar* people; in DP you are tied to *random* people. This is the mechanism that directly feeds the selective-influence filter σ(MSP − d_ij).
3. **Higher-order topology also degrades (B), but less uniformly.** Clustering −31%, triangles −31%, modularity −47%, degree assortativity 0.29 → 0. Path length barely moves (−3%) — DP is still a small world.

## What this can and cannot conclude

- **CAN (for the talk / Paper 4.1):** "breaking the empirical network" is overwhelmingly a **destruction of homophily** (the alignment of attributes with structure) — every axis goes to zero, with the largest effect sizes in the panel, and the model-relevant edge distance jumps 42%. Higher-order topology (clustering/modularity) also drops, so both families are candidate culprits — but the homophily collapse is the most complete and the only one that touches the model's d_ij directly.
- **CANNOT (honest caveat):** because DP destroys B and C *simultaneously*, this panel CANNOT causally separate "homophily" from "clustering" as the driver of the −57% adoption premium. That separation requires the **attribute-shuffle null (#3)** — keep the topology byte-identical, permute node attributes → B unchanged, only C destroyed. Deferred to post-Sunbelt (re-simulation needed). This panel is the *descriptive* "map of suspects"; the shuffle is the *causal* test.

## Outputs
- `data/structural_diagnosis_summary.csv`, `data/structural_diagnosis_raw.csv`
- `plots/structural_diagnosis_B.pdf` (3×3 family B), `plots/structural_diagnosis_C.pdf` (1×2 family C)
