# Networks Decide for Us (NDFU)

Do empirically calibrated homophilic network topologies act as an autonomous causal force
on diffusion, holding individual propensities fixed?

The pipeline estimates homophily strengths from ego-network data, imputes plausible whole
networks onto survey respondents who never reported relational data, assigns each node an
empirically measured adoption propensity, and simulates a hybrid social contagion over
those networks against degree-preserving randomized counterfactuals.

**Authors:** Aníbal Olivera Morales, Jorge Fábrega Lacoa (CICS, Universidad del
Desarrollo), Thomas W. Valente (USC).

---

## The adoption rule (unified, bottom-up)

An agent adopts when the utility the situation offers meets its own measured requirement:

```
adopt_i(t)   <=>   Γ + λ · Ẽ_i(t)  ≥  q_i
```

| Symbol | Meaning | Source |
|---|---|---|
| `Γ` (IUL) | intrinsic utility of the behavior | swept, 41 levels |
| `q_i` (MUR) | the agent's aversion / minimum utility requirement | **empirical** (GSS/ATP items) |
| `Ẽ_i` | selective effective exposure: only adopters that pass the social-proximity gate count | dynamic |
| `h` (MSP) | maximum social proximity — the population's openness | swept, 25 levels |
| `λ` | social-proof coupling: what a fully adopted, fully similar neighborhood is worth | calibrated, ≈ 0.6–0.7 |

The network threshold is **derived**, not postulated: `τ*_i = max(0, (q_i − Γ)/λ)`. Agents
with `q_i ≤ Γ` adopt at the first influential contact (the old "rational" channel);
everyone else needs selective exposure proportional to their utility deficit (the old
"social" channel). All heterogeneity comes from the measured `q_i`.

The design note **[docs/MAKING-IT-BOTTOM-UP.pdf](docs/MAKING-IT-BOTTOM-UP.pdf)** derives
the rule, proves it nests Tur et al. (2024) exactly, compares the two operationalizations
(linear-fractional vs. the multiplicative γ-form), and records the calibration and the
structural-premium results.

---

## Pipeline

Run from the repository root (`Rscript scripts/<file>`), in order:

| Script | What it does | Writes to |
|---|---|---|
| `00_install_dependencies.R` | package setup | — |
| `01_ATP_GSS_imputation.R` | survey harmonization, Blau-space attributes | `data/01_*`, `plots/01_*` |
| `02_GSS_network_ergm.R` / `02_ATP_network_ergm.R` | ERGM-imputed plausible networks (homophily strengths from Smith et al. 2014) | `data/02_*_network*` |
| `02_*_network_ergm_tests.R`, `02b_network_plots_EN.R` | topology validation vs. ER / WS / BA benchmarks | `output/02_*`, `plots/02_*` |
| `03_GSS_MUR_calculation.R` / `03_ATP_MUR_calculation.R` | propensity scores `q_i` from survey items | `data/`, `plots/03_*` |
| `04_make_DP_networks.R` | **DP null model**: degree-preserving twins (`rewire(keeping_degseq)`) | `data/02_GSS_DP_network/`, `output/04_null_models/` |
| `04b_make_SH_networks.R` | **SH null model**: topology byte-identical, node attributes permuted | `data/02_GSS_SH_network/`, `output/04_null_models/` |
| `05_unified_diffusion_sweep.R` | **the engine**: unified-rule diffusion, GSS vs DP, + premium BAM | `output/05_unified_diffusion/`, `plots/05_*` |
| `05_unified_diffusion_sweep_main.R` | runner: the engine across all five seeding strategies | (as above, one set per strategy) |
| `06_premium_sensitivity.R` | how the premium moves with λ (odds vs. probability scale) | `output/06_premium_sensitivity/` |
| `07_structural_diagnosis.R` | what randomization destroys (clustering, assortativity, tie distance) | `output/07_structural_diagnosis/` |
| `08_cronbach_alpha.R` | internal consistency of the propensity constructs | `output/08_cronbach_alpha/` |

### Running the simulations

The engine runs **one seeding strategy per invocation**, each writing to its own
checkpoint directory and its own result files:

```bash
Rscript scripts/05_unified_diffusion_sweep.R                    # random (default)
Rscript scripts/05_unified_diffusion_sweep.R --seeding=central
Rscript scripts/05_unified_diffusion_sweep.R --test             # 30-second smoke run
Rscript scripts/05_unified_diffusion_sweep_main.R               # all five, in sequence
```

Strategies: `random`, `central` (highest degree), `closeness`, `eigen`, `marginal`
(bottom 10% by degree). One strategy is 1,180,800 simulations, ≈ 50 min on 8 cores.

Everything **checkpoints per (topology, λ)**, so a run can be interrupted and relaunched:
completed combinations are skipped. `--test` writes to its own sandbox and can never
overwrite real results.

### `tools/` — off-pipeline

`tools/lambda_calibration.R` sweeps λ ∈ {0.3, …, 1.5} against the legacy engine's surface
and is how λ ≈ 0.7 was identified. It is **not part of the study**: it is a one-off
methodological calibration, and it only runs where `legacy/` is present.

### Requirements

`netdiffuseR` (stochastic-transmission branch) is required by the legacy engine:

```r
devtools::install_github("USCCANA/netdiffuseR", ref = "stochastic-transmission")
```

The unified engine in `05_unified_diffusion_sweep.R` implements the update loop directly
(sparse edge arrays + a per-step stochastic gate) and does not depend on it.

---

## Null models

| Label | Construction | Degree | Clustering | Homophily |
|---|---|---|---|---|
| **GSS** | ERGM-imputed, empirical homophily | — | 0.065 | present (assort. 0.10–0.46) |
| **GSS-DP** | degree-preserving double-edge swaps | preserved | 0.044 (−31%) | destroyed (→ 0) |
| **GSS-SH** | same graph, node attributes permuted | preserved | identical | destroyed (→ 0) |

Mean social distance on ties rises 0.229 → 0.325 under DP: in the plausible networks your
neighbors are similar to you; after randomization they are strangers. That is the quantity
the selective-influence gate reads.

> **Naming note.** Earlier versions of this project called the degree-preserving null "ER".
> They were always degree-preserving; only the label was wrong. Active code now uses **DP**
> throughout. Genuine Erdős–Rényi references survive only in `02_*_network_ergm_tests.R`,
> where ER is used as a *theoretical benchmark* for topology validation.

---

## Repository layout

```
scripts/     active pipeline (see table above)
tools/       off-pipeline utilities (lambda calibration)
data/        survey inputs and imputed networks (ERGM output — expensive, tracked)
output/      analysis results per pipeline stage
plots/       figures per pipeline stage
docs/        design notes; MAKING-IT-BOTTOM-UP.{tex,pdf} is the reference document
docs/notes/  working notes (journal fit, manuscript edits, null-model diagnoses)
paper-NDFU/  manuscript, extended abstract, presentations, reports
legacy/      previous engine and its outputs — NOT tracked (see below)
playground/  scratch and exploratory work — NOT tracked
```

### The legacy engine

The previous engine (two separate channels joined by an OR, with an exogenous network
threshold `τ ~ N(τ_mean, τ_sd)` swept over 16 configurations) produced the results in the
current manuscript draft. It now lives in `legacy/` and is **not tracked by git**:

- `legacy/scripts/` — the old `04_*` simulation, `05_*` plotting, `06_regression_analysis.R`,
  `07_phase_transition_deeper.R`, and the EUSN figure script that reads its outputs;
- `legacy/output/`, `legacy/plots/` — its simulation results and figures;
- `legacy/data/02_GSS_HP_network/` — the abandoned homophily-preserving null model
  (see `docs/notes/HOMOPHILY_TOPOLOGY_INSEPARABILITY.md`).

Everything there remains recoverable from git history (it was tracked until the
`making-bottom-up` refactor). `tools/lambda_calibration.R` reads one legacy file to
calibrate λ against the old surface, so it only runs where `legacy/` is present.
