# Exploratory Notes: The Field Direction (Camino-H)

**Status:** Exploratory. Not for the manuscript as-is. Records what the H-path shows and what it would take to make it rigorous.

**Script:** `03_field_direction_sweep_H.R`
**Data:** GSS and GSS-DP, `mean_0.40`, `sd=0.12`, random seeding (same as main-text Figure 1).

---

## Background: two paths to criticality

Borrowing the ferromagnetic mapping from `playground/physics_analogy.tex`:

| Physics | NDFU | Role |
|---|---|---|
| Temperature $T$ | MSP ($h$) | internal / state parameter |
| External field $H$ | IUL ($\Gamma$) | external bias |
| Magnetization $M$ | Adoption $\Phi$ | order parameter |

- **Camino-T** (main analysis, script 07): fix IUL, sweep MSP across $MSP_c$. This is the $H\approx0$, vary-$T$ protocol that yields a **second-order** transition and the susceptibility exponent $\gamma$. **Result: $\gamma_{GSS}\approx 1$, $\gamma_{DP}\approx 4$.**
- **Camino-H** (this exploration): fix MSP, sweep IUL. The classic first-order field-reversal ($H$ changing sign at $T<T_c$) is **not accessible** because IUL $<0$ is meaningless socially. Instead we sit on the critical isotherm ($MSP=MSP_c$) and sweep IUL $>0$, targeting the **critical-isotherm exponent $\delta$**: $\;\Phi-\Phi_c \sim (IUL-IUL_c)^{1/\delta}$. Mean-field predicts $\delta=3$.

---

## Key finding 1 — On the critical isotherm, $\delta\approx 3$ for BOTH topologies

| Topology | $MSP_c$ (from T-path) | $IUL_c$ (inflection) | $\delta$ | $R^2$ |
|---|---|---|---|---|
| GSS | 0.25 | 0.275 | **2.85** | 0.97 |
| GSS-DP | 0.417 | 0.275 | **2.91** | 0.95 |

Both land near the mean-field value $\delta=3$. **Unlike $\gamma$, the exponent $\delta$ does NOT discriminate between empirical and randomized topologies.**

## Key finding 2 — Why $\delta$ is topology-blind (the interpretation worth keeping)

This is the conceptually valuable outcome of the exploration:

- The field **IUL drives the rational channel** ($q_i \le \Gamma$). This is a **node-local** condition — an agent flips on its own intrinsic utility, *independently of the network*. So the field-direction response is governed by the **distribution of MUR ($q_i$)**, not by topology. Both GSS and DP share the same MUR distribution $\Rightarrow$ same $\delta$.
- The temperature **MSP drives the social channel** (selective exposure $\tilde E_i$), which **is** topology-dependent. So only the $T$-path exponent $\gamma$ "sees" the difference between clustered and randomized structure.

**Takeaway:** the two control parameters probe two different mechanisms. IUL/$\delta$ = topology-blind rational channel; MSP/$\gamma$ = topology-sensitive social channel. The universality *discriminator* is therefore necessarily the $T$-path, not the $H$-path. This is a clean answer to "is it worth analyzing both paths?": **yes, but only the T-path discriminates — the H-path explains *why* by isolating the topology-blind channel.**

## Key finding 3 (inconclusive) — Widom scaling relation does not close yet

The mean-field scaling relation $\gamma = \beta(\delta-1)$ would, if satisfied for GSS and violated for DP, give a powerful *internal-consistency* argument that GSS is a genuine critical point and DP is not. Current status:

- Crude $\beta$ fit (order parameter on the ordered side of the T-path, 13-point MSP grid): $\beta_{GSS}\approx 0.09$, $\beta_{DP}\approx 0.12$, both with **poor $R^2$ (0.6–0.7)** and far from mean-field $\beta=0.5$.
- Widom check fails numerically for both ($0.09\times1.85=0.17 \ne 0.98$).

**This is almost certainly a measurement artifact, not a physical result:** $\beta$ is the hardest exponent to estimate and the MSP axis has only **13 levels** (Δ = 1/12), far too coarse to resolve the order-parameter power law near $MSP_c$. The isotherm-family plot also shows visible **discretization steps** (artifactual stepping from binned IUL/MUR), which biases low-lying fits.

**To make this rigorous would require:**
1. A denser MSP sweep near $MSP_c$ (e.g. Δ = 1/48 in a window $\pm0.1$ around $MSP_c$).
2. A denser IUL sweep near $IUL_c$ for a clean $\delta$.
3. Bootstrap CIs on $\beta, \gamma, \delta$.
4. Then test $\gamma = \beta(\delta-1)$ and $\alpha + 2\beta + \gamma = 2$ (Rushbrooke) as joint consistency checks.

---

## Recommendation for the manuscript / presentation

- **Do NOT** put $\delta$ or the Widom test in the paper yet — too preliminary.
- **DO** keep the conceptual point (Finding 2) as framing: it justifies *why* the MSP sweep is the right universality probe and reassures physics-literate reviewers that the choice of control parameter was principled, not arbitrary.
- For the **CPIN presentation**, the isotherm-family figure (`camino_H_isotherm_family.pdf`) is a nice visual: it shows $\Phi(IUL)$ shifting smoothly as MSP increases — the "field response at different temperatures" — and can motivate the T-path as the criticality axis.
- A denser-grid re-run (item 1–2 above) is the natural next computational step if we want the Widom relation as a headline robustness result.

---

**Generated:** exploratory run, May 2026. See `data/camino_H_exponents.csv` and `plots/camino_H_*.pdf`.
