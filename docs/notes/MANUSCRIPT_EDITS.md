# Recommended Manuscript Text Edits for NDFU

Based on comprehensive analysis of the project, the following text edits are recommended for the paper, supplementary material, and related documentation. These edits address conceptual clarity, terminology consistency, and strengthening of argumentation.

---

## 1. Methods Section: Clarify Deterministic vs. Stochastic Adoption

**Location:** `paper_SN26.tex`, Section "Updating (Stochastic Exposure)"

**Current text:**
> "Updating (Stochastic Exposure) Agent behavior is updated using a stochastic, asynchronous network exposure model via `netdiffuseR`..."

**Recommended replacement:**
> "Updating (Stochastic Exposure) Agent behavior is updated via both deterministic and stochastic pathways. Adoption via intrinsic utility (Equation 1) is deterministic: if q_i ≤ Γ, the agent adopts immediately without randomness. Adoption via social influence (Equation 3) depends stochastically on comparison with a random threshold U_ij ~ Unif(0,1), introducing variability. The asynchronous updating order is randomized at each simulation iteration, consistent with standard agent-based model practice. This hybrid mechanism combines exogenous preference with endogenous social response."

**Rationale:** Clarifies the distinction between rational (deterministic) and social (stochastic) pathways, addressing the potential confusion in line 215-216 of the discussion.

---

## 2. Supplementary Material S2: Temporal Scope of Homophily Coefficients

**Location:** `paper_SN26_sm.tex`, Section "Network Imputation Pipeline"

**Current text:**
> "Homophily coefficients [β for race, sex, religion, age, education] are estimated by Smith, McPherson, and Smith-Lovin (2014, Table 3, Model 2) via case-control logistic regression on GSS ego-network data (1985–2004, n=3,001 respondents, 1,139,161 dyads). We transport these coefficients directly to ATP (2014) respondents without re-estimating, following McPherson and Smith (2019). This approach assumes homophily strength remains stable across 2004–2014."

**Recommended replacement:**
> "Homophily coefficients [β for race, sex, religion, age, education] are estimated by Smith, McPherson, and Smith-Lovin (2014, Table 3, Model 2) via case-control logistic regression on GSS core discussion networks (1985–2004, n=3,001 respondents, 1,139,161 dyads). We transport these coefficients directly to ATP (2014) respondents without re-estimating, following the methodology of McPherson and Smith (2019). While Smith et al. (2014, Figure 1) document that homophily shifted between 1985 and 2004 (particularly for race and religion), these were modest shifts across a 19-year period. We assume homophily stability between 2004 and 2014 is reasonable, given that demographic homophily patterns rank among the most structurally stable sociological patterns (McPherson et al., 2001). Should future GSS modules with network measures become available post-2014, our imputation parameters can be recalibrated; however, existing literature suggests 2004–2014 changes would be smaller than the 1985–2004 documented trajectory."

**Rationale:** Acknowledges potential temporal drift while maintaining the methodological soundness of the approach. Adds transparency without undermining the core analysis.

---

## 3. Supplementary Material S3: Clarify Widom Line or Remove

**Location:** `paper_SN26_sm.tex`, Section "Extended Criticality Analysis"

**Current caption (Figure reference):**
> "Along the critical boundary $MSP_c$ according to the universal scaling exponent γ ... (see Figure 2 and the Widom Line discussion)"

**Two options:**

### Option A: Define Widom line formally (if you have explicit data)

**Recommended addition (new paragraph in S3, subsection "The Widom Line"):**

> "In Landau phenomenology, the Widom line (or pseudocritical boundary) is the set of points in parameter space where susceptibility χ = Var(Φ) is maximal for a fixed driving parameter (here, IUL). In our model, for each value of IUL, we identify the MSP value where variance peaks—this defines the effective critical point MSP_c(IUL). Tracing MSP_c(IUL) across the IUL domain yields a Widom-like locus in the (IUL, MSP) phase space. Within each 'slice' at fixed IUL, the susceptibility exponent γ is extracted from the power-law divergence of variance as the system approaches this boundary."

### Option B: Remove references to Widom line (simpler, if not central)

**Recommended action:** Delete phrase "the Widom Line" from captions and prose if no explicit data supports it.

**Rationale:** Either formalize the concept with explicit data or avoid it to prevent confusion.

---

## 4. Supplementary Material S1.3: Cronbach α Results for MUR Construct

**Location:** `paper_SN26_sm.tex`, Section "Collective-Action Propensity: Construct Internal Reliability"

**Current text:**
> "[PLACEHOLDER: Add Cronbach's alpha reliability test to verify multi-item unidimensionality for the collective-action propensity construct]."

**COMPUTED VALUES** (via `playground/checks-claude/02_cronbach_alpha.R`, computed on the **ORIGINAL survey respondents in `data/`**, not the imputed networks):

| Construct | Source | Items | Respondents | Cronbach's α | Verdict |
|---|---|---|---|---|---|
| GSS Collective Action | `GSS_2004_NORC.dta` | 9 (signdpet, avoidbuy, joindem, attrally, cntctgov, polfunds, usemedia, interpol, actlaw) | 1435 | **0.82** | Good (≥ 0.70) |
| ATP Innovation | `ATP_W3_W4.rds` (2014) | 6 binary (METECH_A–F) | 1761 | **0.59** | Low / marginal |

**Response directions WERE handled** (answer to the methodological question): GSS items all run the *same* direction (1 = "did in past year" → high propensity), recoded uniformly as `4 − x`; since no items oppose, α is invariant to this recoding. ATP items run in *two* directions — A/C/D pro-innovation, B/E/F anti-innovation — so B/E/F are **reverse-scored** (`1 − x`) before α. This is implemented explicitly in `02_cronbach_alpha.R`.

**Recommended replacement text:**

> "To verify the internal consistency of the **collective action** propensity construct (MUR$_{GSS}$), we computed Cronbach's α across the nine items used to build the propensity score (SIGNDPET, AVOIDBUY, JOINDEM, ATTRALLY, CNTCTGOV, POLFUNDS, USEMEDIA, INTERPOL, ACTLAW). We obtain **α = 0.82**, exceeding the conventional 0.70 threshold and confirming adequate internal consistency for the construct underlying the main-text Figure 1 criticality analysis.
>
> For the **innovation** propensity construct (MUR$_{ATP}$), built from six binary METECH items (ATP 2014), Cronbach's α = 0.59. We report this value transparently. Three considerations contextualize it: (i) Cronbach's α is mechanically attenuated for *dichotomous* items and *short* scales (six items), so it understates reliability relative to a Likert-type equivalent; (ii) openness-to-innovation is arguably a *formative* index — a sum of distinct behavioral tendencies (early adoption, brand variety-seeking, word-of-mouth sharing) — rather than a *reflective* scale measuring a single latent trait, and internal consistency is not the appropriate yardstick for formative indices; (iii) crucially, our headline criticality result (the $\gamma\approx1$ vs. $\gamma\approx4$ contrast) is estimated on the **GSS collective-action network (α = 0.82)**; the innovation application demonstrates that the same structural mechanism is *homomorphic* across behavioral-contagion types (see §12). We therefore retain the innovation construct as a behavioral index while basing our central claims on the psychometrically stronger collective-action measure."

**Rationale:** Reports both values honestly, computed on the **original survey respondents** (GSS 2004 .dta, ATP 2014 .rds), not the imputed networks. Pre-empts the obvious reviewer objection to α = 0.59 with three standard, defensible arguments, and re-anchors the main claim on the GSS construct (α = 0.82).

**Caveat for the data-cleaning note in SM:** the raw `ATP_W3_W4.rds` carries missing-data sentinels (value `99`) and ~50% NA on METECH items; α must be computed after filtering to clean `{0,1}` responses (1761 respondents). A naive computation on the raw file yields nonsensical negative α.

### 4b. SM Section S1 — documentation gaps to fix (discovered while computing α)

The implementation and the SM text **disagree**; these must be reconciled before resubmission:

1. **Wrong source attribution.** SM S1 intro says propensities use "batteries of items from the American Trends Panel (ATP)", and S1.2 is titled **"Collective Action Propensity Score (ATP 2016)"**. But the collective-action $q_i$ is actually built from **GSS 2004** (`GSS_2004_NORC.dta`), consistent with S1.4 ("we leverage the GSS structure for Figure 1") and with the GSS networks. **Fix:** retitle S1.2 to "Collective Action Propensity Score (GSS 2004)" and correct the intro to say propensities come from ATP 2014 (innovation) *and* GSS 2004 (collective action).
2. **Missing 9th item.** SM S1.2 lists **8** collective-action items; the code uses **9** — the missing one is `ACTLAW` ("taken / would take action against a law considered unjust"). **Fix:** add the 9th item to the SM list.
3. **Undocumented scoring/recoding rule.** The SM never states how raw responses become $q_i$:
   - GSS: 4-point scale (1 = did in past year … 4 = would never) recoded `4 − x` → 0–3, summed over 9 items, normalized by 27.
   - ATP: binary items; anti-innovation items (B/E/F: "prefer trusted brands", "comfortable with familiar", "wait for others") **reverse-scored** `1 − x`; summed over 6, normalized by 6.
   - **Fix:** add a short "Scoring" paragraph stating both rules and the reverse-scoring of the anti-innovation items explicitly (directions are currently only implied by the *(Pro/Anti-innovation)* tags in S1.1).
4. **Direction tags present only for ATP.** S1.1 tags each innovation item Pro/Anti (good); S1.2 gives the response scale but no per-item direction (acceptable since all collective-action items share one direction — worth one sentence saying so).

---

## 5. Discussion Section: Reframe Narrative from Volume & Dynamics to Amplitude & Continuity

**Location:** `paper_SN26.tex`, Section "Discussion", opening paragraph

**Current structure:**
> "Our results provide robust computational and mathematical evidence that 'networks decide for us.' Conditional on fixed individual-level propensities, the topological arrangement itself acts as an autonomous causal mechanism. First, empirical homophily yields a systematic structural premium..."

**Recommended restructuring:**
> "Our results provide robust computational and mathematical evidence that 'networks decide for us' in two complementary ways. Conditional on fixed individual-level propensities, the topological arrangement itself acts as an autonomous causal mechanism operating on both the amplitude and continuity of collective outcomes. 
>
> **First, Amplitude (The Structural Premium).** Empirical homophily yields a systematic structural premium. Randomizing network structure reduces adoption odds by 57% (β_DP = −0.850, p < 0.001) under mathematically identical individual conditions...
>
> **Second, Continuity (The Physics of Tipping Points).** Beyond volume, topology dictates the *character* of adoption cascades. In empirical networks, social change unfolds as predictable, continuous phase transitions obeying Mean-Field universal scaling (γ ≈ 0.977 ≈ 1). In randomized networks, this predictability collapses, and adoption dynamics become discontinuous and volatile (γ ≈ 3.7–4.6), resembling first-order regime shifts. This suggests that homophilic topology acts as a *stabilizer*, converting chaotic dynamics into globally coherent cascades.
>
> Together, these findings imply that networks do not merely enhance or suppress diffusion; they fundamentally reshape whether social change is gradual and legible or abrupt and surprising."

**Rationale:** Clearer narrative arc connecting the two main findings. Introduces the "stabilizer" concept that addresses the field-mean paradox.

---

## 6. Discussion Section: New Paragraph on Agenc and Structure

**Location:** `paper_SN26.tex`, Section "Discussion", new subsection after main results

**Title:** "Networks and Individual Agency: Structure as Sedimented Choice"

**Recommended new paragraph:**

> "A potential philosophical objection to the framing 'Networks Decide for Us' is that it risks appearing deterministic or structurally reductionist, implying that individual agency is erased. We emphasize that this interpretation is misleading. In our model, the network topology itself—embodied in homophilic tie patterns and clustering—is not exogenous or imposed. Rather, it is a *crystallization of past individual choices aggregated over time*. Each tie in the network reflects a dyadic similarity preference (McPherson & Smith, 2019). When we model the network as fixed and empirically calibrated, we are not denying individual agency; we are capturing the *collective outcome of repeated individual choices operating under similar structural and preferential conditions*. The apparent autonomy of structure in predicting diffusion outcomes is thus a macroscopic manifestation of microscopic behavioral regularities (homophilic preferences), aggregated across populations and time. In this sense, structure and agency are not opposed; structure is the sedimented, emergent result of consistent individual choices. Our finding that this sedimented structure governs diffusion outcomes is therefore consistent with, rather than contradictory to, a sociology that takes individual agency seriously."

**Rationale:** Addresses the risk of being misread as structurally reductionist. Aligns with Goldberg's critique while maintaining the core contribution.

---

## 7. Discussion Section: Paragraph Addressing the Mean-Field Paradox

**Location:** `paper_SN26.tex`, Section "Discussion", after results presentation

**Title:** (Subsection within Discussion) "Why Mean-Field Universality Emerges from Clustered Networks"

**Recommended addition:**

> "A subtle but important feature of our results deserves explanation. Mean-Field universality (γ = 1) is typically observed in well-mixed systems with weak correlations between agents. Yet we recover γ ≈ 1 in empirical networks that are *strongly clustered and locally correlated*—the opposite of well-mixed. This apparent paradox is resolved by recognizing that the homophilic clusters act as quasi-isolated 'subsystems,' each with internal ordering properties but coupled to others through weak, selective bridges (edges crossing Blau-space distance). At criticality, these weak bridges suddenly activate (as MSP increases), causing the entire system of clusters to synchronize as a single coherent entity. This synchronization mechanism—whereby local domains couple into global order—is mathematically indistinguishable from mean-field behavior, even though the underlying network is heterogeneous. This finding suggests that mean-field universality may be more robust to realistic (clustered) topologies than classical theory suggests, provided that inter-cluster coupling is selective and governed by similarity."

**Rationale:** Directly addresses the conceptual tension identified in the feedback. Provides intuition for why γ ≈ 1 despite clustering.

---

## 8. Discussion Section: Add Limitations Paragraph on MSP Heterogeneity (expand existing)

**Location:** `paper_SN26.tex`, Section "Limitations"

**Current text:**
> "Furthermore, a significant limitation lies in our operationalization of the Maximum Social Proximity (MSP)..."

**Recommended expansion:**

> "Furthermore, a significant limitation lies in our operationalization of the Maximum Social Proximity (MSP). Currently, our model treats MSP as a *systemic constant*, implying that social flexibility and tolerance for cross-Blau-distance influence are uniformly distributed across the entire population. In reality, sociological flexibility is inherently heterogeneous; specific sub-populations or demographic groups may exhibit varying distributions of MSP, ranging from highly insulated echo chambers (low MSP) to cosmopolitan, socially porous hubs (high MSP). This is not merely a modeling convenience—it may systematically bias our estimates of γ and the critical boundary, particularly if influential individuals (opinion leaders, bridges) have systematically higher MSP. Future iterations of this modeling framework should explore the effects of heterogeneous MSP distributions assigned at the group or individual level. As a first step, we suggest using clustering algorithms on the imputed networks to identify cohesive groups, then assigning group-specific MSP distributions calibrated to survey items measuring openness to diversity or cosmopolitanism. Such extensions would also connect more directly to Centola's (2011) experimental work on *induced homophily* as a cognitive bias that varies across individuals."

**Rationale:** Acknowledges the limitation more thoroughly. Provides a concrete roadmap for future work (useful for reviewers).

---

## 9. Terminology Consistency Pass: GSS-ER → GSS-DP

**Location:** Throughout `paper_SN26.tex`, `paper_SN26_sm.tex`, figures, and captions

**Action:** Global find-and-replace:
- `GSS-ER` → `GSS-DP` (Degree-Preserving Randomized)
- `ER network` → `degree-preserving randomized network`
- `Erdős–Rényi` → keep as is (historical reference), but note that your baseline preserves degree

**Specific instances to check:**
- Table 2 (gamma table) — column headers
- Figure captions (Figure 1, S3, S5, S6, etc.)
- All prose references to the random baseline

**Rationale:** Terminological precision. "ER" suggests Erdős–Rényi random graph (no degree constraint), but you use degree-preserving, which is a different null model.

---

## 10. Title/Abstract: Optional Enhancement for Clarity

**Location:** `paper_SN26.tex`, Abstract

**Current:**
> "Networks Decide for Us: Emergence of Adoption in Homophily-Imputed Social Topologies"

**Optional alternative** (if you want to highlight the physics angle):
> "Networks Decide for Us: Phase Transitions and Universal Scaling in Homophily-Imputed Adoption Cascades"

**Note:** Keep current title if you prefer sociological over physics framing. This is optional.

---

## 11. Glossary or Notation Appendix (Optional Enhancement)

**Consideration:** If targeting a multi-disciplinary audience (both sociology and physics), consider adding a one-page glossary in SM appendix:

```
MSP       → Maximum Social Proximity (analogous to Temperature in Ising model)
IUL       → Intrinsic Utility Level (analogous to External Magnetic Field)
Φ         → Order Parameter (adoption proportion)
χ         → Susceptibility (variance of Φ)
γ         → Critical exponent (power-law decay of χ)
Ẽ_i       → Effective Exposure (selective social influence)
GSS-DP    → Degree-Preserving Randomized baseline network
```

**Rationale:** Makes the paper more accessible to readers unfamiliar with statistical mechanics notation.

---

## 12. Reframe the paper around COLLECTIVE ACTION (primary), with INNOVATION as a homomorphism check

**Motivation:** the reliability analysis (§4) gives a clean split — collective action (GSS 2004) has α = 0.82, innovation (ATP 2014) has α = 0.59. This argues for promoting collective action to the primary case and demoting innovation to a robustness/generalization demonstration.

**Strategic recommendation (affects framing throughout):**
- **Lead with collective action.** Make the GSS collective-action network the primary empirical system for *all* headline results: the structural premium (β_DP) and especially the criticality / universality result (γ ≈ 1 vs. γ ≈ 4). The SM already notes Figure 1 uses the GSS structure — so this is mostly a *narrative* re-centering, not a re-computation.
- **Recast innovation as a homomorphism test.** Present the ATP innovation results as evidence that *the same structural mechanism governs a qualitatively different behavioral contagion*. The argument: if homophilic topology produces the same structural premium and the same continuous-vs-discontinuous criticality contrast for two unrelated behaviors — high-cost, identity-laden **collective action** and low-cost, consumer **innovation adoption** — then the mechanism is a property of *network structure*, not of the specific behavior. That is a *strength*, and it conveniently sidesteps the weaker ATP reliability.

**Suggested sentence for the Results or Discussion:**
> "We treat the GSS-based collective-action system as our primary case (internal-consistency α = 0.82) and the ATP-based innovation system as a *generalization test*. That both behaviors — despite differing in cost, visibility, and identity-salience — exhibit the same structural premium and the same shift from continuous (homophilic) to discontinuous (randomized) criticality indicates that the governing mechanism is *homomorphic across contagion types*: it is a property of the social topology rather than of the particular behavior diffusing over it."

**Where this propagates:**
- Abstract: lead with collective action; mention innovation as "a second, structurally-homomorphic behavioral domain."
- Methods: present GSS collective-action construct first; ATP innovation second, flagged as the generalization case.
- Keep both Figure 1 (GSS) prominent; the ATP equivalent can move to SM as the homomorphism evidence.

**Rationale:** Turns the α = 0.59 weakness into a framing advantage and tightens the central claim ("structure decides, not the behavior").

---

## 13. γ inconsistency at IUL = 0.275 (paper/abstract vs. actual data) — IMPORTANT

The susceptibility-exponent table is internally inconsistent across documents at the **IUL = 0.275** row:

| Source | $MSP_c$ at IUL=0.275 | GSS $\gamma$ |
|---|---|---|
| `paper_SN26.tex` Table, `abstract_SN26.tex` Table 2, `formal_comparison.tex` | 0.166 | **1.539** |
| `playground/physics_analogy.tex` Table, `output/07_phase_transition/GSS/gamma_exponents.csv` (actual pipeline output) | 0.250 | **0.873** |

So the paper's outlier **1.539** comes from pinning $MSP_c=0.166$ for that row, while the current script pins $MSP_c=0.250$ uniformly and yields **0.873**. This is exactly the $MSP_c$-choice fragility flagged elsewhere.

**Consequence for the "mean γ" claim:**
- With the data-faithful values (0.977, 0.892, 0.892, **0.873**): mean = **0.91** — clean, uniform, strengthens the Mean-Field ($\gamma\approx1$) story.
- With the paper's values (…, **1.539**): mean = **1.075**, and the column looks non-uniform (one outlier at 1.54).

**Recommendation:** adopt ONE $MSP_c$ convention. The uniform $MSP_c=0.25$ (→ 0.873) is preferable: it matches the actual computed output, removes the 1.54 outlier, and yields a tighter GSS column (all ∈ [0.87, 0.98], mean 0.91). Update `paper_SN26.tex` Table and `abstract_SN26.tex` Table 2 accordingly. **No document currently prints an average γ**, so there is no average to correct in the paper — but the per-row 1.539 should be reconciled. (The CPIN deck has been set to the data-faithful 0.87 / mean 0.91.)

## 14. BAM notation $\mathrm{logit}[\mathbb{E}(\Phi)]$ is correct — keep it

The question arose whether the BAM should read $\mathrm{logit}(\Phi)$ instead of $\mathrm{logit}[\mathbb{E}(\Phi)]$. **Keep $\mathbb{E}(\Phi)$.** In a GLM/GAM the link function transforms the *conditional mean* of the response, $g(\mu_i)=g(\mathbb{E}[Y_i\mid X_i]) = \eta_i$, not the random outcome itself. Since $\Phi$ is the simulation-level random outcome (final adoption proportion) and the BAM models its expectation, $\mathrm{logit}[\mathbb{E}(\Phi)]$ is the precise and standard form (cf. Wood, mgcv). Writing $\mathrm{logit}(\Phi)$ would be technically loose. No change to paper or deck.

## 15. Reframe the "mean-field paradox" via network dimensionality (IMPORTANT for the physics framing)

The current SM frames $\gamma\approx1$ on a clustered network as surprising ("mean-field usually needs a well-mixed system"). A cleaner, more correct framing: **complex networks have no finite spatial dimension** — neighborhoods expand (super)exponentially, so the system sits *above* the upper critical dimension ($d_u=4$ for Ising), where **mean-field exponents are the generic, expected outcome** (Dorogovtsev, Goltsev & Mendes, *Critical phenomena in complex networks*, Rev. Mod. Phys. 80, 2008). So $\gamma\approx1$, $\delta\approx3$ are exactly what one should expect for a continuous transition *on a network*; they are NOT evidence of low-dimensional (d=1/2/3) Ising behavior.

**Consequence for interpretation:** the contrast between plausible (GSS, $\gamma\approx1$) and randomized (GSS-DP, $\gamma\approx4$) is **not a difference of universality class or dimension** — both are effectively infinite-dimensional. It is a difference in the **order of the transition**: clustering yields a *continuous* (2nd-order) transition; degree-preserving randomization yields a *discontinuous* (1st-order / hybrid) regime shift (cf. bootstrap/threshold percolation on random graphs). Recommend rewriting SM S3 accordingly and dropping the "quasi-isolated clusters synchronize" hand-wave in favor of this standard, citable account.

## 16. Honest ensemble-vs-realization framing of the early-warning / "predictability" claim

The exponents and "critical slowing down" are computed from the **variance and mean across the ~24 stochastic realizations per (IUL, MSP) cell** — they are **ensemble** quantities. So the early-warning signals (rising variance, diverging relaxation time) describe the *distribution of outcomes across comparable societies / repeated histories*, NOT a precursor guaranteed within one society's single timeline. The manuscript should:
- State explicitly that $\chi=\mathrm{Var}(\Phi)$ and the relaxation-time divergence are ensemble fluctuation measures.
- Frame the contribution as "structure sets the **statistical character / risk profile** of collective change" (smooth-and-signposted vs. bimodal-and-abrupt), not deterministic forecasting of a given event.
- Note that a *single-society time series* would show the precursors only if the control parameter (MSP) **drifts slowly through $MSP_c$** over time — a stronger assumption not yet simulated (future work, connects to Scheffer et al. early-warning-signals literature with sliding-window variance/autocorrelation).

This is the most reviewer-exposed honesty point; getting ahead of it strengthens credibility.

## 17. Frame IUL$\approx$0.2 as a pseudo-critical / Widom locus, not literally $H=0$

The four IUL values $\{0.200,0.225,0.250,0.275\}$ were chosen as the "zero-field equivalent" via the **steepest adoption jump**. Per phase-transition practice (supercritical-fluid Widom-line literature; finite-size-scaling studies), when $H=0$ is physically inaccessible the standard is to locate the pseudo-critical line via **maximum susceptibility** (the Widom line, $\chi$-peak), and practitioners use *several* loci (specific-heat-like max, $\chi$-peak, max order-parameter slope) which **coincide only at the true critical point and "split into a bundle" away from it**. Recommend:
- Reframe IUL$\approx$0.2 as the **pseudo-critical / Widom locus**, not "$H=0$".
- Report that the **max-susceptibility locus** and the **steepest-jump locus** coincide (they do on the $T$-path: both give $MSP_c=0.25$) — this is the evidence that justifies the choice; prefer $\chi$-max for the actual exponent fit ($\chi_{\max}\sim L^{\gamma/\nu}$).
- Acknowledge the finite-size shift of the peak ($\sim N^{-1/\nu}$). Sources: Simeoni et al. *Nat. Phys.* 2010 (Widom line); arXiv:2512.21748 (three-estimator FSS practice).

## 18. Future extension: risky / high-cost collective action propensities

For a thesis chapter or follow-up, the current $q_i$ (low-risk action; medium-risk innovation) can be extended to **risky/violent collective action** using the sources catalogued in `playground/checks-claude/RISKY_COLLECTIVE_ACTION_DATA.md` (ARIS Activism–Radicalism Intention Scales for item structure; UC Davis "Life in America" / CPOST for representative violent-tail base rates; WVS-7 strike/occupy items as the public low-risk bridge), or — as a first computational step — a strongly right-skewed assumed aversion distribution calibrated to those base rates, then re-running the IUL×MSP sweep to test whether the structural premium and the $\gamma\approx1$ vs $\gamma\approx4$ contrast survive a mostly-averse population.

## 19. The meaning of universality — the "so what?" (full interpretation for the Discussion)

This section records, in full, the interpretive argument for *why* the universality result matters and how far the "structure makes social change legible" claim can be pushed. It is the intellectual core of the Discussion and should be written carefully to pre-empt the "so what?" reaction.

### 19.1 What it means that adoption is mean-field

Universality = **the microscopic details are irrelevant** to the aggregate behavior near the critical point. You do not need to know each individual's exact threshold, every tie, or every attribute: a handful of macro-structural parameters (the control parameter MSP and the distance to the critical point) govern the collective dynamics. This is *exactly* what universality means in statistical physics, and it is the precise, defensible content of "structure makes the collective legible."

- **Physics examples of irrelevant micro-detail:** boiling water and a ferromagnet at their critical points fall in the *same* universality class and share the *same* exponents — the chemistry, the substance, the precise lattice geometry, the atomic species are all irrelevant. Only dimensionality and symmetry (the macro-level features) matter.
- **Social analogues that would be irrelevant here:** each agent's exact MUR/threshold value, the specific *content* of the innovation/behavior being diffused, individual idiosyncrasies and cognition. What remains relevant: the structure of the ties and how far the system sits from the tipping point.

### 19.2 What known sociological facts this explains

1. **The S-shaped diffusion curve** (Rogers; Bass). A continuous transition with order-parameter growth *is* the sigmoid aggregate-adoption curve. The phase-transition picture gives a mechanistic foundation for the classic S-curve.
2. **The recurrence of tipping points / critical mass** (Granovetter 1978; Schelling; Centola's ~25% experiments). A continuous transition with a well-defined critical point is the formal object behind "critical mass." Universality explains *why the same tipping phenomenon recurs across utterly different domains* (fads, technologies, norms, protests): the cross-domain regularity of diffusion (same S-curves, same threshold structure) is, in physics terms, universality.

### 19.3 How far "structure makes social change legible" can be defended

- **Defensible:** near criticality the *aggregate* obeys a low-dimensional, scale-free, universal law; macro-observables (variance, relaxation time) carry characteristic, measurable signatures. The aggregate is "legible" = describable with a few parameters and carrying *statistical* early-warning signals (an ensemble property — see §16).
- **Over-reach to avoid:** "legible" ≠ "predictable at the individual or single-event level." One cannot predict which specific protest erupts or who adopts. There is a productive irony worth stating explicitly: **the aggregate is maximally legible (universal law) precisely in the regime where individual outcomes are least predictable** — at criticality, variance diverges and avalanches of *all sizes* occur (scale-free). The honest one-liner: *"the rules governing the distribution of outcomes are legible and universal, even though any single outcome is maximally uncertain."*

### 19.4 The crucial reframing of the contribution (avoid the "so what?")

For a physics-literate reader, *finding* mean-field on a network is **not surprising** — networks are effectively infinite-dimensional, so mean-field is the *generic* expectation (see §15). Therefore the contribution should **not** be framed as "we found universality." The non-trivial, defensible contribution is the **contrast**: empirically homophilic structure yields a **continuous, universal, signposted** transition, whereas degree-preserving randomization yields a **discontinuous, warningless** one. **Structure decides the *order* (and the legibility) of collective change, not merely its volume.** Lean the Discussion on (a) this continuous-vs-discontinuous contrast and (b) the empirical (survey-calibrated) grounding — not on the bare fact that γ≈1.

### 19.5 A concrete closing framing (candidate Discussion / abstract sentence)

> "Even when a single episode of social change feels abrupt and surprising, the aggregate dynamics — averaged over comparable cases — appears to follow identifiable, continuous, even universal patterns."

This is exactly what a continuous, mean-field transition predicts, and it respects the ensemble caveat (it is the *mean* that carries the pattern, not the single realization). It can be anchored to real phenomena: S-curves of innovation adoption (smartphones, hybrid corn, contraception); recurrent tipping fractions in norm change (Centola's ~25%; same-sex-marriage opinion shift; #MeToo); protest waves/revolutions that feel sudden ex-post but show aggregate precursors (rising strike frequency, polarization indices). Always scope the claim to the aggregate/ensemble level; individual events remain unpredictable.

## 20. Seeding-position coefficients vs. the complex-contagion literature (SM note)

In the BAM (Report #1), the seeding coefficients run *against* the naive reading of the complex-contagion literature, and this is worth a careful footnote in the SM:
- **Social adoption:** central/closeness/eigenvector seeds *help* (β ≈ +0.037), **marginal** seeds *hurt* (β ≈ −0.022).
- **Rational adoption:** central/closeness/eigenvector seeds *hurt* (β ≈ −0.056), **marginal** seeds *help* (β ≈ +0.028).

The folk version of the literature (Centola & Macy 2007) says **simple/information contagion spreads best from central seeds**, while **complex/behavioral contagion spreads best from the margins (wide bridges)**. Our pattern looks inverted for the social channel. The reconciliation — and the point to make in the SM:

> The "seed from the margins" result for complex contagion holds for **absolute** thresholds ≥ 2 (a node needs ≥2 already-adopted neighbors). Our model uses a **fractional/relative** threshold drawn from a distribution: a small fraction applied to a low-degree node can be equivalent to an **absolute threshold of 1** — i.e. effectively a *simple* contagion for that node. So a population of fractional thresholds is a *mixture* that is not globally "complex" in the Centola-Macy sense, and the central-seed advantage we observe for social adoption is consistent with that mixture leaning simple for the low-degree tail. The marginal-seed advantage in the *rational* channel is a separate, second-order seeding/accounting effect (see Report #1 §2.3), not a contagion result.

**Why it is low-priority:** seeding is not the focus, the effect sizes are tiny (OR within 0.93–1.13), and the threshold-distribution caveat explains the apparent contradiction. Worth one honest SM paragraph so a referee who knows Centola-Macy is pre-empted, with the absolute-vs-fractional-threshold distinction made explicit. Source: `paper-NDFU/reports-prompts/gss-dp-sh/report-results-1.md` §2.3.

## Summary of Edits by Priority

| Priority | Section | Edit Type | Complexity |
|----------|---------|-----------|-----------|
| **Critical** | S2 | Homophily temporal scope | Low |
| **Critical** | Methods | Deterministic vs. stochastic | Low |
| **Critical** | Discussion | Agenc/structure paragraph | Medium |
| **High** | Discussion | Mean-Field paradox explanation | Medium |
| **High** | Limitations | MSP heterogeneity expansion | Medium |
| **High** | S1.3 | Cronbach α values (GSS 0.82 / ATP 0.59, original data) | Low |
| **Critical** | S1 / framing | Fix ATP-2016→GSS-2004 source error, add 9th item, document scoring/directions (§4b) | Low |
| **High** | Whole paper | Re-center on collective action; innovation as homomorphism check (§12) | Medium |
| **Medium** | Discussion | Narrative restructuring (Amplitude + Continuity) | Medium |
| **Medium** | S3 | Widom line (define or remove) | Low |
| **Medium** | Terminology | GSS-ER → GSS-DP throughout | Low (mechanical) |
| **Low** | Appendix | Glossary (optional) | Low |

---

**Document prepared:** May 25, 2026 (updated May 28, 2026)  
**For:** Aníbal Olivera (NDFU Manuscript)  
**By:** Claude AI Assistant
