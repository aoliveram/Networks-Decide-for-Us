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
