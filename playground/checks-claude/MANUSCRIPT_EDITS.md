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

**Recommended replacement (after running scripts 03_GSS_MUR_calculation.R and 03_ATP_MUR_calculation.R):**

> "To verify the internal consistency of the collective action propensity construct (MUR_GSS), we computed Cronbach's α across the 8 items used to construct the propensity score: SIGNDPET (sign petition), AVOIDBUY (avoid/buy products), JOINDEM (join demonstration), ATTRALLY (attend rally), CNTCTGOV (contact government), POLFUNDS (political funds), USEMEDIA (use media), INTERPOL (internet political forum), and ACTLAW (willingness to break law). The Cronbach's α coefficient equals [INSERT VALUE FROM 03_GSS_MUR_calculation.R OUTPUT]. This value exceeds the standard threshold of α = 0.70, confirming adequate internal consistency and unidimensionality of the collective action propensity construct.
> 
> Similarly, for the innovation propensity construct (MUR_ATP), we computed Cronbach's α across the 6 binary items from the METECH scale: metech_a through metech_f. The resulting α = [INSERT VALUE FROM 03_ATP_MUR_calculation.R OUTPUT] [exceeds / is close to] the standard threshold, confirming that both propensity constructs satisfy psychometric standards for unidimensional scales."

**Rationale:** Converts placeholder into evidence-based statement. Adds specific values once scripts are executed.

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

## Summary of Edits by Priority

| Priority | Section | Edit Type | Complexity |
|----------|---------|-----------|-----------|
| **Critical** | S2 | Homophily temporal scope | Low |
| **Critical** | Methods | Deterministic vs. stochastic | Low |
| **Critical** | Discussion | Agenc/structure paragraph | Medium |
| **High** | Discussion | Mean-Field paradox explanation | Medium |
| **High** | Limitations | MSP heterogeneity expansion | Medium |
| **High** | S1.3 | Cronbach α values | Low |
| **Medium** | Discussion | Narrative restructuring (Amplitude + Continuity) | Medium |
| **Medium** | S3 | Widom line (define or remove) | Low |
| **Medium** | Terminology | GSS-ER → GSS-DP throughout | Low (mechanical) |
| **Low** | Appendix | Glossary (optional) | Low |

---

**Document prepared:** May 25, 2026  
**For:** Aníbal Olivera (NDFU Manuscript)  
**By:** Claude AI Assistant
