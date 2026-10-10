# NDFU — Journal Fit Analysis (for resubmission + Sunbelt 2026 framing)

**Date:** 2026-06-18. Web research on SCImago metrics, aims & scope, representative work, and fit.
**Strategic anchor:** the *real* reason the Social Networks ABM special issue rejected the abstract (Federico Bianchi's 2nd email) was **NOT volume** — it was *"the underdevelopment of the behavioural model… the behaviour of the agents is relatively abstract, and the experimental design appears to be governed by macro-level manipulation rather than bottom-up dynamics."* This critique recurs at ABM-centric venues and must drive journal choice.

---

## 1. Summary metrics table

| Journal | SJR (2024) | Best quartile | H-index | CiteScore | JIF (approx) | OA / cost to publish | Fit score |
|---|---|---|---|---|---|---|---|
| **Social Networks** (Elsevier/INSNA) | ~1.17 | **Q1** Sociology & Political Science | **120** | 6.3 | ~2.9 | Hybrid; **free** on subscription track (Gold OA = $3,710) | **8.5** |
| **Sociological Science** (Soc. for Soc. Science) | 1.197 | Q1 Soc. Sciences (misc.) | 45 | 6 | ~2.0 (Scopus proxy) | Full OA, **rank-based fee** (~$575; pushed to professor tier by senior co-authors) | 7 |
| **JASSS** (U. Surrey) | 0.586 | Q1 Soc. Sciences (misc.) | 67 | 5.6 | ~2.1 | **Diamond OA — free** | 7 |
| **Rationality and Society** (SAGE) | 0.602 | Q1 Soc. Sci (misc.) / Q2 Soc. & Pol. Sci. | 54 | 1.6 | ~1.5 | Hybrid; **free** on subscription track | 7 |

*(SCImago's own pages 403'd automated fetch; values triangulated from SCImago-sourced mirrors — verify on scimagojr.com before quoting in a CV.)*

## 2. THE decisive variable — editorial overlap with the people who rejected you

The four editors who signed the rejection (Bianchi, **Flache**, **Squazzoni**, Takács) sit on several of these boards. This determines how hard Federico's critique bites:

| Journal | Editor(s) | Same editorial mind as the rejection? | Does the "macro-manipulation / abstract agent" critique bite? |
|---|---|---|---|
| **Social Networks** | **Thomas Valente** (your co-author!) + Ulrik Brandes | No (Valente is on your side; COI must be declared, handled by Brandes/AE) | **Neutralized** — SNA journal, not ABM journal; Tur-Zeppini-Frenken (2024) published an abstract utility-threshold diffusion model here |
| **Sociological Science** | Independent editorial team | **No** — fully outside the four | **Low** — Grow-Flache-Wittek (2017) is itself a macro-sweep ABM with stylized agents, published here |
| **JASSS** | **Flaminio Squazzoni** (signed the rejection) | **YES** | **HIGHEST** — bottom-up emergence is the journal's *house epistemology*; the same gatekeeper, same standard |
| **Rationality & Society** | **Andreas Flache** (signed the rejection) | **YES** | **Partial** — same desk, but R&S accepts spare/abstract rational-action agents as a feature; the "abstract agent" half does not bite, the "macro-manipulation" half still does |

**Implication:** the two journals where the critique is weakest (Social Networks, Sociological Science) are *also* the two outside the rejecting editorial circle. The two where the same editors preside (JASSS, R&S) are where you must most aggressively re-spine the paper bottom-up.

## 3. Per-journal: lead result + reframing

**Universal reframing across all four:** LEAD with **Result 1 (the structural premium, −57%)** + the **homophily-imputation pipeline**; keep Result 2 (criticality) as a short sociological "tipping points: continuous vs. abrupt" section; **demote Result 3 (universality / γ,δ / Ising / Widom) to a supplement/appendix** — for none of these four audiences is the statistical-physics framing the hook (it ranges from curiosity to liability).

- **Social Networks (8.5 — top pick).** Mission ("how networks emerge, evolve, and have consequences for behaviour") maps 1:1 onto the structural-premium result; ERGM + threshold/exposure contagion + seeding *are* its methods. Title angle: *"Networks Decide for Us: imputing empirically homophilic whole-networks onto survey respondents and isolating the causal premium of social structure on diffusion."* Sell the pipeline as a reusable measurement contribution + the DP-randomized baseline as a clean counterfactual. Caveats: (i) Valente EIC ⇒ declare COI; (ii) referees will likely demand the imputed networks reproduce held-out tie features (validation).
- **Sociological Science (7, behaves like 8 with reframing).** General-sociology, Q1, **30-day decisions**, no length limit, evaluative (non-developmental) review. Pre-empt two referee asks: give the MSP/threshold rule a behavioral micro-justification (homophilous attention/trust); foreground the DP baseline as the identification device. Title angle: *"How Much of Diffusion Is the Network? Quantifying the Structural Premium of Homophily with Imputed Whole-Networks."* Cost: rank-based fee (~hundreds USD, professor tier via senior co-authors).
- **JASSS (7 — highest risk).** Topically perfect (diffusion, seeding, homophily, threshold ABMs, empirical grounding), diamond OA (free). BUT editor = Squazzoni (signed rejection) and bottom-up emergence is the gatekeeping criterion. Only submit with a deliberate generative reframing: recast the IUL×MSP sweep as exploring a *micro-rule's* parameter space; trace an explicit micro→macro emergence path (proximity-gated exposure → structural premium); lean on the empirical homophily imputation to rebut "abstract agents." Submitted as-is it would likely draw the same objection again.
- **Rationality & Society (7 — modest impact).** Scope invites computational/complexity contributions to the rational-action paradigm; the dual rational-choice/social-influence agent is in its tradition (Braun 1995; Esser frame selection). Editor = Flache (signed rejection), but R&S *welcomes* abstract agents. Re-spine around the micro decision rule, present IUL×MSP as comparative statics, cut phase-transition language. Lowest JIF of the four (~1.5).

## 4. Bottom line recommendation

1. **Primary target: Social Networks (regular manuscript).** Best fit, highest prestige (H-120, JIF ~2.9), critique neutralized, co-author is EIC, Federico himself suggested it. Lead with pipeline + structural premium; demote physics.
2. **Strong second / faster: Sociological Science.** Independent of the rejecting circle, fast OA, accepts macro-sweep ABMs. Excellent if you want a quick, clean home and to avoid the COI dance.
3. **JASSS only after a genuine bottom-up re-spin** (and the MSP-heterogeneity extension would directly answer Squazzoni).
4. **R&S** as a fallback if you want the analytical-sociology/rational-action framing and are okay with modest impact.

**For Sunbelt 2026 (classic SNA audience, "Contagion and Diffusion processes through Social Networks", 15+5 min):** present the *Social Networks* framing — pipeline + structural premium as the spine, dynamics as a short coda, physics essentially dropped or one teaser slide. This both fits the audience and doubles as a dry-run of the SN submission.

## 4b. Mathematical intensity — can SN / Sociological Science take a math-heavy paper? (verified)

Direct comparison of two papers from the user's own library, one per journal:
- **Social Networks — Tur, Zeppini & Frenken (2024):** same genre as NDFU (homophily + social reinforcement + threshold/utility diffusion). Comfortably formal: explicit "second-order critical transition" (Stauffer & Aharony 1994), percolation threshold, low/high-diffusion phases, numbered equations. **SN tolerates phase-transition language and moderate explicit formalism in the body** (a few equations per section, with prose).
- **Sociological Science — Macy & Evtushenko (2020) "Threshold Models of Collective Behavior II":** the genre-twin of NDFU in that journal (thresholds, cascades, "phase transition", "critical mass", predictability, Salganik-Dodds-Watts). **Prose-forward**: the cascade condition is stated in italic prose, not as a display equation; very few equations in the body.

Verdict (from submission guidelines + corpus search):
- **Both publish the SUBSTANCE** (phase transitions, critical mass, cascades). Neither rejects the physical *concept*.
- **SN accepts more explicit formalism** (percolation, exponents, equations) in the body → the ERGM pipeline + continuous/discontinuous transition can live in the main text.
- **Sociological Science wants the physics in PROSE**, heavy derivations (Ising/Landau/Widom, γ/δ) in the **supplement** (unlimited supplementary material; "encourages but does not require formalization"; brevity prized; ~15–30 pp). **No Sociological Science article uses Ising/percolation/temperature language overtly** — so translate physics vocabulary into analytical-sociology terms in the body.
- **Neither** abstract should open with γ/δ/Widom.

**Must-cite (Sociological Science precedents that de-risk the framing):**
- Macy & Evtushenko (2020), *Soc. Sci.* 7:628–648 — tipping/threshold/predictability twin; its thesis ("individual idiosyncrasy matters; don't infer individual dispositions from aggregate outcomes") dialogues directly with NDFU. Editors Jesper Sørensen / Olav Sorenson handled it.
- Grow, Flache & Wittek (2017), *Soc. Sci.* 4:611–640 — macro-sweep ABM on clustered networks with stylized agents → structural twin of NDFU's design, accepted without the macro-manipulation objection.
- Arvidsson, Hedström & Keuschnigg (2025), *Soc. Sci.* 12:715–742 — ABM + diffusion, preferences-vs-outcomes decoupling (analytical sociology register).
- Lee, Lazer & Riedl (2025), *Soc. Sci.* 12:685–714 — complex contagion / social reinforcement, country-scale field experiment.
- Berry, Sirianni, Weber, An & Macy (2021), *Soc. Sci.* 8:285 — estimating homophily from dyadic predictions (methodological homophily).

## 4c. Haiko Lietz (CPIN organizer who praised NDFU & suggested Sociological Science)

GESIS computational social scientist; sociologist + engineer; dissertation "Scale-free identity: the emergence of social network science" (arXiv:2403.07595). **His central project IS mapping statistical-physics concepts onto society**: percolation, scaling laws, self-similarity/fractals, universality class, criticality, phase transitions, order/disorder control parameters, emergence — the CPIN agenda is his. He theorizes identity as residing "at the phase transition between stability and change of meaning"; the NetSci-2025 CPIN report (Lietz et al.) argues different social networks "belong to the same universality class, as physicists say."
- He does NOT use the literal word "temperature" or "free energy" as a named social variable, nor a Bourdieu field-theory formalization under his own name → NDFU's "temperature" metaphor is a natural *extension* of his vocabulary, not a verbatim match. His enthusiasm fits because order/disorder control-parameter framing is exactly his program.
- Publishes in: *International Journal of Sociology*, *Quantitative Science Studies*, *Scientometrics*, *Advances in Complex Systems*, ICWSM, arXiv. **Not found** in Social Networks / Sociological Science / EPJ Data Science — so his Sociological Science suggestion reflects taste/fit, not his own publishing track there.
- **Strategic read:** Lietz is a natural ally/champion and a potential reviewer-type for the criticality angle; but he is not a gatekeeper at any of the four target journals. His endorsement matters most for the *CPIN community* and for keeping a criticality-forward version alive (e.g., *Advances in Complex Systems*, *EPJ Data Science*) as a parallel track — NOT for the SNA-classic resubmission, where the physics is demoted.

## 4d. IUL micro-foundation via pairwise comparison (the user's note — verified literature)

The user wants 1 line justifying "IUL = the option that wins the most pairwise comparisons." This is licensed by **random-utility / comparative-judgment theory**:
- **Thurstone (1927)** law of comparative judgment; **Luce (1959)** choice axiom; **Bradley & Terry (1952)** BTL model; **McFadden (1974)** conditional logit — a mutually-equivalent apparatus (unified by **Yellott 1977**) in which options sit on ONE latent utility continuum and P(i≻j) is a monotone function of the IUL difference, so the highest-IUL option is by construction the one that wins the most pairwise comparisons.
- **Salganik, Dodds & Watts (2006), Science 311:854–856** — canonical diffusion-context separation of *intrinsic quality/appeal* (no-social-influence condition) from *socially-driven success*; exactly the role IUL plays vs. adoption.
- **Drop-in sentence:** "We treat each innovation's IUL ∈ [0,1] as a latent univariate quality parameter in the random-utility / comparative-judgment sense (Thurstone 1927; Luce 1959; Bradley & Terry 1952): assuming options lie on a single utility continuum, the highest-IUL option is by construction the one that wins the most pairwise comparisons, with intrinsic quality deliberately distinguished from socially-driven adoption success (Salganik, Dodds & Watts 2006)."
- **Caveats to state:** a scalar IUL exists only under single-dimension / transitivity (else Condorcet cycles, Arrow 1951); and Luce's IIA breaks under similarity/decoy effects (Tversky 1972) — acceptable for a fixed exogenous IUL on distinct non-substitutable innovations, but flag it. This is the "1 valuable line" the user wanted — best placed where IUL is first defined in the manuscript.

## 5. Representative precedents to cite/emulate (de-risk the framing)
- **Social Networks:** Tur, Zeppini & Frenken (2024) — abstract utility-threshold diffusion model, published as "a theoretical model" → kills the "too abstract" objection.
- **Sociological Science:** Grow, Flache & Wittek (2017), *Global Diversity and Local Consensus in Status Beliefs* — macro-sweep (clustering × resistance) ABM with stylized agents → structural twin of NDFU's design, accepted without the macro-manipulation objection.
- **Sociological Science:** González-Bailón et al. (2024), diffusion/reach of (mis)information on Facebook → journal publishes real-network diffusion.
- **Adjacent taste-setters (AJS/ASR):** DellaPosta-Shi-Macy "Why Do Liberals Drink Lattes?"; Goldberg-Stein associative diffusion — the homophily+influence sub-community NDFU speaks to.
