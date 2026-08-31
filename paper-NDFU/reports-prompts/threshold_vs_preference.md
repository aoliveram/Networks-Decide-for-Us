# Report — Adoption thresholds are NOT individual preferences (and why the rational channel ≠ "a threshold of 1/k")

**Date:** 2026-06-22
**Purpose:** Defend, with citations, the talk/paper line that diffusion threshold models "quietly assume *every* adoption needs peer pressure," and that NDFU's **rational/intrinsic channel** is a genuinely *separate* mechanism — not a relabeling of "a node with 3 neighbors and a fractional threshold of 1/3."

The skeptic's objection (to pre-empt):
> "What's the difference between adding a 'rational choice' channel and just saying that someone with 3 neighbors and a threshold of 1/3 adopts? Isn't your rational adoption just a low fractional threshold?"

**The one-sentence answer:** A fractional threshold $\tau_i$ still operates on **social exposure** $\tilde{E}_i$ (the fraction of *already-adopted* neighbors), so it requires **at least one adopted neighbor** — exposure $> 0$ — to ever fire. NDFU's rational channel ($q_i \le \Gamma$) fires at **exposure $= 0$**: the agent adopts from the innovation's own intrinsic value, with no adopted neighbor required. They are different mechanisms, and the threshold literature itself separates them.

---

## 1. Granovetter (1978): a threshold is a *cost–benefit crossing point*, not a preference

Mark Granovetter, "Threshold Models of Collective Behavior," *American Journal of Sociology* 83(6):1420–1443.

**1.1 The definition (the anchor quote).** (p. 1422)
> "The threshold is simply that point where the perceived benefits to an individual of doing the thing in question exceed the perceived costs."

So a threshold encodes *when the cost–benefit calculation turns positive in a given situation* — it is not a stable "preference." The same person can have different thresholds for different issues, because the cost/benefit structure differs by issue.

**1.2 Thresholds are a combination of costs and benefits, so preference ≠ threshold.** (p. 1422)
> "Since a threshold is the result of some (possibly complex) combination of costs and benefits, two individuals whose thresholds are the same may be politically identical, as reflected in the popular expression, 'strange bedfellows.'"

Two people with *identical preferences/politics* can have *different* thresholds, and two people with *different* preferences can share a threshold. Threshold and preference are **not** the same variable.

**1.3 Knowing preferences is necessary but NOT sufficient.** (p. 1420, abstract + p. 1421)
> "Groups with similar average preferences may generate very different results; hence it is hazardous to infer individual dispositions from aggregate outcomes, or to assume that behavior was directed by ultimately agreed-upon norms."

> "knowing the norms, preferences, motives, and beliefs of participants in collective behavior can, in most cases, only provide a necessary but not a sufficient condition for the explanation of outcomes."

**1.4 The radical / instigator: a threshold of 0% — adoption with NO social input.** (p. 1422)
> "A 'radical' will have a low threshold... Some would be sufficiently radical to have a threshold of 0% — people who will riot even when no one else does. These are the 'instigators.'"

This is the crucial precedent: the threshold framework *already* contains agents who act with **zero adopted neighbors** — purely from their own valuation. NDFU's rational channel is the principled, parameterized version of exactly this "threshold-0 / instigator" idea, driven by the innovation's intrinsic utility $\Gamma$ vs the agent's requirement $q_i$.

**1.5 Behavior can run AGAINST one's own preferences/norms (the car-theft example).** (curated in `theoretical_reviewed.md`)
> "Most did not think it 'right' to commit illegal acts or even particularly want to do so... in our terms, their threshold for stealing cars is low... The boys act because their threshold is exceeded, and their utility is maximized, given the situation, by joining in the criminal activity. But in so doing they act contrary to norms they actually hold."

→ Adoption is a *situated cost–benefit crossing*, not a readout of a fixed preference. Someone averse to a behavior can still adopt it when, for that specific case, the benefit exceeds the cost. (User's own examples: a person generally averse to risky collective action who joins *this* protest because the issue is personally salient; an innovation-averse person who buys *this* phone because it's small and they wanted a small phone.)

**1.6 Average preferences ≠ outcome.** (p. 1421)
> "two crowds whose average preferences are nearly identical could generate entirely different results."

---

## 2. NDFU already commits to this distinction (thesis document, internal consistency)

From *Proyecto de Tesis — Aníbal O.* (CICS, 2025), §3.3 Métodos (p. 10) — the project's own text already defines the rational channel as preference-coherent intrinsic adoption, citing Granovetter 1978:

> "un modelo de influencia interna que hace uso tanto de la **elección racional, entendida como la capacidad de adoptar una innovación o comportamiento que es coherente con las preferencias individuales en el marco de la microeconomía**, como de los hallazgos en influencia social (Granovetter, 1978; Centola, 2013)."

> "Si un agente se informa de una innovación con $q_i \le \Gamma$, la **elección racional se lleva a cabo sin necesidad de influencia social**... En este sentido, la innovación se difunde como un patógeno viral, necesitando solo una interacción para ser adoptado. **Solo cuando $q_i > \Gamma$, la influencia social comienza a jugar un rol.**"

This is the exposure-0 vs exposure->0 distinction in the project's own words: the rational channel = adoption with no adopted neighbor required; the social channel (Eq. 4, with $\tau_i \le \tilde{E}_i \wedge d_{ij}\le h$) only activates when the rational one fails.

**Note on $q_i$ vs preference (subtlety to keep honest):** in the thesis, $q_i$ (the Minimum Utility Requirement) is itself called an individual-level "umbral/requisito mínimo." So NDFU has TWO node-level scalars: $q_i$ (intrinsic requirement, compared to $\Gamma$ → rational channel) and $\tau_i$ (network threshold, compared to social exposure $\tilde{E}_i$ → social channel). The Granovetter point is that the *network threshold* $\tau_i$ is a cost–benefit crossing on the SOCIAL axis, distinct from the *intrinsic requirement* $q_i$ on the PRIVATE axis. The rational channel is the private axis; it is not a fractional value of the social threshold.

---

## 3. Why "rational channel" ≠ "a fractional threshold of 1/k" — the exposure=0 argument

The selective-exposure rule in NDFU (Eq. 4) is
$$\tilde{E}_i \equiv \frac{\sum_{j\neq i}\mathbf{X}_{ij}\tilde{a}_j}{\sum_{j\neq i}\mathbf{X}_{ij}}, \qquad a_i = 1 \iff \tau_i \le \tilde{E}_i \;\wedge\; d_{ij}\le h.$$

- A **fractional threshold** $\tau_i = 1/k$ means "adopt once the fraction of adopted neighbors reaches $1/k$." With $k$ neighbors this needs **$\ge 1$ adopted neighbor**: exposure must be **strictly positive**. If *no* neighbor has adopted, $\tilde{E}_i = 0 < \tau_i$ and the node never fires through this channel.
- The **rational channel** $q_i \le \Gamma$ has **no $\tilde{E}_i$ term at all**. It fires at $\tilde{E}_i = 0$. The trigger is the *intrinsic* value of the item, not any neighbor.

So the two channels are distinguishable by a clean, testable property: **does the node adopt when surrounded entirely by non-adopters?** Fractional threshold: no. Rational channel: yes (if $\Gamma \ge q_i$). This is precisely Granovetter's "instigator / threshold-0" agent (§1.4) elevated to a parameterized intrinsic-utility mechanism.

This also maps onto the **simple vs complex contagion** boundary: the rational channel behaves like *simple* contagion / external influence (one source — here, zero sources — suffices), while the social channel is the *complex* threshold mechanism. They are not the same object folded together.

---

## 4. External literature: thresholds ≠ preferences, and the two-channel separation (web-verified)

### 4.1 More Granovetter quotes that sharpen the point (page-exact, cross-corroborated)

- **Situation-specificity** (the strongest support for "averse in general, but positive for THIS issue"), p. 1436:
  > "Thresholds are situation-specific. An individual's riot threshold is not a number that he carries with him from one riot to another but rather results from the configuration of costs and benefits, to him, of different behaviors in one particular riot situation."
- **Threshold is purely behavioral, not normative**, p. 1435:
  > "The concept of threshold, then, is purely behavioral, connoting nothing about what the actor thinks is the 'right' thing to do."
- **Preferences known ⇒ still not enough**, p. 1420:
  > "outcomes cannot be determined by any simple counting of preferences."
- **The action-not-taken case** (directly relevant to NDFU low-IUL story), p. 1425:
  > "Threshold models may be of particular value in understanding situations where the average level of preferences clearly runs strongly in favor of some action, but the action is not taken."

### 4.2 Later threshold work explicitly separating thresholds from preferences

- **Macy (1991), "Chains of Cooperation: Threshold Effects in Collective Action," *American Sociological Review* 56(6):730–747** — the single best secondary statement (p. 731):
  > "Sudden changes in the level of production of a public good need not reflect similar changes in the overall preferences of the actors. Attention should focus instead on the distribution of thresholds and the social ties through which members learn about the actions of others."
- **Young (2009), "Innovation Diffusion in Heterogeneous Populations: Contagion, Social Influence, and Social Learning," *American Economic Review* 99(5):1899–1924** — distinguishes adoption driven by *inertia/own payoff* from adoption driven by *conformity/contagion*; a node can adopt because the innovation's own payoff is high, separately from how many neighbors adopted.
- **Centola & Macy (2007), "Complex Contagions and the Weakness of Long Ties," *AJS* 113(3):702–734** + review **Guilbeault, Becker & Centola (2018), "Complex Contagions: A Decade in Review"** (Springer; arXiv:1710.07606): the simple-vs-complex distinction is about *number of required sources of exposure*; a one-source (or zero-source) trigger is categorically different from a fractional-threshold trigger.

### 4.3 The decisive formalization of the exposure=0 vs fractional-threshold distinction

**Karsai, Iñiguez, Kaski & Kertész (2016), "Local cascades induced global contagion," *Scientific Reports* 6:27272 (arXiv:1601.07995)** — states the user's exact distinction in the Watts-threshold formalism:
  > "[Watts] identified **innovator nodes that spontaneously change state to 1**... Such nodes have **a trivial threshold φ = 0.** Then there are nodes with threshold **0 < φ ≤ 1/k, called vulnerable, which need one adopting neighbour before their own adoption.**"

This is the cleanest citation for the slide line: a "3 neighbors, threshold 1/3" node is *exactly* the **vulnerable** case (1/3 = 1/k, k=3) — by definition it **needs one adopting neighbour** (exposure > 0). The intrinsic/rational adopter is the **φ = 0** case: it acts at exposure 0. Same paper, empirically: spontaneous (φ=0) adopters "govern the global emergence of social spreading," while "vulnerable adoptions — induced by a single adopting neighbour — appear to be important only locally." So the two are not interchangeable: moving an exposure-0 driver into an exposure-positive threshold changes the network dynamics.

### 4.4 Two-channel adoption models (intrinsic/private term + social term) — the constructive precedent

The NDFU rule (rational $q_i\le\Gamma$ OR social $\tau_i\le\tilde E_i$) is in a well-established family where the private channel is a **separate additive component**, not a neighbor count:

- **Bass (1969), "A New Product Growth for Model Consumer Durables," *Management Science* 15(5):215–227** — the foundational split: hazard $= p + q\,F(t)$. The **innovation coefficient $p$ (external/intrinsic)** drives adoption *regardless of how many adopted* (fires at $F=0$); the **imitation coefficient $q$ (internal/social)** scales with the adopter fraction. NDFU's rational channel is the $p$ term; the social channel is the $q$ term.
- **Brock & Durlauf (2001), "Discrete Choice with Social Interactions," *Review of Economic Studies* 68(2):235–260** — the rigorous random-utility decomposition: utility $V(s_i)$ = **private payoff $u_{s_i}$** + **social/conformity term $\sum_j J_{ij}\,s_i\,\mathbb{E}[s_j]$** + noise. This is the formal proof that the intrinsic channel is an **additive private utility**, not a $1/k$ threshold. (Random-utility logic ← McFadden.)
- **Mas Tur, Zeppini & Frenken (2018, *Phys. Rev. E* 97:022302; 2024, *Social Networks* 76:12–21)** — the same intrinsic-preference ("minimum utility requirement") + social-reinforcement split, **on networks, and published in the NDFU target journal**. 2018 abstract: "individual preferences cannot be overlooked." This is the closest sibling to NDFU's hybrid rule.
- **Salganik, Dodds & Watts (2006), "Experimental Study of Inequality and Unpredictability in an Artificial Cultural Market," *Science* 311:854–856** — the *experimental* warrant that an intrinsic channel is real and measurable: "appeal" = market share in the **independent (no-social-influence) condition**, separated from "success" = share in the social-influence worlds. NDFU's IUL is that intrinsic appeal.

### 4.5 The referee's objection, answered

A referee may say: *"any rational rule can be folded into a threshold distribution — observational equivalence."* Two answers:
1. **Mechanistic, not cosmetic.** The two channels produce *different network dynamics* (Karsai 2016): φ=0 drivers are structure-agnostic and seed global cascades; φ>0 vulnerable adoptions matter only locally. NDFU's low-IUL result — rational adopters igniting cascades a purely-social model would call impossible — is exactly this difference made visible.
2. **It is the established decomposition.** Bass, Brock-Durlauf, Tur et al. all model the private channel as a *separate additive term* (utility / hazard / requirement), and recent threshold-utility reinterpretations decompose the threshold into **"attribute" (intrinsic) + "social" utilities** (e.g. arXiv:2405.13224). The field favors the two-channel decomposition over a bare neighbor count; NDFU follows that, it does not invent it.

---

## 5. Drop-in defenses (for the talk / for a referee)

- **One-liner for the slide ("The gap"):** *"In threshold models, a node only adopts once enough neighbors already did — adoption is gated on social exposure. We add a separate, intrinsic channel: a node can adopt on its own when the thing is good enough for it, with no adopted neighbor required — what Granovetter (1978) called the 'instigator' with a zero threshold."*
- **Answer to "isn't rational adoption just a 1/k threshold?":** *"No — a 1/k threshold still needs at least one adopted neighbor; it's gated on exposure. Our rational channel fires at exposure zero, from the item's intrinsic value. Granovetter's own threshold of 0% — the instigator who acts when no one else does — is exactly this, and he's explicit that a threshold is a cost-benefit crossing point, not a preference."*
- **Honest caveat to volunteer:** thresholds and intrinsic requirements are both node-level scalars, so a referee may say "any rule can be folded into a threshold distribution." The reply: the *mechanistic* content is the exposure=0 vs exposure>0 distinction — they produce different dynamics (the rational channel seeds cascades a purely-social model would call impossible), which is empirically the point NDFU makes at low IUL.

---

### Sources

**Local (read directly):**
- Granovetter (1978) PDF — `A - Trabajo 1/A - Papers - 1/Granovetter (1978)...pdf` (pp. 1420–1423 read here; pp. 1424–1436 quotes web-verified against two independent scans)
- Thesis — `A - Proyecto de Tesis - Documento/Proyecto de Tesis - Aníbal O.pdf` (§3.3 Métodos, p. 10)
- Curated quotes — `paper-NDFU/theoretical_reviewed.md` (entry 4)

**External (web-verified, cross-checked ≥2 sources):**
- Granovetter, M. (1978). Threshold Models of Collective Behavior. *AJS* 83(6):1420–1443.
- Macy, M. (1991). Chains of Cooperation: Threshold Effects in Collective Action. *ASR* 56(6):730–747.
- Young, H. P. (2009). Innovation Diffusion in Heterogeneous Populations. *AER* 99(5):1899–1924.
- Centola, D. & Macy, M. (2007). Complex Contagions and the Weakness of Long Ties. *AJS* 113(3):702–734.
- Guilbeault, Becker & Centola (2018). Complex Contagions: A Decade in Review. arXiv:1710.07606.
- Karsai, Iñiguez, Kaski & Kertész (2016). Local cascades induced global contagion. *Sci. Rep.* 6:27272 (arXiv:1601.07995). ← the φ=0 vs φ∈(0,1/k] formalization.
- Bass, F. (1969). A New Product Growth for Model Consumer Durables. *Manage. Sci.* 15(5):215–227.
- Brock, W. & Durlauf, S. (2001). Discrete Choice with Social Interactions. *Rev. Econ. Stud.* 68(2):235–260.
- Mas Tur, Zeppini & Frenken (2018). *Phys. Rev. E* 97:022302; (2024) *Social Networks* 76:12–21.
- Salganik, Dodds & Watts (2006). Experimental Study of Inequality and Unpredictability... *Science* 311:854–856.

**Verification note:** Granovetter/Macy/Young/Watts-taxonomy quotes were extracted verbatim from primary PDFs and page-anchored (Granovetter cross-corroborated against gwern.net + CUHK scans). Brock–Durlauf and Salganik operationalizations were confirmed from secondary academic sources after the primary PDFs returned paywalled/non-extractable content. Centola & Macy (2007) page numbers given at article level (pre-publication manuscript carries no journal pagination).
