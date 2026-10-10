# Data sources for RISKY / high-cost / violent collective action (for q_i)

**Motivation (point 9).** The current propensity scores $q_i$ come from GSS/ATP items that cover only **low-risk** collective action (sign petition, boycott, peaceful demonstration) and **medium-risk** innovation adoption. To model *risky* collective action (violent protest, militant activism), we need either a survey that measures willingness for high-cost action, or a principled assumed distribution.

## 1) Promising US data sources found (web research, May 2026)

**Best instrument (item template for $q_i$):**
- **ARIS — Activism–Radicalism Intention Scales** (Moskalenko & McCauley 2009). Graded behavioral-*intention* scale: AIS = low/medium-risk activism; **RIS = violent/illegal high-cost action**. Spans the full risk gradient in one tool; items freely available (Table 2). Limitation: validated on non-representative (student) samples → use for item wording / risk-ordering, calibrate the distribution against a representative anchor.

**Best nationally-representative dataset (to calibrate the violent tail):**
- **UC Davis "Life in America" Nationwide Survey on Political Violence** (Wintemute et al.). True *graded personal-willingness* battery (damage property → threaten → injure → kill; firearm-use expectation) on a large probability panel weighted to the CPS. Microdata require a data-use agreement.

**Public cross-checks / bridges:**
- **CPOST** (Chicago Project on Security & Threats) — quarterly nationally-representative toplines on support for political violence; reports public.
- **Bright Line Watch** — political-violence items with list-experiment correction for over-statement; reports public.
- **WVS Wave 7 (US, 2017–2022)** — political-action battery incl. *unofficial strikes* and *occupy buildings* (highest-risk it reaches: illegal/disruptive but non-violent), coded *have done / might do / would never*. **Fully public** — the natural low-/medium-risk bridge to our existing GSS/ATP items.
- **ANES / CES (CCES) / Nationscape** — carry violence-*justification* attitudes (often Kalmoe–Mason modules) but **no graded personal-willingness battery**; fully public.
- **Tausch et al. (2011)** — normative vs. non-normative action measurement template (already cited in the project); UK/India/student, not US public data → measurement model, not data source.

**Recommended combination:** ARIS/RIS for item structure & risk ordering of $q_i$; calibrate the violent-tail prevalence with UC Davis "Life in America" (CPOST toplines as public cross-check); WVS-7 strike/occupy items as the public low-risk floor that connects to the current GSS/ATP scores.

## 2) Fallback / immediate computational step: assumed skewed aversion

Independent of acquiring new data, the natural first experiment is to **assume** an aversion distribution for risky action that is strongly right-skewed — most agents are highly averse (high $q_i$), a small tail willing. Concretely, draw $q_i$ for "risky collective action" from e.g. a Beta$(\alpha,\beta)$ with mean $\sim0.8$–$0.9$ and a thin low-$q$ tail (e.g. Beta(5,1.2) or a truncated/exponential-tail shape), assign to nodes, and re-run the full IUL×MSP sweep. Questions to ask:
- Does the structural premium survive / grow when the population is mostly averse?
- Does the critical line $MSP_c$ shift (need much higher openness to ignite)?
- Does the universality contrast (plausible $\gamma\approx1$ vs random $\gamma\approx4$) persist under a skewed, high-aversion population?

This is **deferred** (not now). It is a clean robustness/extension experiment and a good thesis-chapter hook: "what does it take, structurally, to ignite *risky* collective action when almost nobody wants to?"

**Empirically-grounded variant:** instead of an arbitrary Beta, fit the skew to the ARIS-RIS / CPOST / UC-Davis prevalence estimates so the assumed distribution matches measured US willingness-for-violence base rates.
