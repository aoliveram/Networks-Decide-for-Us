# Using SOCIAL adoption as the order parameter

**Script:** `04_social_order_parameter.R` · **Status:** diagnostic / exploratory
**Data:** GSS and GSS-DP, `mean_0.40`, `sd=0.12`, random seeding.

## The idea (well-motivated)

Total adoption decomposes as $\Phi_{\text{total}} = \Phi_{\text{rational}} + \Phi_{\text{social}} + \text{seeds}$.

- $\Phi_{\text{rational}}$ ($q_i \le \Gamma$): monotone in IUL, **flat in MSP**, topology-blind. A non-critical *background*.
- $\Phi_{\text{social}}$: carries all the network physics (selective, MSP-dependent).

So $\Phi_{\text{social}}$ *should* be the physically correct order parameter, and isolating it *should* sharpen the transition. Confirmed empirically: at IUL=0.2 the rational part is ~0.15 and flat across MSP, while $\Phi_{\text{social}}$ jumps from 0.07 (MSP=0.25) to 0.79 (MSP=0.33).

## What we found

### 1. The susceptibility RIDGE is sharp and clean (visual win)
`plots/social_order_heatmaps.pdf` shows a crisp diagonal susceptibility ridge (a Widom-line) for both topologies — arguably cleaner than the total-adoption version. The social order parameter **locates the critical line very well**.

### 2. The critical point is well-defined and DISCRIMINATES by location
Peak of $\mathrm{SD}(\Phi_{\text{social}})$ sits exactly at the critical point for both:
- GSS: $MSP_c = 0.25$
- GSS-DP: $MSP_c = 0.42$

DP needs *more* "temperature" (openness) to ignite social cascades — a genuine, interpretable difference.

### 3. But the super-critical EXPONENT does NOT cleanly extract (the surprise)
Replicating the published Method-4 (fit $\log \mathrm{SD}$ vs $\log|MSP-MSP_c|$ on the super-critical side) with $\Phi_{\text{social}}$:

| | $\gamma_{\text{social}}$ (SD) | $R^2$ | $\gamma_{\text{social}}$ (Var) | $\gamma_{\text{total}}$ (SD, published) |
|---|---|---|---|---|
| GSS | 0.32 | 0.55–0.73 | 0.63 | **1.08** |
| GSS-DP | 0.39 | 0.55–0.69 | 0.79 | **5.42** |

The social exponent is small, poorly fit, and **does not reproduce the clean 1-vs-4 separation**.

### 4. Why: a noise-floor plateau
On the super-critical (high-MSP) side, $\mathrm{SD}(\Phi_{\text{social}})$ drops then **plateaus at a nonzero floor** (~0.0067 for GSS): deep in the ordered phase social adoption saturates near 0.83 but keeps a constant residual fluctuation. A plateau is not a power law $\Rightarrow$ small slope, poor $R^2$.

In contrast, $\Phi_{\text{total}}$ saturates at ~1.0 and its variance decays cleanly toward ~0.0009, producing the steep, clean power-law tail that yields $\gamma_{\text{total}}$. **So the dramatic 1-vs-4 separation rests substantially on how *total* adoption saturates, not purely on the social channel's critical scaling.**

### 5. A revealing mismatch for DP
For GSS-DP total adoption, the *susceptibility peak* is at MSP≈0.08 while the *jump* in adoption is near MSP≈0.42 — peak and jump are in different places. That mismatch is itself a signature that the DP "transition" is not a clean continuous critical point (consistent with the paper's first-order reading). The social order parameter, by contrast, puts the DP peak cleanly at its $MSP_c$=0.42.

## Honest conclusion

- The user's intuition is **correct in spirit**: $\Phi_{\text{social}}$ is the physically meaningful order parameter and it produces a **sharper susceptibility ridge** and a **well-defined critical line** that discriminates topologies by *location* ($MSP_c^{GSS}=0.25 < MSP_c^{DP}=0.42$).
- BUT the **super-critical SD-tail exponent is not a robust observable for $\Phi_{\text{social}}$** (noise-floor plateau). The headline $\gamma\approx1$ vs $\gamma\approx4$ contrast is a property of the *total* order parameter's saturation dynamics. This is a caveat worth knowing, not a refutation.

## Recommendation
Same root cause as the $\beta$/Widom inconclusiveness: the **13-level MSP grid is too coarse** and the **single-side tail fit is fragile**. To extract a trustworthy social-channel exponent:
1. Dense near-critical MSP sweep ($\Delta \approx 1/48$ in a window $\pm0.1$ around $MSP_c$).
2. Fit the **full peak** (both sides) or use **finite-size data collapse**, not just the super-critical tail.
3. Subtract / model the ordered-phase noise floor before fitting.

For now: use the social susceptibility **heatmap** (sharp ridge) as a strong *visual* in talks; keep $\gamma_{\text{total}}$ as the reported exponent; do NOT claim a social-channel $\gamma$ until the dense sweep is run.
