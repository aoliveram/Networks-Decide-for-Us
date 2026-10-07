# Homophily and higher-order topology are (almost) inseparable in these networks

**Date:** 2026-06-22 · **Status:** finding to revisit later (the 4th cell of the 2×2 design)
**Scripts:** `07_make_HP_optionB_ergm.R`, `08_make_HP_optionA_rewire.R`

## The question we tried to answer

The premium decomposition (GSS → GSS-SH → GSS-DP) is a **sequential** decomposition, not a full 2×2 factorial. We have three cells:

| | topology kept | topology destroyed |
|---|---|---|
| homophily kept | **GSS** | **(missing 4th cell)** |
| homophily destroyed | **GSS-SH** | **GSS-DP** |

The missing cell ("homophily kept, topology destroyed") is what we'd need to claim a clean **topology main effect**. So we tried to build it two ways.

## Option B — re-run the homophily ERGM → FAILS (and proves the point)

The original GSS ERGM (`scripts/02_GSS_network_ergm.R`) uses ONLY homophily terms — `~ edges + nodematch(race,sex,relig) + absdiff(age,educ_num)`. There is **no** `triangle`/`gwesp` term. So clustering is never imposed; it is **emergent from homophily**. Re-simulating the ERGM therefore brings the clustering back for free:

| | density | clustering | modularity | triangles | assort_race |
|---|---|---|---|---|---|
| ER (no homophily) | 0.029 | 0.029 | 0.165 | 3,922 | 0.000 |
| **GSS** | 0.029 | 0.063 | 0.297 | 10,628 | 0.435 |
| **Option B (re-ERGM)** | 0.029 | **0.064** | **0.293** | **10,913** | **0.439** |

→ Re-running the homophily ERGM reproduces homophily AND clustering together. It **cannot** make "homophily without clustering." (Confirms the user's hypothesis: clustering emerges from homophily.)

## Option A — homophily-preserving rewiring → PARTIAL, and exposes the real issue

Degree-preserving double-edge swaps accepted only if **distance-neutral** (|Δ summed edge Gower| ≤ 0.02), so the model's own tie distance d_ij is held while the local wiring is reorganized. Result (GSS vs GSS-HP, 15 nets, ~12% accept rate):

| | clustering | modularity | triangles | **mean edge Gower (model d_ij)** | assort_race | assort_age |
|---|---|---|---|---|---|---|
| GSS | 0.064 | 0.292 | 11,088 | **0.230** | 0.440 | 0.461 |
| GSS-HP | 0.053 | 0.199 | 9,160 | **0.232** | 0.280 | 0.310 |

- **What the model actually uses — the social distance d_ij on the ties — IS preserved** (0.230 → 0.232). Modularity and clustering drop somewhat.
- BUT the per-attribute assortativities fall more than wanted (race 0.44 → 0.28), and clustering only drops partially (not all the way to ~0.029). An earlier variant that *minimized* edge distance did the opposite — it **raised** clustering to 0.18, because packing similar people together creates triangles.

## The deeper lesson (the part worth keeping)

**In empirical homophilic networks, homophily and higher-order topology are nearly the same thing seen two ways.** Putting similar people together *creates* triangles and communities; breaking triangles *loosens* the similarity. The two cannot be cleanly separated — not because the algorithm is weak, but because it is a property of the world. This is arguably a **publishable result in itself**, and it *reinforces* the main claim: you cannot reproduce the realistic topology without the homophily, nor vice versa — they come as a package.

It also means "homophily preserved" has two non-equivalent definitions: (1) the model's tie distance d_ij (Option A preserves this), vs (2) the per-attribute assortativity (sociological standard). They diverge under rewiring.

## Decision pending (post-Sunbelt)

To finish the 4th cell we must choose which "homophily" to pin:
- **Path 1 (fast, mechanism-faithful):** accept Option A as-is — d_ij preserved, topology partially randomized. Defensible if we frame "homophily = the tie distance the model uses".
- **Path 2 (cleaner sociologically, ~1 day):** change acceptance to preserve per-attribute assortativities and raise the swap factor to randomize clustering more; may still hit the homophily↔clustering wall.

Then run the IUL×MSP sweep on the chosen GSS-HP and complete the factorial. **Not needed for Sunbelt** — there we soften the claim to "homophily reproduces ~91% of the penalty on its own."
