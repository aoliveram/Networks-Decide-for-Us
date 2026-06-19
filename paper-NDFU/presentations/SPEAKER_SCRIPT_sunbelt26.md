# Speaker Script — Sunbelt 2026 (≈15 minutes + 5 Q&A)

*Session: "Contagion and Diffusion processes through Social Networks." Classic-SNA audience. Everyday language, like explaining the study to a friend. Lead with the pipeline + the structural premium; the physics is just one teaser slide. Times are cumulative.*

---

## Slide 1 — Title (0:00–0:30)

> Hi everyone, thanks for being here. I'm Aníbal Olivera, a PhD candidate at the Center for Social Complexity in Chile. This is joint work with Jorge Fábrega and Tom Valente.
> The talk is about a simple but stubborn question: how much of where a behavior spreads is about *who people are*, and how much is about *how they're wired together*?

---

## Slide 2 — Collective behavior spreads on networks (0:30–1:15)

> Here's the setup. Behaviors — adopting a new product, joining a protest — spread person to person, like a contagion. And real human networks have a very particular texture: we tend to cluster with people like us. Same age, same education, same background. That's homophily.
>
> So a natural question hangs over all of this: how much of *where* a behavior ends up is about *who* people are, and how much is about *how* they're wired together?

---

## Slide 3 — The gap (1:15–2:15)

> Let me put three standing problems on the table. One: diffusion models usually use synthetic networks — small-world, scale-free — not real ones. Two: the surveys that *do* measure propensities, like the ATP or the GSS, have no network data at all.
>
> And three — this is the one I want you to hold onto. A lot of threshold models quietly assume that *every* adoption needs peer influence. And there's a beautiful paper by Macy and Evtushenko in Sociological Science: when a group with weak interest somehow still tips, the explanation they reach for is *stochastic noise* — a random spark that spontaneously instigates the cascade.
>
> Our bet is different. We want to see whether *measured preference* and *real structure* — not noise, not made-up graphs — can do that work.

---

## Slide 4 — Where we go (2:15–3:00)

> So our paper makes two moves. First, we *impute* a plausible whole-network onto real survey respondents — we build them a realistic social network out of empirical homophily. We're building on McPherson and Smith's idea of imputing social context from surveys, but we push it from a single latent "context" number all the way to full, simulable networks.
>
> Second, we measure how much that structure actually *causes* diffusion — and then we dig into *what about the structure* matters.
>
> The precise question, in the box: holding individual propensities fixed, does the *shape of the network alone* — not who its members are — decide how much spreads?

---

## Slide 5 — Adoption, part 1: Rational Choice (3:00–3:50)

> The model has two ways a person adopts. The first is simple and selfish: if the thing is attractive enough to me personally, I adopt it on my own. We call the attractiveness of the behavior its "Intrinsic Utility Level" — think of it as a quality score. And each person has a personal bar it has to clear — their minimum requirement.
>
> The key thing about this first channel — look at the little networks — is that it completely ignores the network. People adopt on their own, wherever they happen to sit.

*(Point at the three snapshots: black = adopted.)*

---

## Slide 5 — Adoption, part 2: Selective Social Influence (3:30–4:45)

> The second way is social: I adopt because enough of my contacts already did. Classic threshold idea. But here's the one twist that matters: **I don't count everybody.** There's a limit on how *socially different* a neighbor can be before I stop letting them influence me. If someone is too far from me in the demographic sense, their adoption just doesn't register.
>
> We call that limit MSP — maximum social proximity. When it's low, I only listen to people just like me — echo chambers. When it's high, I'm open to influence from very different people. So influence is filtered by social distance — and that's where homophily enters the dynamics.

---

## Slide 6 — The imputation pipeline (4:45–6:45)  ← the methodological heart

> Okay, this is the engine of the paper, and the part I most want you to take home. Three steps.
>
> Step one: we borrow *empirical homophily strengths* — how strongly people actually sort by race, religion, age, education, sex — from large US studies where that's been carefully measured. Step two: we feed those into an ERGM to *grow* whole synthetic-but-realistic networks. Step three: we attach real survey answers — innovativeness from the ATP, collective-action propensity from the GSS — to each node.
>
> The tie probability is just a function of distance in Blau space — closer people are more likely to be connected. And the payoff is on the right: these imputed networks reproduce the lumpy, right-skewed degree structure of real social networks. The textbook random models simply don't do that.
>
> So now we have something rare: realistic networks *and* real individual propensities, for people who never gave us their actual ties.

---

## Slide 7 — Breaking structure: what it destroys, part 1 (topology) (6:45–7:45)

> Before I run anything, I want to ask a sharper question. Our counterfactual is going to be a "scrambled" network — same number of ties per person, but rewired at random. That scrambling breaks *a lot* of things at once. So: what exactly does it break? We measured it, across a hundred real networks and a hundred scrambled twins.
>
> This first plot is the *topology*. Scrambling damages it — clustering drops about a third, triangles drop a third, the community structure, modularity, drops by half. So topology is a suspect. But notice the path length barely moves, and — not shown here — degree, density and the giant component are *identical*. So this isn't a size artifact; it's a structural one.

---

## Slide 8 — Breaking structure: what it destroys, part 2 (homophily) (7:45–8:45)  ← the key one

> But here's the dramatic one — and it's the cleaner culprit. On the left: homophily by each demographic axis — age, education, race, religion, sex. In the real networks it's strongly positive: similar people are connected. In the scrambled ones, every single one collapses to *zero*. Random mixing.
>
> And on the right is the cleanest number in the whole study: the average social distance *along the ties* — the exact quantity my influence rule uses — jumps 42 percent. In a real network your neighbors are *similar* to you; scramble it, and your neighbors become *random strangers*. That's the mechanism, in one number. So when we destroy "structure," what we mostly destroy is homophily. Keep that in mind for the next slide.

---

## Slide 9 — Simulation design (8:45–9:30)

> Okay, now we run it — a lot. We sweep the attractiveness of the behavior, the openness parameter MSP, the thresholds, and different ways of seeding the first adopters. About four million simulations, on the real networks versus their scrambled twins. Then we fit a big additive model to pull out the effect of *just* the structure.

---

## Slide 10 — Result 1: The structural premium (9:30–11:00)  ← headline finding

> And here's the headline. When we scramble the network — destroy the homophily but keep everything else — adoption odds drop by about **57 percent**.
>
> Let that sink in: same people, same preferences, same thresholds. The *only* thing that changed is the wiring — and, as we just saw, mostly the homophily. And it cuts adoption by more than half.
>
> Why? Because homophilic clusters are echo chambers, and echo chambers turn out to be *useful* here — they let social reinforcement build up safely in a corner until it's strong enough to ignite. So homophily isn't a barrier to spread. It's a catalyst.

---

## Slide 11 — Conclusions (11:00–12:15)

> So, to wrap up — three takeaways.
> A *method*: you can impute realistic, *simulable* whole-networks onto survey respondents from empirical homophily. That's reusable anywhere you've got propensities but no ties.
> A *finding*: holding individual preferences fixed, real structure gives you a 57-percent adoption premium.
> And a *mechanism*: what "breaking structure" destroys, above all, is homophily — your neighbors stop being similar to you.
>
> And let me close the loop with Macy and Evtushenko, who I mentioned at the start. They explain surprising cascades with *stochastic noise*. Here, the thing that drives the outcome isn't noise — it's *measured preference* and real structure. The same agents, on a different topology, produce a different society.
>
> The big message: where a behavior ends up isn't the sum of individual preferences. It's governed by the structure of our ties. Networks decide for us.

---

## Slide 12 — Thanks / References (12:15–12:45)
*(Stay on Thanks; take questions. QR links to the preprint. Finishing a bit early leaves room for Q&A — and for the backup slide if there's appetite.)*

> Thank you — the QR there links to the draft, and I'm very happy to take questions.

---

## Backup slide (only if time remains, or in Q&A) — the shape of change

> If there's time, one teaser. Beyond *how much* spreads, the structure also governs *how* the change arrives. On a realistic network the tipping is smooth and predictable — a continuous transition. On a scrambled one, it's abrupt — all-or-nothing.
>
> And there's a deeper signature: near the tipping point, the realistic networks show what physicists call *universality* — mean-field scaling, the same math as a phase transition in matter. In the table, the susceptibility exponent gamma is about 1 for the real networks and about 4 for the scrambled ones — so gamma *tells them apart*. The field exponent delta is about 3 for both, so it doesn't. That's a whole companion paper, and I'd love to get into it.

---

### Timing cheat-sheet
| Block | Target |
|---|---|
| Motivation + gap + where-we-go | 0:00–3:00 |
| Model (2 slides) | 3:00–4:45 |
| Pipeline | 4:45–6:45 |
| Diagnosis: topology (B) + homophily (C) | 6:45–8:45 |
| Simulation design | 8:45–9:30 |
| Result 1 (premium) | 9:30–11:00 |
| Conclusions (incl. Macy loop) | 11:00–12:15 |
| Thanks / Q&A | 12:15–15:00 |
| *(Backup: shape-of-change + γ/δ — only if asked)* | — |

**Note:** Macy is now *front-loaded* in "The gap" (as motivation) and *closed* in Conclusions — no standalone middle slide. The physics is a single **backup slide after Thanks**, used only if time/appetite remains or in Q&A. This deliberately finishes ~12–13 min, leaving generous Q&A.
**Release valve if running long:** compress the topology slide (B) to two sentences — "scrambling also hurts clustering and modularity, but the controls are identical, and homophily is the most complete collapse" — and move straight to the homophily slide (C). The diagnosis *sets up* the premium, so it should feel like one continuous build.

### Anticipated questions
- *"Are the imputed networks validated?"* → they reproduce the right-skewed degree distribution and realistic clustering of real networks (the degree-dist figure); held-out tie-feature validation is the natural next robustness check.
- *"Isn't MSP a macro knob you're just turning?"* → yes, in this version it's a global control parameter; we read it as a population's *social porosity*. Making it heterogeneous and individual-level is the next step (and the honest answer to the bottom-up critique).
- *"Why degree-preserving randomization and not something gentler?"* → it's the cleanest counterfactual that holds degree fixed; the attribute-shuffle null (identical topology, permuted attributes) is the even-cleaner test we're running next.
- *"What's this universality thing?"* → near the tipping point the aggregate follows mean-field scaling; the realistic network gives a continuous transition (exponent γ≈1), the scrambled one an abrupt regime shift. Companion paper — happy to go deeper.
- *"Homophily coefficients are from 2004 data applied to later surveys?"* → they come from the same studies that validated the method; homophily is among the most stable social patterns, and we transport them without re-estimating, following McPherson & Smith.
