# Speaker Script — Sunbelt 2026 (≈15 minutes + 5 Q&A)

*Session: "Contagion and Diffusion processes through Social Networks." Classic-SNA audience. Everyday language, like explaining the study to a friend. Lead with the pipeline + the structural premium; the physics is just one teaser slide. Times are cumulative.*

---

## Slide 1 — Title (0:00–0:30)

> Hi everyone, thanks for being here. I'm Aníbal Olivera, a PhD candidate at the Center for Social Complexity in Chile. This is joint work with Jorge Fábrega and Tom Valente.
> The talk is about a simple but stubborn question: how much of where a behavior spreads is about *who people are*, and how much is about *how they're wired together*?

---

## Slide 2 — The question (0:30–1:30)

> Here's the setup. Behaviors — adopting a new product, joining a protest — spread person to person, like a contagion. And real human networks have a very particular texture: we tend to cluster with people like us. Same age, same education, same background. That's homophily.
>
> But here's the frustrating part. Most diffusion models run on *made-up* networks. And most of our great surveys measure *people* — their attitudes, their propensities — but they don't tell us who's connected to whom.
>
> So the question I want to answer is: if we *freeze* everyone's individual propensity, does the *shape of the network alone* — not who the people are — decide how much spreads?

---

## Slide 3 — The gap, and where we go (1:30–2:30)

> Let me put three standing problems on the table. One: diffusion models usually use synthetic networks — small-world, scale-free — not real ones. Two: the surveys that *do* measure propensities, like the ATP or the GSS, have no network data at all. And three: a lot of threshold models quietly assume that *every* adoption needs peer influence.
>
> Our paper makes two moves against this. First, we *impute* a plausible whole-network onto real survey respondents — we build them a realistic social network out of empirical homophily. Second, we measure how much that structure actually causes diffusion, and we dig into *what about the structure* matters.
>
> We're building on McPherson and Smith's idea of imputing social context from surveys — but we push it from a single latent "context" number all the way to full, simulable networks.

---

## Slide 4 — Adoption, part 1: Rational Choice (2:30–3:30)

> The model has two ways a person adopts. The first is simple and selfish: if the thing is attractive enough to me personally, I adopt it on my own. We call the attractiveness of the behavior its "Intrinsic Utility Level" — think of it as a quality score; the highest-quality option is the one that would win the most head-to-head comparisons. And each person has a personal bar it has to clear — their minimum requirement.
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

## Slide 7 — Simulation design (6:45–7:30)

> Then we just run it — a lot. We sweep the attractiveness of the behavior, the openness parameter MSP, the thresholds, and different ways of seeding the first adopters. About four million simulations.
>
> And here's the crucial comparison. For every realistic network, we build a "scrambled" twin: same number of connections per person, but we rewire *whom* they connect to at random. That destroys the homophily but keeps the degree. That scrambled version — we call it degree-preserving randomization — is our counterfactual for "the same people, a different structure."

---

## Slide 8 — Result 1: The structural premium (7:30–9:00)  ← headline finding

> First result, and it's the headline. We fit a big additive model to isolate the effect of *just* the structure. When we scramble the network — destroy the homophily but keep everything else — adoption odds drop by about **57 percent**.
>
> Let that sink in: same people, same preferences, same thresholds. The *only* thing that changed is the wiring. And it cuts adoption by more than half.
>
> Why? Because homophilic clusters are echo chambers, and echo chambers turn out to be *useful* here — they let social reinforcement build up safely in a corner until it's strong enough to ignite. So homophily isn't a barrier to spread. It's a catalyst.

---

## Slide 9 — Result 2a: What does breaking structure destroy? (9:00–10:30)  ← the new diagnosis

> Okay — so structure matters. But *what about* the structure? When we scramble the network, we break a lot of things at once. So we measured, across a hundred real networks and a hundred scrambled twins, exactly what changes.
>
> This plot is the answer, and it's dramatic. On the left: homophily by each demographic axis — age, education, race, religion, sex. In the real networks it's strongly positive: similar people are connected. In the scrambled ones, every single one collapses to *zero*. Random mixing.
>
> And on the right is the cleanest number in the paper: the average social distance *along the ties* — the exact quantity my influence rule uses — jumps 42 percent. In a real network your neighbors are *similar* to you; scramble it, and your neighbors become *random strangers*. That's the mechanism, in one number. All of this is wildly significant — p less than ten to the minus thirty.

---

## Slide 10 — Result 2b: Topology degrades too (10:30–11:15)

> Now, to be honest with you: scrambling also damages the *topology* itself. Clustering drops about a third, triangles drop a third, the community structure — modularity — drops by half. So those are suspects too.
>
> But notice two things. The controls — degree, density, the giant component — are *identical*, so this isn't a size artifact. And the most complete collapse, the one that touches my model directly, is the homophily. So homophily is the prime suspect. Cleanly separating it from clustering is what we're doing next — I'll come back to that.

---

## Slide 11 — Where this sits (11:15–12:15)

> Let me place this against a paper many of you know — Macy and Evtushenko, in Sociological Science. They reintroduce individual idiosyncrasy into threshold models as *stochastic noise*, and they show something lovely: a little noise can *spontaneously instigate* a cascade that a weak-interest group would otherwise never get off the ground.
>
> We go one level *behind* that. Our idiosyncrasy isn't noise — it's *measured preference*. At low intrinsic utility, the people who adopt on their own — the rational adopters — *are* the instigators. They're not random sparks; they're real, surveyed preferences. And on top of that, we show it's the *structure*, not just the spread of thresholds, that governs the outcome.
>
> The one-line version: the same agents, on a different topology, produce a different society.

---

## Slide 12 — Teaser: the shape of change (12:15–13:15)

> I'll leave you with one teaser, because I think it's the fun part. Beyond *how much* spreads, the structure also governs *how* the change arrives.
>
> On a realistic network, the tipping is smooth and predictable — a continuous transition. On a scrambled network, it's abrupt — all-or-nothing, with no warning. There's actually a deep signature here that physicists call universality — the aggregate behaves like a phase transition in matter. But that's a whole other talk, and I'd love to get into it in the questions.

---

## Slide 13 — Conclusions (13:15–14:15)

> So, to wrap up — three takeaways.
> A method: you can impute realistic, *simulable* whole-networks onto survey respondents from empirical homophily. That's reusable anywhere you've got propensities but no ties.
> A finding: holding individual preferences fixed, real structure gives you a 57-percent adoption premium.
> And a mechanism: what "breaking structure" destroys, above all, is homophily — your neighbors stop being similar to you.
>
> The big message: where a behavior ends up isn't the sum of individual preferences. It's governed by the structure of our ties. Networks decide for us.

---

## Slide 14 — Limitations & next steps (14:15–14:45)

> Quickly, where we're headed. First, the cleanest test to separate homophily from clustering: keep the network *byte-for-byte identical* and just shuffle the attributes across nodes — that isolates pure homophily. Second, let openness vary by group instead of being one global number. Third, push the propensities to higher-stakes behavior — risky collective action. And fourth, the criticality story from the teaser, as a companion paper.

---

## Slide 15 — Thanks / References (14:45–15:00)
*(Stay on Thanks; take questions. QR links to the preprint.)*

> Thank you — the QR there links to the draft, and I'm very happy to take questions.

---

### Timing cheat-sheet
| Block | Target |
|---|---|
| Title + question + gap | 0:00–2:30 |
| Model (2 slides) | 2:30–4:45 |
| Pipeline + design | 4:45–7:30 |
| Result 1 (premium) | 7:30–9:00 |
| Result 2 (diagnosis, 2 slides) | 9:00–11:15 |
| Macy connection | 11:15–12:15 |
| Teaser (physics) | 12:15–13:15 |
| Conclusions + next | 13:15–14:45 |
| Thanks / Q&A | 14:45–15:00 |

**Release valve if running long:** compress Slide 10 (topology degrades) to two sentences — "scrambling also hurts clustering and modularity, but the controls are identical and homophily is the most complete collapse" — and skip straight to Macy.

### Anticipated questions
- *"Are the imputed networks validated?"* → they reproduce the right-skewed degree distribution and realistic clustering of real networks (the degree-dist figure); held-out tie-feature validation is the natural next robustness check.
- *"Isn't MSP a macro knob you're just turning?"* → yes, in this version it's a global control parameter; we read it as a population's *social porosity*. Making it heterogeneous and individual-level is the next step (and the honest answer to the bottom-up critique).
- *"Why degree-preserving randomization and not something gentler?"* → it's the cleanest counterfactual that holds degree fixed; the attribute-shuffle null (identical topology, permuted attributes) is the even-cleaner test we're running next.
- *"What's this universality thing?"* → near the tipping point the aggregate follows mean-field scaling; the realistic network gives a continuous transition (exponent γ≈1), the scrambled one an abrupt regime shift. Companion paper — happy to go deeper.
- *"Homophily coefficients are from 2004 data applied to later surveys?"* → they come from the same studies that validated the method; homophily is among the most stable social patterns, and we transport them without re-estimating, following McPherson & Smith.
