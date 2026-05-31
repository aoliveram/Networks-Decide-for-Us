# Speaker Script — CPIN 2026 (≈12 minutes)

*Everyday language, with brief technical anchors. One block per slide. Times are cumulative targets.*

---

## Slide 1 — Title (0:00–0:40)

> Hello everyone, thank you for being here. I'm Aníbal Olivera, a PhD candidate in Social Complexity Sciences at Universidad del Desarrollo in Chile. 
> The idea behind this provocative title is that the way our social ties are arranged doesn't just **help** the spread of behavior; it also looks, mathematically, like a phase transition."

---

## Slide 2 — The question (0:40–1:40)

> Here's the setup. Collective behaviors — adopting an innovation, joining a protest — spread person to person, like a contagion on a network.
>
> The problem is that most models of this use *synthetic* networks, where people are simple relays that just pass things along. But real human networks have a very specific texture: we usually cluster with people like us. Same age, same education, same background. Real networks follow **multidimensional homophily**.
>
> So the question is: conditional on *fixed individual propensities*, can the **shape of more realistic networks alone** change the outcome? Can the structure itself be an autonomous force, deciding not just *how much* spreads, but *how* it spreads?"

*(Technical anchor, optional: "We hold individual propensities constant and vary only the topology.")*

---

## Slide 3 — Model, part 1: Rational Choice (1:40–2:40)

> The model presents two ways a person can adopt. 
> 1- The first is simple and selfish: if the thing is attractive enough to *me personally*, I adopt it on my own. 
> Following the literature, we call the attractiveness of the innovation *Intrinsic Utility Level* 'IUL'.
> And each person has a *Minimum Utility Requirement* 'MUR', representing their aversion to adopting the behaviour. 
>
> The key thing about this first channel: it's **deterministic** and it ignores the topology.

*(Point at the node diagram: black = adopted, white = not.)*

---

## Slide 4 — Model, part 2: Selective Social Influence (2:40–4:00)

> The second way is by social influence: several of my contacts already adopted, so the cost of joining is lower and I will adopt even if I didn't by the rational channel.
>
> But one ingredient that makes our model different is that **I don't count everybody**. When an innovation or a collective action is spreading, there's a limit on how *socially different* a neighbor can be before I stop letting them influence me.
> Following the literature, if the demographics are very different, they are socially too far.
> We call the limit of influence **Maximum Social Proximity** 'MSP'.
>
> So, when MSP is low, I only listen to people just like me — echo chambers. When MSP is high, I'm open to influence from very different people.
>
> This channel **is** stochastic, and it **does** depend on the topology.

---

## Slide 5 — Imputing a real network (4:00–5:00)

> Now, where does the network come from? We borrow the *homophily strengths* — how strongly people sort by race, religion, age, education, sex — from McPherson & Smith 2019 and use them to construct synthetic-but-realistic networks that reproduce those homophilic tendencies using ERGMs, and we attach real survey answers that give an *aversion-to-innovation score* and a *collective-action propensity score* for each node.
>
> The resulting networks have the lumpy, right-skewed structure of real social networks — that's something the random models simply don't reproduce.

---

## Slide 6 — Simulation design (5:00–5:45)

> Then we just run it — a lot. We sweep the attractiveness of the innovation, the openness parameter MSP, the thresholds, different ways of seeding the first adopters. Around four million simulations.
>
> And crucially, we compare each realistic network against a 'scrambled' version (that has the same number of connections per person, but the connections are random and so the homophily is destroyed).

*(Technical anchor: "The scramble is a degree-preserving randomization — GSS-DP.")*

---

## Slide 7 — The physics mapping (5:45–7:00)

> This model is, almost line for line, a *ferromagnetic Ising model*.
>
> A spin flipping up or down ↔ a person adopting or not. 
> Magnetization (the order parameter) ↔ overall adoption. 
> The external field ↔ the innovation's intrinsic attractiveness (external to the individual). 
> And — this is *important* — temperature ↔ our openness parameter MSP.
>
> Why is MSP the temperature? Because temperature is what melts order. When the social openness is low, the network is frozen into rigid little clusters; when high, you're 'melting' those social boundaries so influence can flow across the whole system.
>
> So now you can ask: is there a *critical point* on adoption? Does the susceptibility diverge?"

---

## Slide 8 — Two paths through the critical point (7:00–8:00)

> Just as in the ferromagnetic situation, you can test criticallity following two paths.
>
> One path — the green one — you fix the external field and change the temperature. That gives you a *continuous 2nd order* transition. 
> The other — the blue one — you fix the temperature below the critical point and change the field. That one can give you a *1st order* transition.

---

## Slide 9 — Result 1: Criticality signatures (dynamics) (8:00–8:50)

> First, the dynamics. These three plots map adoption over the attractiveness–openness plane, and you can see how *rational adoption* grows as the attractiveness grows, but *social adoption* is more complex, so we have similar total adoption over a large region just because social adoption is relevant even when the attractiveness is low.
> Look at the bottom left: there's a sharp ridge where the *variance* of outcomes explodes (on the right you only have the variance of the number of steps of the simulation) — these are the fingerprints of a critical point.

---

## Slide 10 — Result 2: Structural premium (8:50–9:35)

> Second, the volume. If you scramble the network — destroying homophily but keeping everything else — adoption drops by about **57 percent**. Same people, same preferences, just a different wiring.
>
> Why? Because homophilic clusters are echo chambers, and echo chambers are actually *useful* here: they let social reinforcement build up safely in a corner until it's strong enough to ignite — even for innovations almost nobody finds attractive. So homophily is a *catalyst* for spread, not an obstacle.

---

## Slide 11 — Result 3: Universality, the discriminator (9:35–10:45)

> Third, the critical exponents. Near the critical point, susceptibility follows a power law, and the exponent *gamma* is essentially **1** on the realistic networks. That's the classical *mean-field value*.
>
> In the scrambled network, gamma is about **4** — which isn't a real critical exponent at all. The power law breaks. The system stops being predictable: it either does nothing or explodes, with no warning.
>
> So the punchline: this clean, universal behavior is **not generic** — it shows up *only* on the realistic, homophily-built network. Scramble the homophily and you lose it. The structure is what makes social change legible.

---

## Slide 12 — Why only the MSP axis discriminates (10:45–11:15)

> One subtlety we recently found, and I think it's elegant. Remember the two paths. The temperature path — MSP, the social channel — is the one that tells realistic and scrambled networks apart: gamma 1 versus 4. The field path — changing attractiveness, the individual channel — gives the *same* exponent for both, around 3, because that channel doesn't care about the network at all; it only reflects people's private preferences. So the structure-sensitivity lives exactly where it should: in the social axis.

*(Keep this light; flag as recent/exploratory.)*

---

## Slide 13 — The mean-field paradox (11:15–11:40)

> Quick honest puzzle. Mean-field behavior usually shows up in *well-mixed* systems — and a clustered, homophilic network is the opposite of well-mixed. 
> Our resolution: the clusters behave like little subsystems weakly linked by similarity bridges, and right at the critical point those bridges switch on together and synchronize the whole thing — which *looks* mean-field. Homophily gives the system a kind of 'social viscosity' that turns chaos into a smooth cascade.

---

## Slide 14 — Conclusions (11:40–12:10)

> So, to wrap up. Networks decide for us in two ways. 
> They set the *volume* of collective behavior — the 57 percent premium. 
> And they set its *character* — homophily makes tipping points continuous and predictable; destroy it and they become abrupt and volatile.
>
> The big message for this community: collective social change is not the sum of individual preferences. At the tipping point it behaves like critical matter — and the control parameter is the structure of our ties. Thank you."

---

## Slide 15 — References / Thanks
*(Stay on Thanks; take questions.)*

---

### Timing cheat-sheet
| Block | Target |
|---|---|
| Intro + question | 0:00–1:40 |
| Model (2 slides) | 1:40–4:00 |
| Data + design | 4:00–5:45 |
| Physics mapping + paths | 5:45–8:00 |
| Results (3 slides) | 8:00–10:45 |
| Subtlety + paradox | 10:45–11:40 |
| Conclusions | 11:40–12:10 |

### Anticipated questions
- *"Isn't mean-field trivial / expected for N=1000?"* → finite-size caveat; the **contrast** with the scrambled network (γ≈4) is the result, not the absolute value.
- *"Is δ≈3 robust?"* → exploratory; the MSP grid is coarse, β not yet pinned down, Widom relation not yet closed. Honest about it.
- *"Homophily from 2004 data applied to 2014?"* → coefficients are from the same surveys the method was validated on; homophily is among the most stable social patterns; 2004→2014 drift is smaller than the 1985→2004 drift those papers already absorbed.
- *"Why is adoption high at HIGH MSP if that's the hot phase?"* → in this mapping the high-adoption phase is the melted/open one; the sign convention is opposite to a cold ferromagnet but the criticality is identical.
