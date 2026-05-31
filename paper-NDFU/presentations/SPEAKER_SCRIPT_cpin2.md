# Speaker Script — CPIN 2026 (≈12 minutes)

*Everyday language, with brief technical anchors. One block per slide. Times are cumulative targets.*

---

## Slide 1 — Title (0:00–0:40)

> Hello everyone, thank you for being here. I'm Aníbal Olivera, a PhD candidate in Social Complexity Sciences at Universidad del Desarrollo in Chile. 
> The idea behind this provocative title is that the way our social ties are arranged doesn't just **help** the spread of behavior; it also looks, on the mathematical point of view, like a phase transition."

---

## Slide 2 — The question (0:40–1:40)

> Here's the setup. Collective behaviors — adopting an innovation, joining a protest — spread person to person, like a contagion on a network.
>
> The problem is that most models of this use *synthetic* networks, where people are simple relays that just pass things along. But real human networks have a very specific texture: we usually cluster with people like us. Same age, same education, same background. Real networks follow **multidimensional homophily**.
>
> So the question is: conditional on *fixed individual characteristics*, can the **shape of more realistic networks alone** change the outcome? Can the structure itself be an autonomous force, deciding not just *how much* spreads, but *how* it spreads?"

*(Technical anchor, optional: "We hold individual propensities constant and vary only the topology.")*

---

## Slide 3 — Model, part 1: Rational Choice (1:40–2:40)

> The model presents two ways a person can adopt.
> 1- The first is simple and selfish: if the thing is attractive enough to *me personally*, I adopt it on my own.
> Following the literature, we call the attractiveness of the collective behavior *Intrinsic Utility Level* 'IUL', and has a value between 0 and 1.
> And each person has a *Minimum Utility Requirement* 'MUR', representing their aversion to adopting the behaviour. 
>
> So, an agent will adopt if the collective behavior is attractive enough for his requirement.

*(Point at the node diagram: black = adopted, white = not.)*

---

## Slide 4 — Model, part 2: Selective Social Influence (2:40–4:00)

> The second way is by social influence: several of my contacts already adopted, so the cost of joining is lower over the time, and I will eventually adopt, even if I didn't by the rational channel.
>
> But a key ingredient in the model is that **I don't count everybody**, but there's a limit on how *socially different* a neighbor can be before I stop letting them influence me.
> Following the literature, if the demographics are very different, they are socially too far. We call the limit of influence **Maximum Social Proximity** 'MSP'.
>
> So, when *MSP* is low, I only listen to people just like me — echo chambers. When *MSP* is high, I'm open to influence from very different people.

---

## Slide 5 — Imputing a real network (4:00–5:00)

> Now, where do the networks come from? 
> We borrow the *homophily strengths* — that is how strongly people sort by demographics like age, education, or sex — from McPherson & Smith 2019 and use them to construct realistics, or *plausible*, networks that reproduce those homophilic tendencies using a Logit, and we finally attach real survey answers that give us an *aversion-to-innovation score* and a *collective-action propensity score* for each node.
>
> The resulting networks have clustering and degree distribution of real social networks — that's something the synthetic models simply don't reproduce.

---

## Slide 6 — Simulation design (5:00–5:45)

> Then we just run it — a lot. We sweep the attractiveness of the diffusion, the openness parameter MSP, the thresholds, different ways of seeding the first adopters... In total, we have around four million simulations.
>
> And crucially, we compare each realistic network against a 'scrambled' version (where the number of connections of each person is the same, but the connections are random, and so the homophily is destroyed).

*(Technical anchor: "The scramble is a degree-preserving randomization — GSS-DP.")*

---

## Slide 7 — The physics mapping (5:45–7:00)

> This model is, almost line for line, a *ferromagnetic Ising model*.
>
> A spin flipping up or down ↔ a person adopting or not. 
> Magnetization (the order parameter) ↔ mean adoption. 
> The external field ↔ the collective action's intrinsic attractiveness (external to the society). 
> And temperature ↔ our openness parameter (internal to the society).
>
> Why is MSP the temperature? Because temperature is what melts order. When the social openness is low, the network is frozen into rigid little echo chambers, you're only influenced by very similar people; when high, you're 'melting' those social boundaries so influence can flow across the whole system.
>
> So now you can ask: is there a *critical point* on adoption? 

---

## Slide 8 — Two paths through the critical point (7:00–8:00)

> Just as in the ferromagnetic situation, you can test criticallity following two paths.
>
> One path — the green one — you fix the external field and change the temperature. That gives you a *continuous 2nd order* transition. 
> The other — the blue one — you fix the temperature below the critical point and change the field. That one can give you a *1st order* transition.

---

## Slide 9 — Result 1: Adoption and criticality  (dynamics) (8:00–8:50)

> First, the dynamics. These three plots map adoption over the attractiveness–openness plane, and you can see how *rational adoption* grows as the attractiveness grows, but *social adoption* is more complex, so can we have similar total adoption over a large region just because social adoption is relevant **even when the attractiveness is low**.
> 
> Look at the bottom left: there's a sharp ridge where the *variance* explodes... and on the right you have regions where the number of final steps varies a lot.

---

## Slide 10 — Result 2: Structural premium (8:50–9:35)

> Second, the volume. If you scramble the network — destroying homophily but keeping everything else — adoption drops by about **57 percent**. Same people, same preferences, just a different wiring.
>
> Why? Because homophilic clusters are echo chambers, and echo chambers are actually *useful* here: they let social reinforcement build up safely in a corner until it's strong enough to ignite — even for non-attractive behaviors. So realistic networks help adoption.

---

## Slide 11 — Result 3 (part 1): Field direction, δ ≈ 3 for both (9:35–10:15)

> Now the critical exponents. There are two directions we can probe. First, the *field* direction: we sit right at the critical openness and turn up the attractiveness. Here both networks respond the same way, with an exponent delta of about 3. 
> So along this axis the two networks look identical; the structure makes no difference.

---

## Slide 12 — Result 3 (part 2): Temperature direction, γ discriminates (10:15–10:55)

> Now the *temperature* direction: we fix the attractiveness and turn up the social openness. 
> On the realistic networks the susceptibility diverges as a clean power law with exponent gamma essentially equal to **1** — which is the classic mean-field value.
>  On the scrambled network gamma is about **4**, which isn't a real critical exponent at all: the power law breaks, and the system either does nothing or explodes with no warning.
>
> So the punchline: this clean, universal behavior is **not generic** — it appears *only* on the realistic, homophily-built network. Scramble the homophily and you lose it. The structure is what makes social change legible.
>
> And why does that matter? In physics, *universality* means the microscopic details don't matter near the critical point: boiling water and a magnet losing its magnetism share the *same* exponents — the substance, the lattice, the chemistry are all irrelevant. The same logic here: near the tipping point the aggregate doesn't care about each person's exact threshold, the specific content of the behavior, or individual quirks — only the structure and how far we sit from the critical point. That's the precise sense in which structure makes the collective *legible*: a handful of macro-structural parameters govern it, and the microscopic mess washes out.

---

## Slide 13 — Conclusions (10:55–11:40)

> So, to conclude — let me answer the two questions I opened with.
> Can the shape of the network *alone* change the outcome? 
> Yes: empirical homophily gives a 57-percent premium in adoption general.
> Also, when we have an attractive behaviour and an unattractive, we ended up with almost the same adoption — just because of social influence.
> 
> Is the structure an autonomous force deciding *how* change spreads?
> Yes: on a realistic network the transition is continuous, and across comparable societies it carries statistical early warning — rising volatility, slower recovery — before it tips. Scramble the structure and that warning disappears; change becomes abrupt. This signposted regime is exclusive to homophilic structure. [REVISIT]
>
> This reframes a familiar intuition: even when a single episode feels abrupt and surprising, the *aggregate* — averaged over comparable cases — follows identifiable, even universal patterns. It's the same reason the S-curve, and the roughly one-quarter tipping fraction, keep reappearing across fax machines, hashtags, and protest waves.
> 
> The big message for this community: collective social change is not the sum of individual preferences. At the tipping point it behaves like critical matter — and the control parameter is the structure of our ties. Thank you.

---

## Slide 14 — References / Thanks
*(Stay on Thanks; take questions. The QR links to the GitHub repo.)*

---

## Appendix slide — The mean-field paradox (only if asked)

> A puzzle you might raise: mean-field behavior usually shows up in *well-mixed* systems, and a clustered homophilic network is the opposite of well-mixed. The resolution is that networks have no spatial dimension — your neighbors-of-neighbors grow exponentially, so the system effectively sits *above* the dimension where mean-field becomes exact. So mean-field is actually the *expected* outcome on a network. What really differs between our two networks isn't the dimension, it's the *order* of the transition: the clustered one is continuous, the scrambled one is discontinuous. Homophily gives the system a kind of 'social viscosity' that turns an abrupt shatter into a gradual, signposted cascade.

---

### Timing cheat-sheet (per slide, after universality additions)
| # | Slide | Duration | Cumulative |
|---|---|---|---|
| 1 | Title | 0:35 | 0:35 |
| 2 | The question | 0:55 | 1:30 |
| 3 | Rational Choice | 0:55 | 2:25 |
| 4 | Selective Social Influence | 1:15 | 3:40 |
| 5 | Imputation | 1:00 | 4:40 |
| 6 | Simulation design | 0:40 | 5:20 |
| 7 | Physics mapping | 1:10 | 6:30 |
| 8 | Two paths | 0:55 | 7:25 |
| 9 | Result 1 — dynamics | 0:50 | 8:15 |
| 10 | Result 2 — premium | 0:45 | 9:00 |
| 11 | Result 3-i — δ (both) | 0:35 | 9:35 |
| 12 | Result 3-ii — γ **(+ universality / micro-irrelevance)** | 1:20 | 10:55 |
| 13 | Conclusions **(+ real phenomena)** | 1:05 | 12:00 |
| 14 | Thanks | 0:10 | 12:10 |

**Total ≈ 12:10** — slightly over. The additions made slide 12 the heaviest (~1:20) and slide 13 ~1:05. **Release valve:** compress slide 11 (δ, exploratory) from 0:35 to ~0:15 — "along the field axis both networks look identical, δ≈3; the interesting axis is temperature" — which brings you back to ~11:50. Checkpoint: be at "Two paths" (slide 8) by ~7:25.

*(The mean-field paradox is an appendix slide — only if asked.)*

### Anticipated questions
- *"Isn't mean-field trivial / expected for N=1000?"* → finite-size caveat; the **contrast** with the scrambled network (γ≈4) is the result, not the absolute value.
- *"Is δ≈3 robust?"* → exploratory; the MSP grid is coarse, β not yet pinned down, Widom relation not yet closed. Honest about it.
- *"Homophily from 2004 data applied to 2014?"* → coefficients are from the same surveys the method was validated on; homophily is among the most stable social patterns; 2004→2014 drift is smaller than the 1985→2004 drift those papers already absorbed.
- *"Why is adoption high at HIGH MSP if that's the hot phase?"* → in this mapping the high-adoption phase is the melted/open one; the sign convention is opposite to a cold ferromagnet but the criticality is identical.
