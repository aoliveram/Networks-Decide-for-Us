# Speaker Script — CPIN 2026 (≈12 minutes)

*Everyday language, with brief technical anchors. One block per slide. Times are cumulative targets.*

---

## Slide 1 — Title (0:00–0:40)

> Hello everyone, thank you for being here. I'm Aníbal Olivera, a PhD candidate in Social Complexity Sciences at Universidad del Desarrollo in Chile. 
> The idea behind this provocative title is that the way our social ties are arranged doesn't just **help** the spread of behavior; it also looks, on the mathematical point of view, like a phase transition."

---

## Slide 2 — The question (0:40–1:40)

> Here's the setup. Collective behaviors — adopting an innovation, joining a protest — spread as adoption cascades on social networks, as complex contagion.
>
> The problem is that most models use *synthetic* networks, where people are simple relays that just pass things along.
> But real human networks have a very specific texture: we usually cluster with people like us. Same age, same education, same background. Real networks follow **multidimensional homophily**.
> Also, individuals have different preferences to adopting behaviors, which can be measured in surveys.
>
> So we can ask: conditional on *fixed individual characteristics*, can the **shape of more realistic networks alone** change the outcome? Can the structure itself be an autonomous force, deciding not just *how much* spreads, but *how* it spreads?"

*(Technical anchor, optional: "We hold individual propensities constant and vary only the topology.")*

---

## Slide 3 — Model, part 1: Rational Choice (1:40–2:40)

> The model presents two ways a person can adopt. The first is simple and selfish: if the thing is attractive enough to *me*, I will adopt it.
> Following the literature, we say that each collective behavior has an *Intrinsic Utility Level* 'IUL', that is a number between 0 and 1.
> And also each person has a *Minimum Utility Requirement* 'MUR', that represent their aversion to adopting.
>
> So, an agent will adopt if the collective behavior is attractive enough for his *requirement*.

---

## Slide 4 — Model, part 2: Selective Social Influence (2:40–4:00)

> The second way is by social influence: when several of my contacts already adopted,  the cost of joining is lower over the time, and I will eventually adopt, even if it didn't reach my bare utility requirement.
>
> But a key ingredient of the model is that it is selective: **I don't count everybody** as effective influence. I'll listen only to those individuals closelly enogh to me in terms of *social distance* d_ij.
>
> Following the literature, there is a **Maximum Social Proximity** 'MSP' where the limit of influence sits.
> So, when *MSP* is low, I only listen to people just like me — *echo chambers*. When *MSP* is high, I'm open to be influenced from very different people.

---

## Slide 5 — Imputing a real network (4:00–5:00)

> Now, where do the networks come from? 
> We borrow the *homophily strengths* — that is how *strongly people sort* by demographics like age, education, or sex — from McPherson & Smith 2019 and use them to construct realistics networks using the data from survays to *feed a Logit that reproduce those homophilic tendencies*.
> That data have items that allows us to construct an *aversion-to-innovation score* or a *collective-action propensity score* for each node.
>
> Those networks have the right topology, just as the literature shows the real-social-networks topology are.

---

## Slide 6 — Simulation design (5:00–5:45)

> Then we just run it — a lot. We sweep the attractiveness of the diffusion, the openness parameter, the thresholds, different ways of seeding the first adopters... 
> In total, we have around four million simulations.
>
> And crucially, we compare each realistic network against a 'scrambled' version (where the number of connections of each person remains the same, but the connections are now random, and so the homophily is destroyed).

---

## Slide 7 — The physics mapping (5:45–7:00)

> This model is, almost line for line, a *ferromagnetic Ising model*.
>
> A spin up or down ↔ a person adopting or not. 
> Magnetization ↔ mean adoption. 
> The external field ↔ the collective behaviour's intrinsic attractiveness (external to the society). 
> And temperature ↔ our openness parameter (internal to the society).
>
> Why is MSP the temperature? Because temperature is what melts order. When the social openness is low, you're only influenced by very similar people, so we have several little echo chambers; when high, you're 'melting' those social boundaries so influence can flow across the whole system.
>
> So the questions are: does adoption have a *critical point* — and if so, what kind of transition is it?

---

## Slide 8 — Two paths through the critical point (7:00–8:00)

> Just as in the ferromagnetic situation, you can test criticallity following two paths.
>
> In one path — the blue one — you fix the temperature below a critical point and change the field (attractiveness in our case), and so you can compute the critical exponent δ.
> In the other — the green one — you fix the external field and change the temperature to compute the critical exponent γ. 

---

## Slide 9 — Result 1: Adoption and criticality  (dynamics) (8:00–8:50)

> First, the dynamics. These three plots on the top map adoption over the attractiveness–openness plane, and you can see how *rational adoption* grows as the attractiveness grows, but *social adoption* is more complex, so we can have similar total adoption over a large region just because social adoption is relevant **even when the attractiveness is low**.
> 
> Look at the bottom left: there's a sharp ridge where the *variance* explodes... and on the right you have regions where the number of final steps varies a lot.

---

## Slide 10 — Result 2: Structural premium (8:50–9:35)

> Second, the volume. If you scramble the network — destroying homophily but keeping everything else — adoption drops by about **57 percent**.
>
> Why? Because homophilic clusters are echo chambers that are *useful* for complex contegions: they have redundancy so they can grow safely in a corner until they are strong enough to ignite — even for non-attractive behaviors.

---

## Slide 11 — Result 3 (part 1): Field direction, δ ≈ 3 for both (9:35–10:15)

> Now the critical exponents. If we sweep atractiveness at the critical openness, we will have an exponent delta of about 3 *in both topologies*.
> This value is very close to a mean-field value.

---

## Slide 12 — Result 3 (part 2): Temperature direction, γ discriminates (10:15–10:55)

> On the other direction (fixing attractiveness and moving openness), we will see that on realistic networks the susceptibility diverges as a clean *power law* where the exponent gamma essentially **1** — which is the classic mean-field value.
> On the scrambled network gamma is about **4**. Not a real critical exponent.
>
> And why does this matter? In physics, *universality* means the microscopic details don't matter near the critical point: boiling water and a magnet losing its magnetism share the *same* exponents, regardless they different chemistry, structure, etc.
> The logic here is that the aggregate realizations of diffusion doesn't care about the specific content of the behavior — only the structure. 
> So, if we accept this model as a *stylized reality*, we can say that a couple of macro-structural parameters govern the phenomenon.

---

## Slide 13 — Conclusions (10:55–11:40)

> So, to draw some conclusions — let me answer the 2 questions I opened with.
> 
> Can the shape of the network *alone* change the outcome? 
> Yes: empirical homophily yields a structural premium (a 57% more adoption in general).
> Also, there are cases when highly attractive and unattractive behaviors reach nearly the same adoption — just because of the social influence contribution.
> 
> Is the structure an autonomous force deciding *how* change spreads?
> Yes: Even if we feel that there are particular cases when adoption are bimodal, we still *CAN use the usual Roger's S-curve* to explain social diffusions so different as fax machines, hashtags, and protest waves. --> Universality.
> 
> The big message for this community: Macroscopic social change is not the sum of individual preferences. It’s governed by the structure of our real topology.
> Thank you.

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
