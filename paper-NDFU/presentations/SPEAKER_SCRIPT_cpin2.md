# Speaker Script — CPIN 2026 (≈12 minutes)

*Everyday language, with brief technical anchors. One block per slide. Times are cumulative targets.*

---

## Slide 1 — Title (0:00–0:40)

> "Hello everyone, thank you for being here. I'm Aníbal Olivera, a PhD candidate in Social Complexity Sciences at CICS, Universidad del Desarrollo in Chile. This is joint work with Jorge Fábrega and Thomas Valente.
>
> The title is a little provocative — *Networks Decide for Us* — and the claim I want to make today is that the way our social ties are arranged doesn't just *help or slow down* the spread of behavior; it actually controls whether social change is smooth and predictable, or sudden and explosive. And I'll argue that this looks, mathematically, like a phase transition — which is why I'm bringing it to this audience."

---

## Slide 2 — The question (0:40–1:40)

> "Here's the setup. Collective behaviors — adopting an innovation, joining a protest — spread person to person, like a contagion on a network.
>
> The problem is that most models of this use *made-up* networks, and treat people as simple relays that just pass things along. But real human networks have a very specific texture: we cluster with people like us. Same age, same education, same background. That's called **homophily**.
>
> So the question is: if we *freeze* everyone's personal willingness to adopt — their preferences are fixed — can the *shape of the network alone* change the outcome? Can the structure itself be an active force, deciding not just *how much* spreads, but *how* it spreads?"

*(Technical anchor, optional: "We hold individual propensities constant and vary only the topology.")*

---

## Slide 3 — Model, part 1: Rational Choice (1:40–2:40)

> "The model has two ways a person can adopt. The first is simple and selfish: if the thing is attractive enough to *me personally*, I adopt it on my own — I don't need anyone. We call the attractiveness of the innovation 'IUL', and each person has a personal bar it has to clear.
>
> The key thing about this first channel: it's **deterministic** and it completely ignores the network. It only depends on the person."

*(Point at the node diagram: black = adopted, white = not.)*

---

## Slide 4 — Model, part 2: Selective Social Influence (2:40–4:00)

> "The second way is social: I adopt because enough of my contacts already did. This is the classic threshold idea — peer pressure.
>
> But here's the one ingredient that makes our model different, and it's simple: **I don't count everybody**. There's a limit on how *socially different* a neighbor can be before I stop letting them influence me. If someone is too far from me — very different demographics — their adoption just doesn't register for me. We call that limit **MSP**, the maximum social proximity.
>
> So the whole model really only asks one extra thing: *below some threshold of similarity, I simply don't consider that neighbor's influence.* When MSP is low, I only listen to people just like me — echo chambers. When MSP is high, I'm open to influence from very different people.
>
> And unlike the first channel, this social channel **is** stochastic, and it **does** depend on the network."

---

## Slide 5 — Imputing a real network (4:00–5:00)

> "Now, where does the network come from? We don't invent it. We borrow the *homophily strengths* — how strongly people sort by race, religion, age, education, sex — from large US surveys where those have been carefully measured. Then we grow synthetic-but-realistic networks that reproduce those tendencies, and we attach real survey respondents to the nodes.
>
> The payoff is on the right: these networks have the lumpy, right-skewed structure of real social networks — something the textbook random models simply don't reproduce."

---

## Slide 6 — Simulation design (5:00–5:45)

> "Then we just run it — a lot. We sweep the attractiveness of the innovation, the openness parameter MSP, the thresholds, different ways of seeding the first adopters. Around four million simulations.
>
> And critically, we compare each realistic network against a 'scrambled' version — same number of connections per person, but homophily destroyed. That scrambled baseline is our control."

*(Technical anchor: "The scramble is a degree-preserving randomization — GSS-DP.")*

---

## Slide 7 — The physics mapping (5:45–7:00)

> "Here's where it connects to this room. This model is, almost line for line, a **ferromagnet** — the Ising model.
>
> A spin flipping up or down ↔ a person adopting or not. Overall magnetization ↔ overall adoption. The external magnetic field ↔ the innovation's intrinsic attractiveness. And — this is the important one — **temperature ↔ our openness parameter MSP**.
>
> Why is MSP the temperature? Because temperature is what melts order. When MSP is low, the network is frozen into rigid little clusters — like a cold magnet locked into domains. Turn MSP up, and you 'melt' those social boundaries so influence can flow across the whole system. That's exactly thermal agitation.
>
> And once you have that mapping, you get to ask the physicist's question: is there a *critical point*? Does the susceptibility diverge?"

---

## Slide 8 — Two paths through the critical point (7:00–8:00)

> "In a magnet there are two ways to approach the critical point, and they behave very differently.
>
> One path — the green one — you hold the field at zero and change the temperature. That gives you a *smooth, continuous* transition. The other — the blue one — you change the field. That one can give you a *sudden jump*, a discontinuity.
>
> In our language: changing MSP is the temperature path, and it runs through the *social* channel. Changing the innovation's attractiveness is the field path, and it runs through the *individual* channel. Keep this distinction in mind — it pays off at the end."

---

## Slide 9 — Result 1: Structural premium (8:00–8:45)

> "First result, the simple one. If you scramble the network — destroy homophily but keep everything else — adoption drops by about **57 percent**. Same people, same preferences, just a different wiring.
>
> Why? Because homophilic clusters are echo chambers, and echo chambers are actually *useful* here: they let social reinforcement build up safely in a corner until it's strong enough to ignite — even for innovations almost nobody finds attractive. So homophily is a *catalyst* for spread, not an obstacle."

---

## Slide 10 — Result 2: Criticality signatures (8:45–9:30)

> "Second result. These heatmaps map adoption over the attractiveness–openness plane. Look at the bottom row. On the left, there's a sharp ridge where the *variance* of outcomes explodes — that's susceptibility, and a spike in susceptibility is the fingerprint of a critical point. On the right, the time the system takes to settle also blows up — physicists call that **critical slowing down**. Both are textbook signatures of a continuous phase transition."

---

## Slide 11 — Result 3: Universality, the discriminator (9:30–10:45)

> "Third result, and this is the heart of it. Near the critical point, susceptibility follows a power law, and the exponent — gamma — tells you the *universality class*.
>
> In the realistic homophilic network, gamma is essentially **1**. That's the mean-field value — the cleanest, most classical kind of continuous transition. Social change there is an orderly, rolling avalanche.
>
> In the scrambled network, gamma is about **4** — which isn't a real critical exponent at all. The power law breaks. The system stops being predictable: it either does nothing or explodes, with no warning.
>
> So the punchline: **the smooth, predictable, universal behavior is a property of the *realistic* network. Scramble the homophily and you lose it.** The structure is what makes social change legible."

---

## Slide 12 — Why only the MSP axis discriminates (10:45–11:15)

> "One subtlety we recently found, and I think it's elegant. Remember the two paths. The temperature path — MSP, the social channel — is the one that tells realistic and scrambled networks apart: gamma 1 versus 4. The field path — changing attractiveness, the individual channel — gives the *same* exponent for both, around 3, because that channel doesn't care about the network at all; it only reflects people's private preferences. So the structure-sensitivity lives exactly where it should: in the social axis."

*(Keep this light; flag as recent/exploratory.)*

---

## Slide 13 — The mean-field paradox (11:15–11:40)

> "Quick honest puzzle. Mean-field behavior usually shows up in *well-mixed* systems — and a clustered, homophilic network is the opposite of well-mixed. Our resolution: the clusters behave like little subsystems weakly linked by similarity bridges, and right at the critical point those bridges switch on together and synchronize the whole thing — which *looks* mean-field. Homophily gives the system a kind of 'social viscosity' that turns chaos into a smooth cascade."

---

## Slide 14 — Conclusions (11:40–12:10)

> "So, to wrap up. Networks decide for us in two ways. They set the *volume* of collective behavior — the 57 percent premium. And they set its *character* — homophily makes tipping points continuous and predictable; destroy it and they become abrupt and volatile.
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
