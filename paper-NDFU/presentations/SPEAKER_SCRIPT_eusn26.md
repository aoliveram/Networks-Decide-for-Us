# Speaker Script — EUSN 2026 (15 min + 7 min Q&A)

---

## Slide 1 — Title (0:00–0:30)

> Hello everyone, thank you for being here. I'm Aníbal Olivera, a PhD candidate in Social Complexity Sciences at Universidad del Desarrollo, in Chile.
> And the idea behind this provocative title is that our real social ties do causal work in the diffusion of behaviors that cannot be explained by individual attributes alone.

---

## Slide 2 — The question (0:30–2:00)

> Collective behaviors — adopting an innovation, joining a protest — spread as adoption cascades on social networks.
>
> And we know two things. 
> ONE: We can measure *who people are* (because surveys give us their preferences), and we know *how people cluster* (because we usually connect with people similar to us, in a phenomenon called multidimensional homophily).
>
> Also, most studies usually are run on synthetic networks, and use invented preferences for collective behaviors.
>
> Our experiment we address this gap: we **fix the people** (with real preferences), **and we rewire the ties**.
>
> So we have the same individuals, same measured preferences — only the wiring changes. 
> And the question is: conditional on fixed individual preferences, to what extent does the empirical network structure — the one shaped by homophily — govern the diffusion of collective behaviors?

*(Slide 3 is a section divider — just say: "Let me show you the model first.")*

---

## Slide 4 — Adoption, part 1: Rational Choice (2:05–3:20)

> The model gives a person two ways to adopt. The first one is simple: if the thing is attractive enough *for me*, I adopt it.
>
> Following the literature, each behavior has an *Intrinsic Utility Level* — a number between zero and one.
> And each person has a *Minimum Utility Requirement* — their aversion to adopting.
>
> So an agent will adopt if the behavior is attractive enough for their requirement. As you can see in the plots, this channel works like a simple contagion.

---

## Slide 5 — Adoption, part 2: Selective Social Influence (3:20–5:00)

> The second way is social influence. 
> When several of my contacts have already adopted, the cost of joining goes down over time, so I could adopt (even when the behavior did not reach my original requirement).
>
> But there are some nice papers that show that, with the *same* exposure, a socially similar contact is more influential on collective behavior than a distant one.
>
> So we add that key ingredient making our social influence a *selective* social influence, where I only listen to those that are close enoght to me in terms of social distance.
>
>
> Following the literature, there is a **Maximum Social Proximity** 'MSP' where the limit of influence sits.
> So, when *MSP* is high, the society is more flexible; and when *MSP* is low, I only listen to people very similar to me — *echo chambers* —.

*(Slide 6 is a section divider — just say: "Now, where do the networks come from?")*

---

## Slide 7 — Imputing a network onto a survey (5:05–7:00) ← the methodological heart

> Now, where do the networks come from?
> We borrow the *homophily strengths* —that is how *strongly people sort* by demographics like age, education, or sex— from a nice paper from McPherson & Smith, and use them to construct realistic networks using the data from surveys —that don't have structural information— to *feed a Logit that reproduces those homophilic tendencies*.
> And the cool thing is that that data has items that allow us to construct either an *aversion-to-innovation score* or a *collective-action propensity score* for each node.
>
> Those networks have the right topology, just as the literature shows the real-social-networks topologies are,
> and the nodes have realistic propensities scores.
> 
> Something that synthetic models simply cannot reproduce.
>

---

## Slide 8 — The experiment: fix the people, rewire the ties (7:00–8:30)

> Now the experiment. 
> So first, we will call the realistic networks "GSS", just because we use the General Social Survey respondents' data, which have items to construct a "collective action" propensity score.
> And, for each realistic network, we build a *randomized twin* with the same people, same size — but the ties are rewired randomly.
>
> And that should improve the diffusion. Random rewiring creates **long-range bridges**. And classic theory — the strength of weak ties — says those bridges should *help* diffusion.
>
> Or, if we accept the naive point of view, rewiring the ties should do nothing...

---

## Slide 9 — Simulation design (8:30–9:15)

> OK. With those networks, we just run the model — a lot. We sweep the attractiveness of the behavior, the openness parameter, the threshold distributions, different ways of seeding the initial adopters. In total, around four million simulations.
>
> And crucially we fit a Big Additive Model to compare the results between the realistic networks and the rewired version.

*(Slide 10 is a section divider — just say: "The results.")*

---

## Slide 11 — Result 1: How adoption plays out (9:20–10:45)

> First, the dynamics. These three plots map total adoption (and its decomposition into adoption via rational and social influence channels) over the *attractiveness–openness plane*... and you can see how *rational adoption* grows as the attractiveness grows, but *social adoption* is more complex, so we can have similar total adoption over a large region just because social adoption is relevant **even when the attractiveness is low**.

> Now we will make a single cut around here...

---

## Slide 12 — Result 2: A non-structural blocking of complex contagion (10:45–12:00)

> To show you that there is a huge sensibility on the parameters of the behavior.
> 
> As you can see, the rational and social influence channels are complementary, so the total adoption is flat most of the time,
> but the social influence channel shows an abrupt interruption at a certain critical value... and the macro adoption just falls completely because the *Complex contagion term dies* suddenly.
> 
> And we already know that Centola and Macy showed *structural* barriers to complex diffusion; but This is a *not structural* barrier to diffusion.

---

## Slide 13 — Result 3: And now, what happens when we rewire? (12:00–14:00) ← the headline, slow down here

> Third result — what happens when we do the rewiring?
> Well, here we have the results from the Big Aditive Model, so we can see how, When we swap the realistic network for the randomized twin, the **odds of adoption fall by more than half**.
>
> But they are the same individuals, same preferences; And the rewiring created long-range bridges that *should have helped*.
>
> But the Social Influence fall was just too big.
>
> This shows that, even when we have both Simple and Complex contagion, For socially-gated adoption, the bridges are useless...
>
> Why does the realistic network win? Because homophilic *echo chambers* give social reinforcement a safe place to build up — locally, with redundancy — until it is strong enough to ignite.

---

## Slide 14 — Conclusions (14:00–15:00)

> So, to close.
>
> One *method*: you can build realistic, simulable whole-networks for survey respondents out of empirical homophily. That is reusable anywhere you have attitudes but no network data.
>
> Three *findings*. One: social influence explains adoption even for barely attractive behaviors. Two: complex contagion can be blocked in a way that is *not structural* — a small shift in the social parameters flips it from all to none. Three: rewiring the ties cuts the odds of adoption by fifty-seven percent — even though the new bridges "should" have helped.
>
> So, where a behavior ends up is not the sum of individual preferences. It is governed by the structure of our real ties. That is why **networks decide for us**. Thank you.

*(Slides 15–16: References and Thanks. Go to Thanks; take questions. The QR links to the preprint.)*

---

# Q&A — anticipated questions (backup slides 18–22)

**Q: "Is it the homophily, or is it the clustering, that produces the premium?"** → go to **Backup 18**
> Great question — we built a third network exactly for this. We keep the graph byte-identical — every triangle, every clique — and we only *shuffle who sits where*. So the neighbors stay, but they become social strangers.
> That alone reproduces about ninety percent of the full penalty. So *neighbor similarity* carries the premium.
> And why not a clean "topology-only" test? Because in these networks the clustering *emerges from* the homophily — if we re-simulate the homophily-only model, the triangles come back. They are a package: clustering is crystallized homophily.

**Q: "What do you mean by 'social temperature'?"** → go to **Backup 19**
> The model maps, almost line by line, to a ferromagnetic system. Adopting or not is a spin; the attractiveness is an external field; and the openness MSP plays the role of temperature — because temperature is what melts order. Low MSP: frozen echo chambers. High MSP: the boundaries melt and influence flows. That mapping is a whole companion analysis — happy to talk after.

**Q: "Is the abrupt flip the same in both networks? / how do the transitions compare?"** → go to **Backup 20** (and **21** if they want numbers)
> Good question. Both networks flip, but not in the same way. The realistic network tips *earlier*, and it passes through intermediate states — you can see a middle point around fifty percent. The randomized one stays flat and then jumps in one single step — hit or flop, the same pattern Tur and colleagues found theoretically in *Social Networks*.
> And if you want it quantified: near the tipping point the variance of outcomes grows as a power law on the realistic networks — exponent about one, a continuous transition. On the randomized ones there is no clean power law at all — the fingerprint of a discontinuous jump.
> One honest caveat: these are statements about the *distribution* of outcomes across comparable runs — not a forecast for one single history.

**Q: "Are the imputed networks validated?"**
> They reproduce the right-skewed degree distribution and realistic clustering of real personal networks, and they respect the demographic marginals and the estimated homophily. Checking them against held-out tie features is the natural next step.

**Q: "Isn't this circular? You assumed influence is gated by similarity."**
> The mechanism is not free — it has experimental support: Centola 2011 shows that similar contacts are more influential at the same exposure.
> And the tautology would only predict "something drops". It does not predict the size — half the adoption. It does not predict *where* — the premium is largest when society is rigid, and disappears when openness is high. And it does not predict the channel split — the rational channel barely moves. Those are findings, not assumptions.

**Q: "Your agents look abstract; MSP is a global knob."**
> Fair. Making MSP heterogeneous — per person or per group — is the next step. And we already have a fully bottom-up variant of the model, where each person's threshold is *derived* from their measured preference; it reproduces the same phase structure. Work in progress.

**Q: "What exactly is the randomization?"**
> The baseline keeps the same people, the same network size, and the same density, and rewires the ties at random. Degree-preserving refinements are part of ongoing robustness work.

---

### Timing cheat-sheet (target 15:00)

| # | Slide | Duration | Cumulative |
|---|---|---|---|
| 1 | Title | 0:30 | 0:30 |
| 2 | The question | 1:30 | 2:00 |
| 4 | Rational Choice | 1:15 | 3:20 |
| 5 | Selective Social Influence | 1:40 | 5:00 |
| 7 | Imputation pipeline | 1:55 | 7:00 |
| 8 | The experiment (rewire the ties) | 1:30 | 8:30 |
| 9 | Simulation design | 0:45 | 9:15 |
| 11 | Result 1 — maps | 1:25 | 10:45 |
| 12 | Result 2 — non-structural blocking | 1:15 | 12:00 |
| 13 | Result 3 — the rewiring comparison | 2:00 | 14:00 |
| 14 | Conclusions | 1:00 | 15:00 |

**Delivery notes**
- **Slow down at Slide 13** — it is the headline and carries the weak-ties inversion; that is the part this audience has not seen.
- **Release valve if running long:** compress Slide 9 to one sentence ("about four million simulations, then a regression to isolate the network effect") and Slide 11 to ~1:00 (just point at the low-attractiveness region).
- Numbers discipline: always say "**the odds of adoption** fall by 57%" — not "adoption falls 57%". In most of the parameter space this means roughly half the adoption, and you can say that as an approximation ("roughly, adoption is cut in half").
- The words IUL / MUR / MSP are on the slides; in speech prefer "attractiveness", "requirement / aversion", "openness".
- The word "temperature" no longer appears in the main deck. If someone asks about the tipping or the physics in Q&A, backups 19 (ferromagnet) and 20 (tipping compared) are your moment.
