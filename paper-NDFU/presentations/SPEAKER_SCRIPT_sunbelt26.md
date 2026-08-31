# Speaker Script — Sunbelt 2026

---

## Slide 1 — Title (0:00–0:30)

> Hi everyone, thanks for being here. I'm Aníbal Olivera, a PhD student in Social Complexity Sciences at Universidad del Desarrollo in Chile.
> This is a work that we have been doing with my advisor, Jorge Fábrega,
> And the idea behind this provocative title is that *our real ties* actually plays a crucial role on the diffusion of behaviors.

---

## Slide 2 — Collective behavior spreads on networks (0:30–1:15)

> Here's the setup. We know that Collective Behaviors (like adopting an innovation, or joining a protest) spread as adoption cascades on social networks.
> But there are a couple of gaps...
> for example, most models use *synthetic* networks, where people are simple relays that just pass things along.
> But real human networks have a very specific texture: we usually cluster with people like us. Similar age, similar education, similar background; so we say that real networks follow **multidimensional homophily**.
> Also, even when we know individuals have different preferences for adopting behaviors, which can be measured in surveys,
> still, a lot of studies use peer influence as the only mechanism for adoption.

---

## Slide 4 — Where we go (2:15–3:00)

> So we make two moves. First, we *build* realistic networks for real survey people. More details later...
>
> Second, we measure how much that structure actually *causes* the spreading.
>
> In this way, we can ask ourselves: conditional on *fixed individual characteristics*, does collective behavior depend on the network's *shape*, or on *who sits next to whom*?

---

## Slide 5 — Adoption, part 1: Rational Choice (3:00–3:45)

> The model presents two ways a person can adopt. The first is simple and selfish: if the thing is attractive enough to *me*, I will adopt it.
> Following the literature, we say that each collective behavior has an *Intrinsic Utility Level* 'IUL' that is a number between 0 and 1.
>
> And also each person has a *Minimum Utility Requirement* 'MUR', that represents their aversion to adopting.
>
> So, an agent will adopt if the collective behavior is attractive enough for his *requirement*.
>
> You can see in the plots that this channel is very similar to a simple contagion, but with a rule that have to be met.

---

## Slide 6 — Adoption, part 2: Selective Social Influence (3:45–4:45)

> The second way is by social influence; because when several of my contacts have already adopted, the cost of joining gets lower over time, and I will eventually adopt, even if it didn't reach my minimum utility requirement.
>
> But a key ingredient of the model is that it is selective: **I don't count everybody** as effective influence. I'll listen only to those adopters close enough to me in terms of *social distance* d_ij.
>
> Following the literature, there is a **Maximum Social Proximity** 'MSP' where the limit of influence sits.
> So, when *MSP* is low, I only listen people very similar to me —*echo chambers*-, and when *MSP* is high, the society is more flexible.

---

## Slide 7 — Imputing a network onto a survey (4:45–6:15)  ← the methodological heart

> Now, where do the networks come from?
> We borrow the *homophily strengths* —that is how *strongly people sort* by demographics like age, education, or sex— from a nice paper from McPherson & Smith, and use them to construct realistic networks using the data from surveys -that doesn't have structural information- to *feed a Logit that reproduces those homophilic tendencies*.
> And the cool thing is that that data has items that allow us to construct whether an *aversion-to-innovation score* or a *collective-action propensity score* for each node.
>
> Those networks have the right topology, just as the literature shows the real-social-networks topologies are.
> Something that synthetic models simply cannot reproduce.
>
> So now we have something rare: realistic networks *and* real individual propensities, for people who never gave us their relational data.

---

## Slide 8 — Two ways to "break" the network (6:15–7:15)

> Now, before running anything, I want to set up the experiment carefully.
>
> So first, we will call the plausible networks as "GSS", just because we will use the General Social Survey respondents' data, which have items to construct a "collective action" propensity score.
>
> And we will compare the results against "broken" versions.
> In the first one —the *SH* —, we keep the network *exactly* as it is —every single triangle and clique stays the same— and we just *shuffle the people around between the nodes* -so it's like picking a node and, leaving all the connections the same, sending all the characteristics of that node to another place. That's why we call it SH, from shuffle.
> So the topology is identical; the *only* thing we destroy is homophily.
>
> On the second — called *DP* — we preserve the number of connections of each person, but we cut the edges and rewire them *randomly*.
> So now, *both* the homophily and topology are destroyed.
>
>
> Therefore, comparing the realistic to the shuffle version isolates the effect of that lack of *homophily* in adoption;
> and comparing the shuffle to the rewired version tell us how much the *topology* is adding on top of that.

---

## Slide 9 — What breaking destroys: topology (7:15–8:00)

> So let's have a look at how the topology measures are changed.
>
> In this plot we see a lot of *topology* measures —clustering measures, number of triangles, degree assortativity, average path length, etc.;
> And in all of those the plausible networks and the shuffle versions have *identical* values, because shuffling people characteristics doesn't move a single tie. But, the version with rewire shows very different values.
> So this double-checks that topology only changes when we rewire.. (actually all the measures that are high when homophily is present)

---

## Slide 10 — What breaking destroys: homophily (8:00–8:50)  ← the mirror

> And what about *homophily* ?
> If we see the assortativity by each demographic —age, education, race, religion, sex-, we can see that the highest levels are shown by the realistic networks, and in both the shuffle and the rewired versions, the homophily is destroyed.
> On the right we see the average *social distance along the ties*... And we see how the mean social distance is low in the realistic network because your neighbors are *more or less similar* to you because of homophily... but if you modify that composition, your neighbors become *random strangers* and the social distance jumps up. 
> ...remember that those values are used in our selective social influence mechanism.
>
> Ok, so the "broken" versions are exactly what we expected them to be.

---

## Slide 11 — Simulation design (8:50–9:25)

> And we just need to run the model — a lot. We sweep the attractiveness of the behavior, the openness parameter, the thresholds, different ways of seeding the first adopters...
> 
> In total, we have around four million simulations.
>
> And, crucially, we then fit a Big Additive Model to compare the results of each kind of networks.

---

## Slide 12 — Result 1: How adoption plays out (9:25–10:00)

> First, the dynamics. These three plots map total adoption (and its decomposition into rational and social influence terms) over the *attractiveness–openness plane*... and you can see how *rational adoption* grows as the attractiveness grows, but *social adoption* is more complex, so we can have similar total adoption over a large region just because social adoption is relevant **even when the attractiveness is low**.

---

## Slide 13 — Result 2: Breaking structure cuts adoption (10:00–10:50)  ← headline

> Now, how are the adoption results when we compare with the broken versions?
> Well, when we break the network (keeping the exact same people and preferences, just changing how they are connected) adoption drops by more than half.
> The rewire cuts the odds by 57 percent; the shuffle by 53.
> So they are very similar values. And remember, in *both cases the homophily is destroyed*
>
> Think about this for a second: if they are the same individuals, with the same preferences, why does adoption drop off so drastically?
>
> The thing is that we can think that homophily creates echo chambers, and echo chambers turn out to be *useful* here, because they let social reinforcement build up safely in a corner until ot can ignite (with a lot of redundancy, even those with high threshold can reach their requirements).
>

---

## Slide 14 — Result 2: Homophily alone reproduces most of the drop (10:50–11:50)  ← the new, clean result

> But what's the role of topology?
> 
> The full rewire — that's the green one — gives a drop of 57%. But remember, that mixes *both* broken homophily and broken topology.
> The shuffle — the blue one — breaks *only* homophily, keeping the topology exactly in place, and it already gives a 53% drop. That's about *91 percent* of the full effect.
>
> So homophily, on its own, recovers most of the drop — the realistic topology adds only a little on top of it.
> Then, *Homophily is doing most of the work*. So if you are simulating with networks that sits the wrong people next to each other, you're missing the piece that matters most.

---

## Slide 15 — Conclusions (11:50–13:00)

> So, let's draw some conclusions:
> A *method*: you can build realistic, simulable whole-networks for survey respondents out of empirical homophily. That's reusable anywhere you have people's attitudes but no network.
>
> Two *findings*:
> 1. Social influence explains adoption even for non-attractive behaviors.
> 2. A realistic structure gives you a +57% adoption premium, over the same individuals -but with their ties randomised-.
> 3. In realistic social networks, topology it not all: *who sits next to whom* is very important.
>
> A *mechanism*: what's doing the work is **homophily** — your neighbors being similar to you — not the abstract topology.
>
> And let me close the loop by emphasizing the cases where we have high adoption for non-atractive behaviors: There are studies that report the same phenomenon, finding that *random noise* explains those cases. But here, what drives those outcome isn't stochastic noise — it's *measured preference* along with *realistic structure*.
>
> So the big message is: where a behavior ends up is not the sum of individual preferences, but it's **governed by the structure of our real ties. That's why Networks decide for us**. Thank you.





--------------------------------------------------------------------------------



## Slide 16 — References / Slide 17 — Thanks (13:00–13:20)
*(Go to Thanks; take questions. The QR links to the preprint draft.)*

> Thanks — the QR there links to the draft, and I'm very happy to take questions.

---

## Backup slide (only if time remains, or in Q&A) — the shape of change

> If there's time, one teaser. Beyond *how much* spreads, the structure also governs *how* the change arrives. On a real network the tipping is smooth and predictable — a continuous transition. On the rewired one, it's abrupt — all or nothing.
>
> And there's a deeper twist that fits today's double story. There's a signature physicists call *universality* — basically the same math as a phase transition in matter. The exponent gamma is about 1 for the real networks, about 4 for the rewired ones — so gamma *tells them apart*. But here's the nice part: the *shuffle*, which has no homophily, still keeps gamma near 1. So this smoothness is anchored by the *topology*, not the homophily. Which means: homophily drives *how much* spreads, and topology drives *how* it tips. Two different structural things, two different jobs. That's a whole companion paper — happy to dig in.

---

### Timing cheat-sheet (target ~14 min)
| Block | Slides | Target |
|---|---|---|
| Motivation + gap + where-we-go | 2–4 | 0:30–3:00 |
| Model (rational + social) | 5–6 | 3:00–4:45 |
| Imputation pipeline | 7 | 4:45–6:15 |
| Two ways to break + diagnosis (topology, homophily) | 8–10 | 6:15–8:50 |
| Simulation design | 11 | 8:50–9:25 |
| Result 1 (adoption maps) | 12 | 9:25–10:00 |
| Result 2a (cuts adoption) + 2b (culprit = homophily) | 13–14 | 10:00–11:50 |
| Conclusions | 15 | 11:50–13:00 |
| Thanks / Q&A | 16–17 | 13:00 onward |
| *(Backup: shape-of-change + γ — only if asked)* | — | — |

**Notes for delivery:**
- This is the spine of the whole talk: **homophily → how much** (Result 2), **topology → how it tips** (backup). Slides 8–10 set it up, slide 14 pays it off, the backup completes it.
- Macy appears twice: as *motivation* in "The gap", and as the *closing loop* in Conclusions. No standalone middle slide.
- **Release valve if running long:** on slide 9 (topology), just say "the real network and the shuffle are identical here — only the rewire breaks the topology", and move straight on. The two diagnosis slides should feel like one quick build.
- Slide 14 (the culprit) is the freshest result — slow down a touch there; that's the part this audience hasn't seen before.

### Anticipated questions
- *"Are the imputed networks validated?"* → they reproduce the right-skewed degree distribution and realistic clustering of real networks; checking they match held-out tie features is the natural next step.
- *"Isn't MSP just a macro knob you turn?"* → yes, here it's a global parameter; we read it as a population's *social porosity*. Making it vary by person/group is the next step (and the honest answer to the bottom-up critique).
- *"Why the shuffle on top of the rewire?"* → exactly to separate homophily from clustering. The rewire breaks both; the shuffle breaks homophily *only*, with the topology byte-identical. That's how we get the 91% number.
- *"What's this universality thing?"* → near the tipping point the aggregate follows mean-field scaling; the real network gives a continuous transition (γ≈1), the rewire an abrupt one (γ≈4), and the shuffle stays continuous — so it's the topology, not homophily, behind it. Companion paper.
- *"Homophily coefficients are from 2004 data?"* → they come from the same studies that validated the method; homophily is among the most stable social patterns, and we transport them without re-estimating, following McPherson & Smith.
