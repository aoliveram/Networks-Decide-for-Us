# Theoretical Perspectives - Already Reviewed
*Based on the “Networks Decide for Us” project presentations and reports.*

This is a living repository of sociological/theoretical notes that have already been integrated into the project’s framing. Each entry records what is conceptually “taken” from the work (and, when relevant, how it maps to the project’s mechanisms: individual utility/incentives, social influence/peer pressure, thresholds, and observed adoption as an emergent outcome of the networked system).

## A. Discontinuous Adoption, “Tipping Points”, and Case Narratives

### 1) Gladwell, M. (2002). “The tipping point: how little things can make a big difference”. Back Bay Books, Boston, 1st back bay pbk. edition.
- Role in the project: mainly a *catalogue of cases* where adoption changes appear abrupt (i.e., the “tipping point” intuition), often challenging what one would expect from smooth S-shaped diffusion curves.
- Contribution: useful as a motivation for studying discontinuities/threshold-like dynamics, but it is not a formal modeling framework.
- Illustrative cases described in the book:
  - Between the end of 94’ and the beginning of 95’, Hush Puppies went from selling 30k pairs of shoes to 430k, without paid publicity.
  - East New York and Brownsville crime reduction: In five years, murders fell by 64.3 percent, while total crimes fell by almost half.
  - Abrupt growth of fax users: In only two years, fax sales went from selling a few thousand to more than a million.
  - Book: 'Divine Secrets of the Ya-Ya Sisterhood' was a book that everyone read, even when it doesn't score very well in lists. The key was small, cohesive lecturer groups.

### 2) Farrell, W. (1998). “How hits happen: forecasting predictability in a chaotic Marketplace”. HarperBusiness, New York, 1st edition.
- Role in the project: similar to Gladwell, it collects narratives of “hits”/sudden popularity and market tipping dynamics.
- Illustrative cases described in the book:
  - Movies: The Full Monty
  - Market Share: AT&T (’once the company began to lose market share customers seemed to flee to competitors at a desperate pace’)
  - Literature: Hootie & the Blowfish.

## B. Collective Action, Collective Behavior, and Threshold Discontinuities (Formal Modeling Lineage)

### 3) Oliver, P. E. (1993). Formal models of collective action. Annual Review of Sociology, 19, 271–300.
- Why this matters for NDFU: it provides a clean bridge between “collective behavior” and “collective action” as a modelable object, consistent with interpreting diffusion/adoption as the aggregation of *interdependent* individual decisions.
- Definition used to anchor the GSS framing: “any action which provides a collective good"
- Historical/theoretical positioning (as used in the project’s narrative):
  - Pre-Olson assumptions about “natural” mobilization.
  - Olson (1965), *The Logic of Collective Action*: the free-rider problem and why inaction became the default expectation.
  - Post-Olson formal work: individual decision models → aggregation/interdependence models → collective decision models → dynamic interaction (movement–regime) models.
  - Critical mass theory and threshold/discontinuity results appear as recurring patterns in formal models.

### 4) Granovetter, M. (1978). Threshold models of collective behavior. American Journal of Sociology, 83(6), 1420–1443.
- Foundational threshold definition:
  - “The threshold is simply that point where the perceived benefits to an individual of doing the thing in question exceed the perceived costs.”
- Scope claim (beyond “collective behavior” narrowly defined):
  - "Further, once we abandon the definition of collective behavior situations as those in which people develop new norms or abandon existing ones, the range of situations which can be considered broadens. Thus, these models can be applied to processes not usually called "collective behavior," such as voting, residential segregation, diffusion of innovations, educational attainment, strikes, migration, and markets-as well as the more typical processes of crowd behavior and social movements."
- Binary decision focus:
  - "… the analysis is meant to apply to any appropriate binary decision.
    (…) Rumors, Strikes, Voting, Diffusion of Innovation, etc."
- Key interpretation for NDFU: thresholds are not “preferences” per se; they encode when utility becomes positive in context, and thus can generate behavior that conflicts with stated norms:
  - "Most did not think it 'right' to commit illegal acts or even particularly want to do so. But group interaction was such that none could admit this without loss of status; in our terms, their threshold for stealing cars is low (...). The boys act because their threshold is exceeded, and their utility is maximized, given the situation, by joining in the criminal activity. But in so doing they act contrary to norms they actually hold. That this is so indicates not that norms are irrelevant, but rather that they are only one causal influence on behavior and are not always decisive."
- Distributional sensitivity and discontinuity (micro similarities can yield macro divergence):
  - "There is no obvious sociological way to explain why a slight perturbation of the normal distribution around the critical standard deviation should have a wholly discontinuous, striking qualitative effect. (...) This is particularly the case since sociological theory is not at all oriented to analyzing the effects of changes in exact distributions of properties but concentrates, rather, on the effects imputed to average values. This example shows again how two crowds whose average preferences are nearly identical could generate entirely different results."
- Open note/critique to keep in mind: the paper can be read as treating the participation threshold as the only indicator of “individual preference”, potentially hiding relevant details that may need to be modeled explicitly depending on the application.

## C. Diffusion of Innovations in Networks (Thresholds with Relational Information)

### 5) Valente, T. W. (1996). Social network thresholds in the diffusion of innovations. Social Networks, 18(1), 69–89
- Core shift: from whole “social systems” to *personal networks* as the reference frame for thresholds and interpersonal influence:
  - "The present article measures thresholds with respect to personal networks rather than whole social systems to understand more fully the role of interpersonal influence in adoption behavior."
- Framing: diffusion of innovations as a form of collective/aggregate behavior produced by interdependence in adoption decisions.
- Connection to classic diffusion work (time-of-adoption categories as a baseline tradition):
  - "A major contribution to diffusion research has been the categorization of adopters based on innovativeness as measured by time-of-adoption (Rogers, 1958). Adopters are classified as (1) early adopters, (2) early majority, (3) late majority, and (4) laggards (Ryan and Gross, 1943, 1950; Beal and Bohlen, 1955; Rogers, 1983, pp. 245-247)."
- Field definition used in the project notes:
  - “a theoretical approach used to explain how new ideas and practices spread within and between communities.”

### 6) Centola, D. and Macy, M. (2007). Complex Contagions and the Weakness of Long Ties. American Journal of Sociology, 113(3):702–734.
- Core contribution for NDFU: formalizes *structural* barriers to the diffusion of “complex” behaviors, and explicitly limits the scope of the “strength of weak ties” claim for certain kinds of diffusion.
- How they position their argument relative to Granovetter’s “weak ties” claim:
  - “Nevertheless, the central claim of this study is the need to circumscribe carefully the scope of Granovetter’s claim.”
  - Granovetter (1973) baseline claim (as quoted in this project’s notes): “whatever is to be diffused can reach a larger number of people, and traverse a greater social distance, when passed through weak ties rather than strong.”
- Simple vs. complex diffusion (definition is about sources, not repetition):
  - “The distinction between simple and complex refers to the number of sources of exposure required for activation, not the number of exposures.”
- Project connection (NDFU note): beyond structural barriers, the simulations also show non-structural barriers in the form of phase transitions—small increases in MSP or IUL can generate abrupt changes in aggregate adoption.

### 7) Tur, E. M., Zeppini, P., and Frenken, K. (2018). “Diffusion with social reinforcement: The role of individual preferences”. Physical Review E, 97(2):022302.
- Why this matters for NDFU: it explicitly adds heterogeneous, individual preferences into a diffusion-with-reinforcement (complex contagion / threshold-like) setting:
  - “In all these cases, agents decide solely on the basis of the number of neighborswho have already adopted. In this paper,we argue that one dimension of social diffusion still remains underspecified: individual preferences of agents.”
- Modeling choice noted for comparison: a percolation-based approach where adopters exert uniform influence on ego (no modulation by alter–ego social distance/similarity).
- Reinforcement mechanism (as interpreted in these project notes): social reinforcement lowers the ego’s minimum utility requirement, directly shifting the effective preference/threshold during diffusion.
  - NDFU contrast: preferences (captured by the Minimum Utility Requirement) are not modified during diffusion; selectivity affects whose exposures “count” rather than shifting the preference itself.
- Scope/limitations to track when comparing with empirical settings:
  - Preferences are drawn from artificial distributions (may or may not resemble real populations).
  - Networks are based on Watts–Strogatz small-world topology.
- Additional results emphasized in these project notes:
  - Homophily removes the standard nonlinearity of percolation: when people with similar preferences are clustered, diffusion scales linearly with preferences (“those who want to adopt, adopt”).
  - Social reinforcement changes this: with strong homophily, reinforcement makes diffusion more efficient in clustered networks; without homophily, the opposite can occur.
  - Final diffusion size cannot be understood from individual preferences alone; it depends on the joint configuration of network structure, homophily, and reinforcement.

### 8) Tur, E. M., Zeppini, P., and Frenken, K. (2024). “Diffusion in small worlds with homophily and social reinforcement: A theoretical model”. Social Networks, 76:12–21
- What the model mixes: network structure (small worlds), homophily, and social reinforcement (threshold-like reinforcement).
- Homophily is treated explicitly as a structural property (who is connected to whom), not as a cognitive “influence bias” nor as a differential weighting of peers.
  - In this framing, homophily is the degree to which individuals with similar minimum utility requirements are topologically clustered (the network is rewired before diffusion so that similar agents are connected).
  - The paper motivates two ways to implement “homophily”: (1) as bias in receptivity to influence, and (2) as a mechanism to construct the network from fixed attributes; this work uses the second.
- Parameterization (as captured in these project notes):
  - $\rho$ → homophily (≈0: random wiring; ≈1: “with guys like you”).
  - $\xi$ → intrinsic utility levels, compared with a “minimum utility requirement”.
- Main theoretical implication emphasized in these notes: introducing structural homophily removes the classic percolation “hit/flop” nonlinearity; in the limit of perfect homophily, diffusion becomes linear and fully predictable (everyone who is willing to adopt does adopt, because information circulates within the “right” groups).
- Phase-transition detail flagged for NDFU mapping: Figure 5 is noted as showing examples of second order phase transitions when social reinforcement is present (threshold > 1) and homophily is low.
- Differences vs. NDFU (project note):
  - NDFU uses social distances/similarity to define selective exposure; this model does not.
  - NDFU’s selectivity effectively blends “simple” + “complex” contagion: if intrinsic conditions suffice, a node adopts; otherwise, adoption requires enough similar adopters in view.
  - Percolation-theory interpretation of “complex contagion” differs: social reinforcement operates by lowering the minimum utility requirement.

### Notes / pointers used in the project
- **Percolation Theory**: used as a general language to think about connectivity-dependent diffusion and “giant component” conditions for spread (especially in word-of-mouth style processes).
- To see percolation theory applied to “threshold models,” see Tur (2018).
- To see the original Percolation model, see the book Staufer (2003).
- To see a paper that moderates the influence of network contacts depending on how similar they are, reflecting a confirmation bias, see Konc and Savin (2019).
- Centola (2011) is used in the project notes as a reference point for the idea that social similarity/proximity can confer legitimacy during social contagion (relevant when comparing “uniform influence” assumptions).
- Phase transition mapping noted in the project materials:
  - 1. Tur (2024). ---------------------- yes
  - 2. Konc and Savin (2019)  ------- no, only first transition
  - 3. Tur (2018). ---------------------- yes, but without homophily

## D. Homophily & Social Structure
- **McPherson & Smith-Lovin (1987, 2001)**: Homophily (“Birds of a feather flock together”).
  - **Structural Homophily**: the demographic composition of the network limits who interacts with whom.
  - **Induced Homophily**: homophily arising from psychological preference for similarity.
  - **Imputation Method (as used here)**: using ego-net data (GSS) to calculate homophilic strength (Age, Sex, Education, Race, Religion) and imputing plausible networks to other populations (ATP).
  - Significance for NDFU: this addresses “structural reductionism” by using plausible topologies derived from real data rather than purely synthetic ones (e.g., ER or BA graphs).

### 9) Smith, J. A., McPherson, M., and Smith-Lovin, L. (2014). "Social Distance in the United States: Sex, Race, Religion, Age, and Education Homophily among Confidants, 1985 to 2004." American Sociological Review, 79(3), 432–456.
- Why this matters for NDFU: this is the primary empirical source for the homophily-strength parameters used to build both ATP and GSS plausible networks.
  - In NDFU, the homophily strengths are taken from Table 3, Model 2 (as stated in the project notes).
- What the paper estimates: the strength of homophily in U.S. personal networks, interpreted as a penalty on tie probability as sociodemographic distance increases.
- Data: GSS core discussion networks module (1985–2004), using ego–alter dyads.
- Modeling framing: tie formation as a direct function of distance in Blau space (sex, race, religion, age, education).
- Results emphasized in the project notes:
  - Strong and systematic homophily across all dimensions; strongest for race and religion, followed by education and age; sex is weaker.
  - Penalties persist after controlling for opportunity structure (marginal distributions), so they are not a compositional artifact.
  - Homophily changes little from 1985 to 2004, suggesting structural stability in confidant selection patterns.
- Methodological point (useful for imputation): dyadic ego–alter data + case–control design + bootstrap yields interpretable “homophily strength” parameters that can be transported to surveys without complete relational data.

### 10) McPherson, M. and Smith, J. A. (2019). “Network Effects in Blau Space: Imputing Social Context from Survey Data”. Socius: Sociological Research for a Dynamic World.
- Why this matters for NDFU: it provides the core methodological idea that one can impute social context—and even plausible networks—from standard surveys by transporting homophily parameters estimated from ego-network data.
- NDFU-specific use case: leverage rich individual-level survey measures (e.g., “innovativeness” or “propensity to collective action” scores) in datasets without networks by simulating plausible topologies over the observed sociodemographic composition.
- Motivation: many individual outcomes reflect unobserved relational processes (homophily + influence), so standard regressions can misattribute contextual effects to individual covariates.
- Primary use case in the paper (as captured in the project notes): imputing latent social context
  - Estimate homophily parameters in a survey with ego-network data (GSS) by modeling interaction probability as a function of distance in Blau space.
  - Transport these parameters to another survey from the same population without network measures (ANES).
  - Build a social-weights (latent neighborhood) matrix and compute an imputed contextual term for each person (weighted average of the outcome in their inferred social context).
  - In spatial/SAR-style models, adding the contextual term substantially reduces demographic effects, consistent with “individual effects” being partly contextual.
- Methodological extension (how NDFU uses the logic): generating full plausible networks
  - Homophily coefficients can be interpreted as relative tie preferences and used in an ERGM together with a density term to simulate graphs that respect:
    - empirically grounded homophily intensity,
    - realistic density,
    - and the observed sociodemographic composition of the target survey.
  - This shifts the approach from latent context to explicit network topologies, enabling dynamic simulations (diffusion, influence, contagion) on empirically anchored graphs.

### 11) Centola, D. (2011). An Experimental Study of Homophily in the Adoption of Health Behavior. Science, 334(6060):1269–1272.
- Why this matters for NDFU: it provides causal evidence that homophily operates at two distinct levels in diffusion, which directly motivates models where social influence is modulated by ego–alter similarity (rather than assuming uniform influence from all adopters).
- Two levels of homophily in diffusion (as summarized in the project notes):
  - Structural homophily: networks that connect similar people increase exposure to the behavior.
  - Homophily as an influence bias: similar alters are more influential than dissimilar alters even when exposure is held constant.
- Key empirical implication: conditional on exposure, adoption is more likely when ego and alter are closer in social distance (Blau-space-like attributes; e.g., health-related characteristics).
- Design features emphasized in the project notes:
  - Controlled experiment (not observational).
  - Conducted in an online laboratory (web platform) with artificial social networks.
  - Network topology is held constant (same network, degree, clustering) while homophily is manipulated exogenously via neighbor assignment.
  - 710 real participants recruited from an online fitness program.
  - Adoption is behavioral (use of an online diet diary), not self-reported intention.
- Central conclusion used in NDFU: homophily increases diffusion not only by increasing exposure, but because social influence is stronger through socially similar ties—homophily is both a network-formation mechanism and a mechanism for weighting influence.

## E. Social Influence vs. Rational Choice (Mechanisms in the Project)
- **Rational Choice**: individuals adopt if the innovation’s utility exceeds their personal threshold, regardless of peers (or with minimal peer information).
- **Social Influence / Social Reinforcement**: adoption is driven by exposure and reinforcement from peers.
- **Selective Social Influence (project mechanism)**: influence is not global/uniform; it is weighted by **Social Proximity** (Gower’s distance). Only peers within a “Maximum Social Proximity” (MSP / h) exert effective influence.
  - Key construct: “Effective Exposure” ($\tilde{E}_i$) depends on both the number of adopters *and* their similarity to ego.

## F. Collective Action & Emergence (Project Interpretation)
- **Olson (1965)**: *The Logic of Collective Action* (free-rider problem; in this project’s language, it motivates why inaction can be the baseline and why individual incentives may be insufficient for large-scale adoption/mobilization).
- **Emergence / micro–macro link**: macro adoption is not a linear sum of individual propensities.
  - Claim: “Networks decide for us”. Structure governs outcomes beyond what is implied by the distribution of individual attributes alone.
  - Empirical/simulation hook in NDFU: comparing plausible networks vs. Erdős–Rényi (ER) networks yields significant drops in adoption, implying clustering/homophily/topology are enabling conditions for diffusion.

## G. Sociological Interpretation (How to Read the Results)
- **Micro–meso–macro**: individual attributes (micro) + network structure (meso) generate societal outcomes (macro).
- **Context of decision**: the network supplies the social context; a “flexible” population (high MSP) allows influence to travel across social boundaries, while a “rigid” population (low MSP) produces echo chambers and diffusion barriers.

## H. Homophily + Social Influence in ABMs (Blau Space / Network Segregation)

### 12) DellaPosta et al (2015) - Why Do Liberals Drink Lattes?
- What the model does (as summarized in these project notes): an ABM that combines homophily and social influence in networks to explain why clusters of *dynamic* attributes (opinions; binary) become associated, conditional on *static* attributes (binary). Initial attributes are random.
- Space/geometry: uses Blau spaces to model positions of sociodemographic attributes as well as dynamic opinions; both static and dynamic attributes are used as a social distance.
- Network substrate (as used in these notes): a “connected caveman” graph resembling some properties of social networks; as paraphrased in Goldberg (2018): “connected caveman” small-world network in which individuals are segregated into sparsely interconnected and densely intraconnected clusters”.
- Project takeaway: homophily + influence can reproduce nonrandom “bundles” of opinions (e.g., “liberals” and “drink lattes”) anchored in sociodemographic structure; however, a key caveat repeatedly raised in later work is that cultural differentiation depends on (or is amplified by) an underlying structurally segmented network.

## I. Associative Diffusion, Cognition, and Structural Reductionism

### 13) Goldberg, Stein (2018) - Beyond Social Contagion: Associative Diffusion and the Emergence of Cultural Variation
- Core idea: cultural differentiation can emerge from social cognition without requiring a segregated network or preexisting division into groups (“associative diffusion”).
- Two cognitive ingredients highlighted in these project notes:
  - Cognitive associations between concepts (a network where concepts are linked if they are associated).
  - Constraint satisfaction: agents search for congruence among associations; if exposure would move beliefs in an incongruent direction relative to prior associations, it is dropped.
- Critique of homophily-based diffusion accounts: in models like DellaPosta et al. (2015), cultural differentiation is ultimately epiphenomenal to a structurally segmented world; moreover, these models do not require homophilous interaction if interactions are random.
- Quotes recorded in the project notes:
  - 1. Culture is often measured as the distribution of beliefs in a population. A consistent finding in the sociology of culture is that beliefs are not randomly distributed. (Martin 2002).
  - 2. To explain systemic cultural variation, social contagion theories need to assume the existence of structural boundaries to diffusion. The first mechanism is exposure (social influence, right?): models generates social differentiation as long as network ties are denser within groups than they are between them, so there is a structural boundary (example, Centola and Macy 2007). The second is choice homophily: Such a proclivity leads to the emergence of culturally homogenous clusters either because individuals choose to interact homophilously. But also, scholars have proposed a variety of models exploring the complementary effects of homophily and social influence: some emphasize the formation and dissolution of network ties (e.g., Carley 1991; Centola et al. 2007; Mark 1998); others focus on changes in the strength of social influence as a function of cultural similarity (e.g., Dandekar et al. 2013; Flache and Macy 2011)

### 14) Goldberg, A. (2021). “Associative Diffusion and the Pitfalls of Structural Reductionism”. American Sociological Review, 86(6):1205–1210.
- Short response piece (as used in these notes) that argues for incorporating more substantive characteristics of nodes (e.g., cognitive content) rather than relying on structurally reduced homophily-only mechanisms.
- Key claim emphasized in the project notes: homophily-based models (e.g., DellaPosta et al. 2015) create cultural differentiation only if the underlying network structure is already segmented.

## J. Idea Diffusion in Science (Empirical / Relational Ecologies)

### 15) Cheng et al. (2024) - How New Ideas Diffuse in Science
- Research question (as captured in these notes): what explains ideas that are capable of being maintained “until now” (long-run persistence/career)?
- Method: linear regressions.
- Key finding emphasized in these notes: relational embeddedness in both an ideational ecology and a social ecology predicts long-run diffusion (Table 4; Model 4 adds social and ideational prominence/consistency/embeddedness).
- Variables highlighted in these notes (Table 4, Model 4): Social prominence, Social consistency, Social embeddedness, Ideational prominence, Ideational consistency, Ideational embeddedness.
- Quote recorded in the project notes:
  - “We found that the ways ideas become relationally situated in their ideational ecology is especially decisive for their long-run career. Ideas tend to diffuse more broadly when their authors link them to prominent, foundational ideas within science, consistently interrelate and position them within an established network of ideas, and thickly interrelate them-filling out research topics and networks with greater links. We also found that the successful diffusion of ideas depends on their social ecology, that is, the communities of scholars who take them up. Specifically, new ideas tend to diffuse when focally central researchers publish on them and when those authors span diverse research communities. Table 4 highlights this story.”

## K. Reviews and Baseline Diffusion Models (Bass / Word-of-Mouth)

### 16) Guidolin (2023) - Innovation Diffusion Processes Concepts Models Predictions
- Notes focus: review of diffusion models (including Bass).
- Quote recorded in the project notes (word-of-mouth often dominates):
  - “Here, we rely on the empirical fact that, in most innovation diffusions, word-of-mouth—the internal or epidemiological component—exceeds the external component, which is often negligible or absent (Bass 1969, Guidolin & Mortarino 2010, Rao & Kishore 2010, Bunea et al. 2020).”
- Bass model framing (as summarized in these notes): an ODE with internal and external influence terms; a simple toy model for aggregate adoption.

## L. Adoption, Disadoption, and Temporal Alignment (Bi-threshold Models)

### 17) Alipour (2024) - Enough but not too many: A bi-threshold model for behavioral diffusion
- Core idea (as captured in these notes): a bi-threshold model with disadoption; the disadoption mechanism can be interpreted as (1) negative network externalities, (2) distinction, and (3) saturation.
- Empirical application (as noted here): predicts adoption rates (not always monotonic) of hashtags on Twitter (?) in news-related diffusion; reported to perform well (open question in these notes about “complexity” of this behavior).
- Methods pointer in these notes: uses recent advancements for estimating linear thresholds by Tran et al. (23, 24); “NOT SEEING, I’ll wait to see the methods instead of reading the paper”.
- Two empirical contributions recorded in these notes:
  - CTL method (Tran et al.) for bi-threshold.
  - Transmission framing as temporal alignment:
    - “More generally, T_effect and T_login make explicit the importance of temporal alignment for transmission of a social contagion. In other words, network ties do not automatically conduct all available information; rather, transmission depends on the coincidence of the adopter’s actions and the potential adopter’s observation of those actions.”

## M. Social Percolation, Preference Updating, and Hit/Flop Market Outcomes

### 18) Salomon, Weisbuch (2000) - Social percolation models
- Uses percolation theory (percolation threshold: minimum density of nodes/links to allow a global connection) to model word-of-mouth diffusion and how “quality” evolves when nodes update their inherent requirements.
- Algorithm summary recorded in these notes:
  - Use a network with a percolation threshold $p_c$.
  - Give each node an inherent requirement $p$ between 0 and 1.
  - If $p < q$, the node adopts.
  - After adopting / not adopting, update $p$ by $\pm \Delta q$.
- Note for comparison: in Tur (2024), a related percolation logic appears, but social reinforcement lowers the requirement under reinforcement.
- Market “hits vs failures” pointers recorded in these notes:
  - W. Farrell, How hits happen, HarperCollins, New York, 1998.
  - B.W. Arthur, D.A. Lane, Struct. Changes Econom. Dyn. 4 (1993) 81.
  - G. Weisbuch, G. Boudjema, Adv. Complex Systems 2 (1999) 11.
