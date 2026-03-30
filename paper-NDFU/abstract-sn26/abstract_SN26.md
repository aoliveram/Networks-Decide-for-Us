# Extended Abstract for Social Networks Special Issue 2026
**Agent-Based Modelling for Social Network Research**

## Title
Networks Decide for Us: Emergence of Adoption in Homophily-Imputed Social Topologies

## Author Details
- **Name:** Aníbal Olivera Morales
- **Affiliation:** Centro de Investigación en Complejidad Social, Universidad del Desarrollo, Santiago, Chile
- **Email:** an.oliveram@udd.cl

- **Name:** Jorge Fábraga Lacoa
- **Affiliation:** Centro de Investigación en Complejidad Social, Universidad del Desarrollo, Santiago, Chile
- **Email:** jfabrega@udd.cl

- **Name:** Thomas W. Valente
- **Affiliation:** Department of Population and Public Health Sciences, Keck School of Medicine of University of Southern California, Los Angeles, United States of America
- **Email:** tvalente@usc.edu
---

## Research Question
The diffusion of innovations, social norms, and collective behaviors is frequently modeled as either a result of individual utility maximization or a contagion process. However, empirical reality suggests a hybrid mechanism where individuals weigh both their intrinsic preferences and the influence of their peers. Furthermore, the structure of these peer networks is rarely random; it is shaped by deeply ingrained homophilic forces.

This research asks: *To what extent does the structure of social networks, shaped by empirical homophily, independently govern the diffusion of collective behaviors, conditional on fixed individual-level propensities?* We propose that the level of adoption in a society is an emergent property of the relational structure, not just the aggregate of individual wills. We advance the claim that “networks decide for us”: the topological arrangement itself acts as a non-neutral scaffold that actively enables or constrains social change.

## Data and Methods
We address the methodological challenge of “structural reductionism”—where models either lack cognitive depth or assume arbitrary network topologies—by developing a novel imputation pipeline. This combines representative survey data with empirically grounded network parameters:

1. **Homophily Estimation:** We use the General Social Survey (GSS 2004) core discussion networks module to estimate relative homophilic strengths across key dimensions of Blau space: Age, Sex, Education, Race, and Religion.
2. **Imputting Plausible Networks:** These relational structures are imputed onto respondents from the American Trends Panel (ATP 2016). We generate $N=1000$ networks that closely respect the marginal distributions of demographics while strictly adhering to the estimated homophily tendencies, producing realistic clustering coefficients and degree distributions.
3. **Hybrid Agent-Based Model:** We simulate a hybrid diffusion process over these populations by integrating *Rational Choice* and *Selective Social Influence*. 
   - *Rational Choice* triggers adoption if the intrinsic utility of an innovation exceeds the individual’s minimum requirement.
   - *Selective Social Influence* enables adoption through peer contagion. Crucially, the influence of peers is weighted by social proximity (Gower’s distance). An individual adopts socially only if exposed to a critical mass of adopters who are socially similar to them, constrained by a Maximum Social Proximity ($h$) parameter.

Through extensive simulations ($\approx 4$ million runs), we explore the parameter space of intrinsic utility and social flexibility across both our Plausible networks and Erdős-Rényi (Random) baselines of identical density.

## Preliminary Findings
By comparing the diffusion outcomes on Plausible networks against their Random baselines, we quantify the "structural premium" of realistic network topologies. Our Generalized Additive Models (BAM) regressions reveal several key dynamics:

1. **The Structural Premium:** Randomizing the network structure significantly depresses overall adoption. The specific clustering and homophilic segregation inherent in realistic networks facilitate behavioral cascades that random networks cannot sustain. Under random conditions, the odds of total adoption are drastically reduced (BAM coefficient = -0.850, $p < 0.001$).
2. **Non-Linear Phase Transitions:** We identify distinct interaction regimes between intrinsic utility and social flexibility. When the population is highly rigid (low homophilic tolerance), diffusion is stifled regardless of utility. As social flexibility increases, a “tipping point” is reached where social reinforcement overcomes localized resistance, enabling widespread systemic adoption.
3. **The Role of Selective Influence:** The data shows that "who we are" is inextricably linked to "who we know," but the structural pathway dictates the outcome. Our findings directly challenge models assuming uniform influence, emphasizing that homophily operates as a dynamic gate for diffusion.

Our results demonstrate that applying ABM to empirically calibrated network topologies uncovers macro-social outcomes that statistical aggregations miss. The findings contribute to social network research by illustrating that network structure is a fundamental causal mechanism in the dynamic unfolding of collective behavior.

---

*(Word count: ~615 words, leaving ample room up to the 800-word limit)*

### Target Figures and Tables
*For the submission, we will include the following two items from the main paper as allowed by the guidelines (up to 2 tables and 2 figures):*

**Figure 1: Heatmap of Total Adoption**
*(Reference `plots_m04_random_GSS.pdf` in the existing manuscript. This figure visually demonstrates the non-linear phase transitions where adoption jumps from low to high as social flexibility increases).*

**Table 1: BAM Regression Results - Effect of Network Structure and Seeding**
*(Reference Table 1 in the existing manuscript. This table empirically proves the structural premium by showing the negative effect of the ER baseline network on Rational and Social Adoption).*
