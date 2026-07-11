# Technical Description and Usage Guide

## High-Level Description

The Statistical Power Simulator is an interactive tool designed to evaluate how statistical power scales with the number of human subjects ($N$) and the number of test items or trials ($M$) in repeated-measures experimental designs. Its primary purpose is to expose the non-linear trade-offs between recruiting more participants and asking more questions, particularly in small-$N$ applied research contexts (e.g., UX testing, system evaluations).

The simulation generates synthetic data under an additive mixed-effects structural model. It assumes that observations follow a continuous Gaussian distribution and that individual responses are a linear combination of latent participant ability, latent item difficulty, measurement error, and the experimental treatment effect.

### Estimation Strategy: Clark's $\min F'$ vs. Multilevel Modeling

In modern psychometrics, crossed repeated measures are typically analyzed using linear mixed-effects models (e.g., `lme4::lmer()`). However, running iterative Maximum Likelihood optimization hundreds of times across a parameter grid is computationally prohibitive for a browser-based WebAssembly application. 

To bypass this without sacrificing rigor, the simulator implements Herbert Clark’s (1973) quasi-$F$ statistic ($\min F'$). This closed-form calculation derives a lower-bound $p$-value by mathematically combining a by-subject variance ratio ($F_1$) and a by-item variance ratio ($F_2$). This correctly penalizes the test for both subject-by-condition and item-by-condition interactions, matching the statistical behavior of a crossed random-effects model in a fraction of the computational time.

### Derivation of the Output Table

The output table displays statistical power as a percentage. For each combination of $M$ and $N$, the engine runs a Monte Carlo simulation loop (defaulting to 200 independent experiments). In each iteration, it generates synthetic datasets from the underlying parameters, applies the selected statistical test, and records the resulting $p$-value. 

The displayed power is simply the proportion of those simulated experiments where the result reached statistical significance ($p < 0.05$).

---

## Parameter Definitions and Effects

| Parameter | Function and Effect in Simulation |
| :--- | :--- |
| **Study Design** | *Between-subjects* splits $N$ across conditions; power asymptotes quickly as $M$ increases. *Within-subjects* exposes all $N$ to both conditions (splitting items); ability cancels out, removing the subject-intercept limit on power. |
| **Sample Sizes ($N$, $M$)** | Defines the dimensions of the power grid. |
| **Effect Size** | The target true latent shift caused by the intervention, expressed in Cohen's $d$ standard deviation units. |
| **Reliability** | The ratio of true score variance to total observation variance ($0$ to $1$). Lower reliability introduces measurement error, attenuating the observed effect size. Increasing $M$ mathematically averages out this noise. |
| **N_weight / M_weight** | The standard deviations of the baseline subject abilities and item difficulties. `N_weight` acts as the hard asymptotic ceiling for between-subjects power. |
| **sd_subj_slope** | Variance in how much individuals respond to the treatment. Values $>0$ create a "plateau effect" where adding more questions ($M$) eventually stops improving power, necessitating more subjects. |
| **sd_item_slope** | Variance in how much the treatment impacts different stimuli/tasks. Values $>0$ create a plateau where adding subjects ($N$) stops improving power, demanding a larger item pool. |

---

## Grounding Parameters in Empirical Data

To set realistic baseline variances (`N_weight`, `M_weight`) and treatment heterogeneities (`sd_subj_slope`, `sd_item_slope`), you can extract these directly from pilot data using a linear mixed model. Ensure your pilot data is unaggregated (trial-level).

### 1. Estimating Subject Parameters

Fit a model estimating baseline ability and individual treatment response:

```r
library(lme4)

# Fit model estimating baseline ability and individual treatment response
fit_subj <- lmer(score ~ condition + 
                   (1 + condition | subject) + 
                   (1 | item), 
                 data = pilot_data)

summary(fit_subj)