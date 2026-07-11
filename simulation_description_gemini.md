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

To move beyond hypotheticals, you can extract realistic baseline variances (`N_weight`, `M_weight`) and treatment heterogeneities (`sd_subj_slope`, `sd_item_slope`) directly from empirical pilot data. 

Your pilot dataset must be unaggregated (trial-level format), meaning every row represents one response from one subject on one item. To extract item slopes, your experimental design must be fully crossed (i.e., items must appear in both the control and treatment conditions across your participant pool).

### Step 1: Fit the Maximal Mixed Model

Using the `lme4` package in R, fit a model that estimates random intercepts and random slopes for both subjects and items:

```r
library(lme4)

# Fit the model with crossed random intercepts AND random slopes
fit_max <- lmer(score ~ condition + 
                  (1 + condition | subject) + 
                  (1 + condition | item), 
                data = pilot_data)

# View the variance components
summary(fit_max)
```

### Step 2: Extract the Variance Components

In the model summary, locate the `Random effects` section. It will look similar to this:

```text
Random effects:
 Groups   Name        Variance Std.Dev. Corr
 subject  (Intercept) 1.4400   1.2000       
          condition   0.1225   0.3500   0.15
 item     (Intercept) 0.6400   0.8000       
          condition   0.0900   0.3000  -0.10
 Residual             2.2500   1.5000       
```

You need the values from the **Std.Dev.** column:
* **Subject (Intercept) `1.20`**: Represents baseline individual differences.
* **Subject Condition `0.35`**: Represents variance in human treatment responsiveness.
* **Item (Intercept) `0.80`**: Represents baseline item difficulty variance.
* **Item Condition `0.30`**: Represents variance in task treatment effectiveness.

### Step 3: Translate to Simulator Parameters

You have two options for inputting these values into the simulator.

**Option A: Use Raw Empirical Values (Recommended)**
Enter the extracted standard deviations directly into the simulator UI or CLI script. 
* `N_weight` = 1.20
* `M_weight` = 0.80
* `sd_subj_slope` = 0.35
* `sd_item_slope` = 0.30
* `effect` = Your empirical fixed effect estimate for `condition` divided by the total pooled baseline standard deviation.

**Option B: Standardized Ratios (Leaving Weights at 1.0)**
If you prefer to leave `N_weight` and `M_weight` at their default `1.0` settings to think purely in standardized units, you must convert your empirical slope standard deviations into ratios relative to their respective intercepts:

$$ \text{Simulator sd\_subj\_slope} = \frac{\text{Empirical Subject Slope SD}}{\text{Empirical Subject Intercept SD}} = \frac{0.35}{1.20} = 0.29 $$

$$ \text{Simulator sd\_item\_slope} = \frac{\text{Empirical Item Slope SD}}{\text{Empirical Item Intercept SD}} = \frac{0.30}{0.80} = 0.375 $$

This perfectly calibrates the simulator's treatment heterogeneity to match the exact signal-to-noise ratio observed in your empirical environment.

---

## Benchmarking Against Classical Power Formulas (`pwr`)

Because standard power calculators like the R `pwr` package assume a fixed, deterministic world devoid of heterogeneous slopes or measurement error, comparing the simulator to `pwr` requires temporarily turning off the simulation's psychometric realism.

To successfully replicate `pwr.t.test` outputs:
1. **Remove Treatment Variance:** Set `sd_subj_slope = 0` and `sd_item_slope = 0`.
2. **Remove Measurement Error:** Set `reliability = 0.99` to ensure the true latent effect is equal to the observed effect. 
3. **Account for Design Structure:**
   * For **Within-Subjects**, you must mathematically adjust the Cohen's $d$ passed to `pwr` to account for the split-item variance scaling (the signal-to-noise ratio scales with the square root of $M$).
   * For **Between-Subjects**, you can match the `pwr` output directly by pushing $M$ to a very large number (e.g., $M = 200$), which entirely averages away the item variance and asymptotes perfectly at your target effect size.