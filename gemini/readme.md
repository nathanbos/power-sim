# Statistical Power Simulator: Subjects × Items (Version 2.0)

An interactive, browser-based simulation illustrating how statistical power depends on the number of test subjects ($N$) and the number of test items or trials per subject ($M$). 

Designed for researchers planning studies with repeated measures—such as cognitive tasks, usability testing, surveys, or psycholinguistic experiments—this tool demonstrates the non-linear trade-offs between recruiting more participants versus collecting more observations per participant.

**Run it online (no installation required):**  
https://nathanbos.github.io/power-sim/

**Source code:**  
https://github.com/nathanbos/power-sim

---

## What’s New in Version 2.0

While Version 1 focused on basic signal averaging under an additive model, **Version 2.0** introduces publication-grade psychometric realism and crossed-random-effects estimation without sacrificing real-time browser execution speed:

1. **Clark’s $\min F'$ Estimator (Crossed Random Effects):** Version 1 relied on subject-mean aggregation and standard $t$-tests, which implicitly treat test items as fixed rather than sampled. Version 2.0 implements Herbert Clark’s (1973) quasi-$F$ statistic ($\min F'$), delivering a fully crossed test that simultaneously bounds Type I error against both subject and item random variation—running $400\times$ faster than iterative linear mixed models (`lmer`).
2. **Heterogeneous Treatment Slopes:** Interventions rarely affect every person or every test item uniformly. Version 2.0 allows you to simulate treatment variation across participants (`sd_subj_slope`) and items (`sd_item_slope`). This demonstrates how real-world item heterogeneity degrades statistical power substantially faster than classical textbook formulas predict.
3. **Methodological Comparison Toggle:** Users can now toggle between Clark's $\min F'$ and legacy aggregated $t$-tests side-by-side to directly observe **Clark's Fixed-Effect Fallacy**—showing how aggregating over items artificially inflates statistical power when treatment-by-item interactions exist.

---

## Purpose & Key Takeaways

Classical power calculators (e.g., G*Power, the R `pwr` package) assume one independent observation per participant. In multi-item designs, ignoring the nested data structure wastes information, while treating repeated items as independent observations pseudoreplicates degrees of freedom and produces severely inflated false-positive rates.

By simulating the exact Data Generating Process (DGP), this tool highlights three critical methodological lessons:

* **Between-Subjects Designs hit an Asymptotic Ceiling:** Adding items helps average down observation noise, but only up to a hard ceiling set by between-subject variability. You cannot power your way out of individual differences by simply adding more questions.
* **Within-Subjects Designs eliminate Subject Intercepts:** Because participant ability cancels out when taking condition differences, adding items increases power without an asymptotic ceiling. However, splitting items across conditions introduces item difficulty noise into the comparison.
* **Item Treatment Heterogeneity destroys Power:** If an intervention helps on certain questions but not others ($\sigma^2_{\text{item\_slope}} > 0$), aggregating data masks this variance. Accounting for it properly via crossed random effects requires significantly larger sample sizes than standard power tools suggest.

---

## Running Locally

To run the interactive Shiny application on your local machine:

```r
# Install dependencies (if not already installed)
install.packages(c("shiny", "ggplot2", "DT"))

# Launch the app
shiny::runApp("path/to/power-sim")