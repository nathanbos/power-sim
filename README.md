# PIPS Power Simulation: Subjects × Items

An interactive, browser-based simulation illustrating how statistical power depends on the number of test subjects (N) and the number of test items or trials per subject (M).

Designed for researchers planning studies with repeated measures—such as cognitive tasks, usability testing, surveys, or psycholinguistic experiments—this tool demonstrates the non-linear trade-offs between recruiting more participants versus collecting more observations per participant.

**Run it online (no installation required):**
https://nathanbos.github.io/power-sim/

**Source code:**
https://github.com/nathanbos/power-sim

---

## Key Features & Methodological Realism

This simulator goes beyond basic signal averaging by introducing publication-grade psychometric realism and crossed-random-effects estimation without sacrificing real-time browser execution speed:

1. **Clark's $\min F'$ Estimator (Crossed Random Effects):** Standard aggregated tests implicitly treat test items as fixed rather than randomly sampled. This tool implements Herbert Clark's (1973) quasi-$F$ statistic ($\min F'$), delivering a fully crossed test that simultaneously bounds Type I error against both subject and item random variation—running 400x faster than iterative linear mixed models.
2. **Heterogeneous Treatment Slopes:** Interventions rarely affect every person or every test item uniformly. The simulator allows you to model real-world treatment variation across participants (`sd_subj_slope`) and items (`sd_item_slope`) to see how real-world heterogeneity degrades statistical power substantially faster than classical textbook formulas predict.
3. **Methodological Comparison Toggle:** Users can toggle between Clark's $\min F'$ and legacy aggregated tests side-by-side to directly observe **Clark's Fixed-Effect Fallacy**—demonstrating how aggregating over items artificially inflates statistical power when treatment-by-item interactions exist.

---

## Purpose & Key Takeaways

Classical power calculators (e.g., G*Power, the R `pwr` package) assume one independent observation per participant. In multi-item designs, ignoring the nested data structure wastes information, while treating repeated items as independent observations pseudoreplicates degrees of freedom and produces severely inflated false-positive rates.

By simulating the exact Data Generating Process (DGP), this tool highlights three critical methodological lessons for applied researchers:

* **Between-Subjects Designs hit an Asymptotic Ceiling:** Adding items helps average down measurement noise, but only up to a hard ceiling set by between-subject variability. You cannot power your way out of individual differences by simply adding more questions.
* **Within-Subjects Designs eliminate Subject Intercepts:** Because participant ability cancels out when evaluating condition differences, adding items increases power without an asymptotic ceiling. However, this assumes the treatment affects everyone equally. If individuals respond differently to the intervention (subject treatment heterogeneity), power will plateau, and more humans must be recruited.
* **Item Treatment Heterogeneity destroys Power:** If an intervention helps on certain tasks but not others ($\sigma^2_{\text{item\_slope}} > 0$), aggregating data masks this variance. Accounting for it properly via crossed random effects requires significantly larger item pools and sample sizes than standard tools suggest.

---

## Documentation & Technical Guide

For a deep dive into the underlying mathematics — including the $\min F'$ derivation and R code for extracting `sd_subj_slope` / `sd_item_slope` from your own pilot data using `lme4` — see [`simulation_description_gemini.md`](simulation_description_gemini.md).

---

## Running Locally

```r
# Install dependencies (if not already installed)
install.packages(c("shiny", "ggplot2", "DT"))

# Launch the app
shiny::runApp("path/to/power-sim")
```

Alternatively, source `simulate.R` directly to call `simulate_power()` from the R console without the Shiny interface.

---

## Files

| File | Purpose |
|------|---------|
| `simulate.R` | Core simulation function and command-line entry point |
| `app.R` | Shiny application (sources `simulate.R`) |
| `simulation_description_gemini.md` | Technical description of the estimation strategy and parameters |
| `docs/` | Static shinylive export served by GitHub Pages |
| `.github/workflows/pages.yml` | GitHub Actions workflow that deploys `docs/` to Pages on every push |

---

## Deployment

The app is deployed as a static site using [shinylive](https://posit-dev.github.io/r-shinylive/), which compiles the Shiny application to WebAssembly via [webR](https://docs.r-wasm.org/webr/latest/). The R runtime and all package dependencies run **entirely in the user's browser** — there is no server-side computation.

To regenerate the static site after modifying the app:

```r
# From a staging directory containing only app.R and simulate.R:
shinylive::export(".", "docs")
```

Then commit and push; the GitHub Actions workflow (`.github/workflows/pages.yml`) redeploys Pages automatically.

---

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
