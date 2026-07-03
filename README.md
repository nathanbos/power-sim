# Statistical Power Simulator: Subjects × Items

An interactive simulation illustrating how statistical power depends on the number of test subjects (N) and the number of test items or questions per subject (M). The tool is aimed at researchers designing studies that involve multiple items per participant — cognitive tests, surveys, rating scales, language experiments — where the choice of how many items to include is as consequential as the choice of how many participants to recruit.

**Run it online (no installation required):**
https://nathanbos.github.io/power-sim/

**Source code:**
https://github.com/nathanbos/power-sim

---

## Purpose

Classical power analysis (e.g., G\*Power, the R `pwr` package) treats a study as having one observation per participant. In many research designs, participants respond to multiple items, and each response is a separate (noisy) measurement of the same underlying construct. Ignoring this structure either wastes information or, if items are treated as independent observations, inflates degrees of freedom and produces anti-conservative results.

This simulator makes the multilevel structure explicit. It asks: given a true effect size and a known (or assumed) item reliability, how does power change as you add more participants versus more items? The key findings the simulator is designed to illustrate are:

- In a **between-subjects** design (condition A vs. condition B, different participants), adding more items helps because measurement noise averages down — but only up to a ceiling set by between-subject variability. No matter how many items you add, you cannot eliminate the noise from participant-to-participant differences.
- In a **within-subjects** design (every participant does both conditions), participant ability cancels out of the comparison entirely, so adding items improves power without an asymptotic ceiling. The price is that items must be *divided* between conditions (each participant sees each item only once, to avoid learning effects), so item difficulty also contributes noise.
- The tradeoff between N and M is non-linear and depends on the reliability of items, the design, and the effect size. This simulator makes those relationships visible.

---

## Running Locally

```r
# Install dependencies (if not already installed)
install.packages(c("shiny", "ggplot2", "DT"))

# Run the app
shiny::runApp("path/to/power-sim")
```

Alternatively, source `simulate.R` directly to call `simulate_power()` from the R console without the Shiny interface.

---

## Files

| File | Purpose |
|------|---------|
| `simulate.R` | Core simulation function and command-line entry point |
| `app.R` | Shiny application (sources `simulate.R`) |
| `docs/` | Static shinylive export served by GitHub Pages |
| `.github/workflows/pages.yml` | GitHub Actions workflow that deploys `docs/` to Pages on every push |

---

## Parameters

### Study Design
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `design` | `"between"` or `"within"` | `"between"` | **Between-subjects:** N is split evenly across two conditions; all subjects see all M items. **Within-subjects:** all N subjects experience both conditions; the M items are randomly divided, with each subject seeing M/2 items in condition A and M/2 in condition B. |

### Sample Sizes
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `N_values` | integer vector | | Total number of subjects. For between-subjects, N/2 are assigned to each condition. Sweep multiple values to produce the power table. |
| `M_values` | integer vector | | Number of items per subject. Within-subjects requires M ≥ 2. Between-subjects works with M ≥ 1 but items should be at least a few to see the averaging benefit. |

### Effect and Measurement
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `effect` | numeric | 0.5 | Target effect size in Cohen's *d* units. By construction, `effect` equals the *d* you would recover with a classical one-item design (between: M=1 two-sample t-test; within: M=4 paired t-test given equal subject and item variability). Conventional benchmarks: 0.2 = small, 0.5 = medium, 0.8 = large. |
| `reliability` | numeric (0–1) | 0.5 | Per-item reliability: the proportion of a single item's score variance that reflects true subject ability, as opposed to measurement noise. Equivalent to the intraclass correlation for a single item. Low reliability means items are noisy; high reliability means items are precise. Internally converted to a noise standard deviation via E = N\_weight × √((1 − r) / r). |

### Variability Weights
These parameters scale the spread of subject ability and item difficulty relative to each other. Both default to 1, which gives a symmetric model in which subjects and items are equally variable.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `N_weight` | numeric | 1.0 | Standard deviation of the subject ability distribution. Scales the per-item noise E proportionally (see Score Model below), so that `reliability` remains the signal-to-noise ratio regardless of the value chosen. For the between-subjects design, N\_weight cancels algebraically from Cohen's *d* and has no effect on power. For within-subjects, it matters only relative to M\_weight. |
| `M_weight` | numeric | 1.0 | Standard deviation of the item difficulty distribution. Has **no effect** on between-subjects power (items are shared across groups, so difficulty is a constant offset that cancels). For within-subjects designs, item difficulty varies between conditions (different items in each), so M\_weight contributes noise and does affect power. |

When N\_weight = M\_weight (the default), both weights cancel from all power formulas and can be safely left at 1. They are useful when you want to represent a design where items are substantially more or less variable in difficulty than participants are in ability.

### Simulation Control
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `C` | integer | 200 | Number of independent simulated experiments per (N, M) cell. Each simulated experiment draws a new random sample of subjects, items, and noise, runs the appropriate statistical test, and records whether p < α. Power is estimated as the proportion of significant results. Higher C gives more stable estimates (SE ≈ √(p(1−p)/C)); C = 200 gives ±3–4 percentage points at p = 0.5; C = 1000 gives ±1.5 percentage points. |
| `alpha` | numeric | 0.05 | Significance threshold. Fixed at 0.05 in the Shiny interface; adjustable when calling `simulate_power()` directly. |
| `seed` | integer or NULL | NULL | Optional random seed for reproducibility. NULL uses the current RNG state. |

---

## Score Model

Every simulated score follows the same additive model:

```
score_ijk = ability_i − difficulty_j + noise_ijk + raw_effect × (condition == B)
```

Where:
- `ability_i ~ N(0, N_weight)` — subject *i*'s latent ability
- `difficulty_j ~ N(0, M_weight)` — item *j*'s difficulty (higher = harder, lowers score)
- `noise_ijk ~ N(0, E)` — independent per-observation measurement error
- `E = N_weight × √((1 − reliability) / reliability)` — noise SD derived from reliability
- `raw_effect = effect × √(N_weight² + E²)` — the raw score-scale effect, scaled so that `effect` equals Cohen's *d* at the reference test

The reliability formula comes from the classical test theory definition: reliability = variance(true score) / variance(observed score) = N\_weight² / (N\_weight² + E²). Solving for E gives the formula above.

---

## Statistical Tests

### Between-Subjects Design

1. N subjects are divided alternately: subjects 1, 3, 5, … to condition A; subjects 2, 4, 6, … to condition B (N/2 per group).
2. All N subjects respond to the same M items.
3. Each subject's M item scores are averaged to produce one mean score per subject.
4. A two-sample t-test (equal variances assumed) is run on the N/2 group-A means vs. the N/2 group-B means.

**Why aggregate to subject means rather than use a mixed model?** Because items are shared across both groups, item difficulty is a constant offset that cancels exactly in the group comparison. Aggregating to subject means produces a valid test at df = N − 2, with power that increases monotonically with M (as noise averages down). A naïve mixed model with crossed random effects for subjects and items (lmer) was initially tried but produced Satterthwaite degrees of freedom ≈ M instead of ≈ N, causing power to *decrease* as M increased — a known pathology of between-subjects conditions in crossed designs.

**Effective Cohen's *d* for between-subjects:**
```
d(M) = effect × √(N_weight² + E²) / √(N_weight² + E²/M)
```
At M = 1, d = effect exactly. As M → ∞, d approaches effect / √(reliability), which is the ceiling — the effect size you could recover with perfectly reliable items. N\_weight and M\_weight cancel algebraically and do not affect between-subjects power.

**Correspondence to `pwr`:**
```r
pwr.t.test(n = N/2, d = effect, type = "two.sample")  # M=1 baseline
```

### Within-Subjects Design (Split Items)

1. All N subjects respond to both conditions.
2. For each subject, the M items are randomly assigned: M/2 to condition A, M/2 to condition B. No item is seen twice by the same subject.
3. For each subject *i*, a per-subject mean difference is computed: mean(scores_B_i) − mean(scores_A_i).
4. A one-sample t-test is run on the N differences (testing whether the mean difference is zero).

**Why split items rather than give all items in both conditions?** In most applied settings, presenting the same item to the same person in two different conditions is problematic: learning, order, and fatigue effects confound the condition comparison. Splitting items avoids this at the cost of introducing item-difficulty noise into the within-subject comparison (since different items appear in each condition for each subject).

**Why not use lmer here?** `lmer(score ~ condition + (1|subject))` is mathematically equivalent to the paired t-test on subject mean differences when condition is fully within-subject and items are nested within subject-condition. The paired t-test approach was used because it is 400× faster (avoids the iterative REML optimization in lme4) and produces identical results.

**Effective Cohen's *d* for within-subjects:**
```
d(M) = raw_effect / (2 × √((M_weight² + E²) / M))
```
When N\_weight = M\_weight (the default):
```
d(M) = effect × √M / 2
```
This increases with M without an asymptotic ceiling, and — notably — is **independent of reliability** when N\_weight = M\_weight. The reason: lower reliability increases the noise per item (E), but it also increases the raw\_effect proportionally (to maintain the same Cohen's *d*), and these two effects cancel exactly. Reliability only affects within-subjects power when items are more or less variable than subjects (M\_weight ≠ N\_weight).

**Correspondence to `pwr` (when N\_weight = M\_weight = 1, reliability = 0.5):**
```r
pwr.t.test(n = N, d = effect * sqrt(M) / 2, type = "paired")
```
For other weight/reliability combinations, the d formula above must be used.

---

## Algorithmic Description of `simulate.R`

`simulate_power()` runs a doubly-nested loop over M\_values × N\_values. For each (M, N) cell:

1. Draw N subject abilities from N(0, N\_weight) and M item difficulties from N(0, M\_weight). These are redrawn for every simulation repetition.
2. Compute per-observation scores using the additive model.
3. **Between-subjects:** Assign subjects to conditions, compute per-subject means using `tapply`, run `t.test(..., var.equal = TRUE)`.
4. **Within-subjects:** For each subject, call `sample.int(M, M/2)` to randomly select which items go to condition A; compute mean scores per condition per subject; run `t.test` on the N differences.
5. Record whether p < α. Repeat C times. Power for this cell = (number significant) / C.
6. Report progress to console (`cat`) with timing.

The function returns a matrix with dimnames `M=X` × `N=Y` containing power estimates as proportions.

`print_power_table()` formats this matrix with percentage signs and prints it to the console.

The `if (!interactive())` block at the bottom runs a default scenario when the script is executed via `Rscript simulate.R`, but is skipped when the file is sourced in an interactive R session or from `app.R`.

### `app.R` Extensions

`simulate_power_shiny()` wraps the core simulation loop with two additions:
- A progress callback (`progress_fn`) that calls Shiny's `incProgress()` after each cell, updating the progress bar in the UI.
- Per-run data capture: for the first min(C, 100) repetitions per cell, the per-subject condition means are stored (along with p-value and significance flag). After each cell's loop, `build_cell_plot()` converts this stored data into a `ggplot2` faceted plot (one facet per run) and the raw data is discarded. The function returns a list with `$power` (the matrix) and `$plots` (a named list of ggplot objects, one per cell).

`build_cell_plot()` creates a 10-column `facet_wrap` grid. Each facet shows:
- Two overlapping kernel density curves: condition A in blue, condition B in red (per-subject mean scores in each condition)
- Dashed vertical lines at the group means
- A text label in the top-right corner: p-value and observed Cohen's *d* for that simulated experiment
- Green panel background if the result was significant (p < 0.05); light gray if not

Colors are pre-computed as data frame columns and rendered via `scale_fill_identity()` / `scale_colour_identity()`, avoiding the need for the `ggnewscale` package (which is not available in the WebAssembly environment used by shinylive).

---

## Dependencies

### Runtime (simulation and Shiny app)
| Package | Version tested | Purpose |
|---------|---------------|---------|
| R | 4.5.3 | Language runtime |
| `shiny` | ≥ 1.7 | Web application framework |
| `ggplot2` | ≥ 3.4 | Faceted distribution plots |
| `DT` | ≥ 0.20 | Color-coded interactive power table |

`simulate.R` alone has **no package dependencies** — it uses only base R functions (`rnorm`, `t.test`, `tapply`, `sample.int`).

### Export Only (not needed to run the app)
| Package | Purpose |
|---------|---------|
| `shinylive` | Compiles the Shiny app to WebAssembly for static hosting |

### Packages Considered and Rejected
- **`lme4` / `lmerTest`:** Initially used to fit crossed random-effects models. Removed after discovering that Satterthwaite degrees of freedom for the between-subjects condition were ≈ M instead of ≈ N − 2, causing power to decrease as M increased. The aggregation approach is statistically equivalent and does not have this problem.
- **`ggnewscale`:** Initially used to render two independent `fill` scales in the faceted plot (one for panel background, one for density curves). Removed because it is not available in the webR/WebAssembly package repository. Replaced by pre-computing color columns and using `scale_fill_identity()`.

---

## Deployment

The app is deployed as a static site using [shinylive](https://posit-dev.github.io/r-shinylive/), which compiles the Shiny application to WebAssembly via [webR](https://docs.r-wasm.org/webr/latest/). The R runtime and all package dependencies run **entirely in the user's browser** — there is no server-side computation.

**Practical implications:**
- First load requires downloading approximately 86 MB of WebAssembly binaries (R runtime, packages, fonts). Subsequent visits use the browser cache and load in seconds.
- Computation speed depends on the user's device. The simulation is pure R loops; a 5×4 grid at C = 200 takes roughly 5–20 seconds depending on hardware.
- No data leaves the user's browser.

To regenerate the static site after modifying the app:

```r
# From the project directory in R:
shinylive::export(".", "docs")
```

Then commit and push; the GitHub Actions workflow (`.github/workflows/pages.yml`) redeploys Pages automatically.

---

## Limitations and Assumptions

- The score model is purely additive with Gaussian noise. It does not model ceiling/floor effects, non-normal ability distributions, guessing (for multiple-choice items), or order effects.
- The within-subjects design assumes items are randomly and independently assigned to conditions for each participant. Counterbalancing or Latin square designs are not modelled.
- The between-subjects design assumes equal group sizes (N/2 per condition). Unequal allocation is not currently supported.
- The significance threshold α is fixed at 0.05 in the Shiny interface. It can be changed by calling `simulate_power()` directly with an `alpha` argument.
- Power estimates have Monte Carlo variance proportional to 1/√C. At C = 200, estimates near 50% power have a standard error of approximately ±3.5 percentage points.
- The simulation does not currently model multiple comparisons, covariates, or designs with more than two conditions.
