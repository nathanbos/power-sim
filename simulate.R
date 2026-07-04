# -----------------------------------------------------------------------------
# simulate_power()
#
# Returns a matrix of observed power (% significant) where
#   rows  = item counts (M_values)
#   cols  = subject counts (N_values)
#
# Score model:  score = subject_ability - item_difficulty + noise + effect_raw*(cond==B)
#   subject_ability ~ N(0, N_weight)
#   item_difficulty ~ N(0, M_weight)
#   noise           ~ N(0, E),  E = N_weight * sqrt((1 - reliability) / reliability)
#   effect_raw      = effect * sqrt(N_weight² + E²)
#
# design = "between"
#   N/2 subjects in condition A, N/2 in B.  All subjects see the same M items.
#   Item difficulty cancels in the A-vs-B comparison (shared items); N_weight
#   also cancels from Cohen's d, so only reliability and M affect power.
#   At M=1: d = effect exactly.  Reference: pwr.t.test(n=N/2, d=effect, type="two.sample")
#
# design = "within"
#   All N subjects do both conditions.  The M items are randomly split: each
#   subject sees M/2 items in A and M/2 in B (no item seen twice, avoiding
#   learning effects).  Subject ability cancels in B-A; item difficulty does NOT
#   (different items in each condition), so M_weight affects noise.
#   SD of per-subject diffs = 2 * sqrt((M_weight² + E²) / M)
#   When N_weight = M_weight: d(M) = effect*sqrt(M)/2  (independent of reliability)
#   Reference: pwr.t.test(n=N, d=raw_effect/(2*sqrt((M_weight²+E²)/M)), type="paired")
# -----------------------------------------------------------------------------
simulate_power <- function(
  N_values,              # integer vector — numbers of subjects to sweep
  M_values,              # integer vector — numbers of items to sweep
  N_weight    = 1.0,     # SD of subject ability; scales noise E proportionally
  M_weight    = 1.0,     # SD of item difficulty; affects within-subjects noise only
  reliability = 0.5,     # per-item reliability: signal/(signal+noise), in (0, 1)
  effect      = 0.5,     # Cohen's d (see design notes above)
  C           = 200L,    # simulation repetitions per cell
  alpha       = 0.05,    # significance threshold
  seed        = NULL,    # optional RNG seed; NULL uses current RNG state
  design      = "between" # "between" or "within"
) {
  if (!is.null(seed)) set.seed(seed)

  reliability <- min(reliability, 0.999)
  E          <- N_weight * sqrt((1 - reliability) / reliability)
  # raw_effect <- effect * sqrt(N_weight^2 + E^2)   # = effect * N_weight / sqrt(reliability)
  raw_effect <- effect * N_weight
  power <- matrix(
    NA_real_,
    nrow     = length(M_values),
    ncol     = length(N_values),
    dimnames = list(paste0("M=", M_values), paste0("N=", N_values))
  )

  n_cells <- length(N_values) * length(M_values)
  cell    <- 0L

  for (mi in seq_along(M_values)) {
    M <- M_values[mi]

    for (ni in seq_along(N_values)) {
      N    <- N_values[ni]
      cell <- cell + 1L
      t0   <- proc.time()[["elapsed"]]
      n_sig <- 0L

      for (rep in seq_len(C)) {
        subj_ability <- rnorm(N, sd = N_weight)
        item_diff    <- rnorm(M, sd = M_weight)

        if (design == "between") {
          cond_subj <- rep_len(0:1, N)
          subj_idx  <- rep(seq_len(N), each = M)
          item_idx  <- rep(seq_len(M), times = N)

          score <- (subj_ability[subj_idx]
                    - item_diff[item_idx]
                    + rnorm(N * M, sd = E)
                    + cond_subj[subj_idx] * raw_effect)

          subj_means <- tapply(score, subj_idx, mean)
          p <- t.test(subj_means[cond_subj == 1],
                      subj_means[cond_subj == 0],
                      var.equal = TRUE)$p.value

        } else {
          M_a <- M %/% 2L
          M_b <- M - M_a

          diffs <- numeric(N)
          for (i in seq_len(N)) {
            a_idx    <- sample.int(M, M_a)
            b_idx    <- seq_len(M)[-a_idx]
            score_a  <- subj_ability[i] - item_diff[a_idx] + rnorm(M_a, sd = E)
            score_b  <- subj_ability[i] - item_diff[b_idx] + rnorm(M_b, sd = E) + raw_effect
            diffs[i] <- mean(score_b) - mean(score_a)
          }
          p <- t.test(diffs)$p.value
        }

        if (!is.na(p) && p < alpha) n_sig <- n_sig + 1L
      }

      power[mi, ni] <- n_sig / C
      elapsed <- proc.time()[["elapsed"]] - t0
      cat(sprintf("[%d/%d]  N=%3d  M=%3d  power = %.3f  (%.1f s)\n",
                  cell, n_cells, N, M, power[mi, ni], elapsed))
    }
  }
  power
}


# -----------------------------------------------------------------------------
# print_power_table()  —  pretty-print with % formatting
# -----------------------------------------------------------------------------
print_power_table <- function(mat) {
  pct <- apply(mat, c(1, 2), function(x) sprintf("%.1f%%", x * 100))
  cat("\nPower table  (rows = items M, cols = subjects N)\n")
  cat(strrep("-", 6 + ncol(pct) * 9), "\n")
  print(noquote(pct))
  invisible(mat)
}


# =============================================================================
# Run  (skipped when sourced interactively; executes via Rscript on command line)
# =============================================================================
if (!interactive()) {

params <- list(
  N_values    = c(10, 20, 50, 100, 200),
  M_values    = c(2, 5, 10, 20),
  N_weight    = 1.0,
  M_weight    = 1.0,
  reliability = 0.5,
  effect      = 0.5,
  C           = 200L,
  design      = "within"
)

cat(sprintf("Design: %s\n", params$design))
cat(sprintf("Parameters: N_weight=%.2f  M_weight=%.2f  reliability=%.2f  effect=%.2f  C=%d\n\n",
            params$N_weight, params$M_weight, params$reliability, params$effect, params$C))

result <- do.call(simulate_power, params)
print_power_table(result)

} # end if (!interactive())
