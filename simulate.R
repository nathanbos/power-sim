# =============================================================================
# Statistical Power Simulator: Subjects × Items (v2.0 - Methodological Upgrade)
#
# NEW FEATURES IN v2.0:
# 1. Heterogeneous Treatment Slopes: Models real-world variance where items or 
#    subjects respond differently to the experimental intervention.
# 2. Clark's min F' Estimator: Provides a publication-grade, fully crossed 
#    random-effects test that accounts for both subject and item variance 
#    simultaneously without the browser-freezing latency of lmer().
# 3. Methodological Comparison: Allows toggling between 'min_F' and naive 
#    't_test' to illustrate Clark's Fixed-Effect Fallacy and inflated Type I error.
# =============================================================================

simulate_power <- function(
    N_values,              # integer vector — numbers of subjects to sweep
    M_values,              # integer vector — numbers of items to sweep
    N_weight      = 1.0,   # SD of baseline subject ability (random intercept u_0i)
    M_weight      = 1.0,   # SD of baseline item difficulty (random intercept w_0j)
    sd_subj_slope = 0.0,   # SD of subject treatment responsiveness (random slope u_1i)
    sd_item_slope = 0.0,   # SD of item treatment responsiveness (random slope w_1j)
    reliability   = 0.5,   # per-item reliability: signal/(signal+noise), in (0, 1)
    effect        = 0.5,   # target baseline Cohen's d
    C             = 200L,  # simulation repetitions per cell
    alpha         = 0.05,  # significance threshold
    seed          = NULL,  # optional RNG seed
    design        = "between", # "between" or "within"
    method        = "min_F"    # "min_F" (Clark's min F') or "t_test" (legacy aggregated)
) {
  if (!is.null(seed)) set.seed(seed)
  if (!method %in% c("min_F", "t_test")) stop("method must be 'min_F' or 't_test'")
  
  reliability <- min(reliability, 0.999)
  E          <- N_weight * sqrt((1 - reliability) / reliability)
  raw_effect <- effect * sqrt(N_weight^2 + E^2)
  
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
        # 1. Draw Random Intercepts & Slopes from Data Generating Process (DGP)
        subj_ability <- rnorm(N, sd = N_weight)
        item_diff    <- rnorm(M, sd = M_weight)
        subj_slope   <- if (sd_subj_slope > 0) rnorm(N, sd = sd_subj_slope) else numeric(N)
        item_slope   <- if (sd_item_slope > 0) rnorm(M, sd = sd_item_slope) else numeric(M)
        
        if (design == "between") {
          # Alternate subjects into Condition A (0) and Condition B (1)
          cond_subj <- rep_len(0:1, N)
          
          # Vectorized grid setup: N subjects x M items
          s_grid <- rep(seq_len(N), each = M)
          i_grid <- rep(seq_len(M), times = N)
          c_grid <- cond_subj[s_grid]
          
          # Score Model with Heterogeneous Slopes
          score <- subj_ability[s_grid] - item_diff[i_grid] + rnorm(N * M, sd = E) +
            c_grid * (raw_effect + subj_slope[s_grid] + item_slope[i_grid])
          
          if (method == "t_test" || M < 2) {
            # Legacy By-Subject Aggregating t-test
            subj_means <- tapply(score, s_grid, mean)
            p <- t.test(subj_means[cond_subj == 1], subj_means[cond_subj == 0], var.equal = TRUE)$p.value
          } else {
            # Clark's min F' (Crossed Random Effects)
            # F1: By-Subject analysis (averaging over items)
            subj_means <- tapply(score, s_grid, mean)
            t1 <- t.test(subj_means[cond_subj == 1], subj_means[cond_subj == 0], var.equal = TRUE)
            F1 <- unname(t1$statistic^2)
            df1 <- unname(t1$parameter)
            
            # F2: By-Item analysis (averaging over subjects within each condition)
            # Since all items appear in both conditions across different subjects, compare item means
            item_means_A <- tapply(score[c_grid == 0], i_grid[c_grid == 0], mean)
            item_means_B <- tapply(score[c_grid == 1], i_grid[c_grid == 1], mean)
            t2 <- t.test(item_means_B, item_means_A, paired = TRUE)
            F2 <- unname(t2$statistic^2)
            df2 <- unname(t2$parameter)
            
            # Compute min F' and resulting p-value
            min_F <- (F1 * F2) / (F1 + F2)
            df_den <- ((F1 + F2)^2) / ((F1^2 / df2) + (F2^2 / df1))
            p <- 1 - pf(min_F, df1 = 1, df2 = df_den)
          }
          
        } else {
          # WITHIN-SUBJECTS DESIGN (Split Items)
          M_a <- M %/% 2L
          M_b <- M - M_a
          
          # Track scores and item assignments to compute valid by-item statistics
          score_list_A <- vector("list", M)
          score_list_B <- vector("list", M)
          subj_diffs   <- numeric(N)
          
          for (i in seq_len(N)) {
            a_idx <- sample.int(M, M_a)
            b_idx <- seq_len(M)[-a_idx]
            
            score_a <- subj_ability[i] - item_diff[a_idx] + rnorm(M_a, sd = E)
            # Apply intervention effect + subject slope + item slope for Condition B
            score_b <- subj_ability[i] - item_diff[b_idx] + rnorm(M_b, sd = E) + 
              (raw_effect + subj_slope[i] + item_slope[b_idx])
            
            subj_diffs[i] <- mean(score_b) - mean(score_a)
            
            for (k in seq_along(a_idx)) score_list_A[[a_idx[k]]] <- c(score_list_A[[a_idx[k]]], score_a[k])
            for (k in seq_along(b_idx)) score_list_B[[b_idx[k]]] <- c(score_list_B[[b_idx[k]]], score_b[k])
          }
          
          if (method == "t_test" || M < 2) {
            p <- t.test(subj_diffs)$p.value
          } else {
            # Clark's min F' for Split-Item Within-Subjects
            # F1: By-Subject difference test (H0: mean diff == 0)
            t1 <- t.test(subj_diffs)
            F1 <- unname(t1$statistic^2)
            df1 <- unname(t1$parameter)
            
            # F2: By-Item analysis across conditions
            item_means_A <- sapply(score_list_A, mean)
            item_means_B <- sapply(score_list_B, mean)
            t2 <- t.test(item_means_B, item_means_A, paired = TRUE)
            F2 <- unname(t2$statistic^2)
            df2 <- unname(t2$parameter)
            
            min_F <- (F1 * F2) / (F1 + F2)
            df_den <- ((F1 + F2)^2) / ((F1^2 / df2) + (F2^2 / df1))
            p <- 1 - pf(min_F, df1 = 1, df2 = df_den)
          }
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
# Command-Line Execution / Local Demonstration
# =============================================================================
if (!interactive()) {
  
  cat("=== RUN 1: Baseline Additive DGP with Clark's min F' ===\n")
  params_base <- list(
    N_values      = c(12, 24, 48),
    M_values      = c(4, 10, 20),
    N_weight      = 1.0,
    M_weight      = 1.0,
    sd_subj_slope = 0.0,
    sd_item_slope = 0.0,
    reliability   = 0.5,
    effect        = 0.5,
    C             = 200L,
    design        = "between",
    method        = "min_F",
    seed          = 42
  )
  res_base <- do.call(simulate_power, params_base)
  print_power_table(res_base)
  
  cat("\n=== RUN 2: Heterogeneous Slopes (Item Variance Trap) ===\n")
  cat("Notice how power drops when items respond differently to treatment!\n")
  params_het <- params_base
  params_het$sd_item_slope <- 0.4  # Add item treatment heterogeneity
  res_het <- do.call(simulate_power, params_het)
  print_power_table(res_het)
  
}