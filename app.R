library(shiny)
library(ggplot2)
library(DT)

source("simulate.R")

# ── Helpers ───────────────────────────────────────────────────────────────────

parse_int_vector <- function(x, min_val = 1L, max_val = 500L) {
  parts <- trimws(unlist(strsplit(x, "[,\\s]+")))
  parts <- parts[nchar(parts) > 0]
  if (length(parts) == 0) return(NULL)
  vals <- suppressWarnings(as.integer(parts))
  if (any(is.na(vals)) || any(vals < min_val) || any(vals > max_val)) return(NULL)
  sort(unique(vals))
}

safe_id <- function(x) gsub("[^a-zA-Z0-9]", "_", x)

# ── Per-cell distribution plot ─────────────────────────────────────────────────
# run_data: list of length <= 100, each element:
#   list(scores_a, scores_b, p, sig)
# scores_a/b are per-subject mean scores in each condition

build_cell_plot <- function(run_data, N, M, alpha, design) {
  n_runs <- length(run_data)
  
  score_dfs <- vector("list", n_runs)
  summ_rows <- vector("list", n_runs)
  
  for (r in seq_len(n_runs)) {
    rd <- run_data[[r]]
    na <- length(rd$scores_a)
    nb <- length(rd$scores_b)
    
    pool_sd <- sqrt(
      ((na - 1) * var(rd$scores_a) + (nb - 1) * var(rd$scores_b)) / (na + nb - 2)
    )
    obs_d <- (mean(rd$scores_b) - mean(rd$scores_a)) / max(pool_sd, 1e-8)
    
    score_dfs[[r]] <- data.frame(
      run       = r,
      condition = c(rep("A", na), rep("B", nb)),
      score     = c(rd$scores_a, rd$scores_b)
    )
    summ_rows[[r]] <- data.frame(
      run     = r,
      sig     = rd$sig,
      mean_a  = mean(rd$scores_a),
      mean_b  = mean(rd$scores_b),
      obs_d   = obs_d,
      p_label = if (is.na(rd$p)) "p=NA\nd=NA" else sprintf("p=%.3f\nd=%.2f", rd$p, obs_d)
    )
  }
  
  df   <- do.call(rbind, score_dfs)
  summ <- do.call(rbind, summ_rows)
  
  # Pre-compute colors as data columns so a single identity scale covers everything,
  # avoiding the need for ggnewscale (which may not be available in webR/shinylive).
  summ$bg_fill    <- ifelse(summ$sig, "#A5D6A7", "#E8E8E8")
  df$cond_fill    <- ifelse(df$condition == "A", "#42A5F5", "#EF5350")
  df$cond_colour  <- ifelse(df$condition == "A", "#1565C0", "#B71C1C")
  
  cond_label <- if (design == "between") "Group" else "Condition"
  title_str  <- sprintf(
    "N = %d, M = %d  —  %d simulation runs  (green panel = p < %.2f  |  blue = %s A, red = %s B)",
    N, M, n_runs, alpha, cond_label, cond_label
  )
  
  ggplot() +
    geom_rect(
      data = summ,
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
          fill = bg_fill, group = run),
      alpha = 0.3
    ) +
    geom_density(
      data = df,
      aes(x = score, fill = cond_fill, colour = cond_colour,
          group = interaction(run, condition)),
      alpha = 0.35, linewidth = 0.25
    ) +
    scale_fill_identity() +
    scale_colour_identity() +
    geom_vline(
      data = summ,
      aes(xintercept = mean_a),
      colour = "#1565C0", linetype = "dashed", linewidth = 0.3
    ) +
    geom_vline(
      data = summ,
      aes(xintercept = mean_b),
      colour = "#B71C1C", linetype = "dashed", linewidth = 0.3
    ) +
    geom_text(
      data    = summ,
      aes(label = p_label),
      x = Inf, y = Inf, hjust = 1.08, vjust = 1.5,
      size = 1.7, colour = "gray25"
    ) +
    facet_wrap(~ run, ncol = 10) +
    labs(title = title_str, x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      strip.text      = element_blank(),
      panel.spacing   = unit(0.12, "cm"),
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "gray75", fill = NA, linewidth = 0.3),
      legend.position = "top",
      plot.title      = element_text(size = 10, hjust = 0.5)
    )
}

# ── Simulation wrapper: runs simulation + builds plots ─────────────────────────

simulate_power_shiny <- function(
    N_values, M_values,
    N_weight      = 1.0,
    M_weight      = 1.0,
    sd_subj_slope = 0.0,
    sd_item_slope = 0.0,
    reliability   = 0.5,
    effect        = 0.5,
    C             = 200L,
    alpha         = 0.05,
    seed          = NULL,
    design        = "between",
    method        = "min_F",
    progress_fn   = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  reliability <- min(reliability, 0.999)
  E          <- N_weight * sqrt((1 - reliability) / reliability)
  # raw_effect <- effect * sqrt(N_weight^2 + E^2)
  raw_effect <- effect * N_weight
  n_M     <- length(M_values)
  n_N     <- length(N_values)
  n_cells <- n_M * n_N
  cell    <- 0L
  n_store <- min(C, 100L)
  
  power <- matrix(
    NA_real_,
    nrow     = n_M,
    ncol     = n_N,
    dimnames = list(paste0("M=", M_values), paste0("N=", N_values))
  )
  plots <- list()
  
  for (mi in seq_along(M_values)) {
    M <- M_values[mi]
    
    for (ni in seq_along(N_values)) {
      N    <- N_values[ni]
      cell <- cell + 1L
      n_sig    <- 0L
      run_data <- vector("list", n_store)
      
      for (rep in seq_len(C)) {
        subj_ability <- rnorm(N, sd = N_weight)
        item_diff    <- rnorm(M, sd = M_weight)
        subj_slope   <- if (sd_subj_slope > 0) rnorm(N, sd = sd_subj_slope) else numeric(N)
        item_slope   <- if (sd_item_slope > 0) rnorm(M, sd = sd_item_slope) else numeric(M)
        
        if (design == "between") {
          cond_subj <- rep_len(0:1, N)
          s_grid    <- rep(seq_len(N), each = M)
          i_grid    <- rep(seq_len(M), times = N)
          c_grid    <- cond_subj[s_grid]
          
          score <- (subj_ability[s_grid]
                    - item_diff[i_grid]
                    + rnorm(N * M, sd = E)
                    + c_grid * (raw_effect + subj_slope[s_grid] + item_slope[i_grid]))
          
          subj_means <- tapply(score, s_grid, mean)
          scores_a   <- as.numeric(subj_means[cond_subj == 0])
          scores_b   <- as.numeric(subj_means[cond_subj == 1])
          
          if (method == "t_test" || M < 2) {
            p <- t.test(scores_b, scores_a, var.equal = TRUE)$p.value
          } else {
            t1 <- t.test(scores_b, scores_a, var.equal = TRUE)
            F1 <- unname(t1$statistic^2)
            df1 <- unname(t1$parameter)
            
            item_means_A <- tapply(score[c_grid == 0], i_grid[c_grid == 0], mean)
            item_means_B <- tapply(score[c_grid == 1], i_grid[c_grid == 1], mean)
            t2 <- t.test(item_means_B, item_means_A, paired = TRUE)
            F2 <- unname(t2$statistic^2)
            df2 <- unname(t2$parameter)
            
            min_F <- (F1 * F2) / (F1 + F2)
            df_den <- ((F1 + F2)^2) / ((F1^2 / df2) + (F2^2 / df1))
            p <- 1 - pf(min_F, df1 = 1, df2 = df_den)
          }
          
        } else {
          M_a      <- M %/% 2L
          M_b      <- M - M_a
          scores_a <- numeric(N)
          scores_b <- numeric(N)
          score_list_A <- vector("list", M)
          score_list_B <- vector("list", M)
          subj_diffs   <- numeric(N)
          
          for (i in seq_len(N)) {
            a_idx      <- sample.int(M, M_a)
            b_idx      <- seq_len(M)[-a_idx]
            sa         <- subj_ability[i] - item_diff[a_idx] + rnorm(M_a, sd = E)
            sb         <- subj_ability[i] - item_diff[b_idx] + rnorm(M_b, sd = E) + 
              (raw_effect + subj_slope[i] + item_slope[b_idx])
            
            scores_a[i]   <- mean(sa)
            scores_b[i]   <- mean(sb)
            subj_diffs[i] <- scores_b[i] - scores_a[i]
            
            for (k in seq_along(a_idx)) score_list_A[[a_idx[k]]] <- c(score_list_A[[a_idx[k]]], sa[k])
            for (k in seq_along(b_idx)) score_list_B[[b_idx[k]]] <- c(score_list_B[[b_idx[k]]], sb[k])
          }
          
          item_means_A <- vapply(score_list_A, function(v) if (length(v)) mean(v) else NA_real_, numeric(1))
          item_means_B <- vapply(score_list_B, function(v) if (length(v)) mean(v) else NA_real_, numeric(1))
          both_seen    <- !is.na(item_means_A) & !is.na(item_means_B)

          if (method == "t_test" || M < 2 || sum(both_seen) < 2) {
            # Too few items had data in both conditions this rep (can happen when N is
            # small relative to M, since item-condition assignment is per-subject and
            # random) — fall back to the subject-level test rather than crash.
            p <- t.test(subj_diffs)$p.value
          } else {
            t1 <- t.test(subj_diffs)
            F1 <- unname(t1$statistic^2)
            df1 <- unname(t1$parameter)

            t2 <- t.test(item_means_B[both_seen], item_means_A[both_seen], paired = TRUE)
            F2 <- unname(t2$statistic^2)
            df2 <- unname(t2$parameter)

            min_F <- (F1 * F2) / (F1 + F2)
            df_den <- ((F1 + F2)^2) / ((F1^2 / df2) + (F2^2 / df1))
            p <- 1 - pf(min_F, df1 = 1, df2 = df_den)
          }
        }
        
        sig <- !is.na(p) && p < alpha
        if (sig) n_sig <- n_sig + 1L
        
        if (rep <= n_store) {
          run_data[[rep]] <- list(
            scores_a = scores_a,
            scores_b = scores_b,
            p        = p,
            sig      = sig
          )
        }
      }
      
      power[mi, ni] <- n_sig / C
      
      cell_name         <- paste0("N=", N, ", M=", M)
      plots[[cell_name]] <- build_cell_plot(run_data, N, M, alpha, design)
      
      if (!is.null(progress_fn)) progress_fn(cell, n_cells, N, M, power[mi, ni])
    }
  }
  
  list(power = power, plots = plots)
}

# ── UI ─────────────────────────────────────────────────────────────────────────

ui <- fluidPage(
  tags$head(tags$style(HTML("
    .sidebar-section    { margin-top: 14px; }
    .help-block         { font-size: 11.5px; color: #666; margin-top: 2px; }
    .section-label      { font-weight: 600; font-size: 13px; margin-bottom: 4px; }
    .error-box          { color: #c0392b; background: #fdf0ef;
                          border: 1px solid #e8b4b0; border-radius: 4px;
                          padding: 8px 12px; margin-bottom: 10px; }
  "))),
  
  titlePanel("Statistical Power Simulator"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      div(class = "section-label", "Study Design"),
      radioButtons(
        "design", label = NULL,
        choices  = c("Between-subjects" = "between",
                     "Within-subjects"  = "within"),
        selected = "within"
      ),
      helpText(
        "Between: subjects split across conditions; all see the same items.",
        "Within: all subjects do both conditions; items are randomly divided between conditions."
      ),
      
      div(class = "section-label", style = "margin-top: 10px;", "Statistical Method"),
      radioButtons(
        "method", label = NULL,
        choices  = c("Clark's min F' (Crossed)" = "min_F",
                     "Aggregated t-test (Legacy)" = "t_test"),
        selected = "min_F"
      ),
      helpText(
        "min F': Robust fully crossed test accounting for subject & item variance.",
        "Legacy t-test: Aggregates over items (vulnerable to fixed-effect fallacies)."
      ),
      
      div(class = "sidebar-section", div(class = "section-label", "Sample Sizes")),
      textInput("N_values", "Subject counts (N)", value = "4, 8, 16, 32, 64, 128"),
      helpText("Comma-separated integers (2–500). At most 8 values."),

      textInput("M_values", "Item counts (M)", value = "4, 8, 16, 32, 64, 128"),
      helpText("Comma-separated integers. At most 8 values.",
               "Within-subjects requires M ≥ 2."),
      
      div(class = "sidebar-section", div(class = "section-label", "Effect & Measurement")),
      numericInput("effect", "Effect size (Cohen’s d)",
                   value = 0.5, min = 0.01, max = 3, step = 0.05),
      helpText("0.2 = small · 0.5 = medium · 0.8 = large"),
      
      sliderInput("reliability", "Per-item reliability",
                  min = 0.10, max = 0.99, value = 0.50, step = 0.01),
      helpText(
        "Proportion of score variance that reflects true ability.",
        "Low = noisy items; high = precise items. Typical range: 0.3–0.9."
      ),
      
      div(class = "sidebar-section", div(class = "section-label", "Variability Weights")),
      numericInput("N_weight", "Subject baseline SD",
                   value = 1, min = 0.1, max = 10, step = 0.1),
      helpText("SD of baseline subject ability."),
      
      numericInput("M_weight", "Item baseline SD",
                   value = 1, min = 0.1, max = 10, step = 0.1),
      helpText("SD of baseline item difficulty."),
      
      numericInput("sd_subj_slope", "Subject slope SD (Treatment heterogeneity)",
                   value = 0.25, min = 0, max = 5, step = 0.1),
      helpText("0 = uniform intervention effect across subjects; >0 = variable responsiveness."),
      
      numericInput("sd_item_slope", "Item slope SD (Treatment heterogeneity)",
                   value = 0.0, min = 0, max = 5, step = 0.1),
      helpText("0 = uniform intervention effect across items; >0 = intervention affects certain items more."),
      
      div(class = "sidebar-section", div(class = "section-label", "Simulation")),
      textInput("C", "Simulations per cell", value = "200"),
      helpText("Whole number between 10 and 10 000.",
               "More = stabler estimates; up to 100 runs are shown in the drill-down panels."),
      
      numericInput("seed", "Random seed (optional)",
                   value = NA, min = 1, max = 2147483647, step = 1),
      helpText("Leave blank for a different result each run.",
               "Enter any positive integer for a reproducible result."),
      
      br(),
      actionButton("run_btn", "Run Simulation",
                   class = "btn-primary btn-block", width = "100%")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "outer_tabs",
        
        tabPanel(
          "Power Table",
          br(),
          uiOutput("error_box"),
          uiOutput("table_header"),
          DT::DTOutput("power_table")
        ),
        
        tabPanel(
          "Simulation Runs",
          br(),
          uiOutput("run_tabs_ui")
        )
      )
    )
  )
)

# ── Server ─────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {
  
  sim_results <- reactiveVal(NULL)
  
  observeEvent(input$run_btn, {
    
    # ── Input validation ──────────────────────────────────────────────────────
    N_values <- parse_int_vector(input$N_values, min_val = 2L,  max_val = 500L)
    M_values <- parse_int_vector(input$M_values, min_val = 1L,  max_val = 200L)
    C_val    <- suppressWarnings(as.integer(trimws(input$C)))
    
    errs <- character(0)
    
    if (is.null(N_values) || length(N_values) == 0)
      errs <- c(errs, "Subject counts: enter comma-separated integers between 2 and 500.")
    else if (length(N_values) > 8)
      errs <- c(errs, "Subject counts: maximum 8 values.")
    else if (input$design == "between" && any(N_values < 4))
      errs <- c(errs, "Subject counts: between-subjects design needs N ≥ 4 (at least 2 per group).")
    
    if (is.null(M_values) || length(M_values) == 0)
      errs <- c(errs, "Item counts: enter comma-separated integers between 1 and 200.")
    else if (length(M_values) > 8)
      errs <- c(errs, "Item counts: maximum 8 values.")
    else if (input$design == "within" && any(M_values < 2))
      errs <- c(errs, "Item counts: within-subjects design needs M ≥ 2 (items are split between conditions).")
    
    if (is.na(input$effect) || input$effect < 0.01 || input$effect > 3)
      errs <- c(errs, "Effect size must be a number between 0.01 and 3.0.")
    
    if (is.na(input$N_weight) || input$N_weight < 0.1 || input$N_weight > 10)
      errs <- c(errs, "Subject variability weight must be between 0.1 and 10.")
    
    if (is.na(input$M_weight) || input$M_weight < 0.1 || input$M_weight > 10)
      errs <- c(errs, "Item variability weight must be between 0.1 and 10.")
    
    if (is.na(input$sd_subj_slope) || input$sd_subj_slope < 0 || input$sd_subj_slope > 5)
      errs <- c(errs, "Subject slope SD must be between 0 and 5.")
    
    if (is.na(input$sd_item_slope) || input$sd_item_slope < 0 || input$sd_item_slope > 5)
      errs <- c(errs, "Item slope SD must be between 0 and 5.")
    
    if (is.na(C_val) || C_val < 10L || C_val > 10000L)
      errs <- c(errs, "Simulations per cell: enter a whole number between 10 and 10 000.")
    
    output$error_box <- renderUI({
      if (length(errs) == 0) return(NULL)
      div(class = "error-box", tags$ul(lapply(errs, tags$li)))
    })
    
    if (length(errs) > 0) return()
    
    seed_val <- if (is.na(input$seed)) NULL else as.integer(input$seed)
    
    # ── Run simulation ────────────────────────────────────────────────────────
    result <- withProgress(message = "Running simulation…", value = 0, {
      simulate_power_shiny(
        N_values      = N_values,
        M_values      = M_values,
        N_weight      = input$N_weight,
        M_weight      = input$M_weight,
        sd_subj_slope = input$sd_subj_slope,
        sd_item_slope = input$sd_item_slope,
        reliability   = input$reliability,
        effect        = input$effect,
        C             = C_val,
        alpha         = 0.05,
        seed          = seed_val,
        design        = input$design,
        method        = input$method,
        progress_fn   = function(cell, n_cells, N, M, pwr) {
          incProgress(
            1 / n_cells,
            detail = sprintf("N=%d, M=%d → %.1f%%", N, M, pwr * 100)
          )
        }
      )
    })
    
    sim_results(result)
  })
  
  # ── Power table ───────────────────────────────────────────────────────────
  output$table_header <- renderUI({
    req(sim_results())
    p(style = "color:#555; font-size:13px;",
      "Proportion of simulated experiments reaching p < 0.05.",
      " Rows = items (M), columns = subjects (N).")
  })
  
  output$power_table <- DT::renderDT({
    result <- sim_results()
    req(result)
    
    mat  <- result$power
    df   <- as.data.frame(mat)
    df   <- cbind(Items = rownames(mat), df)
    rownames(df) <- NULL
    n_cols <- names(df)[-1]
    
    DT::datatable(
      df,
      rownames = FALSE,
      options  = list(dom = "t", paging = FALSE, ordering = FALSE),
      class    = "compact stripe hover"
    ) |>
      DT::formatPercentage(n_cols, digits = 1) |>
      DT::formatStyle(
        columns         = n_cols,
        backgroundColor = DT::styleInterval(
          c(0.2, 0.4, 0.6, 0.8),
          c("#FFCDD2", "#FFE0B2", "#FFF9C4", "#C8E6C9", "#A5D6A7")
        ),
        fontWeight = "bold"
      )
  })
  
  # ── Simulation-run drill-down tabs ────────────────────────────────────────
  output$run_tabs_ui <- renderUI({
    result <- sim_results()
    req(result)
    
    tab_list <- lapply(names(result$plots), function(nm) {
      pid <- paste0("plt_", safe_id(nm))
      tabPanel(nm, plotOutput(pid, height = "900px"))
    })
    
    do.call(tabsetPanel, c(tab_list, list(id = "run_tabs")))
  })
  
  observe({
    result <- sim_results()
    req(result)
    
    for (nm in names(result$plots)) {
      local({
        n   <- nm
        pid <- paste0("plt_", safe_id(n))
        output[[pid]] <- renderPlot(result$plots[[n]], res = 120)
      })
    }
  })
}

shinyApp(ui, server)