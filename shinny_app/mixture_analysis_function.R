# --- Internal Function ---
compute_C <- function(alpha, beta, y, n) {
  lchoose(n, y) + lbeta(alpha + y, beta + n - y) - lbeta(alpha, beta)
}

# --- Main Analysis Function ---
mixture_analysis <- function(
    data_child,
    data_adult,
    prior_flat = list(a = 1, b = 1),
    weight_grid = seq(0, 1, 0.04),
    seed = 123,
    nsim = 1e5
) {
  set.seed(seed)
  
  results <- data.frame(
    weight = weight_grid,
    median_rr = NA, lower_95_rr = NA, upper_95_rr = NA,
    ess_total_rr = NA, ess_treat_rr = NA, ess_control_rr = NA
  )
  
  p1_noborrow <- rbeta(nsim, prior_flat$a + data_child$treat$y,
                       prior_flat$b + data_child$treat$n - data_child$treat$y)
  p2_noborrow <- rbeta(nsim, prior_flat$a + data_child$control$y,
                       prior_flat$b + data_child$control$n - data_child$control$y)
  var_p1_nob <- var(p1_noborrow)
  var_p2_nob <- var(p2_noborrow)
  
  prior_params <- list(
    treat = if (data_adult$treat$n > 0) {
      list(a = data_adult$treat$y + 1, b = data_adult$treat$n - data_adult$treat$y + 1)
    } else NULL,
    control = if (data_adult$control$n > 0) {
      list(a = data_adult$control$y + 1, b = data_adult$control$n - data_adult$control$y + 1)
    } else NULL
  )
  
  for (i in seq_along(weight_grid)) {
    w <- weight_grid[i]
    
    k_treat <- if (!is.null(prior_params$treat)) {
      C_info <- compute_C(prior_params$treat$a, prior_params$treat$b, data_child$treat$y, data_child$treat$n)
      C_flat <- compute_C(prior_flat$a, prior_flat$b, data_child$treat$y, data_child$treat$n)
      (w * C_info) / (w * C_info + (1 - w) * C_flat)
    } else 0
    
    k_control <- if (!is.null(prior_params$control)) {
      C_info <- compute_C(prior_params$control$a, prior_params$control$b, data_child$control$y, data_child$control$n)
      C_flat <- compute_C(prior_flat$a, prior_flat$b, data_child$control$y, data_child$control$n)
      (w * C_info) / (w * C_info + (1 - w) * C_flat)
    } else 0
    
    post_p1 <- ifelse(
      runif(nsim) < k_treat,
      rbeta(nsim, prior_params$treat$a + data_child$treat$y,
            prior_params$treat$b + data_child$treat$n - data_child$treat$y),
      rbeta(nsim, prior_flat$a + data_child$treat$y,
            prior_flat$b + data_child$treat$n - data_child$treat$y)
    )
    
    post_p2 <- ifelse(
      runif(nsim) < k_control,
      rbeta(nsim, prior_params$control$a + data_child$control$y,
            prior_params$control$b + data_child$control$n - data_child$control$y),
      rbeta(nsim, prior_flat$a + data_child$control$y,
            prior_flat$b + data_child$control$n - data_child$control$y)
    )
    
    rr_post <- post_p1 - post_p2
    
    ess_treat <- if (!is.null(prior_params$treat)) (var_p1_nob / var(post_p1) - 1) * data_child$treat$n else 0
    ess_control <- if (!is.null(prior_params$control)) (var_p2_nob / var(post_p2) - 1) * data_child$control$n else 0
    
    results[i, ] <- c(w, median(rr_post), quantile(rr_post, 0.025), quantile(rr_post, 0.975),
                      ess_treat + ess_control, ess_treat, ess_control)
  }
  
  tipping_idx <- which(results$lower_95_rr > 0)[1]
  
  if (is.na(tipping_idx)) {
    return(list(
      tipping_point = NA,
      ess_treat = NA,
      ess_control = NA,
      results_df = results
    ))
  }
  
  return(list(
    tipping_point = results$weight[tipping_idx],
    ess_treat = results$ess_treat_rr[tipping_idx],
    ess_control = results$ess_control_rr[tipping_idx],
    results_df = results
  ))
}

# --- Plotting Function ---
run_mixture_plot <- function(data_child, data_adult) {
  library(ggplot2)
  library(plotly)
  
  analysis <- mixture_analysis(data_child, data_adult)
  results <- analysis$results_df
  k_tp <- analysis$tipping_point
  
  p <- ggplot(results, aes(x = weight, y = median_rr,
                           text = paste0("k=", round(weight, 2),
                                         "<br>ESS Trt: ", round(ess_treat_rr, 1),
                                         "<br>ESS Ctl: ", round(ess_control_rr, 1),
                                         "<br>95% CI: [", round(lower_95_rr, 3), ", ", round(upper_95_rr, 3), "]"))) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_95_rr, ymax = upper_95_rr), width = 0.01) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = k_tp, linetype = "dashed", color = "blue") +
    labs(title = "Treatment Effect with Mixture Prior",
         x = "Prior Weight (k)",
         y = "Treatment Effect (Drug - Placebo)") +
    theme_minimal(base_size = 14)
  
  ggplotly(p, tooltip = "text")
}
