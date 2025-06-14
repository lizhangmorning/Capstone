# mixture_analysis.R
# Consolidated functions for mixture prior analysis and plotting

library(ggplot2)
library(plotly)
library(dplyr)
library(scales)

# --- Internal Helper Function ---
compute_C <- function(alpha, beta, y, n) {
  lchoose(n, y) + lbeta(alpha + y, beta + n - y) - lbeta(alpha, beta)
}

# --- Main Analysis Function ---
mixture_analysis <- function(child_data, adult_data, alpha_min = 0.001, alpha_max = 0.01, alpha_steps = 50) {
  
  flat_prior <- list(a = 1, b = 1)
  
  # Generate weight sequence
  weights <- seq(alpha_min, alpha_max, length.out = alpha_steps)
  
  results <- data.frame(
    weight = weights,
    median_or = NA,
    lower_95_or = NA,
    upper_95_or = NA,
    ess_total_or = NA,
    ess_treat_or = NA,
    ess_control_or = NA,
    median_rr = NA,
    lower_95_rr = NA,
    upper_95_rr = NA,
    ess_total_rr = NA,
    ess_treat_rr = NA,
    ess_control_rr = NA,
    k_treat = NA,
    k_control = NA
  )
  
  set.seed(123)
  
  # No-borrowing posteriors
  p1_noborrow <- rbeta(1e5, flat_prior$a + child_data$treat$y,
                       flat_prior$b + child_data$treat$n - child_data$treat$y)
  p2_noborrow <- rbeta(1e5, flat_prior$a + child_data$control$y,
                       flat_prior$b + child_data$control$n - child_data$control$y)
  or_noborrow <- (p1_noborrow / (1 - p1_noborrow)) / (p2_noborrow / (1 - p2_noborrow))
  rr_noborrow <- p1_noborrow - p2_noborrow
  log_or_var_noborrow <- var(log(or_noborrow))
  var_rr_noborrow <- var(rr_noborrow)
  var_p1_noborrow <- var(log(p1_noborrow / (1 - p1_noborrow)))
  var_p2_noborrow <- var(log(p2_noborrow / (1 - p2_noborrow)))
  
  step_summary <- data.frame()  # New: to store step-by-step summary
  
  for (i in seq_along(weights)) {
    w <- weights[i]
    
    # Mixture weights (k values)
    C_treat <- list(
      info = exp(compute_C(adult_data$treat$y + 1, adult_data$treat$n - adult_data$treat$y + 1,
                           child_data$treat$y, child_data$treat$n)),
      flat = exp(compute_C(1, 1, child_data$treat$y, child_data$treat$n))
    )
    k_treat <- (w * C_treat$info) / (w * C_treat$info + (1 - w) * C_treat$flat)
    
    C_control <- list(
      info = exp(compute_C(adult_data$control$y + 1, adult_data$control$n - adult_data$control$y + 1,
                           child_data$control$y, child_data$control$n)),
      flat = exp(compute_C(1, 1, child_data$control$y, child_data$control$n))
    )
    k_control <- (w * C_control$info) / (w * C_control$info + (1 - w) * C_control$flat)
    
    post_p1 <- ifelse(
      runif(1e5) < k_treat,
      rbeta(1e5, adult_data$treat$y + 1 + child_data$treat$y,
            adult_data$treat$n - adult_data$treat$y + 1 + child_data$treat$n - child_data$treat$y),
      rbeta(1e5, 1 + child_data$treat$y,
            1 + child_data$treat$n - child_data$treat$y)
    )
    
    post_p2 <- ifelse(
      runif(1e5) < k_control,
      rbeta(1e5, adult_data$control$y + 1 + child_data$control$y,
            adult_data$control$n - adult_data$control$y + 1 + child_data$control$n - child_data$control$y),
      rbeta(1e5, 1 + child_data$control$y,
            1 + child_data$control$n - child_data$control$y)
    )
    
    or_post <- (post_p1 / (1 - post_p1)) / (post_p2 / (1 - post_p2))
    rr_post <- post_p1 - post_p2
    
    log_or_var <- var(log(or_post))
    var_rr <- var(rr_post)
    
    ess_total_or <- (log_or_var_noborrow / log_or_var - 1) * (child_data$treat$n + child_data$control$n)
    ess_treat_or <- (var_p1_noborrow / var(log(post_p1 / (1 - post_p1))) - 1) * child_data$treat$n
    ess_control_or <- (var_p2_noborrow / var(log(post_p2 / (1 - post_p2))) - 1) * child_data$control$n
    
    ess_total_rr <- (var_rr_noborrow / var_rr - 1) * (child_data$treat$n + child_data$control$n)
    ess_treat_rr <- (var(p1_noborrow) / var(post_p1) - 1) * child_data$treat$n
    ess_control_rr <- (var(p2_noborrow) / var(post_p2) - 1) * child_data$control$n
    
    # Store results
    results[i, ] <- list(
      w,
      median(or_post),
      quantile(or_post, 0.025),
      quantile(or_post, 0.975),
      ess_total_or,
      ess_treat_or,
      ess_control_or,
      median(rr_post),
      quantile(rr_post, 0.025),
      quantile(rr_post, 0.975),
      ess_total_rr,
      ess_treat_rr,
      ess_control_rr,
      k_treat,
      k_control
    )
    
    # New: store step summary
    step_summary <- rbind(step_summary, data.frame(
      weight = w,
      median_rr = median(rr_post),
      lower_95 = quantile(rr_post, 0.025),
      upper_95 = quantile(rr_post, 0.975),
      Favorable = quantile(rr_post, 0.025) > 0,
      delta_rr = median(rr_post),
      ess_treat_rr = ess_treat_rr,
      ess_control_rr = ess_control_rr,
      ess_total_rr = ess_total_rr,
      k_treat = k_treat,
      k_control = k_control
    ))
  }
  
  # Find tipping points
  rr_tip_idx <- which(results$lower_95_rr > 0)[1]
  
  tipping_point_summary <- list(
    RR = if (!is.na(rr_tip_idx)) {
      rr_row <- results[rr_tip_idx, ]
      list(
        weight = rr_row$weight,
        ess_treat = rr_row$ess_treat_rr,
        ess_control = rr_row$ess_control_rr,
        ess_total = rr_row$ess_total_rr,
        lower_95 = rr_row$lower_95_rr,
        upper_95 = rr_row$upper_95_rr
      )
    } else {
      list(weight = NA_real_, ess_treat = NA_real_, ess_control = NA_real_, ess_total = NA_real_, 
           lower_95 = NA_real_, upper_95 = NA_real_)
    }
  )
  
  return(list(
    results = results,
    ess_results = results,  # for compatibility
    step_summary = step_summary,  # ✅ key addition
    tipping_point = tipping_point_summary$RR$weight,
    tipping_point_summary = tipping_point_summary
  ))
}


# Plotting function with step size control
plot_mixture_results <- function(analysis_results, alpha_min, alpha_max, alpha_steps) {
  results <- analysis_results$results
  tipping_point <- analysis_results$tipping_point
  
  plot_data <- results %>%
    mutate(
      Favorable = lower_95_rr > 0,
      text = paste0(
        "α = ", round(weight, 4),
        "<br>ESS Trt: ", round(ess_treat_rr, 1),
        "<br>ESS Ctl: ", round(ess_control_rr, 1),
        "<br>95% CI: [", round(lower_95_rr, 3), ", ", round(upper_95_rr, 3), "]"
      )
    )
  
  error_bar_width <- (alpha_max - alpha_min) / alpha_steps
  
  p <- ggplot(plot_data, aes(
    x = weight, 
    y = median_rr,
    color = Favorable,
    text = text
  )) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lower_95_rr, ymax = upper_95_rr), width = error_bar_width) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray")) +
    labs(
      title = "Treatment Effect with Mixture Prior",
      x = "Prior Weight (α)",
      y = "Difference in Response Rates"
    ) +
    theme_minimal(base_size = 15) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  if (!is.na(tipping_point)) {
    tipping_row <- plot_data[which.min(abs(plot_data$weight - tipping_point)), ]
    p <- p +
      geom_vline(xintercept = tipping_point, linetype = "dotted", color = "red") +
      annotate("text",
               x = tipping_point * 0.6,
               y = tipping_row$median_rr + max(plot_data$median_rr) * 0.1,
               label = paste0("Tipping Point = ", round(tipping_point, 4)),
               hjust = 0, vjust = 0, size = 4.5, fontface = "italic", color = "red")
  }
  
  ggplotly(p, tooltip = "text")
}