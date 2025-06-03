library(here)

run_bayesian_tipping_analysis <- function(
    ped_ctrl_resp, ped_ctrl_total,
    ped_treat_resp, ped_treat_total,
    adult_ctrl_resp, adult_ctrl_total,
    adult_treat_resp, adult_treat_total,
    n_adult = 494,
    p_adult = 0.091,
    sigma_values = seq(2, 1, by = -0.1)) {
  
  # Combine data
  pediatric_clean <- data.frame(
    Group = c("Placebo", "DrugA"),
    n = c(ped_ctrl_resp, ped_treat_resp),
    N = c(ped_ctrl_total, ped_treat_total),
    Trial = "Pooled"
  ) |> mutate(Source = "Pediatric")
  
  adult_clean <- data.frame(
    Group = c("Placebo", "DrugA"),
    n = c(adult_ctrl_resp, adult_treat_resp),
    N = c(adult_ctrl_total, adult_treat_total),
    Trial = "Pooled"
  ) |> mutate(Source = "Adult")
  
  combined_data <- bind_rows(pediatric_clean, adult_clean) |> 
    mutate(
      group = ifelse(Group == "DrugA", 1, 0),
      trial_str = paste(Source, Trial),
      trial = as.integer(factor(trial_str)),
      y = n
    )
  
  info_per_subject <- n_adult * p_adult * (1 - p_adult)
  
  n_ped <- sum(combined_data$N[combined_data$group == 0 & grepl("Pediatric", combined_data$trial_str)])
  n_ped_trt <- sum(combined_data$N[combined_data$group == 1 & grepl("Pediatric", combined_data$trial_str)])
  
  tipping_results <- data.frame(
    fixed_sigma_alpha = numeric(),
    OR_median = numeric(),
    OR_lower = numeric(),
    OR_upper = numeric(),
    ESS_prior = numeric(),
    ESS_FDA = numeric(),
    Borrowed_FDA = numeric(),
    ESS_FDA_trt = numeric(),
    Borrowed_FDA_trt = numeric(),
    V1 = numeric(),
    V2 = numeric(),
    V1_trt = numeric(),
    V2_trt = numeric()
  )
  
  ### 加入进度条
  withProgress(message = 'Running Bayesian Hierarchical Model...', value = 0, {
    
    total_steps <- length(sigma_values)
    
    for (i in seq_along(sigma_values)) {
      s <- sigma_values[i]
      
      # 每一轮更新进度
      incProgress(1 / total_steps)
      
      jags_data <- list(
        y = combined_data$y,
        n = combined_data$N,
        group = combined_data$group,
        trial = combined_data$trial,
        N = nrow(combined_data),
        J = length(unique(combined_data$trial)),
        mu_alpha_mean = 0,
        mu_alpha_prec = 0.01,
        fixed_sigma_alpha = s
      )
      
      params <- c("mu_alpha", "alpha", "beta")
      model_file <- here::here("R", "models", "Bayesian_hierarchical_model", "rjags", "baseline_pooled_tipping.txt")
      
      jags_mod <- rjags::jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000,
                                    inits = list(.RNG.name = "base::Super-Duper", .RNG.seed = 123))
      update(jags_mod, 1000)
      samples <- rjags::coda.samples(jags_mod, variable.names = params, n.iter = 50000)
      post <- as.matrix(samples)
      
      beta_draws <- post[, grep("^beta\\[", colnames(post))]
      alpha_draws <- post[, grep("^alpha\\[", colnames(post))]
      or_draws <- exp(beta_draws)
      or_ped <- or_draws[, "beta[2]"]
      
      alpha2_samples <- alpha_draws[, "alpha[2]"]
      alpha1_samples <- alpha_draws[, "alpha[1]"]
      var_post <- var(alpha2_samples)
      var_post_trt <- var(alpha1_samples)
      ess_prior <- (1 / var_post) / info_per_subject
      
      # === Vague model ===
      jags_data_vague <- jags_data
      jags_data_vague$fixed_sigma_alpha <- 1e32
      
      jags_mod_vague <- rjags::jags.model(model_file, data = jags_data_vague, n.chains = 3, n.adapt = 1000,
                                          inits = list(.RNG.name = "base::Super-Duper", .RNG.seed = 123))
      update(jags_mod_vague, 1000)
      samples_vague <- rjags::coda.samples(jags_mod_vague, variable.names = params, n.iter = 5000)
      post_vague <- as.matrix(samples_vague)
      
      alpha2_vague <- post_vague[, "alpha[2]"]
      alpha1_vague <- post_vague[, "alpha[1]"]
      V1 <- var(alpha2_vague)
      V2 <- var_post
      V1_trt <- var(alpha1_vague)
      V2_trt <- var_post_trt
      
      ess_fda <- n_ped * V1 / V2
      ess_borrowed <- ess_fda - n_ped
      ess_fda_trt <- n_ped_trt * V1_trt / V2_trt
      ess_borrowed_trt <- ess_fda_trt - n_ped_trt
      
      tipping_results <- rbind(tipping_results, data.frame(
        fixed_sigma_alpha = s,
        OR_median = median(or_ped),
        OR_lower = quantile(or_ped, 0.025),
        OR_upper = quantile(or_ped, 0.975),
        ESS_prior = ess_prior,
        ESS_FDA = ess_fda,
        Borrowed_FDA = ess_borrowed,
        ESS_FDA_trt = ess_fda_trt,
        Borrowed_FDA_trt = ess_borrowed_trt,
        V1 = V1,
        V2 = V2,
        V1_trt = V1_trt,
        V2_trt = V2_trt
      ))
    }
  })  # <- 进度条括号结束
  
  return(tipping_results)
}



plot_hierarchical_tipping <- function(tipping_results, tipping_row_fda) {
  
  p <- ggplot(tipping_results, aes(
    x = fixed_sigma_alpha,
    y = log(OR_median),
    color = Significant_FDA,
    text = paste0(
      "Sigma: ", round(fixed_sigma_alpha, 2), "<br>",
      "OR: ", round(OR_median, 2), " [", round(OR_lower, 2), ", ", round(OR_upper, 2), "]<br>",
      "Placebo ESS: ", round(Borrowed_FDA, 1), "<br>",
      "Treatment ESS: ", round(Borrowed_FDA_trt, 1)
    )
  )) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = log(OR_lower), ymax = log(OR_upper)), width = 0.05) +  # 同步 log
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +                # log(1) = 0
    {if(nrow(tipping_row_fda) > 0) 
      geom_vline(xintercept = tipping_row_fda$fixed_sigma_alpha, linetype = "dotted", color = "red")} +
    {if(nrow(tipping_row_fda) > 0)
      annotate("text",
               x = tipping_row_fda$fixed_sigma_alpha * 0.6,
               y = log(tipping_row_fda$OR_median) + 0.5,
               label = paste0("Tipping Point = ", round(tipping_row_fda$fixed_sigma_alpha, 2)),
               hjust = 0, vjust = 0, size = 4.5, fontface = "italic", color = "red")} +
    scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray"))+
    labs(
      title = "Treatment Effect with Hierarchical Model",
      x = "Sigma Alpha",
      y = "log(Odds Ratio)",
      color = "CI Exclude 0" 
    ) +
    theme_minimal(base_size = 15)
  
  
  return(p)
}
