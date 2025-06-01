# ===== Power Prior Analysis Functions =====
# 这个文件包含了Power Prior分析的封装函数，可以直接整合到主Shiny应用中

library(rjags)
library(coda)
library(dplyr)
library(ggplot2)
library(scales)

# PowerPrior核心函数
power_prior_rjags <- function(ped_data, adult_data, alpha, n_iter = 10000) {
  # 定义JAGS模型
  model_string <- "
  model {
    # Likelihood for pediatric data
    y_ped_placebo ~ dbin(p_placebo, n_ped_placebo)
    y_ped_drug ~ dbin(p_drug, n_ped_drug)
    
    # Power Prior (weighted adult data)
    a_placebo <- 1 + alpha * y_adult_placebo
    b_placebo <- 1 + alpha * (n_adult_placebo - y_adult_placebo)
    
    a_drug <- 1 + alpha * y_adult_drug
    b_drug <- 1 + alpha * (n_adult_drug - y_adult_drug)
    
    # Prior distribution
    p_placebo ~ dbeta(a_placebo, b_placebo)
    p_drug ~ dbeta(a_drug, b_drug)
    
    # Response rate difference
    rr_diff <- p_drug - p_placebo
  }
  "
  
  # 准备JAGS数据
  jags_data <- list(
    y_ped_placebo = ped_data$y_placebo,
    n_ped_placebo = ped_data$n_placebo,
    y_ped_drug = ped_data$y_drug,
    n_ped_drug = ped_data$n_drug,
    y_adult_placebo = adult_data$y_placebo,
    n_adult_placebo = adult_data$n_placebo,
    y_adult_drug = adult_data$y_drug,
    n_adult_drug = adult_data$n_drug,
    alpha = alpha
  )
  
  # 初始化JAGS模型
  jags_model <- jags.model(textConnection(model_string), 
                           data = jags_data, 
                           n.chains = 3, 
                           quiet = TRUE)
  
  # 预热
  update(jags_model, 5000, progress.bar = "none")
  
  # 抽取后验样本
  samples <- coda.samples(jags_model, 
                          variable.names = c("p_placebo", "p_drug", "rr_diff"), 
                          n.iter = n_iter, 
                          progress.bar = "none")
  
  # 转换为矩阵处理
  samples_matrix <- as.matrix(samples)
  
  rr_diff <- samples_matrix[, "rr_diff"]
  p_placebo <- samples_matrix[, "p_placebo"]  
  p_drug <- samples_matrix[, "p_drug"] 
  median_rr_diff <- median(rr_diff)
  ci_95 <- quantile(rr_diff, c(0.025, 0.975))
  prob_gt_0 <- mean(rr_diff > 0)
  
  return(list(
    median_rr_diff = median_rr_diff,
    ci_95 = ci_95,
    prob_gt_0 = prob_gt_0,
    samples = rr_diff,
    samples_placebo = p_placebo,
    samples_drug = p_drug
  ))
}

# 计算基于方差的ESS
calculate_ess_from_variance_rjags <- function(ped_data, adult_data, alpha_seq, n_samples = 10000) {
  ess_results <- data.frame(alpha = alpha_seq)
  
  # 计算无借用时的后验方差
  no_borrow_res <- power_prior_rjags(ped_data, adult_data, 0)
  
  V1_placebo <- var(no_borrow_res$samples_placebo)
  V1_drug <- var(no_borrow_res$samples_drug)
  V1_diff <- var(no_borrow_res$samples)  
  
  ess_values <- numeric(length(alpha_seq))
  borrowed_values <- numeric(length(alpha_seq))
  variance_ratio_values <- numeric(length(alpha_seq))
  
  placebo_ess_values <- numeric(length(alpha_seq))
  drug_ess_values <- numeric(length(alpha_seq))
  placebo_borrowed_values <- numeric(length(alpha_seq))
  drug_borrowed_values <- numeric(length(alpha_seq))
  placebo_variance_ratio <- numeric(length(alpha_seq))
  drug_variance_ratio <- numeric(length(alpha_seq))
  
  for (i in 1:length(alpha_seq)) {
    alpha <- alpha_seq[i]
    
    # 运行RJAGS模型
    borrow_res <- power_prior_rjags(ped_data, adult_data, alpha)
    
    V2_placebo <- var(borrow_res$samples_placebo)
    V2_drug <- var(borrow_res$samples_drug)
    V2_diff <- var(borrow_res$samples) 
    
    variance_ratio_placebo <- V1_placebo / V2_placebo
    variance_ratio_drug <- V1_drug / V2_drug
    variance_ratio_diff <- V1_diff / V2_diff
    
    # 计算ESS
    n_placebo <- ped_data$n_placebo
    n_drug <- ped_data$n_drug
    n_total <- n_placebo + n_drug
    
    # 计算借用的样本量
    placebo_ess <- n_placebo * (variance_ratio_placebo - 1)
    drug_ess <- n_drug * (variance_ratio_drug - 1)
    total_ess <- n_total * (variance_ratio_diff - 1)
    
    placebo_borrowed <- placebo_ess
    drug_borrowed <- drug_ess
    total_borrowed <- total_ess
    
    ess_values[i] <- total_ess
    borrowed_values[i] <- total_borrowed
    variance_ratio_values[i] <- variance_ratio_diff
    
    placebo_ess_values[i] <- placebo_ess
    drug_ess_values[i] <- drug_ess
    placebo_borrowed_values[i] <- placebo_borrowed
    drug_borrowed_values[i] <- drug_borrowed
    placebo_variance_ratio[i] <- variance_ratio_placebo
    drug_variance_ratio[i] <- variance_ratio_drug
  }
  
  ess_results$ess_total <- ess_values
  ess_results$borrowed <- borrowed_values
  ess_results$variance_ratio <- variance_ratio_values
  ess_results$variance_noborrow <- V1_diff
  ess_results$variance_borrow <- numeric(length(alpha_seq))  
  
  ess_results$placebo_ess <- placebo_ess_values
  ess_results$placebo_borrowed <- placebo_borrowed_values
  ess_results$placebo_variance_ratio <- placebo_variance_ratio
  ess_results$placebo_variance_noborrow <- V1_placebo
  
  ess_results$drug_ess <- drug_ess_values
  ess_results$drug_borrowed <- drug_borrowed_values
  ess_results$drug_variance_ratio <- drug_variance_ratio
  ess_results$drug_variance_noborrow <- V1_drug
  
  for (i in 1:length(alpha_seq)) {
    alpha <- alpha_seq[i]
    borrow_res <- power_prior_rjags(ped_data, adult_data, alpha)
    ess_results$variance_borrow[i] <- var(borrow_res$samples)
  }
  
  return(ess_results)
}

# 查找临界点
find_tipping_point <- function(results) {
  tp_idx <- which(results$lower_ci > 0)[1]
  tipping_point <- ifelse(is.na(tp_idx), NA, results$alpha[tp_idx])
  return(tipping_point)
}

# PowerPrior临界点分析主函数
run_power_prior_analysis <- function(child_data, adult_data, alpha_min = 0.001, alpha_max = 0.01, steps = 50) {
  # 转换数据格式以匹配原始函数
  ped_data <- list(
    y_placebo = child_data$control$y,
    n_placebo = child_data$control$n,
    y_drug = child_data$treat$y,
    n_drug = child_data$treat$n
  )
  
  adult_data_formatted <- list(
    y_placebo = adult_data$control$y,
    n_placebo = adult_data$control$n,
    y_drug = adult_data$treat$y,
    n_drug = adult_data$treat$n
  )
  
  alpha_seq <- seq(alpha_min, alpha_max, length.out = steps)
  results <- data.frame()
  
  for (i in 1:length(alpha_seq)) {
    alpha <- alpha_seq[i]
    rjags_res <- power_prior_rjags(ped_data, adult_data_formatted, alpha)
    
    results <- rbind(results, data.frame(
      alpha = alpha,
      median_rr_diff = rjags_res$median_rr_diff,
      lower_ci = rjags_res$ci_95[1],
      upper_ci = rjags_res$ci_95[2],
      prob_gt_0 = rjags_res$prob_gt_0
    ))
  }
  
  tipping_point <- find_tipping_point(results)
  ess_results <- calculate_ess_from_variance_rjags(ped_data, adult_data_formatted, alpha_seq)
  
  return(list(
    results = results,
    ess_results = ess_results,
    tipping_point = tipping_point
  ))
}

# 生成Power Prior Tipping Point图表 - 返回plotly对象，与Mixture Prior保持一致
run_power_prior_plot <- function(child_data, adult_data) {
  library(ggplot2)
  library(plotly)
  
  # 运行分析
  analysis_results <- run_power_prior_analysis(child_data, adult_data)
  
  results <- analysis_results$results
  ess_results <- analysis_results$ess_results
  tipping_point <- analysis_results$tipping_point
  
  # 合并结果和ESS数据
  combined_results <- merge(results, ess_results, by = "alpha")
  
  # 创建与Mixture Prior相似的ggplot对象
  p <- ggplot(combined_results, aes(x = alpha, y = median_rr_diff,
                                    text = paste0("α=", round(alpha, 4),
                                                  "<br>ESS Trt: ", round(drug_ess, 1),
                                                  "<br>ESS Ctl: ", round(placebo_ess, 1),
                                                  "<br>95% CI: [", round(lower_ci, 3), ", ", round(upper_ci, 3), "]"))) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.0002) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_vline(xintercept = ifelse(!is.na(tipping_point), tipping_point, NA), 
               linetype = "dashed", color = "blue") +
    labs(title = "Treatment Effect with Power Prior",
         x = "Prior Weight (α)",
         y = "Treatment Effect (Drug - Placebo)") +
    theme_minimal(base_size = 14) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
  
  # 返回plotly对象，与Mixture Prior保持一致
  ggplotly(p, tooltip = "text")
}

# Power Prior分析函数 - 返回摘要信息
power_prior_analysis <- function(child_data, adult_data) {
  # 运行分析
  analysis_results <- run_power_prior_analysis(child_data, adult_data)
  
  tipping_point <- analysis_results$tipping_point
  ess_results <- analysis_results$ess_results
  
  # 构建tipping point摘要
  if (!is.na(tipping_point)) {
    # 找到tipping point对应的ESS
    tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
    
    if (length(tp_idx) > 0) {
      tipping_point_summary <- list(
        RR = list(
          weight = tipping_point,
          ess_treat = ess_results$drug_ess[tp_idx],
          ess_control = ess_results$placebo_ess[tp_idx],
          ess_total = ess_results$ess_total[tp_idx]
        )
      )
    } else {
      tipping_point_summary <- list(RR = list(weight = NA, ess_treat = NA, ess_control = NA, ess_total = NA))
    }
  } else {
    tipping_point_summary <- list(RR = list(weight = NA, ess_treat = NA, ess_control = NA, ess_total = NA))
  }
  
  return(list(
    tipping_point_summary = tipping_point_summary,
    full_results = analysis_results
  ))
}