library(shiny)
library(rjags)
library(here)
library(ggplot2)
library(coda)
library(dplyr)
library(scales)

# PowerPrior函数封装
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
run_power_prior_analysis <- function(ped_data, adult_data, alpha_min, alpha_max, steps) {
  alpha_seq <- seq(alpha_min, alpha_max, length.out = steps)
  results <- data.frame()
  
  withProgress(message = 'Running analysis...', value = 0, {
    for (i in 1:length(alpha_seq)) {
      alpha <- alpha_seq[i]
      rjags_res <- power_prior_rjags(ped_data, adult_data, alpha)
      
      results <- rbind(results, data.frame(
        alpha = alpha,
        median_rr_diff = rjags_res$median_rr_diff,
        lower_ci = rjags_res$ci_95[1],
        upper_ci = rjags_res$ci_95[2],
        prob_gt_0 = rjags_res$prob_gt_0
      ))
      
      incProgress(1/length(alpha_seq), detail = paste("Alpha:", round(alpha, 4)))
    }
  })
  
  tipping_point <- find_tipping_point(results)
  ess_results <- calculate_ess_from_variance_rjags(ped_data, adult_data, alpha_seq)
  
  return(list(
    results = results,
    ess_results = ess_results,
    tipping_point = tipping_point
  ))
}

# UI定义
ui <- fluidPage(
  titlePanel("Pediatric Extrapolation Analytics"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Pediatric Data",
                 numericInput("ped_ctrl_resp", "Placebo Responders", 1),
                 numericInput("ped_ctrl_total", "Placebo Total", 6),
                 numericInput("ped_treat_resp", "Treatment Responders", 15),
                 numericInput("ped_treat_total", "Treatment Total", 30)
        ),
        tabPanel("Adult Borrowing Data",
                 numericInput("adult_ctrl_resp", "Placebo Responders", 45),
                 numericInput("adult_ctrl_total", "Placebo Total", 494),
                 numericInput("adult_treat_resp", "Treatment Responders", 375),
                 numericInput("adult_treat_total", "Treatment Total", 509)
        ),
        tabPanel("Model Settings",
                 selectInput("model_type", "Select Method:",
                             choices = c("Mixture Prior", "Bayesian Hierarchical Model", "Power Prior"),
                             selected = "Power Prior"),
                 
                 # PowerPrior特定设置
                 conditionalPanel(
                   condition = "input.model_type == 'Power Prior'",
                   numericInput("alpha_min", "Minimum Alpha", 0.001, min = 0, max = 1, step = 0.001),
                   numericInput("alpha_max", "Maximum Alpha", 0.01, min = 0, max = 1, step = 0.001),
                   numericInput("alpha_steps", "Number of Steps", 50, min = 10, max = 200)
                 )
        )
      ),
      
      actionButton("analyze", "Run Analysis", class = "btn-primary btn-lg")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tipping Point",
                 plotOutput("tippingPlot", height = "400px"),
                 verbatimTextOutput("tippingConclusion")
        ),
        tabPanel("Probability",
                 plotOutput("probabilityPlot", height = "400px")
        ),
        tabPanel("ESS Analysis",
                 plotOutput("essPlot", height = "400px"),
                 plotOutput("groupEssPlot", height = "400px"),
                 verbatimTextOutput("essDetails")
        ),
        tabPanel("Variance Ratio",
                 plotOutput("varianceRatioPlot", height = "400px"),
                 plotOutput("groupVarianceRatioPlot", height = "400px")
        )
      )
    )
  )
)

# 服务器逻辑
server <- function(input, output, session) {
  # 反应式数据和状态变量
  analysis_results <- reactiveVal(NULL)
  is_analyzing <- reactiveVal(FALSE)
  
  # 处理分析按钮点击
  observeEvent(input$analyze, {
    is_analyzing(TRUE)
    
    # 收集数据
    ped_data <- list(
      y_placebo = input$ped_ctrl_resp,
      n_placebo = input$ped_ctrl_total,
      y_drug = input$ped_treat_resp,
      n_drug = input$ped_treat_total
    )
    
    adult_data <- list(
      y_placebo = input$adult_ctrl_resp,
      n_placebo = input$adult_ctrl_total,
      y_drug = input$adult_treat_resp,
      n_drug = input$adult_treat_total
    )
    
    if (input$model_type == "Power Prior") {
      # 运行PowerPrior分析
      results <- run_power_prior_analysis(
        ped_data, 
        adult_data, 
        input$alpha_min, 
        input$alpha_max, 
        input$alpha_steps
      )
      
      analysis_results(results)
    } else if (input$model_type == "Bayesian Hierarchical Model") {
      # 这里实现您现有的BHM代码
      # ...
    } else if (input$model_type == "Mixture Prior") {
      # 这里实现MixturePrior代码
      # ...
    }
    
    is_analyzing(FALSE)
  })
  
  # 渲染临界点图
  output$tippingPlot <- renderPlot({
    req(analysis_results())
    
    results <- analysis_results()$results
    tipping_point <- analysis_results()$tipping_point
    
    # 将ggplot对象赋值给变量p
    p <- ggplot(results, aes(x = alpha)) +
      geom_line(aes(y = median_rr_diff), size = 1, color = "darkblue") +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = "steelblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = "PowerPrior Tipping Point Analysis",
           subtitle = "Response Rate Difference with 95% Credible Intervals",
           x = "Borrowing Weight (α)",
           y = "Response Rate Difference (Drug - Placebo)") +
      theme_minimal() +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    # 如果找到tipping point，添加额外的图层
    if (!is.na(tipping_point)) {
      p <- p + 
        geom_vline(xintercept = tipping_point, linetype = "dotted", color = "purple", size = 1) +
        annotate("text", 
                 x = tipping_point + (max(results$alpha) - min(results$alpha)) * 0.05, 
                 y = min(results$lower_ci) + (max(results$upper_ci) - min(results$lower_ci)) * 0.9,
                 label = paste0("Tipping Point: ", scales::percent(tipping_point, accuracy = 0.1)),
                 color = "purple", hjust = 0)
    }
    
    # 返回绘制好的图表
    return(p)
  })
  
  # 渲染概率图
  output$probabilityPlot <- renderPlot({
    req(analysis_results())
    
    results <- analysis_results()$results
    tipping_point <- analysis_results()$tipping_point
    
    p <- ggplot(results, aes(x = alpha, y = prob_gt_0)) +
      geom_line(size = 1, color = "darkgreen") +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
      labs(title = "Posterior Probability of Positive Treatment Effect",
           subtitle = "Probability that Response Rate Difference > 0",
           x = "Borrowing Weight (α)",
           y = "P(RR Diff > 0)") +
      theme_minimal() +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    if (!is.na(tipping_point)) {
      p <- p + geom_vline(xintercept = tipping_point, linetype = "dotted", color = "purple", size = 1)
    }
    
    return(p)
  })
  
  # 渲染ESS图
  output$essPlot <- renderPlot({
    req(analysis_results())
    
    ess_results <- analysis_results()$ess_results
    tipping_point <- analysis_results()$tipping_point
    
    p <- ggplot(ess_results, aes(x = alpha)) +
      geom_line(aes(y = ess_total), size = 1, color = "darkorange") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.8) +
      labs(title = "Borrowed Sample Size Analysis",
           subtitle = "ESS = n * (V1/V2 - 1)",
           x = "Borrowing Weight (α)",
           y = "Borrowed Sample Size") +
      theme_minimal() +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    if (!is.na(tipping_point)) {
      tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
      tp_ess <- ess_results$ess_total[tp_idx]
      
      p <- p + 
        geom_vline(xintercept = tipping_point, linetype = "dotted", color = "purple", size = 1) +
        annotate("text", 
                 x = tipping_point + (max(ess_results$alpha) - min(ess_results$alpha)) * 0.05, 
                 y = max(ess_results$ess_total) * 0.5,
                 label = paste0("Borrowed ESS at TP: ", round(tp_ess, 1)),
                 color = "purple", hjust = 0)
    }
    
    return(p)
  })
  
  # 渲染分组ESS图
  output$groupEssPlot <- renderPlot({
    req(analysis_results())
    
    ess_results <- analysis_results()$ess_results
    tipping_point <- analysis_results()$tipping_point
    
    p <- ggplot(ess_results, aes(x = alpha)) +
      geom_line(aes(y = drug_ess, color = "Drug Group"), size = 1) +
      geom_line(aes(y = placebo_ess, color = "Placebo Group"), size = 1) +
      geom_line(aes(y = ess_total, color = "Total"), size = 1, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.6) +
      labs(title = "Group-Specific Borrowed Sample Size Analysis",
           subtitle = "Based on Group-Specific Variance Ratios",
           x = "Borrowing Weight (α)",
           y = "Borrowed Sample Size",
           color = "Group") +
      scale_color_manual(values = c("Drug Group" = "blue", 
                                    "Placebo Group" = "red", 
                                    "Total" = "purple")) +
      theme_minimal() +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      theme(legend.position = "bottom")
    
    if (!is.na(tipping_point)) {
      tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
      tp_drug_ess <- ess_results$drug_ess[tp_idx]
      tp_placebo_ess <- ess_results$placebo_ess[tp_idx]
      tp_total_ess <- ess_results$ess_total[tp_idx]
      
      p <- p + 
        geom_vline(xintercept = tipping_point, linetype = "dotted", color = "black", size = 0.8) +
        annotate("text", 
                 x = tipping_point + (max(ess_results$alpha) - min(ess_results$alpha)) * 0.05, 
                 y = max(c(ess_results$drug_ess, ess_results$placebo_ess, ess_results$ess_total)) * 0.8,
                 label = paste0("At Tipping Point (α = ", 
                                scales::percent(tipping_point, accuracy = 0.1), "):",
                                "\nDrug Borrowed: ", round(tp_drug_ess, 1),
                                "\nPlacebo Borrowed: ", round(tp_placebo_ess, 1),
                                "\nTotal Borrowed: ", round(tp_total_ess, 1)),
                 color = "black", hjust = 0, size = 3)
    }
    
    return(p)
  })
  
  # 渲染方差比图
  output$varianceRatioPlot <- renderPlot({
    req(analysis_results())
    
    ess_results <- analysis_results()$ess_results
    tipping_point <- analysis_results()$tipping_point
    
    p <- ggplot(ess_results, aes(x = alpha)) +
      geom_line(aes(y = variance_ratio), size = 1, color = "purple") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(title = "Variance Ratio Analysis",
           subtitle = "Ratio of posterior variances: V1/V2",
           x = "Borrowing Weight (α)",
           y = "Variance Ratio (V1/V2)") +
      theme_minimal() +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))
    
    if (!is.na(tipping_point)) {
      tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
      tp_var_ratio <- ess_results$variance_ratio[tp_idx]
      
      p <- p + 
        geom_vline(xintercept = tipping_point, linetype = "dotted", color = "purple", size = 1) +
        annotate("text", 
                 x = tipping_point + (max(ess_results$alpha) - min(ess_results$alpha)) * 0.05, 
                 y = max(ess_results$variance_ratio) * 0.8,
                 label = paste0("Var Ratio at TP: ", round(tp_var_ratio, 2)),
                 color = "purple", hjust = 0)
    }
    
    return(p)
  })
  
  # 渲染分组方差比图
  output$groupVarianceRatioPlot <- renderPlot({
    req(analysis_results())
    
    ess_results <- analysis_results()$ess_results
    tipping_point <- analysis_results()$tipping_point
    
    p <- ggplot(ess_results, aes(x = alpha)) +
      geom_line(aes(y = drug_variance_ratio, color = "Drug Group"), size = 1) +
      geom_line(aes(y = placebo_variance_ratio, color = "Placebo Group"), size = 1) +
      geom_line(aes(y = variance_ratio, color = "Total"), size = 1, linetype = "dashed") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(title = "Group-Specific Variance Ratio Analysis",
           subtitle = "Ratio of posterior variances: V1/V2",
           x = "Borrowing Weight (α)",
           y = "Variance Ratio (V1/V2)",
           color = "Group") +
      scale_color_manual(values = c("Drug Group" = "blue", 
                                    "Placebo Group" = "red", 
                                    "Total" = "purple")) +
      theme_minimal() +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      theme(legend.position = "bottom")
    
    if (!is.na(tipping_point)) {
      tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
      tp_drug_var_ratio <- ess_results$drug_variance_ratio[tp_idx]
      tp_placebo_var_ratio <- ess_results$placebo_variance_ratio[tp_idx]
      tp_total_var_ratio <- ess_results$variance_ratio[tp_idx]
      
      p <- p + 
        geom_vline(xintercept = tipping_point, linetype = "dotted", color = "black", size = 0.8) +
        annotate("text", 
                 x = tipping_point + (max(ess_results$alpha) - min(ess_results$alpha)) * 0.05, 
                 y = max(c(ess_results$drug_variance_ratio, ess_results$placebo_variance_ratio, 
                           ess_results$variance_ratio)) * 0.8,
                 label = paste0("At Tipping Point (α = ", 
                                scales::percent(tipping_point, accuracy = 0.1), "):",
                                "\nDrug Var Ratio: ", round(tp_drug_var_ratio, 2),
                                "\nPlacebo Var Ratio: ", round(tp_placebo_var_ratio, 2),
                                "\nTotal Var Ratio: ", round(tp_total_var_ratio, 2)),
                 color = "black", hjust = 0, size = 3)
    }
    
    return(p)
  })
  
  # 输出详细ESS信息
  output$essDetails <- renderPrint({
    req(analysis_results())
    
    ess_results <- analysis_results()$ess_results
    tipping_point <- analysis_results()$tipping_point
    
    if (!is.na(tipping_point)) {
      tp_idx <- which(abs(ess_results$alpha - tipping_point) == min(abs(ess_results$alpha - tipping_point)))
      
      cat("===== PowerPrior Tipping Point Analysis Results =====\n\n")
      cat("Tipping Point Alpha:", scales::percent(tipping_point, accuracy = 0.001), "\n\n")
      
      cat("Total Borrowed ESS:", round(ess_results$ess_total[tp_idx], 2), "\n")
      cat("Total Variance Ratio:", round(ess_results$variance_ratio[tp_idx], 4), "\n\n")
      
      cat("Placebo Group Borrowed ESS:", round(ess_results$placebo_ess[tp_idx], 2), "\n")
      cat("Placebo Group Variance Ratio:", round(ess_results$placebo_variance_ratio[tp_idx], 4), "\n\n")
      
      cat("Drug Group Borrowed ESS:", round(ess_results$drug_ess[tp_idx], 2), "\n")
      cat("Drug Group Variance Ratio:", round(ess_results$drug_variance_ratio[tp_idx], 4), "\n")
    } else {
      cat("Tipping point not found within the specified alpha range.\n")
      cat("Try expanding the alpha range to find a valid tipping point.")
    }
  })
  
  # 输出临界点结论
  output$tippingConclusion <- renderPrint({
    req(analysis_results())
    
    tipping_point <- analysis_results()$tipping_point
    
    if (!is.na(tipping_point)) {
      cat("PowerPrior Tipping Point (RR Diff lower bound > 0):", 
          scales::percent(tipping_point, accuracy = 0.001), "\n")
      
      # 获取临界点处的RR_diff
      results <- analysis_results()$results
      tp_idx <- which(abs(results$alpha - tipping_point) == min(abs(results$alpha - tipping_point)))
      
      if (length(tp_idx) > 0) {
        cat("At tipping point:\n")
        cat("- Median RR Diff:", round(results$median_rr_diff[tp_idx], 4), "\n")
        cat("- 95% CI: [", round(results$lower_ci[tp_idx], 4), ", ", 
            round(results$upper_ci[tp_idx], 4), "]\n")
        cat("- P(RR Diff > 0):", scales::percent(results$prob_gt_0[tp_idx], accuracy = 0.1), "\n")
      }
    } else {
      cat("PowerPrior Tipping Point: Not found in the specified alpha range\n")
      cat("Try increasing the maximum alpha value to find a valid tipping point.")
    }
  })
}

# 运行Shiny应用
shinyApp(ui = ui, server = server)