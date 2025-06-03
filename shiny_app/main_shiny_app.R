
library(shiny)
library(plotly)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(rjags)

# === Load analysis functions ===
source(here::here("shiny_app", "mixture_analysis_function.R"))
source(here::here("shiny_app", "hierarchical_tipping_function.R"))
source(here::here("shiny_app", "power_prior_functions_new.R"))

ui <- fluidPage(
  titlePanel("Pediatric Extrapolation with Multiple Prior Methods"),
  
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
        )
      ),
      
      tags$hr(style = "border-top: 2px solid #aaa;"),
      h4("Settings"),
      strong("Baseline:"),
      actionButton("fisher", "Fisher's Exact Test"),
      selectInput("model_type", "Bayesian Analysis:",
                  choices = c("Mixture Prior", "Bayesian Hierarchical Model", "Power Prior")),
      
      # Power Prior设置（条件显示）
      conditionalPanel(
        condition = "input.model_type == 'Power Prior'",
        hr(),
        h5("Power Prior Settings"),
        numericInput("alpha_min", "Minimum Alpha", 0.001, min = 0, max = 1, step = 0.001),
        numericInput("alpha_max", "Maximum Alpha", 0.01, min = 0, max = 1, step = 0.001),
        numericInput("alpha_steps", "Number of Steps", 50, min = 10, max = 200),
        helpText("Alpha represents the borrowing weight. Higher values mean more borrowing from adult data.")
      ),
      
      actionButton("analyze", "Run Analysis")
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.analysisVisible == true",
        tabsetPanel(
          tabPanel("📊 Tipping Point Plot",
                   withSpinner(plotlyOutput("tippingPlot")),
                   # 添加 Error Bar 信息显示（只对 Power Prior 显示）
                   conditionalPanel(
                     condition = "input.model_type == 'Power Prior'",
                     uiOutput("errorBarInfo")
                   ),
                   verbatimTextOutput("tippingConclusion")
          ),
          tabPanel("📈 ESS Calculation",
                   withSpinner(tableOutput("essTable")),
                   uiOutput("essNote")
          )
        )
      ),
      conditionalPanel(
        condition = "output.showFisherPanel == true",
        tags$hr(),
        h4("Fisher's Exact Test Result"),
        verbatimTextOutput("fisherResult")
      )
    )
  )
)

server <- function(input, output, session) {
  
  values <- reactiveValues(analysisDone = FALSE, showFisher = FALSE)
  
  # 创建 reactive 表达式来存储 Power Prior 分析结果
  powerPriorResults <- eventReactive(input$analyze, {
    if (input$model_type == "Power Prior") {
      
      # 验证输入参数
      if (input$alpha_min >= input$alpha_max) {
        showNotification("Alpha minimum must be less than maximum", type = "error")
        return(NULL)
      }
      
      child_data <- list(
        treat = list(y = input$ped_treat_resp, n = input$ped_treat_total),
        control = list(y = input$ped_ctrl_resp, n = input$ped_ctrl_total)
      )
      adult_data <- list(
        treat = list(y = input$adult_treat_resp, n = input$adult_treat_total),
        control = list(y = input$adult_ctrl_resp, n = input$adult_ctrl_total)
      )
      
      # 运行分析
      withProgress(message = 'Running Power Prior Analysis...', value = 0, {
        incProgress(0.1, detail = "Initializing...")
        
        # 运行分析
        analysis_results <- run_power_prior_analysis(
          child_data, 
          adult_data, 
          input$alpha_min,
          input$alpha_max,
          input$alpha_steps
        )
        
        # 获取详细分析结果
        detailed_results <- get_power_prior_analysis(
          child_data, 
          adult_data,
          input$alpha_min,
          input$alpha_max,
          input$alpha_steps
        )
        
        incProgress(0.9, detail = "Finalizing...")
        
        return(list(
          analysis = analysis_results,
          detailed = detailed_results,
          child_data = child_data,
          adult_data = adult_data,
          alpha_min = input$alpha_min,
          alpha_max = input$alpha_max,
          alpha_steps = input$alpha_steps
        ))
      })
    } else {
      return(NULL)
    }
  })
  
  observeEvent(input$analyze, {
    
    values$showFisher <- FALSE  # 运行分析时隐藏Fisher结果
    values$analysisDone <- TRUE
    
    child_data <- list(
      treat = list(y = input$ped_treat_resp, n = input$ped_treat_total),
      control = list(y = input$ped_ctrl_resp, n = input$ped_ctrl_total)
    )
    adult_data <- list(
      treat = list(y = input$adult_treat_resp, n = input$adult_treat_total),
      control = list(y = input$adult_ctrl_resp, n = input$adult_ctrl_total)
    )
    
    if (input$model_type == "Mixture Prior") {
      
      output$tippingPlot <- renderPlotly({
        run_mixture_plot(child_data, adult_data)
      })
      
      output$tippingConclusion <- renderText({
        res <- mixture_analysis(child_data, adult_data)
        tp <- res$tipping_point_summary$RR
        if (!is.null(tp) && !is.na(tp$weight)) {
          paste0("Tipping Point Weight = ", round(tp$weight, 2), "\n",
                 "95% CI: [", round(tp$lower_95_rr, 2), ",", round(tp$upper_95_rr, 2), "]\n")
        } else {
          "No tipping point found within the tested weight range."
        }
      })
      
      output$essTable <- renderTable({
        res <- mixture_analysis(child_data, adult_data)
        tp <- res$tipping_point_summary$RR
        
        if (!is.null(tp) && !is.na(tp$weight)) {
          tibble::tibble(
            Group = c("Treatment", "Control"),
            `Delta ESS` = c(round(tp$ess_treat, 1), round(tp$ess_control, 1))
          )
        } else {
          tibble::tibble(
            Group = character(0), 
            `Delta ESS` = numeric(0)
          )
        }
      })
      output$essNote <- renderUI({
        res <- mixture_analysis(child_data, adult_data)
        tp <- res$tipping_point_summary$RR
        
        if (!is.null(tp) && !is.na(tp$weight)) {
          HTML(paste0(
            "<small><i>Tipping Point Weight = ", round(tp$weight, 2), "</i></small>"
          ))
        } else {
          HTML("<small><i>No tipping point found.</i></small>")
        }
      })
      
    }
    
    if (input$model_type == "Bayesian Hierarchical Model") {
      
      ## === Step 1: Construct pediatric and adult data from user input ===
      pediatric_clean <- data.frame(
        Group = c("Placebo", "DrugA"),
        n = c(input$ped_ctrl_resp, input$ped_treat_resp),
        N = c(input$ped_ctrl_total, input$ped_treat_total),
        Trial = "Pooled"
      ) |> 
        mutate(Source = "Pediatric")
      
      adult_clean <- data.frame(
        Group = c("Placebo", "DrugA"),
        n = c(input$adult_ctrl_resp, input$adult_treat_resp),
        N = c(input$adult_ctrl_total, input$adult_treat_total),
        Trial = "Pooled"
      ) |> 
        mutate(Source = "Adult")
      
      combined_data <- bind_rows(pediatric_clean, adult_clean) |>
        mutate(
          group = ifelse(Group == "DrugA", 1, 0),
          trial_str = paste(Source, Trial),
          trial = as.integer(factor(trial_str)),
          y = n
        )
      
      ## === Step 2: Bayesian model loop ===
      tipping_results <- run_bayesian_tipping_analysis(combined_data)
      
      ## === Step 3: 输出图与结论 ===
      tipping_results$Significant_FDA <- tipping_results$OR_lower > 1
      tipping_row_fda <- tipping_results |>
        filter(Significant_FDA == TRUE) |>
        slice_min(ESS_FDA)
      
      output$tippingConclusion <- renderText({
        if (nrow(tipping_row_fda) > 0) {
          paste0("Tipping Point: sigma = ", round(tipping_row_fda$fixed_sigma_alpha, 2),
                 "\n95% CI: [", round(tipping_row_fda$OR_lower, 2), ", ", round(tipping_row_fda$OR_upper, 2), "]")
        } else {
          "No tipping point found in the tested sigma range."
        }
      })
      
      output$tippingPlot <- renderPlotly({
        p <- ggplot(tipping_results, aes(
          x = fixed_sigma_alpha,
          y = OR_median,
          color = Significant_FDA,
          text = paste0(
            "Sigma: ", round(fixed_sigma_alpha, 2), "<br>",
            "OR: ", round(OR_median, 2), " [", round(OR_lower, 2), ", ", round(OR_upper, 2), "]<br>",
            "Placebo ESS: ", round(ESS_FDA, 1), "<br>",
            "Treatment ESS: ", round(ESS_FDA_trt, 1)
          )
        )) +
          geom_point(size = 3) +
          geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
          geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
          {if(nrow(tipping_row_fda) > 0) 
            geom_vline(xintercept = tipping_row_fda$fixed_sigma_alpha, linetype = "dotted", color = "red")} +
          {if(nrow(tipping_row_fda) > 0)
            annotate("text",
                     x = tipping_row_fda$fixed_sigma_alpha * 0.6,
                     y = tipping_row_fda$OR_median + 80,
                     label = paste0("Tipping Point = ", round(tipping_row_fda$fixed_sigma_alpha, 2)),
                     hjust = 0, vjust = 0, size = 4.5, fontface = "italic", color = "red")} +
          scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray")) +
          scale_y_continuous(
            trans = scales::log10_trans(),
            breaks = c(1, 10, 50, 100)
          ) +
          labs(
            title = "Pediatric OR vs Borrowing Strength (σ)",
            x = "sigma alpha",
            y = "Odds Ratio (DrugA vs Placebo)"
          ) +
          theme_minimal(base_size = 15)
        
        # 关键一步：转为交互图
        ggplotly(p, tooltip = "text")
      })
      
      output$essTable <- renderTable({
        if (nrow(tipping_row_fda) > 0) {
          tipping_row <- tipping_results %>%
            filter(Significant_FDA == TRUE) %>%
            slice_min(ESS_FDA)
          
          placebo_borrowed <- round(tipping_row$Borrowed_FDA, 1)
          placebo_total <- round(tipping_row$ESS_FDA, 1)
          
          treatment_borrowed <- round(tipping_row$Borrowed_FDA_trt, 1)
          treatment_total <- round(tipping_row$ESS_FDA_trt, 1)
          
          total_borrowed <- round(placebo_borrowed + treatment_borrowed, 1)
          total_total <- round(placebo_total + treatment_total, 1)
          
          tibble::tibble(
            Group = c("Placebo", "Treatment"),
            `Delta ESS` = c(placebo_borrowed, treatment_borrowed)
          )
        } else {
          tibble::tibble(
            Group = character(0),
            `Delta ESS` = numeric(0)
          )
        }
      })
      output$essNote <- renderUI({
        if (nrow(tipping_row_fda) > 0) {
          tipping_row <- tipping_results %>%
            filter(Significant_FDA == TRUE) %>%
            slice_min(ESS_FDA)
          
          HTML(paste0(
            "<small><i>Tipping Point Sigma Alpha = ",
            round(tipping_row$fixed_sigma_alpha, 2),
            "</i></small>"
          ))
        } else {
          HTML("<small><i>No tipping point found under FDA criteria.</i></small>")
        }
      })
      
    }
    
    # Power Prior分析 - 修改为使用 reactive 结果
    if (input$model_type == "Power Prior") {
      
      # 生成图表 - 基于 reactive 结果
      output$tippingPlot <- renderPlotly({
        results <- powerPriorResults()
        if (!is.null(results)) {
          run_power_prior_with_params(
            results$child_data, 
            results$adult_data,
            results$alpha_min,
            results$alpha_max,
            results$alpha_steps
          )
        }
      })
      
      # Error Bar 信息 - 基于 reactive 结果
      output$errorBarInfo <- renderUI({
        results <- powerPriorResults()
        if (!is.null(results)) {
          # 重新计算统计信息
          analysis_results <- results$analysis
          combined_results <- merge(analysis_results$results, analysis_results$ess_results, by = "alpha")
          combined_results$CI_Exclude_0 <- combined_results$lower_ci > 0
          
          total_error_bars <- nrow(combined_results)
          significant_error_bars <- sum(combined_results$CI_Exclude_0)
          non_significant_error_bars <- total_error_bars - significant_error_bars
          
          HTML(paste0(
            "<div style='background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin: 10px 0;'>",
            "<strong>📊 Error Bar Statistics:</strong><br>",
            "📈 Total Error Bars: <span style='color: #007bff;'><strong>", total_error_bars, "</strong></span><br>",
            "✅ Significant (CI excludes 0): <span style='color: #28a745;'><strong>", significant_error_bars, "</strong></span><br>",
            "❌ Non-significant (CI includes 0): <span style='color: #6c757d;'><strong>", non_significant_error_bars, "</strong></span><br>",
            "🎯 Alpha Range: ", scales::percent(results$alpha_min, accuracy = 0.1), " - ", scales::percent(results$alpha_max, accuracy = 0.1),
            "</div>"
          ))
        }
      })
      
      # 生成结论 - 基于 reactive 结果
      output$tippingConclusion <- renderText({
        results <- powerPriorResults()
        if (!is.null(results)) {
          tp <- results$detailed$tipping_point_summary$RR
          if (!is.null(tp) && !is.na(tp$weight)) {
            paste0("Tipping Point Weight = ", scales::percent(tp$weight, accuracy = 0.001), "\n",
                   "ESS (Drug) = ", round(tp$ess_treat, 1), ", ESS (Placebo) = ", round(tp$ess_control, 1))
          } else {
            "No tipping point found within the specified alpha range."
          }
        }
      })
      
      # ESS 表格 - 基于 reactive 结果
      output$essTable <- renderTable({
        results <- powerPriorResults()
        if (!is.null(results)) {
          tp <- results$detailed$tipping_point_summary$RR
          
          if (!is.null(tp) && !is.na(tp$weight)) {
            tibble::tibble(
              Group = c("Treatment", "Control"),
              `Delta ESS` = c(round(tp$ess_treat, 1), round(tp$ess_control, 1))
            )
          } else {
            tibble::tibble(
              Group = character(0), 
              `Delta ESS` = numeric(0)
            )
          }
        }
      })
      
      # ESS 注释 - 基于 reactive 结果
      output$essNote <- renderUI({
        results <- powerPriorResults()
        if (!is.null(results)) {
          tp <- results$detailed$tipping_point_summary$RR
          
          if (!is.null(tp) && !is.na(tp$weight)) {
            HTML(paste0(
              "<small><i>Tipping Point Weight = ", round(tp$weight, 3), "</i></small>"
            ))
          } else {
            HTML("<small><i>No tipping point found in Power Prior analysis.</i></small>")
          }
        }
      })
    }
    
    # 输出控制
    output$analysisVisible <- reactive({
      values$analysisDone
    })
    outputOptions(output, "analysisVisible", suspendWhenHidden = FALSE)
    
  }) # end observeEvent for analyze
  
  # Fisher's Exact Test 事件处理
  observeEvent(input$fisher, {
    # 1. 构造2x2列联表
    values$showFisher <- TRUE 
    table_mat <- matrix(
      c(
        input$ped_treat_resp,
        input$ped_treat_total - input$ped_treat_resp,
        input$ped_ctrl_resp,
        input$ped_ctrl_total - input$ped_ctrl_resp
      ),
      nrow = 2,
      byrow = TRUE
    )
    rownames(table_mat) <- c("DrugA", "Placebo")
    colnames(table_mat) <- c("Success", "Failure")
    
    # 2. 执行Fisher's检验
    fisher_res <- fisher.test(table_mat)
    
    # 3. 输出检验结果
    output$fisherResult <- renderPrint({
      cat("2x2 Table (Pediatric Data):\n\n")
      print(table_mat)
      cat("\nFisher's Exact Test Result:\n")
      cat("P-value:", signif(fisher_res$p.value, 4), "\n")
      cat("Odds Ratio:", signif(fisher_res$estimate, 4), "\n")
      cat("95% CI:", paste0("(", signif(fisher_res$conf.int[1], 4), ", ", signif(fisher_res$conf.int[2], 4), ")"), "\n")
      
      # 判断显著性
      if (fisher_res$p.value > 0.05 || (fisher_res$conf.int[1] < 1 && fisher_res$conf.int[2] > 1)) {
        cat("→ Not statistically significant.\n")
      } else {
        cat("→ Statistically significant.\n")
      }
    })
  })
  
  # Fisher面板显示控制
  output$showFisherPanel <- reactive({
    values$showFisher
  })
  outputOptions(output, "showFisherPanel", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)