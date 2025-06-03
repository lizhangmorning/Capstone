library(shiny)
library(plotly)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinycssloaders)
library(rjags)
library(DT)  # 🔥 新增：用于数据表格

# === Load analysis functions ===
source(here::here("shiny_app", "mixture_analysis_function.R"), echo = TRUE, print.eval = TRUE)
source(here::here("shiny_app", "hierarchical_tipping_function.R"), echo = TRUE, print.eval = TRUE)
source(here::here("shiny_app", "power_prior_functions_new.R"), echo = TRUE, print.eval = TRUE)

ui <- fluidPage(
  titlePanel("Pediatric Extrapolation with Bayesian Methods"),
  
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
      
      conditionalPanel(
        condition = "input.model_type == 'Power Prior'|| input.model_type == 'Mixture Prior'",
        hr(),
        h5("Prior Settings"),
        numericInput("alpha_min", "Minimum Alpha", 0.001, min = 0, max = 1, step = 0.001),
        numericInput("alpha_max", "Maximum Alpha", 0.01, min = 0, max = 1, step = 0.001),
        numericInput("alpha_steps", "Number of Steps", 50, min = 10, max = 200),
        helpText("Alpha represents the borrowing weight. Higher values mean more borrowing from adult data.")
      ),
      
      conditionalPanel(
        condition = "input.model_type == 'Bayesian Hierarchical Model'",
        hr(),
        h5("Prior Settings"),
        numericInput("sigma_min", "Minimum Sigma", 1, min = 0, max = 50, step = 1),
        numericInput("sigma_max", "Maximum Sigma", 2, min = 0, max = 1000000, step = 1),
        numericInput("sigma_steps", "Number of Steps", 10, min = 1, max = 200, step = 1),
        helpText("Sigma controls the borrowing level in the Bayesian hierarchical model. Higher sigma allows less borrowing from adult data.")
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
                     condition = "input.model_type == 'Power Prior' ||  'Bayesian Hierarchical Model' || 'Mixture Prior'",
                     uiOutput("errorBarInfo")
                   ),
                   verbatimTextOutput("tippingConclusion"),
                   
                   # 🔥 添加Power Prior解释框
                   conditionalPanel(
                     condition = "input.model_type == 'Power Prior'",
                     div(style = "margin-top: 15px; padding: 10px; background-color: #e8f4f8; border-left: 4px solid #17a2b8; border-radius: 3px;",
                         h5("💡 Power Prior Methodology", style = "color: #17a2b8; margin-bottom: 8px;"),
                         tags$p(style = "margin-bottom: 5px; font-size: 14px;",
                                strong("Prior ~ Beta(1 + αy, 1 + α(n - y))")),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• α (alpha): weight applied to adult study data"),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• Larger α → more borrowing from adult data"),
                         tags$p(style = "margin-bottom: 0px; font-size: 13px; color: #6c757d;",
                                "• Tipping point: minimum α where CI excludes 0")
                     )
                   ),
                   
                   # 🔶 Hierarchical Model 说明框
                   conditionalPanel(
                     condition = "input.model_type == 'Bayesian Hierarchical Model'",
                     div(style = "margin-top: 15px; padding: 10px; background-color: #e8f4f8; border-left: 4px solid #17a2b8; border-radius: 3px;",
                         h5("💡 Bayesian Hierarchical Model Methodology", style = "color: #17a2b8; margin-bottom: 8px;"),
                         tags$p(style = "margin-bottom: 5px; font-size: 14px;",
                                strong("Model Formulation:")),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                HTML("yi ~ Binomial(pi, ni)")),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                HTML("logit(pi) = &alpha;<sub>trial[i]</sub> + &beta;<sub>trial[i]</sub> * group<sub>i</sub>")),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                HTML("&alpha;<sub>trial</sub> ~ N(&mu;<sub>&alpha;</sub>, &sigma;<sub>&alpha;</sub><sup>2</sup>), &beta;<sub>trial</sub> ~ N(&mu;<sub>&beta;</sub>, &sigma;<sub>&beta;</sub><sup>2</sup>)")),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• group_i = 1 for Drug A, 0 for Placebo"),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• trial[i] = 1 for adult, 2 for pediatric"),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• αtrial: baseline log-odds for each population"),
                         tags$p(style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                                "• βtrial: treatment effect (in log-odds) for each population"),
                         tags$p(style = "margin-bottom: 0px; font-size: 13px; color: #6c757d;",
                                "• Tipping point: minimum σ where CI excludes 1")
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.model_type == 'Mixture Prior'",
                     div(
                       style = "margin-top: 15px; padding: 10px; background-color: #e8f4f8; border-left: 4px solid #17a2b8; border-radius: 3px;",
                       h5("💡 Mixture Prior Methodology", style = "color: #17a2b8; margin-bottom: 8px;"),
                       tags$p(
                         style = "margin-bottom: 5px; font-size: 14px;",
                         strong("Prior ~ (1 - α) · Beta(1, 1) + α · Beta(y, n - y)")
                       ),
                       tags$p(
                         style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                         "• α (alpha): weight applied to adult study data"
                       ),
                       tags$p(
                         style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                         "• Beta(1, 1): skeptical prior; Beta(y, n - y): adult effect estimate distribution"
                       ),
                       tags$p(
                         style = "margin-bottom: 5px; font-size: 13px; color: #6c757d;",
                         "• Larger α → more borrowing from adult data"
                       ),
                       tags$p(
                         style = "margin-bottom: 0px; font-size: 13px; color: #6c757d;",
                         "• Tipping point: minimum α where CI excludes 0"
                       )
                     )
                   )
                   
                   
                   
          ),
          tabPanel("📈 Results",
                   conditionalPanel(
                     condition = "input.model_type == 'Power Prior'|| input.model_type =='Mixture Prior' input.model_type == 'Bayesian Hierarchical Model'",
                     withSpinner(DT::dataTableOutput("essTableDetailed"))
                   ),
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
  
  values <- reactiveValues(
    analysisDone = FALSE, 
    showFisher = FALSE
    )
  
  
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
  
  mixturePriorResults <- eventReactive(input$analyze, {
    if (input$model_type == "Mixture Prior") {
      
      # Validate input parameters
      if (input$alpha_min >= input$alpha_max) {
        showNotification("Alpha minimum must be less than maximum", type = "error")
        return(NULL)
      }
      
      # Prepare data
      child_data <- list(
        treat = list(y = input$ped_treat_resp, n = input$ped_treat_total),
        control = list(y = input$ped_ctrl_resp, n = input$ped_ctrl_total)
      )
      adult_data <- list(
        treat = list(y = input$adult_treat_resp, n = input$adult_treat_total),
        control = list(y = input$adult_ctrl_resp, n = input$adult_ctrl_total)
      )
      
      # Run analysis with progress
      withProgress(message = 'Running Mixture Prior Analysis...', value = 0, {
        incProgress(0.1, detail = "Initializing...")
        
        # Run analysis
        analysis_results <- mixture_analysis(
          child_data = child_data,
          adult_data = adult_data,
          alpha_min = input$alpha_min,
          alpha_max = input$alpha_max,
          alpha_steps = input$alpha_steps
        )
        
        # Get detailed results (using the same function in this case)
        detailed_results <- mixture_analysis(
          child_data = child_data,
          adult_data = adult_data,
          alpha_min = input$alpha_min,
          alpha_max = input$alpha_max,
          alpha_steps = input$alpha_steps
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
      
      # Generate plot
      output$tippingPlot <- renderPlotly({
        results <- mixturePriorResults()
        if (!is.null(results)) {
          plot_mixture_results(
            results$analysis,
            results$alpha_min,
            results$alpha_max,
            results$alpha_steps
          )
        }
      })
      
      # Error Bar info
      output$errorBarInfo <- renderUI({
        results <- mixturePriorResults()
        if (!is.null(results)) {
          analysis_results <- results$analysis
          combined_results <- analysis_results$results
          combined_results$Favorable <- combined_results$lower_95_rr > 0
          
          total_error_bars <- nrow(combined_results)
          significant_error_bars <- sum(combined_results$Favorable)
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
      
      # Generate conclusion
      output$tippingConclusion <- renderText({
        results <- mixturePriorResults()
        if (!is.null(results)) {
          tp <- results$analysis$tipping_point_summary$RR
          if (!is.null(tp) && !is.na(tp$weight)) {
            paste0("Tipping Point Weight = ", scales::percent(tp$weight, accuracy = 0.001), "\n",
                   "ESS (Drug) = ", round(tp$ess_treat, 1), ", ESS (Placebo) = ", round(tp$ess_control, 1))
          } else {
            "No tipping point found within the specified alpha range."
          }
        }
      })
      
      output$essTableDetailed <- DT::renderDataTable({
        results <- mixturePriorResults()
        
        if (!is.null(results)) {
          analysis_results <- results$analysis
          display_data <- analysis_results$results %>%
            dplyr::transmute(
              `Weight (α)` = scales::percent(weight, accuracy = 0.001),
              `CI Lower` = round(lower_95_rr, 3),
              `CI Upper` = round(upper_95_rr, 3),
              `Delta ESS Treatment` = round(ess_treat_rr, 1),
              `Delta ESS Control` = round(ess_control_rr, 1)
            )
          
          DT::datatable(
            display_data,
            options = list(
              pageLength = 25,
              scrollX = TRUE
            ),
            rownames = FALSE
          )
          
        } else {
          DT::datatable(
            data.frame(Message = "Please run Mixture Prior analysis first"),
            options = list(dom = 't'),
            rownames = FALSE
          )
        }
      })
      
      
      # Simplified note
      output$essNote <- renderUI({
        results <- mixturePriorResults()
        if (!is.null(results)) {
          tp <- results$analysis$tipping_point_summary$RR
          if (!is.null(tp) && !is.na(tp$weight)) {
            HTML(paste0("<small><i>Tipping Point Weight = ", scales::percent(tp$weight, accuracy = 0.001), "</i></small>"))
          } else {
            HTML("<small><i>No tipping point found.</i></small>")
          }
        }
      })
    }
    
    if (input$model_type == "Bayesian Hierarchical Model") {
      
      sigma_grid <- seq(input$sigma_min, input$sigma_max, length.out = input$sigma_steps)
      tipping_results <- run_bayesian_tipping_analysis(
        ped_ctrl_resp = input$ped_ctrl_resp,
        ped_ctrl_total = input$ped_ctrl_total,
        ped_treat_resp = input$ped_treat_resp,
        ped_treat_total = input$ped_treat_total,
        adult_ctrl_resp = input$adult_ctrl_resp,
        adult_ctrl_total = input$adult_ctrl_total,
        adult_treat_resp = input$adult_treat_resp,
        adult_treat_total = input$adult_treat_total,
        sigma_values = sigma_grid
      )
      
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
        p <- plot_hierarchical_tipping(tipping_results, tipping_row_fda)
        
        # 关键一步：转为交互图
        ggplotly(p, tooltip = "text")
      })
      
      # error bar
      output$errorBarInfo <- renderUI({
        # 直接使用 tipping_results，因为它在外层已经跑过了
        combined_results <- tipping_results
        combined_results$CI_Exclude_0 <- combined_results$OR_lower > 1
        
        total_error_bars <- nrow(combined_results)
        significant_error_bars <- sum(combined_results$CI_Exclude_0)
        non_significant_error_bars <- total_error_bars - significant_error_bars
        
        HTML(paste0(
          "<div style='background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin: 10px 0;'>",
          "<strong>📊 <span style='color:#007bff;'>Error Bar Statistics:</span></strong><br>",
          "📈 Total Error Bars: <span style='color: #007bff;'><strong>", total_error_bars, "</strong></span><br>",
          "✅ Significant (CI excludes 0): <span style='color: #28a745;'><strong>", significant_error_bars, "</strong></span><br>",
          "❌ Non-significant (CI includes 0): <span style='color: #dc3545;'><strong>", non_significant_error_bars, "</strong></span><br>",
          "🎯 Sigma Range: ", round(min(tipping_results$fixed_sigma_alpha), 3), " - ", round(max(tipping_results$fixed_sigma_alpha), 3),
          "</div>"
        ))
      })
      
      
      output$essTableDetailed <- DT::renderDataTable({
        if (!is.null(tipping_results)) {
          display_data <- tipping_results %>%
            dplyr::select(
              `Sigma (σ)` = fixed_sigma_alpha,
              `CI Lower` = OR_lower,
              `CI Upper` = OR_upper,
              `ESS Treatment` = Borrowed_FDA_trt,
              `ESS Control` = Borrowed_FDA
            ) %>%
            dplyr::mutate(
              `CI Lower` = round(`CI Lower`, 3),
              `CI Upper` = round(`CI Upper`, 3),
              `ESS Treatment` = round(`ESS Treatment`, 1),
              `ESS Control` = round(`ESS Control`, 1)
            )
          
          datatable(
            display_data,
            options = list(
              pageLength = 25,
              scrollX = TRUE
            ),
            rownames = FALSE
          )
        } else {
          datatable(
            data.frame(Message = "Please run Power Prior analysis first"),
            options = list(dom = 't'),
            rownames = FALSE
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
    
    # 🔥 简化的Power Prior分析
    if (input$model_type == "Power Prior") {
      
      # 生成图表
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
      
      # Error Bar 信息
      output$errorBarInfo <- renderUI({
        results <- powerPriorResults()
        if (!is.null(results)) {
          analysis_results <- results$analysis
          combined_results <- merge(analysis_results$results, analysis_results$ess_results, by = "alpha")
          combined_results$Favorable <- combined_results$lower_ci > 0
          
          total_error_bars <- nrow(combined_results)
          significant_error_bars <- sum(combined_results$Favorable)
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
      
      # 生成结论
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
      
      # 🔥 根据实际列名显示重要数据，保留三位小数
      output$essTableDetailed <- DT::renderDataTable({
        results <- powerPriorResults()
        if (!is.null(results)) {
          analysis_results <- results$analysis
          combined_results <- merge(analysis_results$results, analysis_results$ess_results, by = "alpha")
          
          display_data <- combined_results %>%
            dplyr::select(
              `Alpha (α)` = alpha,
              `CI Lower` = lower_ci,
              `CI Upper` = upper_ci,
              `ESS Treatment` = drug_ess,
              `ESS Control` = placebo_ess
            ) %>%
            dplyr::mutate(
              `Alpha (α)` = scales::percent(`Alpha (α)`, accuracy = 0.001),
              `CI Lower` = round(`CI Lower`, 3),
              `CI Upper` = round(`CI Upper`, 3),
              `ESS Treatment` = round(`ESS Treatment`, 1),
              `ESS Control` = round(`ESS Control`, 1)
            )
          
          datatable(
            display_data,
            options = list(
              pageLength = 25,
              scrollX = TRUE
            ),
            rownames = FALSE
          )
        } else {
          datatable(
            data.frame(Message = "Please run Power Prior analysis first"),
            options = list(dom = 't'),
            rownames = FALSE
          )
        }
      })
      
      
      # 简化的注释
      output$essNote <- renderUI({
        results <- powerPriorResults()
        if (!is.null(results)) {
          tp <- results$detailed$tipping_point_summary$RR
          if (!is.null(tp) && !is.na(tp$weight)) {
            HTML(paste0("<small><i>Tipping Point Weight = ", scales::percent(tp$weight, accuracy = 0.001), "</i></small>"))
          } else {
            HTML("<small><i>No tipping point found.</i></small>")
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
  
  # Fisher's Exact Test 事件处理 - 切换模式
  observeEvent(input$fisher, {
    
    # 🔥 切换显示状态
    values$showFisher <- !values$showFisher
    
    if (values$showFisher) {
      # 只有在显示时才计算Fisher检验
      # 1. 构造2x2列联表
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
    }
  })
  
  # Fisher面板显示控制
  output$showFisherPanel <- reactive({
    values$showFisher
  })
  outputOptions(output, "showFisherPanel", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)


