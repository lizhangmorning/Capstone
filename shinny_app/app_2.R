library(shiny)
library(rjags)
library(here)
library(ggplot2)
library(tidyr)
library(dplyr)           # ← 加上这个，才能用 mutate()
library(shinycssloaders)
library(plotly)


source(here::here("shinny_app", "hierarchical_tipping_function.R"))

ui <- fluidPage(
  titlePanel("Pediatric Extrapolation"),
  
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
      ),
        tags$hr(style = "border-top: 2px solid #aaa;"),
        h4("Settings"),
        strong("Baseline:"),
        actionButton("fisher", "Fisher's Exact Test"),
        selectInput("model_type", "Bayesian Analysis:",
                    choices = c("Mixture Prior", "Bayesian Hierarchical Model", "Power Prior")),
        actionButton("analyze", "Run Analysis")
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "output.analysisVisible == true",
        tabsetPanel(
          tabPanel("📊 Tipping Point Plot",
                   plotly::plotlyOutput("tippingPlot"),
                   verbatimTextOutput("tippingConclusion")
          ),
          tabPanel("📈 ESS Calculation",
                   withSpinner(tableOutput("essTable"))
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
  )   # ← 不要忘记这对 sidebarLayout 的闭合
)


server <- function(input, output, session) {
  
  values <- reactiveValues(analysisDone = FALSE, showFisher = FALSE)
  
  observeEvent(input$analyze, {
    
    if (input$model_type == "Bayesian Hierarchical Model") {
      values$showFisher <- FALSE
      values$analysisDone <- TRUE
      
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
        paste0("Tipping Point: sigma = ", round(tipping_row_fda$fixed_sigma_alpha, 2),
               "\n95% CI: [", round(tipping_row_fda$OR_lower, 2), ", ", round(tipping_row_fda$OR_upper, 2), "]")
      })
      
      output$tippingPlot <- plotly::renderPlotly({
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
          geom_vline(xintercept = tipping_row_fda$fixed_sigma_alpha, linetype = "dotted", color = "red") +
          annotate("text",
                   x = tipping_row_fda$fixed_sigma_alpha * 0.6,
                   y = tipping_row_fda$OR_median + 80,
                   label = paste0("Tipping Point = ", round(tipping_row_fda$fixed_sigma_alpha, 2)),
                   hjust = 0, vjust = 0, size = 4.5, fontface = "italic", color = "red") +
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
        plotly::ggplotly(p, tooltip = "text")
      })
      
      
      output$essTable <- renderTable({
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
          Group = c("Placebo", "Treatment", "Total"),
          `Borrowed ESS` = c(placebo_borrowed, treatment_borrowed, total_borrowed),
          `Total ESS` = c(placebo_total, treatment_total, total_total)
        )
      })
      
      output$analysisVisible <- reactive({
        values$analysisDone
      })
      outputOptions(output, "analysisVisible", suspendWhenHidden = FALSE)
      
      
      
      
    } # end if BHM
    
  }) # end observeEvent
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
  
  output$showFisherPanel <- reactive({
    values$showFisher
  })
  outputOptions(output, "showFisherPanel", suspendWhenHidden = FALSE)
  
}



shinyApp(ui = ui, server = server)