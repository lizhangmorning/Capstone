library(shiny)
library(plotly)
library(here)
library(ggplot2)
library(dplyr)
library(shinycssloaders)

# === Load analysis functions ===
source(here::here("shinny_app", "mixture_analysis_function.R"))
source(here::here("shinny_app", "hierarchical_tipping_function.R"))

ui <- fluidPage(
  titlePanel("Pediatric Extrapolation with Mixture Prior"),
  
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
      
      h4("Choose Model & Settings"),
      selectInput("model_type", "Select Method:",
                  choices = c("Mixture Prior", "Bayesian Hierarchical Model", "Power Prior")),
      hr(),
      actionButton("analyze", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("\U0001F4CA Tipping Point Plot",
                 uiOutput("dynamicPlotUI"),
                 verbatimTextOutput("tippingConclusion")
        ),
        tabPanel("\U0001F4C8 ESS Calculation",
                 withSpinner(tableOutput("essTable"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$dynamicPlotUI <- renderUI({
    if (input$model_type == "Mixture Prior") {
      withSpinner(plotlyOutput("tippingPlot"))
    } else if (input$model_type == "Bayesian Hierarchical Model") {
      withSpinner(plotOutput("tippingPlot"))
    } else {
      h4("Plot not available for this method.")
    }
  })
  
  observeEvent(input$analyze, {
    
    data_child <- list(
      treat = list(y = input$ped_treat_resp, n = input$ped_treat_total),
      control = list(y = input$ped_ctrl_resp, n = input$ped_ctrl_total)
    )
    data_adult <- list(
      treat = list(y = input$adult_treat_resp, n = input$adult_treat_total),
      control = list(y = input$adult_ctrl_resp, n = input$adult_ctrl_total)
    )
    
    if (input$model_type == "Mixture Prior") {
      output$tippingPlot <- renderPlotly({
        run_mixture_plot(data_child, data_adult)
      })
      
      output$tippingConclusion <- renderText({
        res <- mixture_analysis(data_child, data_adult)
        paste0("Tipping Point = ", round(res$tipping_point, 2), "\n",
               "ESS (Trt) = ", round(res$ess_treat, 1), ", ESS (Ctl) = ", round(res$ess_control, 1))
      })
      
      output$essTable <- renderTable({
        res <- mixture_analysis(data_child, data_adult)
        tibble::tibble(
          Group = c("Treatment", "Control"),
          ESS = c(round(res$ess_treat, 1), round(res$ess_control, 1))
        )
      })
    }
    
    if (input$model_type == "Bayesian Hierarchical Model") {
      pediatric_clean <- data.frame(
        Group = c("Placebo", "DrugA"),
        n = c(input$ped_ctrl_resp, input$ped_treat_resp),
        N = c(input$ped_ctrl_total, input$ped_treat_total),
        Trial = "Pooled") |>
        mutate(Source = "Pediatric")
      
      adult_clean <- data.frame(
        Group = c("Placebo", "DrugA"),
        n = c(input$adult_ctrl_resp, input$adult_treat_resp),
        N = c(input$adult_ctrl_total, input$adult_treat_total),
        Trial = "Pooled") |>
        mutate(Source = "Adult")
      
      combined_data <- bind_rows(pediatric_clean, adult_clean) |>
        mutate(
          group = ifelse(Group == "DrugA", 1, 0),
          trial_str = paste(Source, Trial),
          trial = as.integer(factor(trial_str)),
          y = n
        )
      
      tipping_results <- run_bayesian_tipping_analysis(combined_data)
      tipping_results$Significant_FDA <- tipping_results$OR_lower > 1
      tipping_row_fda <- tipping_results |>
        filter(Significant_FDA == TRUE) |>
        slice_min(ESS_FDA)
      
      output$tippingConclusion <- renderText({
        paste0("Tipping Point: sigma = ", round(tipping_row_fda$fixed_sigma_alpha, 2),
               "\n95% CI: [", round(tipping_row_fda$OR_lower, 2), ", ", round(tipping_row_fda$OR_upper, 2), "]")
      })
      
      output$tippingPlot <- renderPlot({
        ggplot(tipping_results, aes(x = fixed_sigma_alpha, y = OR_median, color = Significant_FDA)) +
          geom_point(size = 3) +
          geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
          geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
          geom_vline(xintercept = tipping_row_fda$fixed_sigma_alpha, linetype = "dotted", color = "red") +
          annotate("text", x = tipping_row_fda$fixed_sigma_alpha * 0.6,
                   y = tipping_row_fda$OR_median + 80,
                   label = paste0("Tipping Point = ", round(tipping_row_fda$fixed_sigma_alpha, 2)),
                   hjust = 0, vjust = 0, size = 4.5, fontface = "italic", color = "red") +
          scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray")) +
          scale_y_continuous(trans = scales::log10_trans(), breaks = c(1, 10, 50, 100)) +
          labs(title = "Pediatric OR vs Borrowing Strength (σ)",
               x = "sigma alpha", y = "Odds Ratio (DrugA vs Placebo)") +
          theme_minimal(base_size = 15)
      })
      
      output$essTable <- renderTable({
        tipping_row <- tipping_results |>
          filter(Significant_FDA == TRUE) |>
          slice_min(ESS_FDA)
        
        placebo_borrowed <- round(tipping_row$Borrowed_FDA, 1)
        placebo_total <- round(tipping_row$ESS_FDA, 1)
        
        treatment_borrowed <- round(tipping_row$Borrowed_FDA_trt, 1)
        treatment_total <- round(tipping_row$ESS_FDA_trt, 1)
        
        tibble::tibble(
          Group = c("Placebo", "Treatment", "Total"),
          `Borrowed ESS` = c(placebo_borrowed, treatment_borrowed, placebo_borrowed + treatment_borrowed),
          `Total ESS` = c(placebo_total, treatment_total, placebo_total + treatment_total)
        )
      })
    }
  })
}

shinyApp(ui = ui, server = server)
