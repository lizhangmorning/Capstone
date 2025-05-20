library(shiny)
library(rjags)
library(here)
library(ggplot2)

ui <- fluidPage(
##Change the title
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
        tabPanel("📊 Tipping Point Plot",
                 plotOutput("tippingPlot"),
                 verbatimTextOutput("tippingConclusion")
        ),
        tabPanel("📈 ESS Calculation",
                 verbatimTextOutput("essResult")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  tipping_results <- reactiveVal(NULL)
  all_samples_list_alpha <- reactiveVal(NULL)
  tp_info <- reactiveVal(NULL)
  
  observeEvent(input$analyze, {
    if (input$model_type == "Bayesian Hierarchical Model") {
      
      combined_data <- data.frame(
        Source = rep(c("Pediatric", "Adult"), each = 2),
        Group = rep(c("DrugA", "Placebo"), 2),
        n = c(input$ped_treat_resp, input$ped_ctrl_resp,
              input$adult_treat_resp, input$adult_ctrl_resp),
        N = c(input$ped_treat_total, input$ped_ctrl_total,
              input$adult_treat_total, input$adult_ctrl_total)
      )
      
      combined_data <- combined_data |> 
        dplyr::mutate(
          group = ifelse(Group == "DrugA", 1, 0),
          trial_str = paste(Source, Group),
          trial = as.integer(factor(trial_str)),
          y = n
        )
      
      results <- data.frame()
      all_draws <- list()
      sigma_seq <- c(seq(2, 1, by = -0.1), seq(0.9, 0.3, by = -0.2))
      n_adult <- input$adult_ctrl_total + input$adult_treat_total
      p_adult <- input$adult_ctrl_resp / input$adult_ctrl_total
      info_per_subject <- n_adult * p_adult * (1 - p_adult)
      n_ped <- input$ped_ctrl_total
      
      for (s in sigma_seq) {
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
        
        model_file <- here::here("R", "models", "Bayesian_hierarchical_model", "rjags", "pooled_tipping.txt")
        params <- c("mu_alpha", "sigma_alpha", "alpha", "beta")
        jags_mod <- jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000)
        update(jags_mod, 1000)
        samples <- coda.samples(jags_mod, variable.names = params, n.iter = 5000)
        post <- as.matrix(samples)
        
        beta_draws <- post[, grep("^beta\\[", colnames(post))]
        alpha_draws <- post[, grep("^alpha\\[", colnames(post))]
        or_draws <- exp(beta_draws)
        or_ped <- or_draws[, "beta[2]"]
        
        V2 <- var(alpha_draws[, "alpha[2]"])
        
        jags_data_vague <- jags_data
        jags_data_vague$fixed_sigma_alpha <- 1e32
        jags_mod_vague <- jags.model(model_file, data = jags_data_vague, n.chains = 3, n.adapt = 1000)
        update(jags_mod_vague, 1000)
        samples_vague <- coda.samples(jags_mod_vague, variable.names = params, n.iter = 5000)
        post_vague <- as.matrix(samples_vague)
        V1 <- var(post_vague[, "alpha[2]"])
        
        ess_prior <- (1 / V2) / info_per_subject
        ess_fda <- n_ped * V1 / V2
        borrowed <- ess_fda - n_ped
        
        results <- rbind(results, data.frame(
          fixed_sigma_alpha = s,
          OR_median = median(or_ped),
          OR_lower = quantile(or_ped, 0.025),
          OR_upper = quantile(or_ped, 0.975),
          ESS_prior = ess_prior,
          ESS_FDA = ess_fda,
          Borrowed_FDA = borrowed
        ))
        
        all_draws[[as.character(s)]] <- alpha_draws
      }
      
      results$Significant_FDA <- results$OR_lower > 1
      tipping_results(results)
      all_samples_list_alpha(all_draws)
      
      tp <- results |> dplyr::filter(OR_lower <= 1) |> dplyr::slice(1)
      tp_info(tp)
    }
  })
  
  output$tippingPlot <- renderPlot({
    req(tipping_results())
    tipping_row <- tipping_results() |> dplyr::filter(Significant_FDA == TRUE) |> dplyr::slice_min(ESS_FDA)
    
    ggplot(tipping_results(), aes(x = fixed_sigma_alpha, y = OR_median, color = Significant_FDA)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      geom_vline(xintercept = tipping_row$fixed_sigma_alpha, linetype = "dotted", color = "red") +
      annotate("text",
               x = tipping_row$fixed_sigma_alpha * 0.6,
               y = tipping_row$OR_median + 78,
               label = paste0("Tipping Point = ", round(tipping_row$fixed_sigma_alpha, 3)),
               hjust = 0, vjust = 0,
               size = 4.5, fontface = "italic", color = "red") +
      annotate("text",
               x = tipping_row$fixed_sigma_alpha * 0.6,
               y = tipping_row$OR_median + 70,
               label = paste0("95% CI: [", round(tipping_row$OR_lower, 2), ", ", round(tipping_row$OR_upper, 2), "]"),
               hjust = 0, vjust = 0,
               size = 4.5, fontface = "italic", color = "red") +
      scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray")) +
      labs(
        title = "Pediatric OR vs Borrowing Strength (\u03c3)",
        x = "sigma alpha",
        y = "Odds Ratio (DrugA vs Placebo)"
      ) +
      theme_minimal(base_size = 15)
  })
  
  output$tippingConclusion <- renderPrint({
    req(tp_info())
    tp <- tp_info()
    
    if (nrow(tp) == 0) {
      "Tipping point not reached (OR CrI always > 1)"
    } else {
      paste0("Tipping point \u03c3 = ", tp$fixed_sigma_alpha,
             " | 95% CrI = [", round(tp$OR_lower, 2), ", ", round(tp$OR_upper, 2), "]")
    }
  })
  
  output$essResult <- renderPrint({
    req(tp_info())
    tp <- tp_info()
    if (nrow(tp) > 0) {
      paste0("ESS_prior = ", round(tp$ESS_prior, 2),
             ", FDA-style ESS = ", round(tp$ESS_FDA, 2),
             ", Borrowed ESS = ", round(tp$Borrowed_FDA, 2))
    } else {
      "ESS not computed (no tipping point found)"
    }
  })
}


shinyApp(ui = ui, server = server)
