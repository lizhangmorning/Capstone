library(rjags)
library(coda)
library(here)
library(dplyr)

# --- Load Data ---
source(here("R", "data_loading", "load_data.R"))

# --- Data Preparation ---

# Pediatric pooled placebo: Trial1 placebo + pediatric_other_placebo
pediatric_placebo_n <- sum(pediatric_trials |> filter(Trial == "Trial1", Group == "Placebo") |> pull(n),
                           pediatric_other_placebo |> pull(n))
pediatric_placebo_N <- sum(pediatric_trials |> filter(Trial == "Trial1", Group == "Placebo") |> pull(N),
                           pediatric_other_placebo |> pull(N))

# Pediatric pooled DrugA (Parsabiv): Trial1 + Trial2 DrugA
pediatric_druga_n <- sum(pediatric_trials |> filter(Group == "DrugA", !is.na(n)) |> pull(n))
pediatric_druga_N <- sum(pediatric_trials |> filter(Group == "DrugA", !is.na(N)) |> pull(N))

# Adult pooled placebo
adult_placebo_n <- adult_trials |> filter(Trial == "Pooled", Group == "Placebo") |> pull(n)
adult_placebo_N <- adult_trials |> filter(Trial == "Pooled", Group == "Placebo") |> pull(N)

# Adult pooled DrugA (Parsabiv)
adult_druga_n <- adult_trials |> filter(Trial == "Pooled", Group == "DrugA") |> pull(n)
adult_druga_N <- adult_trials |> filter(Trial == "Pooled", Group == "DrugA") |> pull(N)

# --- Organize data for JAGS ---

y <- c(pediatric_placebo_n, pediatric_druga_n, adult_placebo_n, adult_druga_n)
N_total <- c(pediatric_placebo_N, pediatric_druga_N, adult_placebo_N, adult_druga_N)

# population: 0 = pediatric, 1 = adult
population <- c(0, 0, 1, 1)

# treatment group: 0 = placebo, 1 = DrugA
trt_group <- c(0, 1, 0, 1)

jags_base_data <- list(
  y = y,
  N = N_total,
  pop = population,
  trt = trt_group,
  n_groups = length(y)
)

# --- Tipping Point Analysis ---

fixed_sigmas <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20, 10000, 100000, 500000,10000000000000000000)
results <- data.frame(Sigma = fixed_sigmas, OR_median = NA, OR_2.5 = NA, OR_97.5 = NA,
                      RD_median = NA, RD_2.5 = NA, RD_97.5 = NA)

model_file <- here("R", "models", "Bayesian_hierarchical_model", "rjags", "pooled_tipping.txt")

for (i in seq_along(fixed_sigmas)) {
  jags_data <- c(jags_base_data, list(fixed_sigma = fixed_sigmas[i]))
  
  model <- jags.model(file = model_file, data = jags_data, n.chains = 3, quiet = TRUE)
  update(model, 1000)
  samples <- coda.samples(model, variable.names = c("beta", "mu"), n.iter = 5000)
  
  samples_matrix <- as.matrix(samples)
  p_pediatric_placebo <- samples_matrix[, "beta[1]"] + samples_matrix[, "mu[1,1]"]
  p_pediatric_druga <- samples_matrix[, "beta[2]"] + samples_matrix[, "mu[1,2]"]
  
  OR_pediatric <- exp(p_pediatric_druga - p_pediatric_placebo)
  inv_logit <- function(x) exp(x)/(1+exp(x))
  RD_pediatric <- inv_logit(p_pediatric_druga) - inv_logit(p_pediatric_placebo)
  
  results$OR_median[i] <- median(OR_pediatric)
  results$OR_2.5[i] <- quantile(OR_pediatric, 0.025)
  results$OR_97.5[i] <- quantile(OR_pediatric, 0.975)
  
  results$RD_median[i] <- median(RD_pediatric)
  results$RD_2.5[i] <- quantile(RD_pediatric, 0.025)
  results$RD_97.5[i] <- quantile(RD_pediatric, 0.975)
}

print(results)

library(ggplot2)

# --- Prepare Data ---
results_plot <- results |> 
  mutate(
    lower = OR_2.5,
    median = OR_median,
    upper = OR_97.5
  )

# --- Plot ---
ggplot(results_plot, aes(x = Sigma, y = median)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(trans = "log10", breaks = results$Sigma) +
  ylim(0, NA) +
  labs(
    title = "Tipping Point Analysis: Pediatric Odds Ratio vs Sigma",
    x = "Fixed Sigma (log scale)",
    y = "Odds Ratio (Parsabiv vs Placebo)"
  ) +
  theme_minimal(base_size = 14)
