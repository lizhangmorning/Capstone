library(rjags)
library(coda)
library(here)
library(dplyr)

# Load data
source(here("R", "data_loading", "load_data.R"))

# pediatric data cleaning
pediatric_clean <- pediatric_trials |>
  filter(!is.na(n) & !is.na(N)) |>
  filter(Trial %in% c("Pooled")) |>
  mutate(Source = "Pediatric")

# adult data cleaning
adult_clean <- adult_trials |>
  filter(Trial %in% c("Pooled")) |>
  mutate(Source = "Adult")

# combined data and generate trial index
combined_data <- bind_rows(pediatric_clean, adult_clean) |>
  mutate(
    group = ifelse(Group == "DrugA", 1, 0),
    trial_str = paste(Source, Trial),
    trial = as.integer(factor(trial_str)),
    y = n
  )

# JAGS data list
jags_data <- list(
  y = combined_data$y,
  n = combined_data$N,
  group = combined_data$group,
  trial = combined_data$trial,
  N = nrow(combined_data),
  J = length(unique(combined_data$trial)),
  mu_alpha_mean = 0,
  mu_alpha_prec = 0.01,
  sigma_alpha_max = 10
)

# Show data structure
str(jags_data)

# parameters to monitor
params <- c("mu_alpha", "sigma_alpha", "beta")

# Load model
model_file <- here("R", "models", "Bayesian_hierarchical_model", "rjags", "model_hierarchical_baseline.txt")
jags_mod <- jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000)
update(jags_mod, 1000)  # Burn-in

# MCMC sampling
samples <- coda.samples(jags_mod, variable.names = params, n.iter = 5000)

# Posterior summary
summary(samples)

# posterior draws
post <- as.matrix(samples)
beta_draws <- post[, grep("^beta\\[", colnames(post))]

# OR draws
or_draws <- exp(beta_draws)

# OR summary
or_summary <- apply(or_draws, 2, function(x) {
  c(mean = mean(x),
    median = median(x),
    lower = quantile(x, 0.025),
    upper = quantile(x, 0.975),
    prob_superior = mean(x > 1))
})

or_summary <- t(or_summary)
rownames(or_summary) <- c("Adult", "Pediatric")
print(round(or_summary, 3))



