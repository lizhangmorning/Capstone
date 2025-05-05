# === Load libraries ===
library(rjags)
library(coda)
library(here)
library(dplyr)
library(ggplot2)
library(scales)

# === Load data ===
source(here("R", "data_loading", "load_data.R"))

# === Clean data ===
pediatric_clean <- pediatric_trials |>
  filter(!is.na(n) & !is.na(N)) |>
  filter(Trial == "Pooled") |>
  mutate(Source = "Pediatric")

adult_clean <- adult_trials |>
  filter(Trial == "Pooled") |>
  mutate(Source = "Adult")

combined_data <- bind_rows(pediatric_clean, adult_clean) |>
  mutate(
    group = ifelse(Group == "DrugA", 1, 0),
    trial_str = paste(Source, Trial),
    trial = as.integer(factor(trial_str)),
    y = n
  )

# === Set adult placebo info ===
n_adult <- 494
p_adult <- 0.091
info_per_subject <- n_adult * p_adult * (1 - p_adult)  # ~40.94

# === Prepare results container ===
tipping_results <- data.frame(
  fixed_sigma_alpha = numeric(),
  OR_median = numeric(),
  OR_lower = numeric(),
  OR_upper = numeric(),
  ESS_prior = numeric()
)

# === Loop over sigma_alpha values ===
sigma_values <- c(
  seq(2, 1, by = -0.1),
  seq(0.9, 0.3, by = -0.2)
)

for (s in sigma_values) {
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
  
  params <- c("mu_alpha", "sigma_alpha", "alpha", "beta")
  
  model_file <- here("R", "models", "Bayesian_hierarchical_model", "rjags", "baseline_pooled_tipping.txt")
  jags_mod <- jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000, inits = list(.RNG.name = "base::Super-Duper", .RNG.seed = 123))
  update(jags_mod, 1000)
  
  samples <- coda.samples(jags_mod, variable.names = params, n.iter = 50000)
  post <- as.matrix(samples)
  
  # Extract relevant posteriors
  beta_draws <- post[, grep("^beta\\[", colnames(post))]
  alpha_draws <- post[, grep("^alpha\\[", colnames(post))]
  or_draws <- exp(beta_draws)
  or_ped <- or_draws[, "beta[2]"]
  
  # ESS estimation: posterior variance of alpha[2]
  alpha2_samples <- alpha_draws[, "alpha[2]"]
  var_post <- var(alpha2_samples)
  ess_prior <- (1 / var_post) / info_per_subject
  
  # === FDA-style ESS estimation ==================================
  # Step 1: posterior variance with borrowing
  V2 <- var_post
  
  # Step 2: refit model with vague prior (no borrowing)
  jags_data_vague <- jags_data
  jags_data_vague$fixed_sigma_alpha <- 1e32  # vague = no borrowing
  
  jags_mod_vague <- jags.model(model_file, data = jags_data_vague, n.chains = 3, n.adapt = 1000, inits = list(.RNG.name = "base::Super-Duper", .RNG.seed = 123))
  update(jags_mod_vague, 1000)
  
  samples_vague <- coda.samples(jags_mod_vague, variable.names = params, n.iter = 5000)
  post_vague <- as.matrix(samples_vague)
  
  alpha2_vague <- post_vague[, "alpha[2]"]
  V1 <- var(alpha2_vague)
  
  # Step 3: pediatric **placebo** sample size only
  n_ped <- sum(pediatric_clean$N[pediatric_clean$Group == "Placebo"])
  
  # Step 4: FDA ESS calculation
  ess_fda <- n_ped * V1 / V2
  ess_borrowed <- ess_fda - n_ped
  
  # Save results
  tipping_results <- rbind(tipping_results, data.frame(
    fixed_sigma_alpha = s,
    OR_median = median(or_ped),
    OR_lower  = quantile(or_ped, 0.025),
    OR_upper  = quantile(or_ped, 0.975),
    ESS_prior = ess_prior,
    ESS_FDA   = ess_fda,
    Borrowed_FDA = ess_borrowed,
    V1 = V1,
    V2 = V2
  ))
}

# === Plot 1: your original OR vs ESS_prior ===
tipping_results$Significant <- tipping_results$OR_lower > 1

# === Identify tipping point ESS ===
tipping_point_ess <- min(tipping_results$ESS_prior[tipping_results$Significant])
cat("Estimated tipping point ESS:", round(tipping_point_ess, 2), "\n")

# find tipping point row
tipping_row <- tipping_results |>
  filter(Significant == TRUE) |>
  slice_min(ESS_prior)

ggplot(tipping_results, aes(x = ESS_prior, y = OR_median, color = Significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = tipping_row$ESS_prior, linetype = "dotted", color = "black") +
  annotate("text",
           x = tipping_row$ESS_prior * 1.5,  
           y = tipping_row$OR_median + 70,   
           label = paste0("Tipping ESS ≈ ", round(tipping_row$ESS_prior, 3)),
           hjust = 0, vjust = 0,
           size = 4.5, fontface = "italic") +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray")) +
  scale_x_continuous(
    trans = scales::log10_trans(),
    breaks = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
    labels = label_number(accuracy = 0.001)
  ) +
  labs(
    title = "Pediatric OR vs Effective Sample Size (ESS)",
    x = "Effective Sample Size (log scale)",
    y = "Odds Ratio (DrugA vs Placebo)"
  ) +
  theme_minimal(base_size = 15)

# === Plot 2: FDA-style ESS tipping analysis ===
tipping_results$Significant_FDA <- tipping_results$OR_lower > 1

tipping_point_ess_fda <- min(tipping_results$ESS_FDA[tipping_results$Significant_FDA])
cat("Estimated tipping point ESS (FDA):", round(tipping_point_ess_fda, 2), "\n")

tipping_row_fda <- tipping_results |>
  filter(Significant_FDA == TRUE) |>
  slice_min(ESS_FDA)

ggplot(tipping_results, aes(x = ESS_FDA, y = OR_median, color = Significant_FDA)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = tipping_row_fda$ESS_FDA, linetype = "dotted", color = "black") +
  annotate("text",
           x = tipping_row_fda$ESS_FDA * 1.5,  
           y = tipping_row_fda$OR_median + 70,   
           label = paste0("Tipping ESS (FDA) ≈ ", round(tipping_row_fda$ESS_FDA, 3)),
           hjust = 0, vjust = 0,
           size = 4.5, fontface = "italic") +
  scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray")) +
  scale_x_continuous(
    trans = scales::log10_trans(),
    breaks = c(0.1, 1, 2, 5, 10, 20, 50, 100),
    labels = label_number(accuracy = 1)
  ) +
  labs(
    title = "Pediatric OR vs FDA-Style Effective Sample Size",
    x = "FDA-Style Effective Sample Size (log scale)",
    y = "Odds Ratio (DrugA vs Placebo)"
  ) +
  theme_minimal(base_size = 15)