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
  seq(2, 1, by = -0.2),
  seq(0.9, 0.3, by = -0.1),
  seq(0.29, 0.05, by = -0.01)
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
  jags_mod <- jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000)
  update(jags_mod, 1000)
  
  samples <- coda.samples(jags_mod, variable.names = params, n.iter = 5000)
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
  
  # Save results
  tipping_results <- rbind(tipping_results, data.frame(
    fixed_sigma_alpha = s,
    OR_median = median(or_ped),
    OR_lower  = quantile(or_ped, 0.025),
    OR_upper  = quantile(or_ped, 0.975),
    ESS_prior = ess_prior
  ))
}

# === Add significance column ===
tipping_results$Significant <- tipping_results$OR_lower > 1

# === Identify tipping point ESS ===
tipping_point_ess <- min(tipping_results$ESS_prior[tipping_results$Significant])
cat("Estimated tipping point ESS:", round(tipping_point_ess, 2), "\n")

# 找到 tipping point row
tipping_row <- tipping_results |>
  filter(Significant == TRUE) |>
  slice_min(ESS_prior)

ggplot(tipping_results, aes(x = ESS_prior, y = OR_median, color = Significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_vline(xintercept = tipping_row$ESS_prior, linetype = "dotted", color = "black") +
  annotate("text",
           x = tipping_row$ESS_prior * 1.5,  # 稍微右移一点
           y = tipping_row$OR_median + 70,   # 放到线上方
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