library(rjags)
library(coda)
library(here)
library(dplyr)
library(ggplot2)

# 载入数据
source(here("R", "data_loading", "load_data.R"))

# 儿科数据清洗
pediatric_clean <- pediatric_trials |>
  filter(!is.na(n) & !is.na(N)) |>
  filter(Trial %in% c("Pooled")) |>
  mutate(Source = "Pediatric")

# 成人数据清洗
adult_clean <- adult_trials |>
  filter(Trial %in% c("Pooled")) |>
  mutate(Source = "Adult")

# 合并数据，并生成 trial index
combined_data <- bind_rows(pediatric_clean, adult_clean) |>
  mutate(
    group = ifelse(Group == "DrugA", 1, 0),
    trial_str = paste(Source, Trial),
    trial = as.integer(factor(trial_str)),
    y = n
  )

tipping_results <- data.frame(
  fixed_sigma_alpha = numeric(),
  OR_median = numeric(),
  OR_lower = numeric(),
  OR_upper = numeric()
)
all_samples_list_alpha <- list()

for (s in seq(1, 2, by = 0.1)){
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
  
  model_file <- here("R", "models", "Bayesian_hierarchical_model", "rjags", "pooled_tipping.txt")
  jags_mod <- jags.model(model_file, data = jags_data, n.chains = 3, n.adapt = 1000)
  update(jags_mod, 1000)  # Burn-in
  
  samples <- coda.samples(jags_mod, variable.names = params, n.iter = 5000)
  
  post <- as.matrix(samples)
  beta_draws <- post[, grep("^beta\\[", colnames(post))]
  alpha_draws <- post[, grep("^alpha\\[", colnames(post))]
  
  or_draws <- exp(beta_draws)
  or_ped <- or_draws[, "beta[2]"] 
  
  all_samples_list_alpha[[as.character(s)]] <- alpha_draws
  tipping_results <- rbind(tipping_results, data.frame(
    fixed_sigma_alpha = s,
    OR_median = median(or_ped),
    OR_lower  = quantile(or_ped, 0.025),
    OR_upper  = quantile(or_ped, 0.975)
  ))
}

tipping_results <- tipping_results %>%
  dplyr::mutate(
    log_OR_median = log(OR_median),
    log_OR_lower = log(OR_lower),
    log_OR_upper = log(OR_upper)
  )

  
ggplot(tipping_results, aes(x = fixed_sigma_alpha, y = log_OR_median)) +
  geom_point(color = "black", size = 2) +  # 点始终是黑色
  geom_errorbar(
    aes(ymin = log_OR_lower, ymax = log_OR_upper, color = log_OR_lower < 0),
    width = 0.05
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black"),
    labels = c("FALSE" = "log(OR) ≥ 0", "TRUE" = "log(OR) < 0"),
    name = "Error Bar Direction"
  ) +
  labs(
    title = "Tipping Point Analysis",
    x = "Fixed Sigma",
    y = "log(Odds Ratio: Treatment vs Placebo)"
  ) +
  theme_minimal(base_size = 14)


ggplot(tipping_results, aes(x = fixed_sigma_alpha, y = OR_median)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.05) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(trans = "log10", breaks = tipping_results$fixed_sigma_beta) +
  scale_y_log10() +  # ✅ 限制y轴范围
  labs(
    title = "Tipping Point Analysis: Pediatric Odds Ratio vs Sigma",
    x = "Fixed Sigma (log scale)",
    y = "Odds Ratio (Parsabiv vs Placebo, log scale)"
  ) +
  theme_minimal(base_size = 14)




# 找出第一个 OR_lower <= 1 的点（CrI 不再显著优效）
tipping_point <- tipping_results |>
  filter(OR_lower <= 1) |>
  slice(1)  # 取第一个满足条件的

print(tipping_point)

sigma_tp <- tipping_point$fixed_sigma_alpha

# 获取对应 posterior draws
tp_draws <- all_samples_list_alpha[[as.character(sigma_tp)]]

# 方法 1: 计算 posterior variance
post_var <- var(tp_draws)
post_var_alpha2 <- post_var["alpha[2]", "alpha[2]"] 

# 方法 2: 如果 prior 是 N(0, σ²_prior), 则 prior variance：
prior_var <- tipping_point$fixed_sigma_alpha^2 # tau_beta = 1e6 → σ² = 1e-6

# 假设 adult total N:
n_adult <- sum(adult_clean$N)

# ESS based on variance shrinkage
ess <- (prior_var / post_var_alpha2) * n_adult
cat(sprintf("Estimated Effective Sample Size borrowed = %.1f\n", ess))

# 方法 3: MCMC Effective Size
ess_mcmc <- effectiveSize(tp_draws)
cat(sprintf("MCMC Effective Sample Size = %.1f\n", ess_mcmc))



