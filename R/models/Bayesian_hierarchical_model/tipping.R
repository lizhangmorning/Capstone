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


#新增——————————————————————————————————————————
# Adult placebo 信息量
n_adult <- 494
p_adult <- 0.091
info_per_subject <- n_adult * p_adult * (1 - p_adult)

# 新建 ESS 向量
ess_prior_vec <- c()
#新增——————————————————————————————————————————————





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
  
  model_file <- here("R", "models", "Bayesian_hierarchical_model", "rjags", "baseline_pooled_tipping.txt")
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
  

# 添加显著性判断列
tipping_results$Significant <- tipping_results$OR_lower > 1

# 绘图
ggplot(tipping_results, aes(x = fixed_sigma_alpha, y = OR_median, color = Significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper), width = 0.1, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_log10(
    breaks = unique(tipping_results$fixed_sigma_alpha),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_color_manual(
    values = c("TRUE" = "blue", "FALSE" = "gray"),
    labels = c("Significant", "Not Significant"),
    name = "95% CrI excludes 1"
  ) +
  ylim(0, NA) +
  labs(
    title = "Tipping Point Analysis: Pediatric OR vs Borrowing Strength",
    x = expression(Sigma[alpha]~"(log scale)"),
    y = "Odds Ratio (DrugA vs Placebo)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )


ess_beta2 <- effectiveSize(samples)[["beta[2]"]]
print(ess_beta2)


ess_selected <- effectiveSize(samples)[c("beta[1]", "beta[2]", "mu_alpha", "sigma_alpha")]
print(ess_selected)

# ————————————————————————————————————————————————等值算法——————————————————————
# Adult placebo arm
n_adult <- 494
p_adult <- 0.091
var_likelihood <- 1 / (n_adult * p_adult * (1 - p_adult))

# 一系列 sigma_alpha
sigma_alpha_seq <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)

# 计算 prior variance 和 ESS
ess_table <- data.frame(
  sigma_alpha = sigma_alpha_seq,
  prior_variance = sigma_alpha_seq^2,
  ESS_prior = var_likelihood / (sigma_alpha_seq^2)
)

print(ess_table)

#————————————————————————————————————————————————————————

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



