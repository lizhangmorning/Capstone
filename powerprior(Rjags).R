find_tipping_point <- function(alpha_seq = seq(0, 0.02, by = 0.001), seed = 42) {
  library(rjags)
  
  # Pediatric data
  y_ped <- c(1, 15)
  n_ped <- c(6, 30)
  
  # Adult data
  y_adult <- c(45, 375)
  n_adult <- c(494, 509)
  
  # JAGS model (same as你之前的)
  model_string <- "
  model {
    for (i in 1:2) {
      y_ped[i] ~ dbin(p[i], n_ped[i])
      logit(p[i]) <- theta[i]
    }
    for (j in 1:n_placebo_hist) {
      y_placebo_hist[j] ~ dbern(p[1])
    }
    for (j in 1:n_druga_hist) {
      y_druga_hist[j] ~ dbern(p[2])
    }
    for (i in 1:2) {
      theta[i] ~ dnorm(0, 0.01)
    }
    logor <- theta[2] - theta[1]
    or <- exp(logor)
  }"
  
  result <- data.frame(alpha = alpha_seq, lower = NA, median = NA, upper = NA)
  
  for (i in seq_along(alpha_seq)) {
    a <- alpha_seq[i]
    
    # Simulate pseudo data
    set.seed(seed)
    pseudo_placebo <- rep(c(1, 0), times = c(y_adult[1], n_adult[1] - y_adult[1]))
    pseudo_drugA <- rep(c(1, 0), times = c(y_adult[2], n_adult[2] - y_adult[2]))
    
    n_sample <- round(a * length(pseudo_placebo))
    pseudo_placebo <- sample(pseudo_placebo, n_sample)
    pseudo_drugA <- sample(pseudo_drugA, n_sample)
    
    # JAGS data list
    jags_data <- list(
      y_ped = y_ped,
      n_ped = n_ped,
      y_placebo_hist = pseudo_placebo,
      y_druga_hist = pseudo_drugA,
      n_placebo_hist = length(pseudo_placebo),
      n_druga_hist = length(pseudo_drugA)
    )
    
    # Run model
    model <- jags.model(textConnection(model_string), data = jags_data,
                        n.chains = 3, n.adapt = 500, quiet = TRUE)
    update(model, 1000, progress.bar = "none")
    samples <- coda.samples(model, c("or"), n.iter = 3000, progress.bar = "none")
    
    or_q <- quantile(as.matrix(samples)[, "or"], c(0.025, 0.5, 0.975))
    result[i, c("lower", "median", "upper")] <- or_q
  }
  
  print(result)
  
  # 找出第一个 lower bound > 1 的 α
  tipping <- result$alpha[which(result$lower > 1)[1]]
  if (is.na(tipping)) {
    message("⚠️ No tipping point found (OR CrI always includes 1)")
    return(result)
  } else {
    message(sprintf("🎯 Tipping point found at alpha = %.2f", tipping))
    return(list(tipping_point = tipping, all_results = result))
  }
}

res <- find_tipping_point(alpha_seq = seq(0, 0.02, by = 0.001))
df <- res$all_results
tp <- res$tipping_point

tp
# 加载绘图包
library(ggplot2)

ggplot(df, aes(x = alpha, y = median)) +
  geom_point(size = 2, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.03, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = tp, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = tp + 0.05, y = 1.2, label = paste0("Tipping Point (", tp, ")"), color = "red", size = 4, hjust = 0) +
  annotate("text", x = 0.02, y = 1.1, label = "95% CrI includes 1", color = "red", size = 4, hjust = 0.5) +
  theme_gray(base_size = 14) +  # 灰色背景更接近你贴图风格
  labs(
    title = "",
    x = "Weight Parameter",
    y = "Odds Ratio"
  ) +
  coord_cartesian(ylim = c(0, max(df$upper) * 1.05))  # 适当限制 y 轴避免太高

