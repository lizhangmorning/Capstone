library(here)
library(dplyr)
library(epitools)

source(here("R", "data_loading", "load_data.R"))

trial_data <- subset(pediatric_trials, Trial %in% c("Trial1", "Trial2"))
trial_data <- na.omit(trial_data)
summary_table <- trial_data |>
  group_by(Group) |>
  summarise(
    Success = sum(n),
    Failure = sum(N - n)
  )

# Create a summary table with the desired format
table_mat <- as.matrix(summary_table[match(c("DrugA", "Placebo"), summary_table$Group), 2:3])
rownames(table_mat) <- c("DrugA", "Placebo")
colnames(table_mat) <- c("Success", "Failure")

print(table_mat)
##        Success Failure
##DrugA        15      15
##Placebo       1       5

# Fisher's Exact Test
fisher.test(table_mat)
## p-value = 0.1962(> 0.05, not significant)
##CI95% = 0.4568545 to 251.5802684(contains 1, not significant)
##odds ratio = 4.8062(> 1, Drug A may be better than placebo)