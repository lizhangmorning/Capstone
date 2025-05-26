library(here)

pediatric_trials <- read.table(here("data", "pediatric_trials.txt"), header = TRUE)
adult_trials <- read.table(here("data", "adult_trials.txt"), header = TRUE)
pediatric_other_placebo <- read.table(here("data", "pediatric_other.txt"), header = TRUE)
