#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### October 23rd, 2023                                            ###
###                                                               ###
### (Results of previous simulation)                              ###
###   Using Ryan Pittman's data generation methods to simulate    ###
### and investigate functional DFFITS.                            ###
#####################################################################




## Starting steps
#####################################################################
# Packages
library(mvtnorm)
library(tidyverse)
library(fda)
library(ffscb)
library(conformalInference.fd)
library(progress)

# Main directory:
main_dir <- file.path("/Users/creutzml/Library/Mobile Documents",
                      "/com~apple~CloudDocs/Documents/Dissertation",
                      "/functional_data_analysis/code")
source(file.path(main_dir, 
                 "fast_and_fair/make_regression_band_mc.R"))
source(file.path(main_dir, 
                 "fast_and_fair/make_band_mc.R"))
source(file.path(main_dir, 
                 "fast_and_fair/confidence_band_mc.R"))
source(file.path(main_dir, 
                 "fast_and_fair",
                 "make_concurrent_regression_band_mc_10_12_23.R"))
source(file.path(main_dir, 
                 "functional_dffits/fRegress_concurrent.R"))
source(file.path(main_dir, 
                 "functional_dffits/pittman_data_gen.R"))
#####################################################################



### Investigate the results
#####################################################################
load(paste0("/Users/creutzml/Library/Mobile Documents",
            "/com~apple~CloudDocs/Documents/Dissertation",
            "/functional_data_analysis/sim_results",
            "/functional_dffits/dffits_exp_mc_11_2_23.RData"))

# Pivot longer:
sim_results_df_longer <- sim_results_df %>%
  pivot_longer(cols = alpha_0.1000:alpha_0.0015, 
               names_to = "alpha", 
               values_to = "n_out")

# Take a look at some boxplots:
ggplot() +
  geom_boxplot(aes(fill = influential, y = n_out, x = alpha),
               data = sim_results_df_longer) +
  geom_hline(yintercept = 100, linetype = "dashed") +
  facet_grid(rows = vars(n), cols = vars(lambda)) +
  theme(axis.text.x = element_text(angle = 45))

## 99.7%ile cutoff:
sim_results_df_longer_997 <- sim_results_df_longer %>% 
  filter(alpha == "alpha_0.0015")
ggplot() +
  geom_boxplot(aes(x = influential, y = n_out),
               data = sim_results_df_longer_997) +
  # geom_hline(yintercept = 100, linetype = "dashed") +
  facet_grid(rows = vars(n), cols = vars(lambda))

# Summary
sim_results_df_sum_997 <- sim_results_df_longer_997 %>%
  mutate(n_out_over_30 = n_out > 30, 
         n_out_over_100 = n_out > 100) %>%
  pivot_longer(c(n_out_over_30, n_out_over_100), 
               names_to = "rule", values_to = "inf_est") %>%
  group_by(n, lambda, rule) %>%
  mutate(tp = case_when(
    influential & inf_est ~ 1, 
    TRUE ~ 0
  ), 
  tn = case_when(
    !influential & !inf_est ~ 1, 
    TRUE ~ 0
  ), 
  fp = case_when(
    !influential & inf_est ~ 1, 
    TRUE ~ 0
  ), 
  fn = case_when(
    influential & !inf_est ~ 1, 
    TRUE ~ 0
  )) %>%
  summarize(N = n(), 
            ntp = sum(tp), 
            ntn = sum(tn), 
            nfp = sum(fp), 
            nfn = sum(fn))
# log_model <- glm(influential ~ n + lambda + n_out, 
#                  data = sim_results_df_longer_997, 
#                  family = binomial)
# summary(log_model)
# sim_results_df_longer_997$log_preds <- log_model$fitted.values
# sim_results_df_longer_997 %<>% 
#   dplyr::arrange(n_out)
# plot(log_preds ~ n_out, 
#      data = sim_results_df_longer_997,
#      type = "l",
#      lwd = 2)

# Summary stats
sim_results_summary <- sim_results_df_longer %>%
  group_by(n, lambda, influential, alpha) %>%
  summarize(n_2.5 = quantile(n_out, probs = 0.025),
            n_mean = mean(n_out, na.rm = T),
            n_median = median(n_out, na.rm = T),
            n_97.5 = quantile(n_out, probs = 0.975))

# Plot of summary stats
ggplot() +
  geom_line(aes(x = n, y = n_median, color = influential), 
            data = sim_results_summary) + 
  geom_point(aes(x = n, y = n_median, color = influential), 
             data = sim_results_summary) +
  facet_grid(rows = vars(alpha), cols = vars(lambda))

# Differenes in counts boxplot... can I?
sim_results_diff <- sim_results_df_longer %>%
  pivot_wider(names_from = influential)
#####################################################################

