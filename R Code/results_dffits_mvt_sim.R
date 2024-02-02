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
            "/functional_dffits/dffits_mvt_mod1_12_19_23.RData"))

# Pivot longer:
# sim_results_na_count <- sim_results_df %>%
#   dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
#   dplyr::group_by(n, lambda, n_inf, alpha_val) %>%
#   dplyr::summarize_all(~sum(is.na(.)))

sim_results_df_sum <- sim_results_df %>%
  dplyr::mutate(b_accuracy = (sensitivity + specificity)/2, 
                precision = n_tp/(n_tp + n_fp), 
                mcc = (n_tp*n_tn - n_fp*n_fn)/
                  sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                         (n_tn + n_fp)*(n_tn + n_fn))) %>%
  dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
  dplyr::group_by(n, lambda, n_inf, alpha_val) %>%
  dplyr::summarise(avg_sensitivity = mean(sensitivity, na.rm = T),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   na_sensitivity = sum(is.na(sensitivity)),
                   avg_specificity = mean(specificity, na.rm = T),
                   sd_specificity = sd(specificity, na.rm = T),
                   na_specificity = sum(is.na(specificity)),
                   avg_accuracy = mean(accuracy, na.rm = T),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   na_accuracy = sum(is.na(accuracy)),
                   avg_b_accuracy = mean(b_accuracy, na.rm = T),
                   sd_b_accuracy = sd(b_accuracy, na.rm = T), 
                   na_b_accuracy = sum(is.na(b_accuracy)),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T),
                   na_mcc = sum(is.na(mcc)), 
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T), 
                   na_precision = sum(is.na(precision)))

sim_results_df_sum_all <- sim_results_df %>%
  dplyr::mutate(b_accuracy = (sensitivity + specificity)/2, 
                precision = n_tp/(n_tp + n_fp), 
                mcc = (n_tp*n_tn - n_fp*n_fn)/
                  sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                         (n_tn + n_fp)*(n_tn + n_fn))) %>%
  dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
  dplyr::group_by(alpha_val) %>%
  dplyr::summarise(n_iters_tot = n(),
                   avg_sensitivity = mean(sensitivity, na.rm = T),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   na_sensitivity = sum(is.na(sensitivity)),
                   avg_specificity = mean(specificity, na.rm = T),
                   sd_specificity = sd(specificity, na.rm = T),
                   na_specificity = sum(is.na(specificity)),
                   avg_accuracy = mean(accuracy, na.rm = T),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   na_accuracy = sum(is.na(accuracy)),
                   avg_b_accuracy = mean(b_accuracy, na.rm = T),
                   sd_b_accuracy = sd(b_accuracy, na.rm = T), 
                   na_b_accuracy = sum(is.na(b_accuracy)),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T),
                   na_mcc = sum(is.na(mcc)), 
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T), 
                   na_precision = sum(is.na(precision)))

# Colorblind palette for plotting
cb_pallette <- c("#E69F00", "#0072B2", "black", "#006b4e")

# Take a look at some plots:
sim_results_df_sum %>%
  mutate(#n = factor(n), 
         n_inf = factor(n_inf), 
         alpha_val = factor(alpha_val)) %>%
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  filter(lam > 1) %>%
  # mutate(lam = factor(lam)) %>%
ggplot() +
  # geom_point(aes(x = sp, y = avg_accuracy, 
  #                color = alpha_val)) +
  geom_line(aes(x = n, y = avg_precision, color = alpha_val, 
                group = alpha_val)) +
  # geom_errorbar(aes(x = sp, 
  #                   ymin = avg_accuracy - 2*sd_accuracy/sqrt(n), 
  #                   ymax = avg_accuracy + 2*sd_accuracy/sqrt(n),
  #                   color = alpha_val)) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`))) +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(expand = c(0,0), 
                     breaks = c(10, 50, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Avg. Precision", 
       color = expression(alpha*" Quantile")) +
  theme_bw() +
  theme(text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1.5, "lines"))
#####################################################################
