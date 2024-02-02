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
library(magrittr)

# Install ffscbExtra
devtools::install_github("creutzml/ffscbExtra")
library(ffscbExtra)

# Directories
functions_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")


# Source Pittman's functions
source(file.path(main_dir, 
                 "functional_dffits/pittman_data_gen.R"))
#####################################################################



### Investigate the results
#####################################################################
load(file.path(data_path, "dffits_exp_mod1_df_11_23_23.RData"))

# sim_results_df2 %<>%
#   dplyr::mutate(
#     mcc = (n_tp*n_tn - n_fp*n_fn)/
#       sqrt((n_tp + n_fp)*(n_tp + n_fn)*(n_tn + n_fp)*(n_tn + n_fn)))
# 
sim_results_df_sum <- sim_results_df %>%
  # dplyr::bind_rows(sim_results_df2) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity) %>%
  # dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
  dplyr::group_by(n, lambda, n_inf, alpha_val, sp) %>%
  dplyr::summarise(avg_n_tp = mean(n_tp), 
                   avg_n_fp = mean(n_fp), 
                   avg_n_tn = mean(n_tn), 
                   avg_n_fn = mean(n_fn), 
                   avg_fpr = mean(fpr),
                   avg_sensitivity = mean(sensitivity),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   avg_specificity = mean(specificity),
                   sd_specificity = sd(specificity, na.rm = T),
                   avg_accuracy = mean(accuracy),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T), 
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T))



# # Pivot longer:
# sim_results_df_sum <- sim_results_df %>%
#   dplyr::bind_rows(sim_results_df2) %>%
#   dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
#   dplyr::group_by(n, lambda, n_inf, alpha_val, sp) %>%
#   dplyr::summarise(avg_sensitivity = mean(sensitivity), 
#                    sd_sensitivity = sd(sensitivity, na.rm = T),
#                    avg_specificity = mean(specificity),  
#                    sd_specificity = sd(specificity, na.rm = T),
#                    avg_accuracy = mean(accuracy),
#                    sd_accuracy = sd(accuracy, na.rm = T),
#                    avg_mcc = mean(mcc, na.rm = T),
#                    sd_mcc = sd(mcc, na.rm = T))

# Colorblind palette for plotting
cb_pallette <- c("#E69F00", "#0072B2", "black", "#006b4e")

# Take a look at some plots:
sim_results_df_sum %>%
  mutate(n = factor(n), 
         n_inf = factor(n_inf)) %>% 
         # lambda = lambda + 1,
         #alpha_val = factor(alpha_val)) %>%
  filter(sp != 0, alpha_val == 0.1) %>%
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
ggplot() +
  # geom_point(aes(x = sp, y = avg_accuracy, 
  #                color = alpha_val)) +
  geom_line(aes(x = sp, y = avg_mcc, color = n)) +
  # geom_errorbar(aes(x = sp, 
  #                   ymin = avg_accuracy - 2*sd_accuracy/sqrt(n), 
  #                   ymax = avg_accuracy + 2*sd_accuracy/sqrt(n),
  #                   color = alpha_val)) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`))) +
  # facet_grid(cols = vars(n), rows = vars(lam), 
  #            labeller = label_bquote(
  #              rows = lambda*" = "*.(lam), 
  #              cols = "n = "*.(n))) +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Sampling Point (t)", 
       y = "Avg. Matthew's Correlation Coefficient (MCC)", 
       color = expression(alpha*" Quantile")) +
  theme_bw() +
  theme(text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#####################################################################

sim_results_df_sum_r <- sim_results_df %>%
  # dplyr::bind_rows(sim_results_df2) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity,
                n_inf_p = n_inf/n) %>%
  # dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
  dplyr::group_by(n_inf_p, lambda, alpha_val, sp) %>%
  dplyr::summarise(avg_n_tp = mean(n_tp), 
                   avg_n_fp = mean(n_fp), 
                   avg_n_tn = mean(n_tn), 
                   avg_n_fn = mean(n_fn), 
                   avg_fpr = mean(fpr),
                   avg_sensitivity = mean(sensitivity),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   avg_specificity = mean(specificity),
                   sd_specificity = sd(specificity, na.rm = T),
                   avg_accuracy = mean(accuracy),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T), 
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T))

# Take a look at some plots:
sim_results_df_sum_r %>%
  # mutate(n_inf_p = factor(n_inf_p)) %>% 
  # lambda = lambda + 1,
  mutate(alpha_val = factor(alpha_val)) %>%
  filter(sp != 0) %>%
  rename(#"Proportion Influential" = "n_inf", 
         "lam" = "lambda") %>%
  ggplot() +
  # geom_point(aes(x = sp, y = avg_accuracy, 
  #                color = alpha_val)) +
  geom_line(aes(x = sp, y = avg_mcc, color = alpha_val)) +
  # geom_errorbar(aes(x = sp, 
  #                   ymin = avg_accuracy - 2*sd_accuracy/sqrt(n), 
  #                   ymax = avg_accuracy + 2*sd_accuracy/sqrt(n),
  #                   color = alpha_val)) +
  facet_grid(cols = vars(n_inf_p), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "Prop. Inf. = "*.(n_inf_p))) +
  # facet_grid(cols = vars(n), rows = vars(lam), 
  #            labeller = label_bquote(
  #              rows = lambda*" = "*.(lam), 
  #              cols = "n = "*.(n))) +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Sampling Point (t)", 
       y = "Avg. Matthew's Correlation Coefficient (MCC) ", 
       color = expression(alpha*" Quantile")) +
  theme_bw() +
  theme(text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
