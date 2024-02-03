#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### December 13th, 2023                                           ###
###                                                               ###
###   Simulation comparing the results of Pittman's DFFITS with   ###
### my own. Pittman's method is executed on smoothed and non      ###
### data, for a more complete comparison. (RESULTS)               ###
#####################################################################


## Starting steps
#####################################################################
# Packages
library(mvtnorm)
library(tidyverse)
library(fda)
library(ffscb)
# library(progress)
library(R.utils)
options(dplyr.summarise.inform = FALSE)
library(ggrepel)

# Install ffscbExtra package
devtools::install_github("creutzml/ffscbExtra")
library(ffcbExtra)

# Directories
functions_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0_T1000 <- sim_results %>%
  dplyr::mutate(alpha_b = 0, n_sp = 1000)

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0.25.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0.25_T1000 <- sim_results %>%
  dplyr::mutate(alpha_b = 0.25, n_sp = 1000)

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0.5.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0.5_T1000 <- sim_results %>%
  dplyr::mutate(alpha_b = 0.5, n_sp = 1000)

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0_T_100.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0_T100 <- sim_results %>%
  dplyr::mutate(alpha_b = 0, n_sp = 100)

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0.25_T_100.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0.25_T100 <- sim_results %>%
  dplyr::mutate(alpha_b = 0.25, n_sp = 100)

# Load in the simulation data
load(file = file.path(data_path, 
                      "dffits_creutz_v_pitt_sim_a_0.5_T_100.RData"))

# Rename, so that all three can be loaded in
sim_results_a_0.5_T100 <- sim_results %>%
  dplyr::mutate(alpha_b = 0.5, n_sp = 100)

# Clear space
rm(sim_results)
#####################################################################



## Summary tables
#####################################################################
# All simulation results
sim_results_all <- sim_results_a_0_T1000 %>%
  dplyr::bind_rows(sim_results_a_0.25_T1000) %>%
  dplyr::bind_rows(sim_results_a_0.5_T1000) %>%
  dplyr::bind_rows(sim_results_a_0_T100) %>%
  dplyr::bind_rows(sim_results_a_0.25_T100) %>%
  dplyr::bind_rows(sim_results_a_0.5_T100) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity)

# Sim results total summary
sim_results_sum_tot <- sim_results_all %>%
  # dplyr::select(-c(n_tp, n_fp, n_tn, n_fn)) %>%
  dplyr::group_by(method, alpha_val, alpha_b) %>%
  dplyr::summarise(n_sims = n(),
                   avg_runtime = mean(run_time, na.rm = T),
                   sd_runtime = sd(run_time, na.rm = T),
                   # avg_n_tp = mean(n_tp), 
                   # avg_n_fp = mean(n_fp), 
                   # avg_n_tn = mean(n_tn), 
                   # avg_n_fn = mean(n_fn), 
                   # avg_fpr = mean(fpr),
                   avg_sensitivity = mean(sensitivity, na.rm = T),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   na_sensitivity = sum(is.na(sensitivity)),
                   avg_specificity = mean(specificity, na.rm = T),
                   sd_specificity = sd(specificity, na.rm = T),
                   na_specificity = sum(is.na(specificity)),
                   avg_accuracy = mean(accuracy, na.rm = T),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   na_accuracy = sum(is.na(accuracy)),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T), 
                   na_mcc = sum(is.na(mcc)),
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T),
                   na_precision = sum(is.na(precision)))


## Table C.2 in Appendix C, Creutzinger (2024)
# Make a publishable form of the Pittman specific results
sim_results_summary_all_pretty <- sim_results_sum_tot %>%
  dplyr::mutate(alpha_val = factor(alpha_val), 
                method = case_when(
                  method == "Creutz" ~ "Theoretical (raw)", 
                  method == "Creutz_Smooth" ~ "Theoretical (smooth)",
                  method == "Pittman" ~ "Bootstrapped (raw)",
                  method == "Pittman_Smooth" ~ "Bootstrapped (smooth)"
                )) %>%
  # dplyr::filter(method %in% c("Bootstrapped (raw)", 
  #                             "Bootstrapped (smooth)")) %>%
  dplyr::filter(alpha_b == 0.5) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  dplyr::mutate(across(where(is.numeric), format, nsmall = 3)) %>%
  mutate(`Run Times (s)` = paste0(avg_runtime, " (", 
                                  sd_runtime, ") "), 
         `Sensitivity` = paste0(avg_sensitivity, " (", 
                                sd_sensitivity, ") "), 
         `Specificity` = paste0(avg_specificity, " (", 
                                sd_specificity, ") "), 
         `Accuracy` = paste0(avg_accuracy, " (", 
                             sd_accuracy, ") "), 
         `Precision` = paste0(avg_precision, " (", 
                              sd_precision, ") "), 
         `MCC` = paste0(avg_mcc, " (", 
                        sd_mcc, ") ")) %>%
  dplyr::select(c(alpha_val, method, `Run Times (s)`, 
                  Sensitivity, Specificity, Accuracy, 
                  Precision, MCC)) %>%
  arrange(alpha_val) #%>%
  # dplyr::ungroup() %>%
  # dplyr::select(-alpha_val)

# Produce latex code
kableExtra::kbl(sim_results_summary_all_pretty, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))

## Results show that Pittman (smooth) with alpha_b = 0.500 and 
## Creutzinger perform the best for each methodologies in the sim
## Next steps focus in on those results for the main body of the pape

## Table 4.1 in Chapter 4, Creutzinger (2024)
# Filter out the lesser results from the total table:
sim_results_sum_tot_main <- sim_results_sum_tot %>%
  dplyr::filter(method %in% c("Creutz", "Pittman_Smooth"), 
                alpha_b == 0.5)

sim_results_summary_all_main_pretty <- sim_results_sum_tot_main %>%
  dplyr::mutate(alpha_val = factor(alpha_val), 
                method = case_when(
                  method == "Creutz" ~ "Theoretical (raw)", 
                  method == "Pittman_Smooth" ~ "Bootstrapped (smooth)"
                )) %>%
  dplyr::mutate(across(where(is.numeric), round, 3)) %>%
  dplyr::mutate(across(where(is.numeric), format, nsmall = 3)) %>%
  mutate(`Run Times (s)` = paste0(avg_runtime, " (", 
                                  sd_runtime, ") "), 
         `Sensitivity` = paste0(avg_sensitivity, " (", 
                                sd_sensitivity, ") "), 
         `Specificity` = paste0(avg_specificity, " (", 
                                sd_specificity, ") "), 
         `Accuracy` = paste0(avg_accuracy, " (", 
                             sd_accuracy, ") "), 
         `Precision` = paste0(avg_precision, " (", 
                              sd_precision, ") "), 
         `MCC` = paste0(avg_mcc, " (", 
                        sd_mcc, ") ")) %>%
  dplyr::select(c(method, alpha_val, `Run Times (s)`, 
                  Sensitivity, Specificity, Accuracy, 
                  Precision, MCC)) %>%
  arrange(method, alpha_val) %>%
  ungroup() %>%
  dplyr::select(-c(method,`Run Times (s)`))

# Produce latex code
kableExtra::kbl(sim_results_summary_all_main_pretty, 
                booktabs = TRUE, 
                format = "latex") %>%
  kableExtra::kable_styling(latex_options = c("scale_down", 
                                              "hold_position"))


# Sim results summary for top performing methods for Creutzinger and
# Pittman
sim_results_sum <- sim_results_all %>%
  dplyr::filter(method %in% c("Creutz", "Pittman_Smooth"), 
                alpha_b == 0.5, alpha_val == 0.005) %>%
  dplyr::select(-c(alpha_b, alpha_val)) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity) %>%
  dplyr::mutate(method = case_when(
    method == "Creutz" ~ "Theoretical (raw)", 
    method == "Pittman_Smooth" ~ "Bootstrapped (smooth)"
  )) %>%
  dplyr::group_by(method, n, n_sp, lambda, n_inf) %>%
  dplyr::summarise(#avg_n_tp = mean(n_tp), 
                   # avg_n_fp = mean(n_fp), 
                   # avg_n_tn = mean(n_tn), 
                   # avg_n_fn = mean(n_fn), 
                   # avg_fpr = mean(fpr),
                   avg_runtime = mean(run_time, na.rm = T),
                   avg_sensitivity = mean(sensitivity, na.rm = T),
                   sd_sensitivity = sd(sensitivity, na.rm = T),
                   na_sensitivity = sum(is.na(sensitivity)),
                   avg_specificity = mean(specificity, na.rm = T),
                   sd_specificity = sd(specificity, na.rm = T),
                   na_specificity = sum(is.na(specificity)),
                   avg_accuracy = mean(accuracy),
                   sd_accuracy = sd(accuracy, na.rm = T),
                   na_accuracy = sum(is.na(accuracy)),
                   avg_mcc = mean(mcc, na.rm = T),
                   sd_mcc = sd(mcc, na.rm = T), 
                   na_mcc = sum(is.na(mcc)),
                   avg_precision = mean(precision, na.rm = T), 
                   sd_precision = sd(precision, na.rm = T),
                   na_precision = sum(is.na(precision)))
#####################################################################



### Figures 4.3 in Chapter 4, C.10 and C.12 in Appendix C, 
### Creutzinger (2024)
### Summary Figures over the sample size on the x axis
#####################################################################
sim_results_sum_n <- sim_results_all %>%
  dplyr::filter(method %in% c("Creutz", "Pittman_Smooth"), 
                alpha_b == 0.5, alpha_val == 0.005) %>%
  dplyr::select(-c(alpha_b, alpha_val)) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity) %>%
  dplyr::mutate(method = case_when(
    method == "Creutz" ~ "Theoretical (raw)", 
    method == "Pittman_Smooth" ~ "Bootstrapped (smooth)"
  )) %>%
  dplyr::group_by(method, n, lambda, n_inf) %>%
  dplyr::summarise(#avg_n_tp = mean(n_tp), 
    # avg_n_fp = mean(n_fp), 
    # avg_n_tn = mean(n_tn), 
    # avg_n_fn = mean(n_fn), 
    # avg_fpr = mean(fpr),
    n_sims_tot = n(),
    avg_runtime = mean(run_time, na.rm = T),
    avg_sensitivity = mean(sensitivity, na.rm = T),
    sd_sensitivity = sd(sensitivity, na.rm = T),
    na_sensitivity = sum(is.na(sensitivity)),
    avg_specificity = mean(specificity, na.rm = T),
    sd_specificity = sd(specificity, na.rm = T),
    na_specificity = sum(is.na(specificity)),
    avg_accuracy = mean(accuracy, na.rm = T),
    sd_accuracy = sd(accuracy, na.rm = T),
    na_accuracy = sum(is.na(accuracy)),
    avg_mcc = mean(mcc, na.rm = T),
    sd_mcc = sd(mcc, na.rm = T), 
    na_mcc = sum(is.na(mcc)),
    avg_precision = mean(precision, na.rm = T), 
    sd_precision = sd(precision, na.rm = T),
    na_precision = sum(is.na(precision)))

# Colorblind palette for plotting
cb_pallette <- c("#E69F00", "#0072B2", "black", "#006b4e")

## Take a look at some plots:
# Average Specificity
sim_results_sum_n %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  ggplot() +
  geom_line(aes(x = n, 
                y = avg_specificity, 
                color = method, 
                linetype = method), 
            linewidth = 1.5) +
  geom_point(aes(x = n, 
                 y = avg_specificity, 
                 color = method, 
                 shape = method), 
            size = 5) +
  # ggrepel::geom_label_repel(aes(x = n, 
  #                               y = avg_specificity, 
  #                               label = round(avg_specificity, 3)), 
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`)), 
             scales = "free") +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(10, 50, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Avg. Specificity", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 18), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1.5, "lines"), 
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text = element_text(size = 16)) 

# Average Sensitivity
sim_results_sum_n %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  filter(lam != 1) %>%
  ggplot() +
  geom_line(aes(x = n, 
                y = avg_sensitivity, 
                color = method, 
                linetype = method), 
            size = 1.5) +
  geom_point(aes(x = n, 
                 y = avg_sensitivity, 
                 color = method, 
                 shape = method), 
             size = 5) +
  # ggrepel::geom_label_repel(aes(x = n, 
  #                               y = avg_sensitivity, 
  #                               label = round(avg_sensitivity, 3)), 
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`))) +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(10, 50, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Avg. Sensitivity", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1.5, "lines"), 
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size=5)))

# Average MCC
sim_results_sum_n %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  filter(lam != 1, n != 10) %>%
  ggplot() +
  geom_line(aes(x = n, 
                y = avg_mcc, 
                color = method, 
                linetype = method), 
            size = 1.5) +
  geom_point(aes(x = n, 
                 y = avg_mcc, 
                 color = method, 
                 shape = method), 
             size = 5) +
  # ggrepel::geom_label_repel(aes(x = n, 
  #                               y = avg_mcc, 
  #                               label = round(avg_mcc, 3)), 
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`))) +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(10, 50, 100)) +
  labs(x = "Sample Size (n)", 
       y = "Avg. Matthew's Correlation Coefficient (MCC)", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing = unit(1.5, "lines"), 
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        strip.text = element_text(size = 16)) +
  guides(color = guide_legend(override.aes = list(size=5)))
#####################################################################



### Figures 4.4 in Chapter 4, C.11 and C.13 in Appendix C, 
### Creutzinger (2024)
### Looking more specifically at the results over the number of 
### sampling points in the domain
#####################################################################
sim_results_sum_T <- sim_results_all %>%
  dplyr::filter(method %in% c("Creutz", "Pittman_Smooth"), 
                alpha_b == 0.5, alpha_val == 0.005) %>%
  dplyr::select(-c(alpha_b, alpha_val)) %>%
  dplyr::mutate(precision = n_tp/(n_tp + n_fp)) %>%
  dplyr::mutate(fpr = 1 - specificity) %>%
  dplyr::mutate(method = case_when(
    method == "Creutz" ~ "Theoretical (raw)", 
    method == "Pittman_Smooth" ~ "Bootstrapped (smooth)"
  )) %>%
  dplyr::group_by(method, n_sp, lambda, n_inf) %>%
  dplyr::summarise(#avg_n_tp = mean(n_tp), 
    # avg_n_fp = mean(n_fp), 
    # avg_n_tn = mean(n_tn), 
    # avg_n_fn = mean(n_fn), 
    # avg_fpr = mean(fpr),
    n_sims_tot = n(),
    avg_runtime = mean(run_time, na.rm = T),
    avg_sensitivity = mean(sensitivity, na.rm = T),
    sd_sensitivity = sd(sensitivity, na.rm = T),
    na_sensitivity = sum(is.na(sensitivity)),
    avg_specificity = mean(specificity, na.rm = T),
    sd_specificity = sd(specificity, na.rm = T),
    na_specificity = sum(is.na(specificity)),
    avg_accuracy = mean(accuracy, na.rm = T),
    sd_accuracy = sd(accuracy, na.rm = T),
    na_accuracy = sum(is.na(accuracy)),
    avg_mcc = mean(mcc, na.rm = T),
    sd_mcc = sd(mcc, na.rm = T), 
    na_mcc = sum(is.na(mcc)),
    avg_precision = mean(precision, na.rm = T), 
    sd_precision = sd(precision, na.rm = T),
    na_precision = sum(is.na(precision)))


# Colorblind palette for plotting
cb_pallette <- c("#E69F00", "#0072B2", "black", "#006b4e")

## Take a look at some plots:
# Average Specificity
sim_results_sum_T %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  # filter(lam != 1) %>%
  ggplot() +
  geom_line(aes(x = n_sp, 
                y = avg_specificity, 
                color = method, 
                linetype = method), 
            linewidth = 1.5) +
  geom_point(aes(x = n_sp, 
                 y = avg_specificity, 
                 color = method, 
                 shape = method), 
             size = 5) +
  # ggrepel::geom_label_repel(aes(x = n_sp, 
  #                               y = avg_specificity, 
  #                               label = round(avg_specificity, 3)), 
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`)), 
             scales = "free") +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(100, 1000)) +
  labs(x = "Total Sampling Points (T)", 
       y = "Avg. Specificity", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 18), 
    plot.margin = unit(c(.5, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.spacing = unit(1.5, "lines"), 
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(), 
    legend.position = "bottom", 
    strip.text = element_text(size = 16)) 

# Average Sensitivity
sim_results_sum_T %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  filter(lam != 1) %>%
  ggplot() +
  geom_line(aes(x = n_sp, 
                y = avg_sensitivity, 
                color = method, 
                linetype = method), 
            linewidth = 1.5) +
  geom_point(aes(x = n_sp, 
                 y = avg_sensitivity, 
                 color = method, 
                 shape = method), 
             size = 5) +
  # ggrepel::geom_label_repel(aes(x = n_sp, 
  #                               y = avg_sensitivity, 
  #                               label = round(avg_sensitivity, 3)), 
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`)), 
             scales = "free") +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(100, 1000)) +
  labs(x = "Total Sampling Points (T)", 
       y = "Avg. Sensitivity", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 18), 
    plot.margin = unit(c(.5, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.spacing = unit(1.5, "lines"), 
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(), 
    legend.position = "bottom", 
    strip.text = element_text(size = 16)) 

# Average MCC
sim_results_sum_T %>%
  mutate(n_inf = factor(n_inf)) %>% 
  rename("# Influential" = "n_inf", 
         "lam" = "lambda") %>%
  filter(lam != 1) %>%
  ggplot() +
  geom_line(aes(x = n_sp, 
                y = avg_mcc, 
                color = method, 
                linetype = method), 
            linewidth = 1.5) +
  geom_point(aes(x = n_sp, 
                 y = avg_mcc, 
                 color = method, 
                 shape = method), 
             size = 5) +
  # ggrepel::geom_label_repel(aes(x = n_sp, 
  #                               y = avg_mcc, 
  #                               label = round(avg_mcc, 3), 
  #                               color = method),
  #                           size = 7) +
  facet_grid(cols = vars(`# Influential`), rows = vars(lam), 
             labeller = label_bquote(
               rows = lambda*" = "*.(lam), 
               cols = "# Influential = "*.(`# Influential`)), 
             scales = "free") +
  scale_color_manual(values = cb_pallette) +
  scale_x_continuous(breaks = c(100, 1000)) +
  labs(x = "Total Sampling Points (T)", 
       y = "Avg. Matthew's Correlation Coefficient (MCC)", 
       color = "Method", 
       shape = "Method", 
       linetype = "Method") +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 18), 
    plot.margin = unit(c(.5, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.spacing = unit(1.5, "lines"), 
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(), 
    legend.position = "bottom", 
    strip.text = element_text(size = 16)) 
#####################################################################