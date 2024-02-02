#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### December 13th, 2023                                           ###
###                                                               ###
###   Using Ryan Pittman's data generation methods to simulate    ###
### and investigate functional DFFITS. Updated the code (11/15/23)###
### to compare and save the results of TP, FP, TN, FN for each    ###
### sampling point on the domain. Edited the data generation      ###
### function to create a more influential observation across the. ###
### entire sampling domain (11/23/23). Recently realized an error ###
### in the computation of the Student t quantiles. We need to     ###
### a student t process, and thus multivariate student t          ###
### quantiles. Using `qmvt` from `mvtnorm` (12/13/2023).          ###
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

load(paste0("/Users/creutzml/Library/Mobile Documents/com~",
            "apple~CloudDocs/Documents/Dissertation/",
            "functional_data_analysis/sim_results/",
            "functional_dffits/mvt_cutoffs_sim_1_19_23.RData"))
#####################################################################




## Run a simulation for functional dffits
#####################################################################
## Fixed values for the loop
# {err_dist = "t"; n = 50; sigs = 1; dfs = 10;}

# Set a seed, so that we can run this twice, with two different sigma
# and get a direct comparison of the results
# set.seed(101911)

## Parameters to loop over:
n_sps <- c(100, 500, 1000)
n_iters <- 100
n <- c(10, 50, 100)
lambda <-  c(1, 1.5, 2)
n_infs <- c(1, 2, 3)
alphas <- c(0.10, 0.05, 0.025, 0.005)
sim_parms <- expand.grid(n, lambda, n_infs)
colnames(sim_parms) <- c("n", "lambda", "n_inf")

# Calculate the cutoffs ahead of time (no need to compute them 1000
# times for each simulation, when it's the same for each iteration)
# t quantiles
mvt_cutoffs <- expand.grid(n, n_sps, alphas) %>%
  as.data.frame() %>%
  dplyr::rename("n_obs" = Var1, 
                "n_sp" = Var2, 
                "alphas" = Var3) %>%
  dplyr::mutate(t_cuts = 0)

for (i in 1:nrow(mvt_cutoffs)) {
  a_val <- mvt_cutoffs$alphas[i]
  t_df <- mvt_cutoffs$n_obs[i] - 3
  scalar <- sqrt(2/(mvt_cutoffs$n_obs[i] - 2))
  n_sp <- mvt_cutoffs$n_sp[i]
  t_quants <- qmvt(p = 1 - a_val,
                   df = t_df,
                   tail = "lower.tail",
                   sigma = diag(t_df/(t_df - 2),
                                nrow = n_sp))
  cuts <- scalar*t_quants$quantile
  mvt_cutoffs$t_cuts[i] <- cuts
}

save(mvt_cutoffs,
     file = paste0("/Users/creutzml/Library/Mobile Documents/com~",
                   "apple~CloudDocs/Documents/Dissertation/",
                   "functional_data_analysis/sim_results/",
                   "functional_dffits/mvt_cutoffs_sim_1_19_23.RData"))


# Empty list to store simulation results
sim_results_list <- vector("list", nrow(sim_parms))

# Length of data frame for each simulation
df_length <- length(alphas)*n_iters

# Start the simulation
#####################################################################
# Loop over the simulation parameters
pb <- progress_bar$new(total = nrow(sim_parms)*n_iters)
for (p in 1:nrow(sim_parms)) {
  
  # Empty data frame for results
  sim_results_temp <- data.frame(
    n = rep(sim_parms$n[p], df_length),
    lambda = rep(sim_parms$lambda[p], df_length),
    n_inf = rep(sim_parms$n_inf[p], df_length),
    iter = rep(1:n_iters, each = length(alphas)),
    alpha_val = rep(alphas, n_iters),
    n_tp = vector("numeric", df_length), 
    n_fp = vector("numeric", df_length), 
    n_tn = vector("numeric", df_length), 
    n_fn = vector("numeric", df_length)
  )
  
  # Loop over the iterations
  for (iter in 1:n_iters) {
    
    ## Generate data:
    #################################################################
    # Number of observations and sampling points, 
    n_obs <- sim_parms$n[p]
    grid <- make_grid(n_sp, rangevals = c(0, 1))
    
    # Generate new functional data using Pittman's approach:
    newData <- GenerateFunctionalDataOut1(
      N = n_obs, length = n_sp, 
      lambda = sim_parms$lambda[p], 
      n_inf = sim_parms$n_inf[p]
    )
    
    # Which one is influential?
    if (sim_parms$lambda[p] != 1) {
      inf_temp <- (1:n_obs %in% newData$OutlierNumber)
    } else {
      inf_temp <- rep(FALSE, n_obs)
    }
    #################################################################
    
    
    ## Functional DFFITS (Michael)
    #################################################################
    ## And now, if we fit with concurrent regression and our bands
    fReg_list <- fRegress_concurrent(y_mat = newData3$Ydata,
                                     x_array = newData3$Xdata)

    ## Pointwise cutoff:
    # Df and mean scalar quantity
    t_df <- n_obs - 3
    scalar <- sqrt(2/(n_obs - 2))

    # t quantiles
    cutoffs <- mvt_cutoffs$t_cuts[mvt_cutoffs$n_obs == n_obs]

    # compare the dffits to the cutoffs
    fReg_dffits3 <- fReg_list$dffits_mat
    dff_test <- lapply(1:length(alphas), FUN = function(i) {
      
      # Get the alpha value
      alpha_temp <- alphas[i]
      t_cut <- cutoffs[i]
      
      # Check the upper bound and lower bound
      dffits_out <- colSums(abs(fReg_dffits) > t_cut) > 0 
      
      n_tp <- sum(dffits_out & inf_temp)
      n_fp <- sum(dffits_out & !inf_temp)
      n_tn <- sum(!dffits_out & !inf_temp)
      n_fn <- sum(!dffits_out & inf_temp)
      
      # Add alpha value to dffits_diag
      dffits_diag <- c(alpha_temp, n_tp, n_fp, n_tn, n_fn)
      
      return(dffits_diag)
    })

    # Combine the list results
    dff_results <- do.call(rbind, dff_test)
    ################################################################
    
    
    ## Functional extenrally studentized residuals (Michael)
    #################################################################
    # ## And now, if we fit with concurrent regression and our bands
    # fReg_list <- fRegress_concurrent(y_mat = newData$Ydata, 
    #                                  x_array = newData$Xdata)
    # 
    # ## Pointwise cutoff:
    # # Df and mean scalar quantity
    # t_df <- n_obs - 3
    # 
    # # t quantiles
    # cutoffs_h <- qt(1 - alphas, t_df)
    # cutoffs_l <- qt(alphas, t_df)
    # 
    # # compare the dffits to the cutoffs
    # fReg_esr <- fReg_list$ext_stu_res_mat
    # esr_test <- lapply(1:4, FUN = function(i) {
    #   # Check the upper bound and lower bound
    #   esr_high_mat <- fReg_esr > cutoffs_h[i]
    #   esr_high_low <- fReg_esr < cutoffs_l[i]
    #   esr_out <- esr_high_mat + esr_high_low
    #   
    #   # Count the total number of points out
    #   esr_out_tot <- colSums(esr_out)
    #   
    #   return(esr_out_tot)
    # })
    # 
    # # Combine the list results
    # esr_results <- do.call(cbind, esr_test)
    ################################################################
      
      
    ## Update data frame and progress
    #################################################################
    # Indices for rows
    s_idx <- (iter - 1)*length(alphas) + 1
    e_idx <- s_idx + length(alphas) - 1
    
    # Update data frame
    sim_results_temp[s_idx:e_idx, 5:9] <- dff_results
    
    # Update progress bar
    pb$tick()
    #################################################################
  }
  
  # Calculate sensitivity, specificity, accurracy, and mcc
  sim_results_temp2 <- sim_results_temp %>%
    mutate(sensitivity = n_tp/(n_tp + n_fn), 
           specificity = n_tn/(n_tn + n_fp), 
           accuracy = (n_tp + n_tn)/(n), 
           mcc = (n_tp*n_tn - n_fp*n_fn)/
             sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                    (n_tn + n_fp)*(n_tn + n_fn)))
  
  ## Save data frame to list:
  ###################################################################
  sim_results_list[[p]] <- sim_results_temp2
  
  # Remove objects for memory purposes
  rm(sim_results_temp, sim_results_temp2)
  ###################################################################
}

# # Create and save data frame from list
sim_results_df <- do.call(rbind, sim_results_list)
# sim_results_df2 <- do.call(rbind, sim_results_list[51:100])
save(sim_results_df,
     file = paste0("/Users/creutzml/Library/Mobile Documents",
                   "/com~apple~CloudDocs/Documents/Dissertation",
                   "/functional_data_analysis/sim_results",
                   "/functional_dffits/dffits_mvt_mod3_12_19_23.RData"))
# 
# rm(sim_results_list)
# sim_results_df_final <- rbind(sim_results_df, sim_results_df2)
# rm(sim_results_df2)

# sim_results_sum_final <- bind_rows(sim_results_df_sum, sim_results_df_sum2)
# rm(sim_results_df_sum2)
