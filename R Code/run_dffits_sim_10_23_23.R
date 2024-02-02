#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### October 23rd, 2023                                            ###
###                                                               ###
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




## Run a simulation for functional dffits
#####################################################################
## Fixed values for the loop
# {err_dist = "t"; n = 50; sigs = 1; dfs = 10;}

# Set a seed, so that we can run this twice, with two different sigma
# and get a direct comparison of the results
# set.seed(101911)

## Parameters to loop over:
n_iters <- 1000
n <- c(10, 20, 50, 100)
lambda <- c(0.01, 0.5, 1, 1.5, 2) #c(.5,.75,.9,1,1.1,1.25,1.5,1.75,2) # -> use fewer
n_infs <- c(1, 2, 3, 4, 5)
alphas <- c(0.05, 0.025, 0.005)
sim_parms <- expand.grid(n, lambda, n_infs)
colnames(sim_parms) <- c("n", "lambda", "n_inf")

# Empty list to store simulation results
sim_results_list <- vector("list", nrow(sim_parms))

# Start the simulation
#####################################################################
# Loop over the simulation parameters
pb <- progress_bar$new(total = nrow(sim_parms)*n_iters)
for (p in 1:nrow(sim_parms)) {
  
  # Empty data frame for results
  sim_results_temp <- data.frame(
    n = rep(sim_parms$n[p], n_iters*sim_parms$n[p]),
    lambda = rep(sim_parms$lambda[p], n_iters*sim_parms$n[p]),
    n_inf = rep(sim_parms$n_inf[p], n_iters*sim_parms$n[p]),
    alpha_0.1000 = vector("numeric", n_iters*sim_parms$n[p]),
    alpha_0.0500 = vector("numeric", n_iters*sim_parms$n[p]),
    alpha_0.0100 = vector("numeric", n_iters*sim_parms$n[p]),
    alpha_0.0015 = vector("numeric", n_iters*sim_parms$n[p]), 
    influential = vector("logical", n_iters*sim_parms$n[p])
  )
  
  # Loop over the iterations
  for (iter in 1:n_iters) {
    
    ## Generate data:
    #################################################################
    # Number of observations and sampling points, 
    n_obs <- sim_parms$n[p]
    n_sp <- 1000
    grid <- make_grid(n_sp, rangevals = c(0, 1))
    
    # Generate new functional data using Pittman's approach:
    newData <- GenerateFunctionalDataOut(
      N = n_obs, length = n_sp, lambda = sim_parms$lambda[p], 
      n_inf = sim_parms$n_inf[p]
    )
    
    # Which one is influential?
    inf_temp <- (1:n_obs %in% newData$OutlierNumber)
    #################################################################
    
    
    ## Functional DFFITS (Michael)
    #################################################################
    ## And now, if we fit with concurrent regression and our bands
    fReg_list <- fRegress_concurrent(y_mat = newData$Ydata,
                                     x_array = newData$Xdata)

    ## Pointwise cutoff:
    # Df and mean scalar quantity
    t_df <- n_obs - 3
    scalar <- sqrt(2/(n_obs - 2))

    # t quantiles
    cutoffs_h <- scalar*qt(1 - alphas, t_df)
    cutoffs_l <- scalar*qt(alphas, t_df)

    # compare the dffits to the cutoffs
    fReg_dffits <- fReg_list$dffits_mat
    dff_test <- lapply(1:4, FUN = function(i) {
      # Check the upper bound and lower bound
      dffits_high_mat <- fReg_dffits > cutoffs_h[i]
      dffits_high_low <- fReg_dffits < cutoffs_l[i]
      dffits_out <- dffits_high_mat + dffits_high_low

      # Count the total number of points out
      dffits_out_tot <- colSums(dffits_out)

      return(dffits_out_tot)
    })

    # Combine the list results
    dff_results <- do.call(cbind, dff_test)
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
    s_idx <- (iter - 1)*n_obs + 1
    e_idx <- s_idx + n_obs - 1
    
    # Update data frame
    sim_results_temp[s_idx:e_idx, 4:7] <- dff_results
    sim_results_temp[s_idx:e_idx, 8] <- inf_temp
    
    # Update progress bar
    pb$tick()
    #################################################################
  }
  
  ## Save data frame to list:
  ###################################################################
  sim_results_list[[p]] <- sim_results_temp
  ###################################################################
}

# Create and save data frame from list
sim_results_df <- do.call(rbind, sim_results_list)
save(sim_results_df, 
     file = paste0("/Users/creutzml/Library/Mobile Documents",
                   "/com~apple~CloudDocs/Documents/Dissertation",
                   "/functional_data_analysis/sim_results",
                   "/functional_dffits/dffits_exp_mc_11_14_23.RData"))

