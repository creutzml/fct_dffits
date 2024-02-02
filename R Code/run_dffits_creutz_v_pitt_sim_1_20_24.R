#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### December 13th, 2023                                           ###
###                                                               ###
###   Simulation comparing the results of Pittman's DFFITS with   ###
### my own. Pittman's method is executed on smoothed and non      ###
### data, for a more complete comparison.                         ###
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


# Read in command arguments from bash files
args <- commandArgs(trailingOnly=TRUE )
if (length(args) != 6) {
  stop("Looping parameters failed to initialize in 'commandArgs'.")
}

# Make new objects from input
sim_mod <- as.character(args[1])
n_obs <- as.numeric(args[2])
lambdas <- as.numeric(args[3])
n_infs <- as.numeric(args[4])
iter <- as.numeric(args[5])

# Test case
# {sim_mod = "1"; n_obs = 10; lambdas = 1; n_infs = 1; iter = 2}

# Trying to run 8 iterations in one, so that all can be submitted at
# the same time on SLURM alpine
s_iter <- (iter-1)*8 + 1
e_iter <- iter*8
iterations <- s_iter:e_iter


# File name for saved results later
file_name <- paste0(
  "Simulation_Model", sim_mod, "_n_obs=", n_obs, "_lambda=", 
  lambdas, "_n_infs=", n_infs, "_iters=", s_iter,"_to_", e_iter,
  "_a=0_T=100.RData"
)

base::cat(file_name, "\n")

# Directories
scratch_path <- "/scratch/alpine/creutzml@colostate.edu"
dffits_code_path <- file.path(scratch_path, "dffits/code")
dffits_results_path <- file.path(scratch_path, "dffits/sim_results")
#####################################################################


# Main directory:
source(file.path(dffits_code_path, 
                 "make_regression_band_mc.R"))
source(file.path(dffits_code_path, 
                 "make_band_mc.R"))
source(file.path(dffits_code_path, 
                 "confidence_band_mc.R"))
source(file.path(dffits_code_path, 
                 "make_concurrent_regression_band_mc_10_12_23.R"))
source(file.path(dffits_code_path, 
                 "fRegress_concurrent.R"))
source(file.path(dffits_code_path, 
                 "pittman_data_gen_nonsmooth.R"))

# Load in the multivariate cutoffs
load(file.path(dffits_code_path, 
               "mvt_cutoffs_sim_1_19_23.RData"))
#####################################################################




## Run a simulation for functional dffits
#####################################################################
# Match the parameter to the correct data function
generate_data_temp <- match.fun(paste0("GenerateFunctionalDataOut", 
                                       sim_mod))

## Parameters to loop over:
n_sp <- 100
alphas <- c(0.10, 0.05, 0.025, 0.005)
methods <- c("Creutz", "Creutz_Smooth", "Pittman_Smooth", "Pittman")

# Length of data frame for each simulation
df_length <- length(alphas)*length(methods)

# Start the simulation
#####################################################################
# Empty data frame for results
sim_results_temp <- data.frame(
  sim_mod = rep(sim_mod, df_length),
  method = rep(methods, each = length(alphas)),
  n = rep(n_obs, df_length),
  lambda = rep(lambdas, df_length),
  n_inf = rep(n_infs, df_length),
  iter = rep(iterations, each = df_length),
  alpha_val = rep(alphas, length(methods)),
  n_tp = vector("numeric", df_length), 
  n_fp = vector("numeric", df_length), 
  n_tn = vector("numeric", df_length), 
  n_fn = vector("numeric", df_length),
  run_time = vector("numeric", df_length)
)

for (loop_iter in 1:8) {
  ## Generate data:
  #################################################################
  # Number of observations and sampling points, 
  n_obs <- n_obs
  grid <- make_grid(n_sp, rangevals = c(0, 1))
  
  # Generate new functional data using Pittman's approach:
  newData <- generate_data_temp(
    N = n_obs, 
    length = n_sp, 
    lambda = lambdas,
    n_inf = n_infs
  )
  
  # Which one is influential?
  if (lambdas != 1) {
    inf_temp <- (1:n_obs %in% newData$OutlierNumber)
  } else {
    inf_temp <- rep(FALSE, n_obs)
  }
  #################################################################
  
  
  ## Functional DFFITS (Michael)
  #################################################################
  # Start the timer
  creutz_s_time <- Sys.time()
  
  ## And now, if we fit with concurrent regression and our bands
  fReg_list <- fRegress_concurrent(y_mat = newData$Ydata,
                                   x_array = newData$Xdata)
  
  ## Pointwise cutoff:
  # Df and mean scalar quantity
  t_df <- n_obs - 3
  scalar <- sqrt(2/(n_obs - 2))
  
  # t quantiles
  cutoffs <- mvt_cutoffs$t_cuts[mvt_cutoffs$n_obs == n_obs &
                                  mvt_cutoffs$n_sp == n_sp]
  
  # Save the dffits as an object
  fReg_dffits <- fReg_list$dffits_mat
  
  # Clock it
  creutz_f_time <- as.numeric(difftime(Sys.time(), 
                                       creutz_s_time,
                                       units = "secs"))
  
  # compare the dffits to the cutoffs
  dff_test <- lapply(1:length(alphas), FUN = function(i) {
    
    # Get the alpha value and cutoff
    alpha_temp <- alphas[i]
    t_cut <- cutoffs[i]
    
    # Check the upper bound and lower bound
    dffits_out <- colSums(abs(fReg_dffits) > t_cut) > 0 
    
    n_tp <- sum(dffits_out & inf_temp)
    n_fp <- sum(dffits_out & !inf_temp)
    n_tn <- sum(!dffits_out & !inf_temp)
    n_fn <- sum(!dffits_out & inf_temp)
    
    # Add alpha value to dffits_diag
    dffits_diag <- c(alpha_temp, 
                     n_tp, n_fp, n_tn, n_fn, 
                     creutz_f_time)
    
    return(dffits_diag)
  })
  
  # Combine the list results
  dff_results <- do.call(rbind, dff_test)
  ################################################################
  
  
  # Set up a tryCatch to avoid any issues that occur with smoothing
  tryCatch({
    
    ## Functional DFFITS with smoothing (Pittman)
    ###############################################################
    # Start the timer
    pitt_smooth_s_time <- Sys.time()
    
    # Run Pittman's bootstrap as is (with smoothing involved)
    pittman_smooth_boot <- absdffit_bootsamedata(
      alpha = alphas, 
      xdata = newData$Xdata, 
      ydata = newData$Ydata, 
      outlier = newData$OutlierNumber,
      nboot = 100
    )
    
    # Clock it
    pitt_smooth_f_time <- as.numeric(difftime(Sys.time(), 
                                              pitt_smooth_s_time,
                                              units = "secs"))
    
    # Identify which observations should be marked out
    pittman_smooth_out <- lapply(1:length(alphas), FUN = function(i) {
      
      # Get the alpha value
      alpha_temp <- alphas[i]
      
      # Check the upper bound and lower bound
      dffits_out <- pittman_smooth_boot$AllObsMeanAbsDFFITS >= 
        pittman_smooth_boot$BootPercentiles[i]
      
      # Calculate true positives, false positives, etc. 
      n_tp <- sum(dffits_out & inf_temp)
      n_fp <- sum(dffits_out & !inf_temp)
      n_tn <- sum(!dffits_out & !inf_temp)
      n_fn <- sum(!dffits_out & inf_temp)
      
      # Add alpha value to dffits_diag
      dffits_diag <- c(alpha_temp, 
                       n_tp, n_fp, n_tn, n_fn, 
                       pitt_smooth_f_time)
      
      # Return the resulting vector
      return(dffits_diag)
    })
    
    # Combine the list results
    pittman_smooth_results <- do.call(rbind, pittman_smooth_out)
    ###############################################################
    
    
    ## Functional DFFITS after smoothing, tested with Creutz
    ###############################################################
    # Save the dffits as an object
    fReg_smooth_dffits <- pittman_smooth_boot$RawDFFITS$altdffits
    
    # compare the dffits to the cutoffs
    dff_smooth_test <- lapply(1:length(alphas), FUN = function(i) {
      
      # Get the alpha value and cutoff
      alpha_temp <- alphas[i]
      t_cut <- cutoffs[i]
      
      # Check the upper bound and lower bound
      dffits_out <- colSums(abs(fReg_smooth_dffits) > t_cut) > 0 
      
      n_tp <- sum(dffits_out & inf_temp)
      n_fp <- sum(dffits_out & !inf_temp)
      n_tn <- sum(!dffits_out & !inf_temp)
      n_fn <- sum(!dffits_out & inf_temp)
      
      # Add alpha value to dffits_diag
      dffits_diag <- c(alpha_temp, 
                       n_tp, n_fp, n_tn, n_fn, 
                       creutz_f_time)
      
      return(dffits_diag)
    })
    
    # Combine the list results
    dff_smooth_results <- do.call(rbind, dff_smooth_test)
    ###############################################################
    
    # Handling an error if it occurs
  }, error = function(e) {
    
    # Set the values to NA if the smoothing did not work correctly
    pittman_smooth_results <<- matrix(data = c(alphas, rep(NA, 20)), 
                                      nrow = 4)
    dff_smooth_results <<- matrix(data = c(alphas, rep(NA, 20)), 
                                  nrow = 4)
  })
  
  
  ################################################################
  
  
  
  ## Functional DFFITS without smoothing (Pittman)
  #################################################################
  # Start the timer
  pitt_nonsmooth_s_time <- Sys.time()
  
  # Run Pittman's bootstrap without smoothing
  pittman_nonsmooth_boot <- absdffit_bootsamedata_ns(
    alpha = alphas, 
    xdata = newData$Xdata, 
    ydata = newData$Ydata, 
    outlier = newData$OutlierNumber,
    nboot = 100
  )
  
  # Clock it
  pitt_nonsmooth_f_time <- as.numeric(difftime(Sys.time(), 
                                               pitt_nonsmooth_s_time,
                                               units = "secs"))
  
  # Identify which observations should be marked out
  pittman_nonsmooth_out <- lapply(
    1:length(alphas), 
    FUN = function(i) {
      
      # Get the alpha value
      alpha_temp <- alphas[i]
      
      # Check the upper bound and lower bound
      dffits_out <- pittman_nonsmooth_boot$AllObsMeanAbsDFFITS >= 
        pittman_nonsmooth_boot$BootPercentiles[i]
      
      n_tp <- sum(dffits_out & inf_temp)
      n_fp <- sum(dffits_out & !inf_temp)
      n_tn <- sum(!dffits_out & !inf_temp)
      n_fn <- sum(!dffits_out & inf_temp)
      
      # Add alpha value to dffits_diag
      dffits_diag <- c(alpha_temp, 
                       n_tp, n_fp, n_tn, n_fn, 
                       pitt_nonsmooth_f_time)
      
      return(dffits_diag)
    }
  )
  
  # Combine the list results
  pittman_nonsmooth_results <- do.call(rbind, pittman_nonsmooth_out)
  #################################################################
  
  
  
  ## Alter data and save the results
  #################################################################
  # Indices for rows
  s_idx <- (loop_iter - 1)*length(alphas)*length(methods) + 1
  e_idx <- s_idx + length(alphas)*length(methods) - 1
  
  # Update data frame
  sim_results_temp[s_idx:e_idx, 7:12] <- rbind(
    dff_results, 
    dff_smooth_results,
    pittman_smooth_results, 
    pittman_nonsmooth_results
  )
}
# Calculate sensitivity, specificity, accurracy, and mcc
sim_results_temp2 <- sim_results_temp %>%
  mutate(sensitivity = n_tp/(n_tp + n_fn), 
         specificity = n_tn/(n_tn + n_fp), 
         accuracy = (n_tp + n_tn)/(n), 
         mcc = (n_tp*n_tn - n_fp*n_fn)/
           (sqrt((n_tp + n_fp)*(n_tp + n_fn)*
                   (n_tn + n_fp)*(n_tn + n_fn))))

save(sim_results_temp2, 
     file = file.path(dffits_results_path, file_name))
###################################################################





