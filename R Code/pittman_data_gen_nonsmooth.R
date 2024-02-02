#This is the Code for |DFBETAS0| Simulation and does 1 iteration

# For running on own computer
# # Main directory:
# main_dir <- file.path("/Users/creutzml/Library/Mobile Documents",
#                       "/com~apple~CloudDocs/Documents/Dissertation",
#                       "/functional_data_analysis/code")
# 
# # Source the original functions these build on
# source(file.path(main_dir, 
#                  "functional_dffits/pittman_data_gen.R"))

# Running as a bash script
source(paste0("/scratch/alpine/creutzml@colostate.edu/dffits/code",
              "/pittman_data_gen.R"))

##Can manually change things like BasisNum, the B0fun, B1fun, 
##B1funout, functions, and length
library(fda)
library(MASS)
library(corpcor) 

##Calculates dffits for each t and each observation
BestDFFITS_fun_ns<-function(Xtmat,Ytmat){
  # Fit concurrent model without smoothing
  fReg_list <- fRegress_concurrent(y_mat = newData$Ydata,
                                   x_array = newData$Xdata)
  
  # Return the dffits
  return(list(dffits = fReg_list$dffits_mat))
}

#generate data before using this function. 
# This function bootstraps given data
#Give it the outliar number from the data
# (MLC 1/2/24): changed function so that alpha in Pittman is 0.5, one
# that performs well for small sample sizes. Then, alpha will 
# correspond to the percentile, rather than their parameter
absdffit_bootsamedata_ns<-function(alpha, 
                                   xdata = sim_data$Xdata, 
                                   outlier = sim_data$OutlierNumber,
                                   ydata = sim_data$Ydata, 
                                   nboot = 100) {
  
  #Say the number of observations we have
  N = ncol(xdata) 
  
  #Now run the concurrent model and calculate mean |dffits| for each
  # of the observation
  #pointwise DFFITS for each observation
  raw_dffits<-BestDFFITS_fun_ns(Xtmat = xdata,
                                Ytmat = ydata) 
  mean_abs_dffits<-apply(abs(raw_dffits$dffits),2,mean)
  adjusted_obs_absdffits<-mean_abs_dffits[outlier]
  theta<-as.vector((1/mean_abs_dffits)^0/
                     sum((1/mean_abs_dffits)^0))
  #Need to set up the bootstrap next SEE OrdinaryRegressionDistribution 
  # for Details
  #Set blank objects here
  boot_meanabsdffits<-c()
  for(i in 1:nboot){
    #resamples the needed data with the given probability 
    index_samp<-sample(1:N,N,replace=TRUE,prob = theta) 
    temp_dffits<-BestDFFITS_fun_ns(Xtmat = xdata[,index_samp],
                                   Ytmat = ydata[,index_samp])$dffits
    boot_meanabsdffits<-c(boot_meanabsdffits,
                          apply(abs(temp_dffits),2,mean))
  }
  #returns percentiles from the N*nboot |dffits| calculated
  boot_percentiles<-quantile(boot_meanabsdffits,1-alpha)
  percentile <- ecdf(boot_meanabsdffits)
  adjusted_point_percentile<- percentile(mean_abs_dffits[outlier])
  out<-list(boot_percentiles,adjusted_point_percentile,mean_abs_dffits)
  #And returns the percentile of this list that the observed |dffits| 
  #falls into
  names(out)<-c("BootPercentiles","AdjPercentile","AllObsMeanAbsDFFITS") 
  return(out)
}
