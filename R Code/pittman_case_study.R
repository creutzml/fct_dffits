#####################################################################
### Michael Creutzinger                                           ###
### Colorado State University                                     ###
### December 21st, 2023                                           ###
###                                                               ###
###   Testing out the case study example from Pittman.            ###
#####################################################################



### Starting steps:
#####################################################################
# Library packages
library(mvtnorm)
library(tidyverse)
library(ggplot2)
library(fdaconcur)
library(readr)

# Install ffscbExtra
devtools::install_github("creutzml/ffscbExtra")
library(ffscbExtra)

# Directories
functions_path <- file.path(here::here(), "R Code")
data_path <- file.path(here::here(), "Data")


# Load in Pittman's functions
source(file.path(functions_path, "pittman_data_gen_nonsmooth.R"))

# Load the data in
y_floods <- read_csv(file = file.path(data_path, 
                                      "FinalYtL1star.csv")) %>%
  dplyr::select(-`...1`) %>%
  dplyr::rename("August_1995" = V1, 
                "February_1998" = V2, 
                "March_2003" = V3, 
                "May_2003" = V4, 
                "September_2004" = V5, 
                "March_2007" = V6, 
                "February_2010" = V7, 
                "May_2013" = V8, 
                "November_2018" = V9, 
                "February_2020" = V10)

y_floods_mat <- as.matrix(y_floods)

x_floods <- read_csv(file = file.path(data_path, 
                                      "FinalXtL1star.csv")) %>%
  dplyr::rename("August_1995" = V1, 
                "February_1998" = V2, 
                "March_2003" = V3, 
                "May_2003" = V4, 
                "September_2004" = V5, 
                "March_2007" = V6, 
                "February_2010" = V7, 
                "May_2013" = V8, 
                "November_2018" = V9, 
                "February_2020" = V10)

x_floods_mat <- as.matrix(x_floods)

Oct15CongHt <- read_table(file = file.path(
  data_path, 
  "Oct15CongHt.txt"
))

# Plot of realigned river heights
# matplot(x_floods_mat, 
#         type = "l", 
#         main = expression(LAL[1]*" Selected Events"), 
#         ylab = "River Height (ft)", 
#         xlab = "Index")

## Ggplot2 plot instead of matplot
x_floods_gg <- x_floods %>%
  dplyr::mutate(Idx = 1:2000) %>%
  pivot_longer(-Idx, names_to = "Year", values_to = "river_height")

ggplot() + 
  geom_line(aes(x = Idx, y = river_height, color = Year), 
            data = x_floods_gg) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 16), 
        plot.margin = unit(c(.5, 1, 1, 1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Index", 
       y = "River Height (ft)", 
       color = "Year") +
  guides(color = "none")
#####################################################################



### Repeating Pittman's analysis
#####################################################################
# Packages
library(fdaconcur)
library(fda)
library(pracma)
library(lubridate)
library(MASS) #for mvrnorm
library(corpcor) #for make positve definite

# Functions from Pittman
file.sources <- list.files(path = file.path(main_dir, 
                                            "code/functional_dffits",
                                            "pittman_code"), 
                           full.names = T)
sapply(file.sources, source, .GlobalEnv)

#Calculates "B0t","B1t","stderr.B0t","stderr.B1t","Residuals",
# "AdjHatMatrix","MeanHat","FullyHat"
BetaOnly<-function(xData,yData,BasisNum=11, 
                   BetaOnly=FALSE,lambda=10^(-1)) {
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/n)^2, 0), rangeval=gaitrange)
  gaitbasis <- create.fourier.basis(gaitrange, nbasis=BasisNum)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  
  gaitSmooth = smooth.basisPar(gaittime, mygaitExp, gaitbasis, 
                               Lfdobj=harmaccelLfd, lambda=lambda) 
  #Could Edit this^ part Original played around no different
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] 
  # Response Stuff This seems like a more detailed version of what 
  #we get from the other textbook.
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd)
  betafdPar = fdPar(gaitbasis, harmaccelLfd)
  betalist  = list(const=betafdPar, cong=betafdPar)
  
  
  gaitRegress= fRegress(cedfd, xfdlist, betalist)
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine) #B0
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term B1
  
  Intercept<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  Slope<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  
  if(BetaOnly==FALSE){
    # MLC (12/21/23): needs to be FD object
    cedhatfd = gaitRegress$yhatfd 
    #1 to 2000 compare this to yHat.t.mat
    cedhatmat = eval.fd(gaittime, cedhatfd)  
    #The 2 represents cedar (response)
    resmat. = mygaitExp[,,2] - cedhatmat 
    SigmaE = cov(t(resmat.))
    
    cedfinemat   = eval.fd(gaitfine, cedfd)
    cedmeanvec   = eval.fd(gaitfine, mean.fd(cedfd))
    cedhatfinemat= eval.fd(gaitfine, cedhatfd)
    #finer grid than the . needed for computation below
    resmat        = cedfinemat - cedhatfinemat 
    #Says the number of flood events
    ncurve        = dim(mygaitExp)[2] 
    resmat0 = cedfinemat - cedmeanvec %*% matrix(1,1,ncurve)
    SSE0 = apply((resmat0)^2, 1, sum)
    SSE1 = apply(resmat^2, 1, sum)
    ced.R2 = (SSE0-SSE1)/SSE0
    #And this is what I had plotted before exactly.
    ylim2=c(0, max(congCoef, ced.R2)) 
    gaitbasismat = eval.basis(gaitfine, gaitbasis)
    #Needed to get error
    y2cMap = gaitSmooth$y2cMap 
    
    cedhatfd = gaitRegress$yhatfd 
    #1 to 2000 compare this to yHat.t.mat
    cedhatmat = eval.fd(gaittime, cedhatfd)
    #The 2 represents cedar (response)
    resmat. = mygaitExp[,,2] - cedhatmat 
    SigmaE = cov(t(resmat.))
    
    fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                             y2cMap=y2cMap, SigmaE=SigmaE)
    
    fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
    betastderrlist = fRegressList2$betastderrlist
    
    error.B0<-predict(betastderrlist[[1]], gaitfine)
    error.B1<-predict(betastderrlist[[2]], gaitfine)
    #This gives the standard errors for B0t and B1t to use for DFBETAS
    stderr.B0<-approx(x=1:length(error.B0), y=error.B0, n = n)$y
    stderr.B1<-approx(x=1:length(error.B1), y=error.B1, n = n)$y
    
    
    #Adjust to be nRow Ncol
    meanHat_i<-c()
    #adds the column of 1's
    data1<-cbind(rep(1,nPred),t(xData)) 
    Hatiit<-matrix(NA,nrow=nPred,ncol=n)
    for(t in 2:(n+1)){
      X<-data1[,c(1,t)]
      Hat<-X%*%solve(t(X)%*%X)%*%t(X)
      Hatiit[,t-1]<-diag(Hat) 
      # Each row is an event. The corresponding 
      # column is the Hii of that event at each timepoint
    }
    meanHat_i<-c(meanHat_i,apply(Hatiit,1,mean)) #will return an 
    # average Hii for each of the 10 sampled curve
    
    out<-list(Intercept,Slope,stderr.B0,stderr.B1,resmat.,
              Hatiit,meanHat_i,cedhatmat)
    names(out)<-c("B0t","B1t","stderr.B0t","stderr.B1t",
                  "Residuals","AdjHatMatrix","MeanHat","FullyHat")
  }
  else {out<-list(Intercept,Slope)
  names(out)=c("B0t","B1t")
  }
  return(out)
}

# Smooths the river flood data
BetaOnly1000<-function(xData,
                       yData,
                       BasisNum=11, 
                       BetaOnly=FALSE,
                       lambda=10^(-1)) {
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/n)^2, 0), rangeval=gaitrange)
  gaitbasis <- create.fourier.basis(gaitrange, nbasis=BasisNum)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  #Could Edit this part Original played around no different
  gaitSmooth = smooth.basisPar(gaittime, mygaitExp, gaitbasis, 
                               Lfdobj=harmaccelLfd, lambda=lambda) 
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  #Predictor Stuff
  congfd  = gaitfd[,1] 
  congfd_eval <- eval.fd(1:1000, congfd)
  
  # Response Stuff This seems like a more detailed version of what 
  # we get from the other textbook.
  cedfd = gaitfd[,2] 
  cedfd_eval <- eval.fd(1:1000, cedfd)
  return(list(congfd_eval, cedfd_eval))
}

smooth_floods <- BetaOnly1000(xData = x_floods_mat, 
                              yData = y_floods_mat)

FinalXtL1star <- x_floods_mat
FinalYtL1star <- y_floods_mat
FullResults<-BetaOnly(FinalXtL1star,FinalYtL1star,BasisNum=11)
B0_i_<-matrix(NA,
              nrow=nrow(FinalXtL1star),
              ncol=ncol(FinalXtL1star))
B1_i_<-matrix(NA,
              nrow=nrow(FinalXtL1star),
              ncol=ncol(FinalXtL1star))
stderrB1_i_2<-matrix(NA,
                     nrow=nrow(FinalXtL1star),
                     ncol=ncol(FinalXtL1star))
stderrB0_i_2<-matrix(NA,
                     nrow=nrow(FinalXtL1star),
                     ncol=ncol(FinalXtL1star))

#Calculates "B0t","B1t","stderr.B0t","stderr.B1t","Residuals",
# "AdjHatMatrix","MeanHat","FullyHat"
#With each of the i events left out using the same number of Fourier 
# basis functions as full run with all 10
for(i in 1:10){
  BasisToUse<-11
  tempBetas<-BetaOnly(FinalXtL1star[,-i],FinalYtL1star[,-i],
                      BetaOnly=FALSE,BasisNum=BasisToUse)
  B0_i_[,i]<-tempBetas$B0t
  B1_i_[,i]<-tempBetas$B1t
  stderrB0_i_2[,i]<-tempBetas$stderr.B0t
  stderrB1_i_2[,i]<-tempBetas$stderr.B1t
}

DFBETAS011<-matrix(NA,
                   nrow=nrow(FinalXtL1star),
                   ncol=ncol(FinalXtL1star))
DFBETAS111<-matrix(NA,
                   nrow=nrow(FinalXtL1star),
                   ncol=ncol(FinalXtL1star))
DFBETA011<-matrix(NA,
                  nrow=nrow(FinalXtL1star),
                  ncol=ncol(FinalXtL1star))
DFBETA111<-matrix(NA,
                  nrow=nrow(FinalXtL1star),
                  ncol=ncol(FinalXtL1star))
for(i in 1:10){
  DFBETAS011[,i]<-(FullResults$B0t-B0_i_[,i])/stderrB0_i_2[,i]
  DFBETAS111[,i]<-(FullResults$B1t-B1_i_[,i])/stderrB1_i_2[,i]
  DFBETA011[,i]<-(FullResults$B0t-B0_i_[,i])
  DFBETA111[,i]<-(FullResults$B1t-B1_i_[,i])
}


##DFFITS
yHat_mini<-matrix(NA,
                  nrow=nrow(FinalXtL1star),
                  ncol=ncol(FinalXtL1star))
for(i in 1:10){
  yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,i]
}

DFFIT<-FullResults$FullyHat-yHat_mini

#For MSE_mini minus i we want the true-predicted at each t
#for all 9 events when the 10th event is left out
MSE_mini_t<-matrix(NA,
                   nrow=nrow(FinalXtL1star),
                   ncol=ncol(FinalXtL1star))
yHatsNoi<-matrix(NA,
                 nrow=nrow(FinalXtL1star),
                 ncol=ncol(FinalXtL1star)-1)
n=1:10
for(i in n){
  k=1
  for(j in n[-i]){
    #goes through all but i using the betas with i left out
    yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,j] 
    k=k+1
  }
  #Get an MSE across the other events when event i is left out
  MSE_mini_t[,i]<-apply((FinalYtL1star[,-i]-yHatsNoi)^2,1,mean) 
}
DFFITS<-matrix(NA,nrow=nrow(DFFIT),ncol=ncol(DFFIT))
for(t in 1:nrow(DFFIT)) {
  DFFITS[t,] <- DFFIT[t,]/
    sqrt(MSE_mini_t[t,]*FullResults$AdjHatMatrix[,t])
}

# Table of average DFFITS values
knitr::kable(cbind(names(x_floods), round(colMeans(abs(DFFITS)), 2)), 
             format = "latex", 
             booktabs = T)

# Make a plot of the DFFITS estimated by Pittman
pitt_dffits_gg <- DFFITS %>%
  as.data.frame() %>%
  dplyr::mutate(Index = 1:2000) %>%
  pivot_longer(-Index, names_to = "Year", values_to = "river_height")

ggplot() + 
  geom_line(aes(x = Index, y = river_height, color = Year), 
            data = pitt_dffits_gg) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 16), 
    plot.margin = unit(c(.5, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "Index", 
       y = "DFFITS(t)", 
       color = "Year") +
  guides(color = "none")

## Estimate the bootstrapped distribution cutoff for alpha = 0.025
pitt_bootstrap <- absdffit_bootsamedata(alpha = 0.025, 
                                        xdata = x_floods_mat, 
                                        ydata = y_floods_mat, 
                                        outlier = 10, 
                                        nboot = 100)

## Pitt Influential
pitt_inf <- colMeans(abs(DFFITS)) > pitt_bootstrap$BootPercentiles
names(x_floods)[pitt_inf]
#####################################################################



## Creutzinger analysis
#####################################################################
# Run the concurrent regression with X and Y from the floods data
flood_regress <- ffscbExtra::fRegress_concurrent(
  y_mat = y_floods_mat, 
  x_array = x_floods_mat
)

# multivariate Student's t distribution quantiles
load(file = file.path(data_path, "mvt_cutoffs_sim.RData"))

# Plot of the dffits
# matplot(flood_regress$dffits_mat, type = "l", 
#         ylab = "DFFITS(t)", 
#         xlab = "Index")
# abline(h = c(mvt_cutoffs$t_cuts[4], -mvt_cutoffs$t_cuts[4]), 
#        lty = "dashed")
creutz_dffits_gg <- flood_regress$dffits_mat %>%
  as.data.frame() %>%
  dplyr::mutate(Index = 1:2000) %>%
  pivot_longer(-Index, names_to = "Year", values_to = "river_height")

ggplot() + 
  geom_line(aes(x = Index, y = river_height, color = Year), 
            data = creutz_dffits_gg) +
  geom_hline(yintercept = c(mvt_cutoffs$t_cuts[34], 
                            - mvt_cutoffs$t_cuts[34]), 
             linetype = "dashed") +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw(base_size = 18) +
  theme(#text = element_text(size = 16), 
    plot.margin = unit(c(.5, 1, 1, 1), "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  labs(x = "Index", 
       y = "DFFITS(t)", 
       color = "Year") +
  guides(color = "none")

## Which observations are influential?
# Calculate the MVT cutoff
# a_val <- 0.005
# t_df <- 7
# scalar <- sqrt(2/(10 - 2))
# n_sp <- 10
# t_quants <- qmvt(p = 1 - a_val,
#                  df = t_df,
#                  tail = "lower.tail",
#                  sigma = diag(t_df/(t_df - 2),
#                               nrow = n_sp))
# cuts <- scalar*t_quants$quantile

inf_idx <- which(colSums(abs(flood_regress$dffits_mat) > 
                           mvt_cutoffs$t_cuts[4]) > 0)
names(x_floods)[inf_idx]
#####################################################################